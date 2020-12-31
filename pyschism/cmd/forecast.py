from argparse import Namespace
from datetime import datetime, timedelta
import json
import os
import pathlib
import shutil
import subprocess

from dateutil.parser import parse as parse_datetime
import numpy as np  # type: ignore[import]
# import psutil  # type: ignore[import]
import pytz

from pyschism.forcing import GlobalForecastSystem as GFS
from pyschism.forcing import Tides
from pyschism.mesh import Vgrid
from pyschism.server.slurm import SlurmConfig
from pyschism.driver import ModelDriver
from pyschism.domain import ModelDomain
from pyschism.forcing.atmosphere import NWS2
from pyschism.param.schout import SurfaceOutputVars


class Args:

    def __set__(self, obj, args: Namespace):
        obj.__dict__['args'] = args
        if args.action == 'init':
            if obj.config_file.is_file() and args.overwrite is False:
                raise IOError(
                    f'The given directory {str(obj.project_directory)} has '
                    'been initialized previously. Use `pyschism forecast '
                    '--overwrite init [...]` to allow overwrite of previous '
                    'initialization options.')
            with open(obj.config_file, 'w') as fp:
                json.dump(obj.args.__dict__, fp, indent=4)
            obj.generate_forecast()
        elif args.action == 'update':
            raise NotImplementedError('update action')
        else:
            raise Exception(f"Unknown argument action {args.action}.")

    def __get__(self, obj, val):
        return obj.__dict__['args']


class ConfigFile:

    def __get__(self, obj, val):
        config_file = obj.__dict__.get('config_file')
        if config_file is None:
            config_file = obj.project_directory / 'config.json'
            obj.__dict__['config_file'] = config_file
        return config_file


class TargetDatetime:

    def __get__(self, obj, val):
        target_datetime = obj.__dict__.get('target_datetime')
        if target_datetime is None:
            utcnow = pytz.timezone('UTC').localize(datetime.utcnow())
            localnow = utcnow.astimezone(pytz.timezone(obj.args.timezone))
            nowcast_cycle = int(obj.forecast_interval * np.floor(
                localnow.hour/obj.forecast_interval))  # % 24
            target_datetime = datetime(
                localnow.year, localnow.month, localnow.day, nowcast_cycle,
                tzinfo=pytz.timezone(obj.args.timezone))
            obj.__dict__['target_datetime'] = target_datetime
        return target_datetime


class Timezone:

    def __get__(self, obj, val):
        timezone = obj.__dict__.get('timezone')
        if timezone is None:
            timezone = obj.args.timezone
            if timezone != 'UTC':
                raise NotImplementedError(
                    'Only UTC timezone allowed for now.')
            obj.__dict__['timezone'] = obj.args.timezone
        return timezone


class ForecastInterval:

    def __get__(self, obj, val):
        forecast_interval = obj.__dict__.get('forecast_interval')
        if forecast_interval is None:
            if obj.args.forecast_interval != 24:
                raise NotImplementedError('Forecast interval only at 24 hours '
                                          'for now.')
            obj.__dict__['forecast_interval'] = obj.args.forecast_interval
        return obj.__dict__['forecast_interval']


class ProjectDirectory:

    def __get__(self, obj, val):
        project_directory = obj.__dict__.get('project_directory')
        if project_directory is None:
            project_directory = pathlib.Path(obj.args.project_directory)
            project_directory.mkdir(parents=True, exist_ok=True)
            obj.__dict__['project_directory'] = project_directory
        return project_directory


class ForecastsDirectory:

    def __get__(self, obj, val):
        forecasts_directory = obj.__dict__.get('forecasts_directory')
        if forecasts_directory is None:
            forecasts_directory = obj.project_directory / 'forecasts'
            forecasts_directory.mkdir(exist_ok=True)
            obj.__dict__['forecasts_directory'] = forecasts_directory
        return forecasts_directory


class ColdstartDirectory:

    def __get__(self, obj, val):
        coldstart_directory = obj.__dict__.get('coldstart_directory')
        if coldstart_directory is None:
            timestamp = str(obj.target_datetime).replace(' ', 'T')
            coldstart_directory = obj.project_directory / \
                'coldstart' / f'{timestamp}'
            coldstart_directory.mkdir(
                parents=True, exist_ok=True)
            obj.__dict__['coldstart_directory'] = coldstart_directory
        return coldstart_directory


class StaticFilesDirectory:

    def __get__(self, obj, val):
        static_files = obj.__dict__.get('static_files')
        if static_files is None:
            static_files = obj.project_directory / 'static'
            static_files.mkdir(exist_ok=True)
            obj.__dict__['static_files'] = static_files
        return static_files


class TargetOutputDirectory:

    def __get__(self, obj, val):
        target_output_directory = obj.__dict__.get('target_output_directory')
        if target_output_directory is None:
            timestamp = str(obj.target_datetime).replace(' ', 'T')
            target_output_directory = obj.forecasts_directory / f'{timestamp}'
            target_output_directory.parent.mkdir(exist_ok=True)
            target_output_directory.mkdir(exist_ok=True)
            obj.__dict__['target_output_directory'] = target_output_directory
        return target_output_directory


class ColdstartDomainDescriptor:

    def __get__(self, obj, val):
        model_domain = obj.__dict__.get('coldstart_domain')
        if model_domain is None:
            model_domain = ModelDomain.open(
                obj.hgrid_path, obj.fgrid_path,
                # TODO: Vgrid.open() raises NotImplementedError
                # obj.vgrid_path,
                hgrid_crs=obj.args.hgrid_crs,
                fgrid_crs=obj.args.fgrid_crs
            )
            if obj.tides is not None:
                model_domain.add_boundary_condition(obj.tides)
            obj.__dict__['coldstart_domain'] = model_domain
        return model_domain


class ColdstartModelDriver:

    def __get__(self, obj, val):
        coldstart = obj.__dict__.get('coldstart_driver')
        if coldstart is None:
            start_date = obj.target_datetime - timedelta(
                days=obj.args.spinup_days)
            rnday = obj.target_datetime - start_date
            coldstart = ModelDriver(
                    model_domain=obj.coldstart_domain,
                    dt=obj.args.timestep,
                    rnday=rnday,
                    dramp=(2./3.)*rnday,
                    start_date=start_date,
                    nhot_write=True,
                    # ibc=obj.ibc,
                    # drampbc=obj.drampbc,
                    server_config=obj.server_config,
                    )
            obj.__dict__['coldstart_driver'] = coldstart
        return coldstart


class HotstartDomainDescriptor:

    def __get__(self, obj, val):
        model_domain = obj.__dict__.get('hotstart_domain')
        if model_domain is None:
            model_domain = ModelDomain.open(
                obj.hgrid_path,
                obj.fgrid_path,
                # TODO: Vgrid.open() raises NotImplementedError
                # obj.vgrid_path,
                hgrid_crs=obj.args.hgrid_crs,
                fgrid_crs=obj.args.fgrid_crs)
            if obj.tides is not None:
                model_domain.add_boundary_condition(obj.tides)
            if obj.nws2 is not None:
                model_domain.set_atmospheric_forcing(obj.nws2)
            obj.__dict__['hotstart_domain'] = model_domain
        return model_domain


class HotstartModelDriver:

    def __get__(self, obj, val):
        hotstart_driver = obj.__dict__.get('hotstart_driver')
        if hotstart_driver is None:
            hotstart_driver = ModelDriver(
                model_domain=obj.hotstart_domain,
                dt=obj.args.timestep,
                rnday=obj.args.forecast_days,
                # ihfskip=obj.args.ihfskip,
                start_date=obj.target_datetime,
                # ibc=obj.ibc,
                # stations=obj.stations,
                nspool=obj.nspool,
                combine_hotstart=obj.previous_run_directory / 'outputs',
                server_config=obj.server_config,
                **obj.surface_outputs
                )
            obj.__dict__['hotstart_driver'] = hotstart_driver
        return hotstart_driver


class TidesDescriptor:

    def __get__(self, obj, val):
        tides = obj.__dict__.get('tides')
        if tides is None:
            if not obj.args.all_constituents \
                    and not obj.args.major_constituents \
                    and not obj.args.constituents:
                return
            else:
                tides = Tides(velocity=obj.args.bnd_vel)
                if obj.args.all_constituents:
                    tides.use_all()
                if obj.args.major_constituents:
                    tides.use_major()
                if obj.args.constituents:
                    for constituent in obj.args.constituents:
                        tides.use_constituent(constituent)
            obj.__dict__['tides'] = tides
        return tides


class SrcSnkDescriptor:

    def __get__(self, obj, val):
        return obj.__dict__.get('srcsnk')


class Sflux1:

    def __get__(self, obj, val):
        sflux_1 = obj.__dict__.get('sflux_1')
        if sflux_1 is None:
            for arg, value in obj.args.__dict__.items():
                if 'gdas' in arg:
                    if value is True:
                        raise NotImplementedError(
                            'GDAS product not implemented.')
                elif 'gfs' in arg:
                    if value is True:
                        if arg == 'gfs':
                            sflux_1 = GFS()
                        else:
                            sflux_1 = GFS(product=arg)
            obj.__dict__['sflux_1'] = sflux_1
        return sflux_1


class Nws2Descriptor:

    def __get__(self, obj, val):
        nws2 = obj.__dict__.get('nws2')
        if nws2 is None:
            if obj.sflux_1 is not None:
                nws2 = NWS2(
                    obj.sflux_1,
                    # self.sflux_2
                    )
            obj.__dict__['nws2'] = nws2
        return nws2


class PreviousRunDirectory:

    def __get__(self, obj, val):
        previous_run_directory = obj.__dict__.get('previous_run_directory')
        if previous_run_directory is None:

            if obj.args.action == 'init':
                # TODO SKIPPING COLDSTART ON INIT
                obj.generate_coldstart()
                previous_run_directory = obj.coldstart_directory

            else:
                dirs = os.listdir(obj.forecasts_directory)
                candidate_runs = list(filter(lambda x: x < obj.target_datetime,
                                      map(parse_datetime, dirs)))
                raise NotImplementedError('There is a candidate run.')

            obj.__dict__['previous_run_directory'] = previous_run_directory

        if previous_run_directory is None:
            raise NotADirectoryError

        return previous_run_directory


class HgridPath:

    def __get__(self, obj, val):
        hgrid_path = obj.__dict__.get('hgrid_path')
        if hgrid_path is None:
            hgrid_path = obj.static_files_directory / 'hgrid.gr3'
            if hgrid_path.is_file() and obj.args.overwrite is True:
                hgrid_path.unlink()
            if not hgrid_path.is_file():
                shutil.copy2(obj.args.hgrid, hgrid_path,
                             follow_symlinks=True)
            obj.__dict__['hgrid_path'] = hgrid_path
        return hgrid_path


class VgridPath:

    def __get__(self, obj, val):
        vgrid_path = obj.__dict__.get('vgrid_path')
        if vgrid_path is None:
            vgrid_path = obj.static_files_directory / 'vgrid.in'
            if vgrid_path.is_file():
                if obj.args.overwrite is True:
                    vgrid_path.unlink()
            if obj.args.vgrid is None:
                Vgrid().write(vgrid_path)
            else:
                shutil.copy2(
                    obj.args.vgrid, vgrid_path, follow_symlinks=True)
            obj.__dict__['vgrid_path'] = vgrid_path
        return vgrid_path


class FgridPath:

    def __get__(self, obj, val):
        fgrid_path = obj.__dict__.get('fgrid_path')
        if fgrid_path is None:
            if obj.args.fgrid_type == 'auto':
                fgrid = pathlib.Path(obj.args.fgrid)
                fgrid_path = obj.static_files_directory / fgrid.name
            else:
                fgrid_path = obj.static_files_directory / \
                    (obj.args.fgrid_type + '.gr3')
            if obj.args.overwrite is True:
                if fgrid_path.is_file():
                    os.remove(fgrid_path)
            shutil.copy2(obj.args.fgrid, fgrid_path, follow_symlinks=True)
            obj.__dict__['fgrid_path'] = fgrid_path
        return fgrid_path


class WindrotPath:

    def __get__(self, obj, val):
        windrot_path = obj.__dict__.get('windrot_path')
        if windrot_path is None:
            windrot_path = obj.static_files_directory / 'windrot_geo2proj.gr3'
            if obj.args.overwrite is True and windrot_path.exists():
                windrot_path.unlink()
            if not windrot_path.exists():
                obj.hotstart.nws._windrot = obj.hotstart_domain.hgrid
                obj.hotstart.nws._windrot.write(
                    windrot_path, overwrite=obj.args.overwrite)
            obj.__dict__['windrot_path'] = windrot_path
        return windrot_path


class ServerConfigDescriptor:

    def __init__(self):
        self.server_config = None

    def __get__(self, obj, val):

        if self.server_config is None:
            if obj.args.server_config == "slurm":
                kwargs = {
                    "account": obj.args.account,
                    "ntasks": obj.args.nproc,
                    "partition": obj.args.partition,
                    "walltime": obj.args.walltime,
                    "mail_type": obj.args.mail_type,
                    "mail_user": obj.args.mail_user,
                    "log_filename": obj.args.log_filename,
                    "modules": obj.args.modules,
                    # "schism_binary": obj.args.schism_binary,
                    "extra_commands": obj.args.extra_commands,
                    "launcher": obj.args.slurm_launcher,
                    "nodes": obj.args.slurm_nodes}
                # if obj.args.slurm_filename is not None:
                #     kwargs.update({"filename": obj.args.ntasks})
                # if obj.args.slurm_rundir is not None:
                #     kwargs.update({"run_directory": obj.args.slurm_rundir})
                # if obj.args.run_name is not None:
                #     kwargs.update({"run_name": obj.args.run_name})
                self.server_config = SlurmConfig(**kwargs)
        return self.server_config


class SurfaceOutputsDescriptor:

    surface_output_vars = SurfaceOutputVars()

    def __get__(self, obj, val):
        surface_outputs = obj.__dict__.get('surface_outputs')
        if surface_outputs is None:
            surface_outputs = {}
            outvars = []
            for vardata in self.surface_output_vars.values():
                for varname, _ in vardata:
                    outvars.append(varname)
            for key, val in obj.args.__dict__.items():
                if key in outvars and val is True:
                    surface_outputs[key] = val
            obj.__dict__['surface_outputs'] = surface_outputs
        return surface_outputs


class Nspool:

    def __get__(self, obj, val):
        nspool = obj.__dict__.get('nspool')
        if nspool is None:
            if obj.args.nspool is not None:
                if '.' in obj.args.nspool:
                    nspool = timedelta(hours=float(obj.args.nspool))
                else:
                    nspool = int(obj.args.nspool)
                obj.__dict__['nspool'] = nspool
        return nspool


class GenerateForecastCli:

    args = Args()
    config_file = ConfigFile()
    project_directory = ProjectDirectory()
    forecasts_directory = ForecastsDirectory()
    coldstart_directory = ColdstartDirectory()
    static_files_directory = StaticFilesDirectory()
    target_output_directory = TargetOutputDirectory()
    target_datetime = TargetDatetime()
    previous_run_directory = PreviousRunDirectory()
    srcsnk = SrcSnkDescriptor()
    coldstart_domain = ColdstartDomainDescriptor()
    hotstart_domain = HotstartDomainDescriptor()
    coldstart = ColdstartModelDriver()
    hotstart = HotstartModelDriver()
    tides = TidesDescriptor()
    nws2 = Nws2Descriptor()
    hgrid_path = HgridPath()
    fgrid_path = FgridPath()
    vgrid_path = VgridPath()
    windrot_path = WindrotPath()
    server_config = ServerConfigDescriptor()
    surface_outputs = SurfaceOutputsDescriptor()
    sflux_1 = Sflux1()
    nspool = Nspool()
    forecast_interval = ForecastInterval()

    def __init__(self, args: Namespace):
        self.args = args
        # self.run()

    def generate_coldstart(self):
        self._symlink_files(self.coldstart_directory,
                            self.coldstart_domain.ics)
        self.coldstart.write(self.coldstart_directory, hgrid=False,
                             vgrid=False, fgrid=False,
                             overwrite=self.args.overwrite)
        if self.args.skip_run is False:
            subprocess.check_call(
                ["make", "run"], cwd=self.coldstart_directory)

    def generate_forecast(self):
        self._symlink_files(self.target_output_directory,
                            self.hotstart_domain.ics,
                            windrot=True
                            )
        self.hotstart.write(self.target_output_directory, hgrid=False,
                            vgrid=False, fgrid=False, wind_rot=False,
                            overwrite=self.args.overwrite)

        if self.args.skip_run is False:
            subprocess.check_call(
                ["make", "run"], cwd=self.target_output_directory)

    def _symlink_files(self, target_directory, ics, windrot=False):
        hgrid_lnk = target_directory / 'hgrid.gr3'
        vgrid_lnk = target_directory / 'vgrid.in'
        fgrid_lnk = target_directory / f'{self.fgrid_path.name}'
        if self.args.overwrite is True:
            if hgrid_lnk.exists():
                hgrid_lnk.unlink()
            if vgrid_lnk.exists():
                vgrid_lnk.unlink()
            if fgrid_lnk.exists():
                fgrid_lnk.unlink()
        os.symlink(
            os.path.relpath(
                self.hgrid_path, target_directory), hgrid_lnk)
        os.symlink(
            os.path.relpath(
                self.vgrid_path, target_directory), vgrid_lnk)
        os.symlink(
            os.path.relpath(
                self.fgrid_path, target_directory), fgrid_lnk)
        if ics == 2:
            self._symlink_hgridll(target_directory)
        if windrot is True:
            self._symlink_windrot(target_directory)

    def _symlink_hgridll(self, target_directory):
        hgridll_lnk = target_directory / 'hgrid.ll'
        if hgridll_lnk.exists():
            hgridll_lnk.unlink()
        if self.args.overwrite is True:
            if self.args.overwrite is True:
                if hgridll_lnk.exists():
                    hgridll_lnk.unlink()
        os.symlink(os.path.relpath(
            self.hgrid_path, target_directory), hgridll_lnk)

    def _symlink_windrot(self, target_directory):
        windrot_lnk = target_directory / 'windrot_geo2proj.gr3'
        if self.args.overwrite is True:
            if windrot_lnk.exists():
                windrot_lnk.unlink()
        os.symlink(
            os.path.relpath(
                self.windrot_path, target_directory), windrot_lnk)
