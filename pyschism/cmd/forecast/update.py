from argparse import Namespace
from collections import namedtuple
from datetime import timedelta
import json
import logging
import subprocess

from pyschism.cmd.forecast.init import ForecastInit
from pyschism.driver import ModelDriver
from pyschism.domain import ModelDomain
from pyschism.forcing.atmosphere.nws.nws2 import NWS2
from pyschism.forcing.atmosphere import GlobalForecastSystem as GFS
from pyschism.forcing.hydrology import NationalWaterModel as NWM
from pyschism.param.schout import SurfaceOutputVars

logger = logging.getLogger()


class HotstartDirectory:

    def __get__(self, obj, val):
        hotstart_directory = obj.__dict__.get('hotstart_directory')
        if hotstart_directory is None:
            timestamp = str(obj.target_datetime).replace(' ', 'T')
            hotstart_directory = obj.forecasts_directory / f'{timestamp}'
            hotstart_directory.parent.mkdir(exist_ok=True)
            hotstart_directory.mkdir(exist_ok=True)
            obj.__dict__['hotstart_directory'] = hotstart_directory
        return hotstart_directory


class HotstartDomain:

    def __get__(self, obj, val):
        hotstart_domain = obj.__dict__.get('hotstart_domain')
        if hotstart_domain is None:
            # obj.logger.debug('Generating hotstart ModelDomain.')
            hotstart_domain = ModelDomain.open(
                obj.hgrid_path,
                obj.fgrid_path,
                # TODO: Vgrid.open() raises NotImplementedError
                # obj.vgrid_path,
                hgrid_crs=obj.args.hgrid_crs,
                fgrid_crs=obj.args.fgrid_crs
            )
            if obj.tides is not None:
                # obj.logger.debug('Adding tides to hotstart_domain.')
                hotstart_domain.add_boundary_condition(obj.tides)
            if obj.nws2 is not None:
                # obj.logger.debug('Adding NWS2 object to hotstart_domain.')
                hotstart_domain.set_atmospheric_forcing(obj.nws2)
            if obj.hydrology is not None:
                for hydrology in obj.hydrology:
                    # obj.logger.debug(
                    #     'Adding hydrology object to hotstart_domain.')
                    hotstart_domain.add_hydrology(hydrology)
            obj.__dict__['hotstart_domain'] = hotstart_domain
        return hotstart_domain


class HotstartDriver:

    # reusing a descriptor from the schout module.
    surface_output_vars = SurfaceOutputVars()

    def __get__(self, obj, val):
        hotstart_driver = obj.__dict__.get('hotstart_driver')
        if hotstart_driver is None:
            logger.debug('Creating hotstart driver.')
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
                **self.surface_outputs(obj)
                )
            obj.__dict__['hotstart_driver'] = hotstart_driver
        return hotstart_driver

    def surface_outputs(self, obj):
        surface_outputs = {}
        outvars = []
        for vardata in self.surface_output_vars.values():
            for varname, _ in vardata:
                outvars.append(varname)
        for key, val in obj.args.__dict__.items():
            if key in outvars and val is True:
                surface_outputs[key] = val
        return surface_outputs


class ForecastsDirectory:

    def __get__(self, obj, val):
        forecasts_directory = obj.__dict__.get('forecasts_directory')
        if forecasts_directory is None:
            forecasts_directory = obj.project_directory / 'forecasts'
            forecasts_directory.mkdir(exist_ok=True)
            obj.__dict__['forecasts_directory'] = forecasts_directory
        return forecasts_directory


class Forcings:

    def __get__(self, obj, val):
        forcings = obj.__dict__.get('forcings')
        if forcings is None:
            forcings = namedtuple(
                'forcings', ['tides', 'atmosphere', 'hydrology'])
            forcings = forcings(
                tides=obj.tides,
                atmosphere=obj.nws2,
                hydrology=obj.hydrology
                )
            obj.__dict__['forcings'] = forcings
        return forcings


class NWS2Descriptor:

    def __get__(self, obj, val):
        nws2 = obj.__dict__.get('nws2')
        if nws2 is None:
            sflux_1 = self.sflux_1(obj)
            if sflux_1 is not None:
                nws2 = NWS2(
                    sflux_1,
                    self.sflux_2(obj)
                    )
                obj.__dict__['nws2'] = nws2
        return nws2

    def sflux_1(self, obj):
        for arg, value in obj.args.__dict__.items():
            if 'gdas' in arg:
                if value is True:
                    raise NotImplementedError(
                        'GDAS product not implemented.')
            elif 'gfs' in arg:
                if value is True:
                    if arg == 'gfs':
                        return GFS()
                    else:
                        return GFS(product=arg)

    def sflux_2(self, obj):
        return


class Sflux2:

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


class HydrologyDescriptor:

    def __get__(self, obj, val):
        hydrology = obj.__dict__.get('hydrology')
        if hydrology is None:
            hydrology = []
            for hydro in obj.args.hydrology:
                if hydro == "NWM":
                    obj.logger.debug('Append NWM object.')
                    hydrology.append(NWM())
            obj.__dict__['hydrology'] = hydrology
        return hydrology


class ForecastUpdate(ForecastInit):

    forecasts_directory = ForecastsDirectory()
    hotstart_directory = HotstartDirectory()
    hotstart_domain = HotstartDomain()
    hotstart_driver = HotstartDriver()
    forcings = Forcings()
    nws2 = NWS2Descriptor()
    hydrology = HydrologyDescriptor()

    def __init__(self, args: Namespace):
        self._args = args
        logger.info("Starting forecast update sequence.")
        self._load_config()

        """
        TODO: Generate slurm-aware pre-processing of inputs.
        This means we need to create "from_file" constructors for each of the
        inputs.
        Essentially, if we are running this script on the "bare metal" all
        processes have to be run sequentially because we don't have enough
        CPU's. In contrast, with slurm and other process managers, we can
        launch multiple concurrent parallel pools.
        The current code base assumes bare metal processing.
        """
        self._symlink_files(
            self.hotstart_directory,
            self.hotstart_domain.ics,
            windrot=(True if isinstance(self.forcings.atmosphere, NWS2)
                     else False),
            )
        logger.info("Calling hotstart write sequence.")
        self.hotstart_driver.write(
            self.hotstart_directory,
            hgrid=False,
            vgrid=False,
            fgrid=False,
            wind_rot=False,
            overwrite=self.args.overwrite
        )
        if self.args.skip_run is False:
            self.logger.info("Executing SCHISM.")
            subprocess.check_call(
                ["make", "run"],
                cwd=self.hotstart_directory
            )

    def _load_config(self):
        config = self.project_directory / 'config.json'
        with open(config) as json_file:
            config = json.load(json_file)
        config.update(self.args.__dict__)
        self._args = Namespace(**config)

    @property
    def nspool(self):
        if self.args.nspool is not None:
            if '.' in self.args.nspool:
                return timedelta(hours=float(self.args.nspool))
            else:
                return int(self.args.nspool)

    @property
    def previous_run_directory(self):
        if (self.coldstart_directory / 'outputs').exists():
            return self.coldstart_directory
        else:
            raise NotImplementedError(
                'Must find previous_run_directory from the hotstarts or '
                'regenerate init')
