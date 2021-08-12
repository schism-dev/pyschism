import argparse
from datetime import timedelta
from enum import Enum

# import json
import logging

# import os
import pathlib

# import shutil

# import geopandas as gpd
from psutil import cpu_count

from pyschism import dates
from pyschism.driver import ModelConfig

# from pyschism.enums import ForecastProduct
# from pyschism.forcing.nws import NWS2, GFS, HRRR
# from pyschism.forcing.source_sink import NWM
# from pyschism.forcing.bctides import Tides
from pyschism.cmd import common
from pyschism.forcing.bctides import Bctides
from pyschism.mesh import gridgr3, prop

# from pyschism.mesh.fgrid import Fgrid, ManningsN, DragCoefficient
from pyschism.param.schout import SurfaceOutputVars

logger = logging.getLogger(__name__)

CONFIG_FILE_NAME = "config.json"
STATIC_DIRECTORY = "static"
FORECAST_DIRECTORY = "forecast"


class GridGr3Type(Enum):
    ALBEDO = gridgr3.Albedo
    DIFFMAX = gridgr3.Diffmax
    DIFFMIN = gridgr3.Diffmin
    WATERTYPE = gridgr3.Watertype
    SHAPIRO = gridgr3.Shapiro
    WINDROT = gridgr3.Windrot
    ELEV_IC = gridgr3.ElevIc


class PropType(Enum):
    FLUXFLAG = prop.Fluxflag
    TVDFLAG = prop.Tvdflag


class GridGr3Descriptor:
    def __init__(self, gridgr3_type):
        self.type = gridgr3_type

    def __set__(self, obj, val):
        if val is not None:
            assert isinstance(val, self.type)
        self.gridgr3 = val

    def __get__(self, obj, val):
        if not hasattr(self, "gridgr3"):
            if obj.args.vgrid.is3D():
                return self.type.default(obj.args.hgrid)
        else:
            return self.gridgr3


class ForecastCliMeta(type):
    def __new__(meta, name, bases, attrs):
        attrs["surface_output_vars"] = SurfaceOutputVars()
        for gr3type in GridGr3Type:
            attrs[gr3type.name.lower()] = GridGr3Descriptor(gr3type.value)

        return type(name, bases, attrs)


class ForecastCli(metaclass=ForecastCliMeta):

    start_date = dates.StartDate()
    end_date = dates.EndDate()

    def __init__(self, args: argparse.Namespace):
        self.start_date = dates.nearest_cycle()
        self.args = args

        if self.args.skip_run is True:
            if self.coldstart is not None:
                if self.args.spinup_days is not None:
                    self.coldstart.write(
                        self.coldstart_directory,
                        overwrite=self.args.overwrite,

                    )
                else:
                    self.coldstart.write(
                        self.hotstart_directory,
                        overwrite=self.args.overwrite

                    )
        else:
            if self.coldstart is not None:
                if self.args.spinup_days is not None:
                    self.coldstart.run(
                        self.coldstart_directory,
                        overwrite=self.args.overwrite
                    )
                else:
                    self.coldstart.run(
                        self.hotstart_directory,
                        overwrite=self.args.overwrite
                    )



        #         hotstart.write(self.hotstart_directory,
        #                        overwrite=self.args.overwrite)

    #         hotstart.run(self.hotstart_directory,
    #                      overwrite=self.args.overwrite)

    # def load_user_arguments(self):
    #     raise NotImplementedError

    # def get_drivers(self):
    #     return self.get_coldstart(), self.get_hotstart()

    @property
    def coldstart(self):

        if self.args.vgrid.is2D() is True:
            if self.args.spinup_days is not None:
                return self.config.coldstart(
                    timestep=self.args.timestep,
                    start_date=self.start_date - self.args.run_days,
                    end_date=self.start_date,
                    dramp=self.args.spinup_days,
                    drampbc=self.args.spinup_days,
                    dramp_ss=self.args.spinup_days,
                    drampwafo=self.args.spinup_days,
                    drampwind=self.args.spinup_days,
                    elev_ic=self.args.elev_ic,
                    temp_ic=self.args.temp_ic,
                    salt_ic=self.args.salt_ic,
                    nspool=None,
                    ihfskip=None,
                    nhot_write=None,
                    stations=None,
                    server_config=None,
                    param_template=self.args.use_param_template,
                    # **self.user_requested_surface_outputs,
                )
            else:
                return self.config.coldstart(
                    timestep=self.args.timestep,
                    start_date=self.start_date,
                    end_date=self.args.run_days,
                    elev_ic=self.args.elev_ic,
                    temp_ic=self.args.temp_ic,
                    salt_ic=self.args.salt_ic,
                    nspool=None,
                    ihfskip=None,
                    nhot_write=None,
                    stations=None,
                    server_config=None,
                    param_template=self.args.use_param_template,
                    # **self.user_requested_surface_outputs, 
                )

        else:
            raise NotImplementedError("Model is 3D.")

    @property
    def hotstart(self):

        raise NotImplementedError("hotstart")

        def get_surface_outputs():
            surface_outputs = {}
            outvars = []
            for vardata in self.surface_output_vars.values():
                for varname, _ in vardata:
                    outvars.append(varname)
            for key, val in self.args.__dict__.items():
                if key in outvars and val is True:
                    surface_outputs[key] = val
            return surface_outputs

        return self.config.hotstart(
            self.get_hotstart_driver(),
            timestep=self.args.timestep,
            end_date=timedelta(days=2) - timedelta(hours=2),
            nspool=self.args.nspool,
            **get_surface_outputs(),
        )

    @property
    def config(self):
        if not hasattr(self, "_config"):
            self._config = ModelConfig(
                self.args.hgrid,
                vgrid=self.args.vgrid,
                fgrid=self.args.fgrid,
                albedo=self.args.albedo,
                diffmin=self.args.diffmin,
                diffmax=self.args.diffmax,
                watertype=self.args.watertype,
                fluxflag=self.args.fluxflag,
                tvdflag=self.args.tvdflag,
                bctides=self.bctides,
                nws=self.args.nws,
                source_sink=self.args.source_sink,
                # waves=self.args.waves,
            )
        return self._config

    @property
    def bctides(self):
        if not hasattr(self, "_bctides"):
            self._bctides = Bctides(
                self.args.hgrid,
                vgrid=self.args.vgrid,
                iettype=self.args.iettype,
                ifltype=self.args.ifltype,
                isatype=self.args.isatype,
                itetype=self.args.itetype,
                # itrtype=self.args.itrtype,
                cutoff_depth=self.args.cutoff_depth,
            )
        return self._bctides

    # @property
    # def static_directory(self):
    #     if not hasattr(self, '_static_directory'):
    #         self._static_directory = self.project_directory / STATIC_DIRECTORY
    #         self._static_directory.mkdir(exist_ok=self.args.overwrite)
    #     return self._static_directory

    @property
    def forecasts_directory(self):
        return self.args.project_directory / FORECAST_DIRECTORY

    @property
    def coldstart_directory(self):
        return self.hotstart_directory / "coldstart"

    @property
    def hotstart_directory(self):
        if not hasattr(self, "_hotstart_directory"):
            timestamp = str(self.start_date).replace(" ", "T")
            self._hotstart_directory = self.forecasts_directory / f"{timestamp}"
            self._hotstart_directory.parent.mkdir(exist_ok=True)
            self._hotstart_directory.mkdir(exist_ok=True)
        return self._hotstart_directory

    # @property
    # def hgrid_path(self):
    #     if not hasattr(self, '_hgrid_path'):
    #         hgrid_path = self.static_directory / 'hgrid.gr3'
    #         if not hgrid_path.exists() or self.args.overwrite is True:
    #             logger.info(
    #                 f'Copying hgrid file from {self.args.hgrid} to '
    #                 f'{hgrid_path}.'
    #                 )
    #             if self.args.overwrite is True:
    #                 if hgrid_path.is_file():
    #                     os.remove(hgrid_path)
    #             shutil.copy2(self.args.hgrid, hgrid_path,
    #                          follow_symlinks=True)
    #         self._hgrid_path = hgrid_path
    #     return self._hgrid_path

    # @property
    # def vgrid_path(self):
    #     if not hasattr(self, '_vgrid_path'):
    #         vgrid_path = self.static_directory / 'vgrid.in'
    #         if not vgrid_path.exists() or self.args.overwrite is True:
    #             if self.args.vgrid_bin is not None:
    #                 logger.info(
    #                     f'Calling vgrid binary from {self.args.vgrid_bin}.')
    #                 Vgrid.from_binary(
    #                     self.hgrid,
    #                     binary=self.args.vgrid_bin
    #                 ).write(
    #                     vgrid_path,
    #                     overwrite=self.args.overwrite
    #                 )
    #             else:
    #                 if self.args.vgrid is None:
    #                     logger.info(
    #                         f'Writing default vgrid file to {vgrid_path}.')
    #                     Vgrid.default().write(
    #                         vgrid_path,
    #                         overwrite=self.args.overwrite)
    #                 else:
    #                     logger.info(
    #                         f'Copying vgrid file from {self.args.vgrid} to '
    #                         f'path {vgrid_path}.')
    #                     if self.args.overwrite is True:
    #                         if vgrid_path.is_file():
    #                             os.remove(vgrid_path)
    #                     shutil.copy2(self.args.vgrid, vgrid_path,
    #                                  follow_symlinks=True)
    #         self._vgrid_path = vgrid_path
    #     return self._vgrid_path

    # @property
    # def fgrid_path(self):
    #     if not hasattr(self, '_fgrid_path'):

    #         # user did not specify an fgrid
    #         if self.args.fgrid is None:

    #             if self.vgrid.is2D():
    #                 fgrid_path = self.static_directory / 'manning.gr3'
    #                 logger.info(
    #                     'Initializing default mannings file for 2D model to '
    #                     f'{fgrid_path}.')
    #                 ManningsN.linear_with_depth(
    #                     self.hgrid).write(fgrid_path,
    #                                       overwrite=self.args.overwrite)
    #             else:
    #                 fgrid_path = self.static_directory / 'drag.gr3'
    #                 logger.info(
    #                     'Initializing default drag file for 3D model to '
    #                     f'{fgrid_path}')
    #                 DragCoefficient.linear_with_depth(
    #                     self.hgrid).write(fgrid_path,
    #                                       overwrite=self.args.overwrite)

    #         else:
    #             # user can override fgrid_type for files with arbitrary names
    #             # auto means that we derive the type of fgrid from the filename
    #             if self.args.fgrid_type == 'auto':
    #                 # just opens the Fgrid to make sure it's valid
    #                 fgrid = Fgrid.open(
    #                     self.args.fgrid,
    #                     crs=self.args.fgrid_crs
    #                 )
    #                 fgrid_path = self.static_directory / fgrid.name
    #                 if len(fgrid.nodes) != len(self.hgrid.nodes):
    #                     raise Exception(
    #                         'Nodes array mismatch between fgrid '
    #                         f'{self.args.fgrid} and hgrid {self.args.hgrid}.')

    #                 if len(fgrid.elements) != len(self.hgrid.elements):
    #                     raise Exception(
    #                         'Elements array mismatch between fgrid '
    #                         f'{self.args.fgrid} and hgrid {self.args.hgrid}.')

    #             # user specified an fgrid_type
    #             else:
    #                 fgrid_path = self.static_directory / \
    #                     (self.args.fgrid_type + '.gr3')

    #             if self.args.overwrite is True:
    #                 if fgrid_path.is_file():
    #                     os.remove(fgrid_path)

    #             logger.info(
    #                 f'Copying friction file {self.args.fgrid} to '
    #                 f'{fgrid_path}.')

    #             shutil.copy2(self.args.fgrid, fgrid_path,
    #                          follow_symlinks=True)

    #         self._fgrid_path = fgrid_path
    #     return self._fgrid_path

    # @property
    # def hgrid(self):
    #     if not hasattr(self, '_hgrid'):
    #         self._hgrid = Hgrid.open(self.hgrid_path, crs=self.args.hgrid_crs)
    #     return self._hgrid

    # @property
    # def vgrid(self):
    #     if not hasattr(self, '_vgrid'):
    #         self._vgrid = Vgrid.open(self.vgrid_path)
    #     return self._vgrid

    # @property
    # def fgrid(self):
    #     if not hasattr(self, '_fgrid'):
    #         self._fgrid = Fgrid.open(self.fgrid_path, crs=self.args.fgrid_crs)
    #     return self._fgrid

    # @property
    # def fluxflag(self):
    #     if not hasattr(self, '_fluxflag'):
    #         if self.args.fluxflag is None:
    #             if self.vgrid.is3D() is True:
    #                 fluxflag = prop.Fluxflag.default(self.hgrid)
    #             else:
    #                 fluxflag = None
    #         else:
    #             fluxflag = self.args.fluxflag
    #         self._fluxflag = fluxflag
    #     return self._fluxflag

    # @property
    # def tvdflag(self):
    #     if not hasattr(self, '_tvdflag'):
    #         if self.args.tvdflag is None:
    #             if self.vgrid.is3D() is True:
    #                 tvdflag = prop.Tvdflag.default(self.hgrid)
    #             else:
    #                 tvdflag = None
    #         else:
    #             tvdflag = self.args.tvdflag
    #         self._tvdflag = tvdflag
    #     return self._tvdflag

    @staticmethod
    def add_subparser_action(subparsers):
        add_forecast_options_to_parser(subparsers.add_parser("forecast"))


def add_forecast_options_to_parser(parser):
    actions = parser.add_subparsers(dest="action")
    add_forecast_init_to_parser(
        actions.add_parser(
            "init",
            help="Initializes a directory for a sequential SCHISM forecast model "
            "deployment.",
        )
    )
    add_forecast_update_to_parser(
        actions.add_parser(
            "update",
            help="Updates a directory  that has been previously initialized for "
            "SCHISM forecast model deployment.",
        )
    )


class ProjectDirectoryAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        tmp_parser = argparse.ArgumentParser(add_help=False)
        tmp_parser.add_argument("--overwrite", action="store_true")
        tmp_args = tmp_parser.parse_known_args()[0]
        project_directory = pathlib.Path(values)
        project_directory.mkdir(exist_ok=tmp_args.overwrite)
        setattr(namespace, self.dest, project_directory)


def add_forecast_init_to_parser(parser):
    parser.add_argument(
        "project_directory",
        action=ProjectDirectoryAction,
        help="System path where the project will be staged.",
    )
    common.add_hgrid_to_parser(parser)
    common.add_vgrid_to_parser(parser)
    common.add_fgrid_to_parser(parser)
    common.add_stratification_to_parser(parser)
    parser.add_argument("--timestep", type=float, required=True)
    parser.add_argument(
        "--forecast-days",
        type=lambda x: timedelta(days=float(x)),
        required=True,
        dest="run_days",
    )
    parser.add_argument(
        "--spinup-days",
        help="Number of days used for model initialization. "
        "Defaults to 15 days spinup.",
        type=lambda x: timedelta(days=float(x)),
    )
    parser.add_argument(
        "--skip-run", action="store_true", help="Skips running the model."
    )
    parser.add_argument("--nproc", type=int, default=cpu_count(logical=False))
    parser.add_argument(
        "--use-param-template",
        # TODO: We need to optionally take a user-provided template
        action="store_true",
    )
    common.add_gridgr3_to_parser(parser)
    common.add_ic_to_parser(parser)
    common.add_prop_to_parser(parser)
    common.add_ibctype_to_parser(parser)
    common.add_bctides_options_to_parser(parser)
    common.add_nws_to_parser(parser)
    common.add_source_sink_to_parser(parser)
    # common.add_waves_to_parser(parser)
    common.add_surface_outputs_to_parser(parser)
    common.add_stations_outputs_to_parser(parser)
    common.add_log_level_to_parser(parser)
    parser.add_argument(
        "--overwrite", action="store_true", help="Allow overwrite of output directory."
    )


def add_forecast_update_to_parser(parser):
    parser.add_argument(
        "project_directory", help="System path of the staged directory to update."
    )


# ------ drafts ----------


































#     init.add_argument('hgrid', help='Horizontal grid file.')
#     init.add_argument('--fgrid', help='Friction grid file.')

#     mesh_options = init.add_argument_group('mesh_options')
#     mesh_options.add_argument('--hgrid-crs')
#     mesh_options.add_argument('--fgrid-crs')
#     mesh_options.add_argument(
#         '--fgrid-type',
#         choices=['auto', 'manning', 'drag', 'rough'],
#         default='auto')
#     init.add_argument("--timestep", type=float, required=True)
#     init.add_argument("--forecast-days", type=float, required=True)
#     init.add_argument(
#         "--spinup-days",
#         help="Number of days used for model initialization. "
#         "Defaults to 15 days spinup.",
#         type=lambda x: timedelta(days=float(x)),
#     )
#     add_prop_options_to_parser(init)
#     init.add_argument(
#         "--skip-run", action="store_true",
#         help="Skips running the model.")
#     init.add_argument('--nproc', type=int, default=cpu_count(logical=False))
#     init.add_argument('--use-param-template', action='store_true')
#     init.add_argument(
#         "--log-level",
#         choices=['info', 'warning', 'debug'],
#         default='warning'
#     )
#     init.add_argument(
#         "--overwrite", action="store_true",
#         help="Allow overwrite of output directory.")
#     _add_tidal_constituents(init)
#     _add_atmospheric_forcing(init)
#     _add_hydrologic_forcing(init)
#     init.add_argument('--shapiro', action='store_true')
#     init.add_argument('--shapiro-bin', action='store_true')

#     # TODO: Additional forcings.
#     # _add_wave_forcing(forecast)
#     model_outputs = init.add_argument_group('model_outputs')
#     # TODO: Stations outputs.
#     _add_stations_outputs(model_outputs)
#     _add_surface_outputs(model_outputs)
#     server_config = init.add_subparsers(dest="server_config")
#     slurm = server_config.add_parser(
#         'slurm', help="Add options for slurm run configuration.")
#     slurm.add_argument('--account')
#     slurm.add_argument('--partition')
#     slurm.add_argument('--walltime', type=float, help="In hours, float.")
#     slurm.add_argument('--slurm-filename')
#     slurm.add_argument('--slurm-rundir')
#     slurm.add_argument('--run-name')
#     slurm.add_argument('--mail-type')
#     slurm.add_argument('--mail-user')
#     slurm.add_argument('--log-filename')
#     slurm.add_argument('--path-prefix')
#     slurm.add_argument('--slurm-nodes')
#     slurm.add_argument('--slurm-launcher', default='srun')
#     slurm.add_argument('--extra-commands', action='append')
#     slurm.add_argument(
#         '--module',
#         default=list(),
#         action='append',
#         dest='modules')


# def add_forecast_update(actions):
#     update = actions.add_parser("update")
#     update.add_argument("project_directory")

# def _add_tidal_constituents(parser):
#     tides = parser.add_argument_group('tides')
#     options = tides.add_mutually_exclusive_group()
#     options.required = True
#     options.add_argument("--all-constituents", action="store_true")
#     options.add_argument("--major-constituents", action="store_true")
#     options.add_argument(
#         "-c", "--constituents",
#         action='append',
#         # choices=["K1", "O1", "P1", "Q1", "MM", "Mf", "M4", "MN4", "MS4",
#         #         "2N2", "S1"],
#         choices=["Z0", "K2", "S2", "M2", "N2", "K1", "P1", "O1", "Q1"],
#         dest='constituents',
#         default=False,
#         help="Tidal constituent to be forced in the model (case-insensitive).")
#     tides.add_argument("--tidal-database", choices=['tpxo', 'hamtide'],
#                        default='tpxo')
#     tides.add_argument("--include-tidal-velocity",
#                        "--bnd-vel", action="store_true", dest='bnd_vel')
#     tides.add_argument("--Z0", type=float)


# def _add_atmospheric_forcing(parser):
#     atmospheric_forcing_1 = parser.add_argument_group(
#         'atmospheric forcing level 1')
#     data_source_1 = atmospheric_forcing_1.add_mutually_exclusive_group()
#     for product in ForecastProduct:
#         data_source_1.add_argument(
#             f"--{product.value.lower().replace('_', '-')}",
#             action="store_true",
#             help=f'Use {product.value} as data source.')
#     data_source_1.add_argument(
#         "--gfs", action="store_true",
#         help="Alias for --gfs-0p25-1hr")
#     data_source_1.add_argument(
#         "--gdas", action='store_true',
#         help="Alias for --gdas-0p25",
#     )
#     atmospheric_forcing_1.add_argument(
#         "--no-air-1", action="store_true",
#         help="Disables air mass forcing option.")
#     atmospheric_forcing_1.add_argument(
#         "--no-prc-1", action='store_true',
#         help="Disables precipitation option.")
#     atmospheric_forcing_1.add_argument(
#         "--no-rad-1", action='store_true',
#         help="Disables solar radiation flux option.")

#     atmospheric_forcing_2 = parser.add_argument_group(
#         'atmospheric forcing level 2')
#     data_source_2 = atmospheric_forcing_2.add_mutually_exclusive_group()
#     for product in ["hrrr"]:
#         data_source_2.add_argument(
#             f"--{product}",
#             help=f'Use {product.upper()} as data source.',
#             action='store_true'
#         )
#     atmospheric_forcing_2.add_argument(
#         "--no-air-2", action="store_true",
#         help="Disables air mass forcing option.")
#     atmospheric_forcing_2.add_argument(
#         "--no-prc-2", action='store_true',
#         help="Disables precipitation option.")
#     atmospheric_forcing_2.add_argument(
#         "--no-rad-2", action='store_true',
#         help="Disables solar radiation flux option.")


# def _add_hydrologic_forcing(parser):
#     src_snk = parser.add_argument_group('Sources and sinks')
#     src_snk.add_argument("--hydrology", action="append",
#                          help="Add source and sink flow.",
#                          choices=['NWM'],
#                          default=[])


# def _add_wave_forcing(parser):
#     wave_forcing = parser.add_argument_group('Waves')
#     data_source = wave_forcing.add_mutually_exclusive_group()
#     for product in ["WWIII"]:
#         data_source.add_argument(f'--{product}')
