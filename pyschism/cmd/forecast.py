import argparse
from datetime import timedelta
from enum import Enum
import logging
import pathlib
from pyschism.server.slurm import SlurmConfig

from psutil import cpu_count

from pyschism import dates
from pyschism.cmd import common
from pyschism.driver import ModelConfig
from pyschism.mesh import gridgr3, prop
from pyschism.param.schout import SurfaceOutputVars
from pyschism.server import ServerConfig, SlurmConfig

logger = logging.getLogger(__name__)

CONFIG_FILE_NAME = "config.json"
STATIC_DIRECTORY = "static"


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
                        self.hotstart_directory, overwrite=self.args.overwrite
                    )

        else:
            if self.coldstart is not None:
                if self.args.spinup_days is not None:
                    self.coldstart.run(
                        self.coldstart_directory, overwrite=self.args.overwrite
                    )
                else:
                    self.coldstart.run(
                        self.hotstart_directory, overwrite=self.args.overwrite
                    )

        if self.args.skip_run is True:
            self.hotstart.write(self.hotstart_directory, overwrite=self.args.overwrite)
        else:
            self.hotstart.run(self.hotstart_directory, overwrite=self.args.overwrite)

    @property
    def coldstart(self):
        if not hasattr(self, "_coldstart"):
            if self.args.vgrid.is2D() is True:
                if self.args.spinup_days is not None:
                    self._coldstart = self.config.coldstart(
                        timestep=self.args.timestep,
                        start_date=self.start_date - self.args.spinup_days,
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
                        server_config=self.server_config,
                        param_template=self.args.use_param_template,
                        # **self.user_requested_surface_outputs,
                    )
                else:
                    self._coldstart = self.config.coldstart(
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
                        server_config=self.server_config,
                        param_template=self.args.use_param_template,
                        # **self.user_requested_surface_outputs,
                    )

            else:
                raise NotImplementedError("Model is 3D.")
        return self._coldstart

    @property
    def hotstart(self):
        return self.config.hotstart(
            self.coldstart,
            timestep=self.args.timestep,
            end_date=self.args.run_days,
            nspool=None,
            ihfskip=None,
            nhot_write=None,
            stations=None,
            server_config=self.server_config,
            param_template=self.args.use_param_template,
            # **self.user_requested_surface_outputs,
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
                iettype=self.args.iettype,
                ifltype=self.args.ifltype,
                itetype=self.args.itetype,
                isatype=self.args.isatype,
                nws=self.args.nws,
                source_sink=self.args.source_sink,
                # waves=self.args.waves,
            )
        return self._config

    @property
    def forecasts_directory(self):
        return self.args.project_directory / "runs"

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

    @property
    def server_config(self):
        if not hasattr(self, "_server_config"):
            if self.args.workload_manager is None:
                self._server_config = ServerConfig(
                    nproc=self.args.ntasks,
                    # symlink_outputs: str = None,
                    schism_binary=self.args.schism_binary,
                    mpi_launcher=None,
                    modules=self.args.modules,
                    modulepath=self.args.modulepath,
                    modules_init=self.args.modules_init,
                )
            elif self.args.workload_manager == "slurm":
                self._server_config = SlurmConfig(
                    account=self.args.account,
                    ntasks=self.args.ntasks,
                    partition=self.args.partition,
                    walltime=self.args.walltime,
                    filename=self.args.filename,
                    run_directory=None,
                    run_name=None,
                    mail_type=None,
                    mail_user=None,
                    log_filename=None,
                    modules=self.args.modules,
                    modulepath=self.args.modulepath,
                    modules_init=self.args.modules_init,
                    schism_binary=self.args.schism_binary,
                    extra_commands=None,
                    launcher=self.args.launcher,
                    nodes=None,
                    # symlink_outputs=None,
                )
        return self._server_config

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
    common.add_tidal_constituents_to_parser(parser)
    common.add_tidal_database_to_parser(parser)
    common.add_baroclinic_database_to_parser(parser)
    common.add_bctides_options_to_parser(parser)
    common.add_nws_to_parser(parser)
    common.add_source_sink_to_parser(parser)
    # common.add_waves_to_parser(parser)
    common.add_surface_outputs_to_parser(parser)
    common.add_stations_outputs_to_parser(parser)
    common.add_log_level_to_parser(parser)
    common.add_schism_binary_to_parser(parser)
    common.add_modules_to_parser(parser)
    common.add_workload_manager_options_to_parser(parser)
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
