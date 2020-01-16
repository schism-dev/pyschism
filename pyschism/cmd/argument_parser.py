import logging
import sys


def add_mesh_options(parser):
    # mesh
    parser.add_argument('hgrid')

    parser.add_argument('--vgrid')

    # Mesh spatial reference
    msg = "Mesh projection information. Defaults to EPSG:4326 which "
    msg += "corresponds to WGS84. For Cartesian meshes use 3395 which "
    msg += "corresponds to Mercator projection."
    parser.add_argument('--crs', default=4326)


def add_general_options(parser, runtype=None):
    add_mesh_options(parser)

    if runtype is not None:
        if runtype == "tidal":
            add_tidal_run_options(parser)

        elif runtype == 'best_track':
            add_best_track_options(parser)

    # output directory
    msg = "Directory to which SCHISM input files will be written to. "
    parser.add_argument('--output-directory', '--outdir', "-o", help=msg)

    msg = 'Allows overwrite of output directory.'
    parser.add_argument('--overwrite', help=msg, action="store_true")

    msg = "Generates and saves input files to the output directory but does "
    msg += "not deploy the SCHISM run."
    parser.add_argument(
        '--generate-only', "--no-run", "--skip-run",
        action="store_true",
        help=msg
        )

    add_log_level_options(parser)
    add_server_options(parser)

    # add tidal constituents
    add_tidal_constituents_options(parser)

    # add surface output requests
    add_surface_output_request('elevation', parser)
    add_surface_output_request('velocity', parser)
    add_surface_output_request('meteorological', parser)
    add_surface_output_request('concentration', parser)

    # parse stations from a file
    msg = "File containing list of stations for outputs. It will parse "
    msg += "the stations below the NOUTE, NOUTV and NOUTM keywords for "
    msg += "their respective stations list."
    parser.add_argument(
        "--stations-file",
        help=msg
        )

    # stations output requests
    add_stations_output_request('elevation', parser)
    add_stations_output_request('velocity', parser)
    add_stations_output_request('meteorological', parser)
    add_stations_output_request('concentration', parser)

    parser.add_argument(
        '--ascii',
        dest='netcdf',
        action='store_false',
        default=True,
        help="Request outputs in ASCII format. NetCDF is the default."
    )


def add_log_level_options(parser):
    log_level = parser.add_mutually_exclusive_group()
    log_level.add_argument(
        '--log-level-info',
        nargs='?',
        const=logging.INFO,
        dest="log_level")
    log_level.add_argument(
        '--log-level-debug',
        nargs='?',
        const=logging.DEBUG,
        dest="log_level")
    log_level.add_argument(
        '--log-level-warning',
        nargs='?',
        const=logging.WARNING,
        dest="log_level")


def add_server_options(parser):

    # flag some options as required when a resource manager is enabled
    _required = "--use-torque" in sys.argv
    _required = _required | ("--use-pbs" in sys.argv)
    _required = _required | ("--use-slurm" in sys.argv)

    # add server options
    parser.add_argument('--hostname')
    parser.add_argument('--port', type=int)
    parser.add_argument("--wdir", required=_required)
    parser.add_argument("--keep-wdir", action="store_true")
    parser.add_argument(
        "--binaries-path", "--binaries-prefix",
        dest="binaries_prefix")
    parser.add_argument("--source-script")
    parser.add_argument("--additional-mpi-options")

    # make nproc required when using ssh
    args = parser.parse_known_args()[0]
    if args.hostname is not None:
        parser.add_argument("--nproc", "--ncpu", type=int, required=True)
    else:
        parser.add_argument("--nproc", "--ncpu", type=int, default=-1)

    # add resource manager option
    manager = parser.add_mutually_exclusive_group()
    manager.add_argument('--use-torque', action="store_true")
    manager.add_argument('--use-pbs', action="store_true")
    manager.add_argument('--use-slurm', action="store_true")

    # resource manager specific options
    parser.add_argument('--account', required=_required)
    parser.add_argument('--walltime', required=_required)
    parser.add_argument(
        '--module',
        default=list(),
        action='append',
        dest='modules'
        )


def add_tidal_constituents_options(parser):
    # tidal constituents
    msg = "Tidal constituent to be forced in the model. Pass "
    msg += "--use-constituent='all' to use all available constituents "
    msg += "(K1, O1, P1, Q1, MM, Mf, M4, MN4, MS4, 2N2, S1) "
    msg += "or use --use-constituent='major'. "
    msg += "For a custom list of forcing constituents, pass -c= for each "
    msg += "individual constituent to use (case-insensitive). "
    msg += "Use None for no tidal forcing. Defaults to 'all'."
    parser.add_argument(
        "--use-constituent", "-c",
        action='append',
        choices=["K1", "O1", "P1", "Q1", "MM", "Mf", "M4", "MN4", "MS4",
                 "2N2", "S1", "all", "major"],
        dest='constituents',
        default=[],
        help=msg
        )


def add_surface_output_request(physical_var, parser):

    # surface output requests
    msg = f"{physical_var.capitalize()} surface output sampling frequency "
    msg += f"(in minutes). When this number is greater than 0, {physical_var} "
    msg += " surface outputs are written to disk during hotstart phase."
    parser.add_argument(
        f"--{physical_var}-surface-sampling-frequency",
        f"--{physical_var[:4]}",
        type=float,
        help=msg
        )

    # surface harmonic analysis
    msg = f"Enables {physical_var} surface harmonic analysis."
    parser.add_argument(
        f'--{physical_var}-surface-harmonic-analysis',
        f"--{physical_var[:4]}-harm",
        action='store_true',
        default=False,
        help=msg
        )


def add_stations_output_request(physical_var, parser):

    # stations output requests
    msg = f"{physical_var.capitalize()} stations sampling frequency in "
    msg += f"minutes. When this number is greater than 0, {physical_var} "
    msg += "stations output is turned on during hotstart phase."
    parser.add_argument(
        f"--{physical_var}-stations-sampling-frequency",
        f"--{physical_var[:4]}-s",
        type=float,
        help=msg
    )

    # stations harmonic analysis
    msg = f"Enables {physical_var} stations harmonic analysis."
    parser.add_argument(
        f'--{physical_var}-stations-harmonic-analysis',
        f"--{physical_var[:4]}-s-harm",
        action='store_true',
        default=False,
        help=msg
        )


def add_best_track_options(parser):
    # storm_id
    msg = "National Hurricane Center (NHC) storm id. "
    msg += " Examples: AL132012 for Sandy2012 or AL152017 for Maria2017."
    parser.add_argument('storm_id', help=msg)
    parser.add_argument('--start_date')
    parser.add_argument('--end_date')


def add_tidal_run_options(parser):
    # start_date
    msg = "Start date is relative to hotstart, that is, this is the "
    msg += "true start date of the model (in UTC time). Use format "
    msg += "%%Y-%%m-%%dT%%H:%%M to specify date. For example, for August "
    msg += "1, 2013, 00:00 hours, write \"2018-08-01T00:00\" (can be used "
    msg += "with or without the quotes). Time zone is assumed to be UTC "
    msg += "and there is no way to change this on this release. In the "
    msg += "future the code might be expanded to include timezone "
    msg += "specification. See https://docs.python.org/3/library/datetime"
    msg += ".html#strftime-strptime-behavior for more details about "
    msg += "datetime formats."
    parser.add_argument('start_date', help=msg)

    # end_date
    parser.add_argument('end_date')

    # spinup_days
    parser.add_argument('--spinup-days', type=float, required=True)
