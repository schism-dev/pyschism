import argparse
import logging
import sys


def mesh(parser):

    # mesh
    parser.add_argument('hgrid')

    parser.add_argument('--vgrid')

    # Mesh spatial reference
    msg = "Mesh projection information. Defaults to EPSG:4326 which "
    msg += "corresponds to WGS84. For Cartesian meshes use 3395 which "
    msg += "corresponds to Mercator projection."
    parser.add_argument('--crs')


def output_directory(parser):
    # output directory
    msg = "Directory to which SCHISM input files will be written to. "
    parser.add_argument('--output-directory', '--outdir', "-o", help=msg)


def allow_overwrite(parser):
    msg = 'Allows overwrite of output directory.'
    parser.add_argument('--overwrite', help=msg, action="store_true")


def generate_only(parser):
    msg = "Generates and saves input files to the output directory but does "
    msg += "not deploy the SCHISM run."
    parser.add_argument(
        '--generate-only', "--no-run", "--skip-run",
        action="store_true",
        help=msg
        )


def log_level(parser):
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


def server(parser):

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


def tidal_constituents(parser):
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


def best_track(parser):
    # storm_id
    msg = "National Hurricane Center (NHC) storm id. "
    msg += " Examples: AL132012 for Sandy2012 or AL152017 for Maria2017."
    parser.add_argument('storm_id', help=msg)
    parser.add_argument('--start_date')
    parser.add_argument('--end_date')


def tidal_run(parser):
    # start_date
    msg = "Start date is relative to hotstart, that is, this is the "
    msg += "true start date of the model (in UTC time). Use format "
    msg += "%%Y-%%m-%%dT%%H:%%M to specify date. For example, for August "
    msg += "1, 2013, 00:00 hours, write \"2018-08-01T00:00\" (can be used "
    msg += "with or without the quotes)."
    parser.add_argument('start_date', help=msg)
    # end_date
    parser.add_argument('end_date')
    # spinup_days
    parser.add_argument('--spinup-days', type=float, required=True)


def hydro():
    """
    hydro output options
    """
    return {
        1: ("elev", "0: off; 1: on - elev. [m]"),
        2: ("air_pressure", "air pressure [Pa]"),
        3: ("air_temperature", "air temperature [C]"),
        4: ("specific_humidity", "Specific humidity [-]"),
        5: ("solar_radiation", "solar (shortwave) radiation [W/m/m]"),
        6: ("sensible_flux", "sensible flux (positive upward) [W/m/m]"),
        7: ("latent_heat", "latent heat flux (positive upward) [W/m/m]"),
        8: ("upward_longwave", "upward longwave radiation (positive upward) [W/m/m]"),
        9: ("downward_longwave", "downward longwave radiation (positive downward) [W/m/m]"),
        10: ("total_heat_flux", "total flux=-flsu-fllu-(radu-radd) [W/m/m]"),
        11: ("evaporation", "evaporation rate [kg/m/m/s]"),
        12: ("precipitation", "precipitation rate [kg/m/m/s]"),
        13: ("bottom_stress", "Bottom stress vector [kg/m/s^2(Pa)]"),
        14: ("wind_speed", "wind velocity vector [m/s]"),
        15: ("wind_stress", "wind stress vector [m^2/s/s]"),
        16: ("dahv", "depth-averaged vel vector [m/s]"),
        17: ("vertical_velocity", "vertical velocity [m/s]"),
        18: ("temp", "water temperature [C]"),
        19: ("salt", "water salinity [PSU]"),
        20: ("water_density", "water density [kg/m^3]"),
        21: ("diffusivity", "eddy diffusivity [m^2/s]"),
        22: ("viscosity", "eddy viscosity [m^2/s]"),
        23: ("TKE", "turbulent kinetic energy"),
        24: ("mixing-lenght", "turbulent mixing length [m]"),
        25: ("hvel", "horizontal vel vector [m/s]"),
        26: ("hvel_side", "horizontal vel vector defined @side [m/s]"),
        27: ("wvel_elem", "vertical vel. @elem [m/s]"),
        28: ("temp_elem", "T @prism centers [C]"),
        29: ("salt_elem", "S @prism centers [PSU]"),
        30: ("pressure_gradient", "Barotropic pressure gradient force vector (m.s-2) @side centers"),

    }


def wwm():
    """
    WWM output options
    """
    return {
        1: ("WWM_1", "sig. height (m)"),
        2: ("WWM_2", "Mean average period (sec) - TM01"),
        3: ("WWM_3", "Zero down crossing period for comparison with buoy (s) - TM02"),
        4: ("WWM_4", "Average period of wave runup/overtopping - TM10"),
        5: ("WWM_5", "Mean wave number (1/m)"),
        6: ("WWM_6", "Mean wave length (m)"),
        7: ("WWM_9", "Mean average energy transport direction (degr) - MWD in NDBC?"),
        8: ("WWM_10", "Mean directional spreading (degr)"),
        9: ("WWM_11", "Discrete peak period (sec) - Tp"),
        10: ("WWM_12", "Continuous peak period based on higher order moments (sec)"),
        11: ("WWM_13", "Peak phase vel. (m/s)"),
        12: ("WWM_14", "Peak n-factor."),
        13: ("WWM_15", "Peak group vel. (m/s)"),
        14: ("WWM_16", "Peak wave number"),
        15: ("WWM_17", "Peak wave length"),
        16: ("WWM_18", "Peak (dominant) direction (degr)"),
        17: ("WWM_19", "Peak directional spreading"),
        18: ("WWM_20", "Discrete peak direction (radian?) "),
        19: ("WWM_21", "Orbital vel. (m/s) "),
        20: ("WWM_22", "RMS Orbital vel. (m/s) "),
        21: ("WWM_23", "Bottom excursion period (sec?) "),
        22: ("WWM_24", "Bottom wave period (sec) "),
        23: ("WWM_25", "Uresell number based on peak period "),
        24: ("WWM_26", "Friction velocity (m/s?) "),
        25: ("WWM_27", "Charnock coefficient "),
        26: ("WWM_28", "Rougness length "),
        27: ("WWM_energy_dir", "WWM_energy vector"),
        28: ("wave-force", "Wave force vector (m.s-2) computed by wwm @side centers and whole levels"),
    }


def gen():
    """
    gen output options
    """
    return {
        1: ("GEN_1", "1st tracer"),
        2: ("GEN_2", "2nd tracer"),
    }


def age():
    """
    age output options
    """
    return {
        1: ("AGE_1", "Indices from \"1\" to \"ntracer_age/2\"; [days]"),
        2: ("AGE_2", "Indices from \"1\" to \"ntracer_age/2\"; [days]"),
    }


def sed():
    """
    sed output options
    """
    return {
        1: ("SED_depth_change", "bottom depth _change_ from init. condition (m)"),
        2: ("SED_D50", " Bed median grain size in the active layer (mm)"),
        3: ("SED_bed_stress", " Bottom shear stress (Pa)"),
        4: ("SED_bed_roughness", " Bottom roughness lenghth (mm)"),
        5: ("SED_TSC", "total suspended concentration (g/L)"),
        6: ("bed_thickness", " total bed thickness @elem (m)"),
        7: ("bed_age", " total bed age over all layers @elem (sec)"),
        8: ("z0st", " Sediment transport roughness length @elem (m) (z0st_elem)"),
        9: ("z0cr", "current-ripples roughness length @elem (m) (z0cr_elem)"),
        10: ("z0sw", "sand-waves roughness length (m) @elem (z0sw_elem)"),
        11: ("z0wr", "wave-ripples roughness length @elem (m) (z0wr_elem)"),
        12: ("SED3D_1", "conc. of 1st class (one output need by each class) [g/L]"),
        13: ("SED_bdld_1", "Bedload transport rate vector (kg.m-1.s-1) for 1st tracer (one output need by tracer)"),
        14: ("SED_bedfrac_1", "Bed fraction 1st tracer (one output need by each class) [-]"),
        15: ("SED3D_2", "conc. of 2nd class"),
        16: ("SED_bdld_2", "Bedload transport of 2nd class"),
        17: ("SED_bedfrac_3", "Bed fraction of 2nd class"),
    }


def eco():
    """
    EcoSim output options
    """
    return {
        1: ("ECO_1", "EcoSim outputs")
    }


def icm():
    """
    ICM output options
    """
    return {
        1: ("ICM_Chl", "Chlorophyll"),
        2: ("ICM_pH", "PH values (ICM_PH on)"),
        3: ("ICM_PrmPrdt", "ICM primary production @elem [gC/m^3/day]"),
        4: ("ICM_DIN", "ICM totoal inorganic nitrogen (DIN) @elem [gN/m^3]"),
        5: ("ICM_PON", "ICM paticulate organic nitrogen (PON) @elem [gN/m^3]"),
        6: ("ICM_SED_BENDOC", "ICM bed sediment flux arrays: SED_BENDOC (output name:ICM_SED_BENDOC) @elem [gC/(m^2 day)]"),
        7: ("ICM_SED_BENNH4", "ICM bed sediment flux arrays: SED_BENNH4 (output name:ICM_SED_BENNH4) @elem [gC/(m^2 day)]"),
        8: ("ICM_SED_BENNO3", "ICM bed sediment flux arrays: SED_BENNO3 (output name:ICM_SED_BENNO3)@elem [gC/(m^2 day)]"),
        9: ("ICM_SED_BENPO4", "ICM bed sediment flux arrays: SED_BENPO4 (output name:ICM_SED_BENPO4) @elem [gC/(m^2 day)]"),
        10: ("ICM_SED_BENCOD", "ICM bed sediment flux arrays: SED_BENCOD (output name:ICM_SED_BENCOD) @elem [gC/(m^2 day)]"),
        11: ("ICM_SED_BENDO", "ICM bed sediment flux arrays: SED_BENDO (output name:ICM_SED_BENDO) @elem [gC/(m^2 day)]"),
        12: ("ICM_SED_BENSA", "ICM bed sediment flux arrays: SED_BENSA (output name:ICM_SED_BENSA) @elem [gC/(m^2 day)]"),
        13: ("ICM_lfsav", "ICM SAV leaf biomass @elem [gC/m^3] (k=1 is surface)"),
        14: ("ICM_stsav", "ICM SAV stem biomass @elem [gC/m^3]"),
        15: ("ICM_rtsav", "ICM SAV root biomass @elem [gC/m^3]"),
        16: ("ICM_tlfsav", "ICM SAV total leaf biomass @elem [gC/m^2]"),
        17: ("ICM_tstsav", "ICM SAV total stem biomass @elem [gC/m^2]"),
        18: ("ICM_trtsav", "ICM SAV total root biomass @elem [gC/m^2]"),
        19: ("ICM_hcansav", "ICM SAV canopy height @elem [m]"),
        20: ("ICM_CNH4", "bottom NH4 conc"),
        21: ("ICM_CNH3", "bottom NO3 conc"),
        22: ("ICM_CPIP", "bottom P conc"),
        23: ("ICM_CPOS", "bottom Si conc"),
        24: ("ICM_CCH4", "bottom CH4 conc"),
        25: ("ICM_CSO4", "bottom SO4 conc"),
        26: ("ICM_CH2S", "bottom H2S conc"),
        27: ("ICM_SEDPON1", "bottom PON g1 conc"),
        28: ("ICM_SEDPON2", "bottom PON g2 conc"),
        29: ("ICM_SEDPON3", "bottom PON g3 conc"),
        30: ("ICM_SEDPOP1", "bottom POP g1 conc"),
        31: ("ICM_SEDPOP2", "bottom POP g2 conc"),
        32: ("ICM_SEDPOP3", "bottom POP g3 conc"),
        33: ("ICM_SEDPOC1", "bottom POC g1 conc"),
        34: ("ICM_SEDPOC2", "bottom POC g2 conc"),
        35: ("ICM_SEDPOC3", "bottom POC g3 conc"),
        36: ("ICM_EROH2S", "erosion flux H2S"),
        37: ("ICM_EROLPOC", "ersoion flux LPOC"),
        38: ("ICM_ERORPOC", "ersoion flux RPOC"),
        39: ("ICM_DO_consumption", "DO consumption"),
        40: ("ICM_GP1", "PB growth #1"),
        41: ("ICM_GP2", "PB growth #2"),
        42: ("ICM_GP3", "PB growth #3"),
        43: ("ICM_1", "Zoo. #1"),
        44: ("ICM_2", "Zoo. #2"),
        45: ("ICM_3", "phyto #1"),
        46: ("ICM_4", "phyto #2"),
        47: ("ICM_5", "phyto #3"),
        48: ("ICM_6", "RPOC"),
        49: ("ICM_7", "LPOC"),
        50: ("ICM_8", "DOC"),
        51: ("ICM_9", "RPON"),
        52: ("ICM_10", "LPON"),
        53: ("ICM_11", "DON"),
        54: ("ICM_12", "NH4"),
        55: ("ICM_13", "NO3"),
        56: ("ICM_14", "RPOP"),
        57: ("ICM_15", "LPOP"),
        58: ("ICM_16", "DOP"),
        59: ("ICM_17", "PO4t"),
        60: ("ICM_18", "Si- biogenic"),
        61: ("ICM_19", "available Si"),
        62: ("ICM_20", "COD: Chemical oxygen demand"),
        63: ("ICM_21", "DO"),
        64: ("ICM_22", "TIC"),
        65: ("ICM_23", "ALK"),
        66: ("ICM_24", "CA"),
        67: ("ICM_25", "CACO3"),
    }


def cos():
    """
    CoSINE output options
    """
    return {
        1: ("COS_1", "COS_1"),
        2: ("COS_2", "COS_2"),
        3: ("COS_3", "COS_3"),
        4: ("COS_4", "COS_4"),
        5: ("COS_5", "COS_5"),
        6: ("COS_6", "COS_6"),
        7: ("COS_7", "COS_7"),
        8: ("COS_8", "COS_8"),
        9: ("COS_9", "COS_9"),
        10: ("COS_10", "COS_10"),
        11: ("COS_11", "COS_11"),
        12: ("COS_12", "COS_12"),
        13: ("COS_13", "COS_13"),
    }


def fib():
    """
    Fecal indicating bacteria output options
    """
    return {
        1: ("FIB_1", "FIB_1")
    }


def sed2d():
    """
    SED2D output options
    """
    return {
        1: ("SED2D_depth_change", "bottom depth _change_ from init. condition (m)"),
        2: ("SED2D_Cd", "drag coefficient used in transport formulae"),
        3: ("SED2D_cflsed", "Courant number (b.qtot.dt / h.dx)"),
        4: ("SED2D_d50", "Top layer d50 (m)"),
        5: ("SED2D_total_transport", "total transport rate vector (kg/m/s)"),
        6: ("SED2D_susp_load", "suspended tranport rate vector (kg/m/s)"),
        7: ("SED2D_bed_load", "bedload transport rate vector (kg/m/s)"),
        8: ("SED2D_average_transport", "time averaged total transport rate vector (kg/m/s)"),
        9: ("SED2D_bottom_slope", "bottom slope vector (m/m); negative uphill"),
        10: ("z0eq2d", "Total roughness length @elem (m) (z0eq)"),
        11: ("z0cr2d", "current-ripples roughness length @elem (m) (z0cr)"),
        12: ("z0sw2d", "sand-waves roughness length @elem (m) (z0sw)"),
        13: ("z0wr2d", "wave-ripples roughness length @elem (m) (z0wr)"),
    }


def mar():
    """
    marsh output options
    """
    return {
        1: ("marsh_flag", "marsh_flag"),
    }


def ice():
    """
    ice output options
    """
    return {
        1: ("ICE_velocity", "ice advective velcoity vector [m/s]"),
        2: ("ICE_strain_rate", "strain rate @ elem [1/sec]"),
        3: ("ICE_net_heat_flux", "net heat flux to ocean (>0 warm up SST) [W/m/m]"),
        4: ("ICE_fresh_water_flux", "net fresh water flux to ocean (>0 freshens up SSS) [kg/s/m/m]"),
        5: ("ICE_top_T", "ice temperature [C] at air-ice interface"),
        6: ("ICE_tracer_1", "ice volume [m]"),
        7: ("ICE_tracer_2", "ice concentration [-]"),
        8: ("ICE_tracer_3", "snow volume [m]"),
    }


def ana():
    return {
        1: ("ANA_air_pres_grad_x", "x-component of 𝛁air_pres/ρ0 [m/s/s]"),
        2: ("ANA_air_pres_grad_y", "y-component of 𝛁air_pres/ρ0 [m/s/s]"),
        3: ("ANA_tide_pot_grad_x", "α*g*𝛁Ψ [m/s/s] (gradient of tidal potential)"),
        4: ("ANA_tide_pot_grad_y", "α*g*𝛁Ψ [m/s/s]"),
        5: ("ANA_hor_viscosity_x", "𝛁·(μ𝛁u) [m/s/s] (horizontal viscosity)"),
        6: ("ANA_hor_viscosity_y", "𝛁·(μ𝛁u) [m/s/s]"),
        7: ("ANA_bclinic_force_x", "-g/rho0* ∫_z^η dr_dx dz  [m/s/s] (b-clinic gradient)"),
        8: ("ANA_bclinic_force_y", "-g/rho0* ∫_z^η dr_dy dz  [m/s/s]"),
        9: ("ANA_vert_viscosity_x", "d (ν du/dz)/dz [m/s/s] - no vegetation effects (vertical viscosity)"),
        10: ("ANA_vert_viscosity_y", "d (ν dv/dz)/dz [m/s/s] - no vegetation effects"),
        11: ("ANA_mom_advection_x", "(u·𝛁) u [m/s/s] (momentum advection)"),
        12: ("ANA_mom_advection_y", "(u·𝛁) u [m/s/s]"),
        13: ("ANA_Richardson", "gradient Richardson number [-]"),
        14: ("ANA_transport_min_dt_elem", "min time step at each element over all subcycles in horizontal transport solver [s]  "),
    }


def outputs(parser):

    outputs = {
        'hyd': hydro(),
        'wwm': wwm(),
        'gen': gen(),
        'age': age(),
        'sed': sed(),
        'eco': eco(),
        'icm': icm(),
        'cos': cos(),
        'fib': fib(),
        'sed2d': sed2d(),
        'mar': mar(),
        'ice': ice(),
        'ana': ana(),
    }
    for short_name, output in outputs.items():
        for id, (long_name, help_msg) in output.items():
            parser.add_argument(
                f"--{long_name.lower().replace('_', '-')}",
                f"-{short_name}{id}",
                dest=f"{short_name}{id}",
                help=help_msg,
                action='store_true'
                )


def timezone(parser):
    parser.add_argument("--timezone")


def get_parser(runtype=None, description=None):
    parser = argparse.ArgumentParser(description=description)
    mesh(parser)
    if runtype is not None:
        if runtype == "tidal":
            tidal_run(parser)
        elif runtype == 'best_track':
            best_track(parser)
    tidal_constituents(parser)
    timezone(parser)
    output_directory(parser)
    allow_overwrite(parser)
    generate_only(parser)
    # parse stations from a file
    msg = "File containing list of stations for outputs. It will parse "
    msg += "the stations below the NOUTE, NOUTV and NOUTM keywords for "
    msg += "their respective stations list."
    parser.add_argument(
        "--stations-file",
        help=msg
        )
    outputs(parser)
    log_level(parser)
    server(parser)
    return parser
