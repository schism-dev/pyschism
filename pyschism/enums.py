from enum import Enum


class Stratification(Enum):
    BAROCLINIC = 0
    BAROTROPIC = 1
    @staticmethod
    def keys():
        return list(map(lambda c: c.name, Stratification))
    @classmethod
    def _missing_(cls, name):
        if isinstance(name, str):
            if name.upper() in cls.keys():
                return cls[name.upper()]
        raise ValueError(f'Argument {name} is not a valid Stratification type.')


class Coriolis(Enum):
    AUTO = 1
    CORICOEFF = 0
    RLATITUDE = -1


class IofWetdryVariables(Enum):
    wetdry_node = "wetdry_node"
    wetdry_elem = "wetdry_elem"
    wetdry_side = "wetdry_side"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid WETDRY output variable.")


class IofZcorVariables(Enum):
    zcor = "zcor"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid ZCOR output variable.")


class IofHydroVariables(Enum):
    elev = "elev"
    air_pressure = "air_pressure"
    air_temperature = "air_temperature"
    specific_humidity = "specific_humidity"
    solar_radiation = "solar_radiation"
    sensible_flux = "sensible_flux"
    latent_heat = "latent_heat"
    upward_longwave = "upward_longwave"
    downward_longwave = "downward_longwave"
    total_heat_flux = "total_heat_flux"
    evaporation = "evaporation"
    precipitation = "precipitation"
    bottom_stress = "bottom_stress"
    wind_speed = "wind_speed"
    wind_stress = "wind_stress"
    dahv = "dahv"
    vertical_velocity = "vertical_velocity"
    temp = "temp"
    salt = "salt"
    water_density = "water_density"
    diffusivity = "diffusivity"
    viscosity = "viscosity"
    TKE = "TKE"
    mixing_length = "mixing_length"
    hvel = "hvel"
    hvel_side = "hvel_side"
    wvel_elem = "wvel_elem"
    temp_elem = "temp_elem"
    salt_elem = "salt_elem"
    pressure_gradient = "pressure_gradient"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid HYDRO output variable.")


class IofDvdVariables(Enum):
    DVD_1 = "DVD_1"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid DVD output variable.")


class IofWwmVariables(Enum):
    WWM_1 = "WWM_1"
    WWM_2 = "WWM_2"
    WWM_3 = "WWM_3"
    WWM_4 = "WWM_4"
    WWM_5 = "WWM_5"
    WWM_6 = "WWM_6"
    WWM_9 = "WWM_9"
    WWM_10 = "WWM_10"
    WWM_11 = "WWM_11"
    WWM_12 = "WWM_12"
    WWM_13 = "WWM_13"
    WWM_14 = "WWM_14"
    WWM_15 = "WWM_15"
    WWM_16 = "WWM_16"
    WWM_17 = "WWM_17"
    WWM_18 = "WWM_18"
    WWM_19 = "WWM_19"
    WWM_20 = "WWM_20"
    WWM_21 = "WWM_21"
    WWM_22 = "WWM_22"
    WWM_23 = "WWM_23"
    WWM_24 = "WWM_24"
    WWM_25 = "WWM_25"
    WWM_26 = "WWM_26"
    WWM_27 = "WWM_27"
    WWM_28 = "WWM_28"
    WWM_energy_dir = "WWM_energy_dir"
    wave_force = "wave_force"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid WWM output variable.")


class IofGenVariables(Enum):
    GEN_1 = "GEN_1"
    GEN_2 = "GEN_2"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid GEN output variable.")


class IofAgeVariables(Enum):
    AGE_1 = "AGE_1"
    AGE_2 = "AGE_2"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid AGE output variable.")


class IofSedVariables(Enum):
    SED_depth_change = "SED_depth_change"
    SED_D50 = "SED_D50"
    SED_bed_stress = "SED_bed_stress"
    SED_bed_roughness = "SED_bed_roughness"
    SED_TSC = "SED_TSC"
    bed_thickness = "bed_thickness"
    bed_age = "bed_age"
    z0st = "z0st"
    z0cr = "z0cr"
    z0sw = "z0sw"
    z0wr = "z0wr"
    SED3D_1 = "SED3D_1"
    SED_bdld_1 = "SED_bdld_1"
    SED_bedfrac_1 = "SED_bedfrac_1"
    SED3D_2 = "SED3D_2"
    SED_bdld_2 = "SED_bdld_2"
    SED_bedfrac_3 = "SED_bedfrac_3"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid SED output variable.")


class IofEcoVariables(Enum):
    ECO_1 = "ECO_1"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid ECO output variable.")


class IofIcmVariables(Enum):
    ICM_Chl = "ICM_Chl"
    ICM_pH = "ICM_pH"
    ICM_PrmPrdt = "ICM_PrmPrdt"
    ICM_DIN = "ICM_DIN"
    ICM_PON = "ICM_PON"
    ICM_SED_BENDOC = "ICM_SED_BENDOC"
    ICM_SED_BENNH4 = "ICM_SED_BENNH4"
    ICM_SED_BENNO3 = "ICM_SED_BENNO3"
    ICM_SED_BENPO4 = "ICM_SED_BENPO4"
    ICM_SED_BENCOD = "ICM_SED_BENCOD"
    ICM_SED_BENDO = "ICM_SED_BENDO"
    ICM_SED_BENSA = "ICM_SED_BENSA"
    ICM_lfsav = "ICM_lfsav"
    ICM_stsav = "ICM_stsav"
    ICM_rtsav = "ICM_rtsav"
    ICM_tlfsav = "ICM_tlfsav"
    ICM_tstsav = "ICM_tstsav"
    ICM_trtsav = "ICM_trtsav"
    ICM_hcansav = "ICM_hcansav"
    ICM_CNH4 = "ICM_CNH4"
    ICM_CNH3 = "ICM_CNH3"
    ICM_CPIP = "ICM_CPIP"
    ICM_CPOS = "ICM_CPOS"
    ICM_CCH4 = "ICM_CCH4"
    ICM_CSO4 = "ICM_CSO4"
    ICM_CH2S = "ICM_CH2S"
    ICM_SEDPON1 = "ICM_SEDPON1"
    ICM_SEDPON2 = "ICM_SEDPON2"
    ICM_SEDPON3 = "ICM_SEDPON3"
    ICM_SEDPOP1 = "ICM_SEDPOP1"
    ICM_SEDPOP2 = "ICM_SEDPOP2"
    ICM_SEDPOP3 = "ICM_SEDPOP3"
    ICM_SEDPOC1 = "ICM_SEDPOC1"
    ICM_SEDPOC2 = "ICM_SEDPOC2"
    ICM_SEDPOC3 = "ICM_SEDPOC3"
    ICM_EROH2S = "ICM_EROH2S"
    ICM_EROLPOC = "ICM_EROLPOC"
    ICM_ERORPOC = "ICM_ERORPOC"
    ICM_DO_consumption = "ICM_DO_consumption"
    ICM_GP1 = "ICM_GP1"
    ICM_GP2 = "ICM_GP2"
    ICM_GP3 = "ICM_GP3"
    ICM_1 = "ICM_1"
    ICM_2 = "ICM_2"
    ICM_3 = "ICM_3"
    ICM_4 = "ICM_4"
    ICM_5 = "ICM_5"
    ICM_6 = "ICM_6"
    ICM_7 = "ICM_7"
    ICM_8 = "ICM_8"
    ICM_9 = "ICM_9"
    ICM_10 = "ICM_10"
    ICM_11 = "ICM_11"
    ICM_12 = "ICM_12"
    ICM_13 = "ICM_13"
    ICM_14 = "ICM_14"
    ICM_15 = "ICM_15"
    ICM_16 = "ICM_16"
    ICM_17 = "ICM_17"
    ICM_18 = "ICM_18"
    ICM_19 = "ICM_19"
    ICM_20 = "ICM_20"
    ICM_21 = "ICM_21"
    ICM_22 = "ICM_22"
    ICM_23 = "ICM_23"
    ICM_24 = "ICM_24"
    ICM_25 = "I6CM_25"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid ICM output variable.")


class IofCosVariables(Enum):
    COS_1 = "COS_1"
    COS_2 = "COS_2"
    COS_3 = "COS_3"
    COS_4 = "COS_4"
    COS_5 = "COS_5"
    COS_6 = "COS_6"
    COS_7 = "COS_7"
    COS_8 = "COS_8"
    COS_9 = "COS_9"
    COS_10 = "COS_10"
    COS_11 = "COS_11"
    COS_12 = "COS_12"
    COS_13 = "COS_13"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid COS output variable.")


class IofFibVariables(Enum):
    FIB_1 = "FIB_1"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid FIB output variable.")


class IofSed2dVariables(Enum):
    SED2D_depth_change = "SED2D_depth_change"
    SED2D_drag_coefficient = "SED2D_drag_coefficient"
    SED2D_cflsed = "SED2D_cflsed"
    SED2D_d50 = "SED2D_d50"
    SED2D_total_transport = "SED2D_total_transport"
    SED2D_susp_load = "SED2D_susp_load"
    SED2D_bed_load = "SED2D_bed_load"
    SED2D_average_transport = "SED2D_average_transport"
    SED2D_bottom_slope = "SED2D_bottom_slope"
    z0eq = "z0eq"
    z0cr2d = "z0cr2d"
    z0sw2d = "z0sw2d"
    z0wr2d = "z0wr2d"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid SED2D output variable.")


class IofMarshVariables(Enum):
    marsh_flag = "marsh_flag"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid MARSH output variable.")


class IofIceVariables(Enum):
    ICE_velocity = "ICE_velocity"
    ICE_strain_rate = "ICE_strain_rate"
    ICE_net_heat_flux = "ICE_net_heat_flux"
    ICE_fresh_water_flux = "ICE_fresh_water_flux"
    ICE_top_T = "ICE_top_T"
    ICE_tracer_1 = "ICE_tracer_1"
    ICE_tracer_2 = "ICE_tracer_2"
    ICE_tracer_3 = "ICE_tracer_3"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid ICE output variable.")


class IofAnaVariables(Enum):
    ANA_air_pres_grad_x = "ANA_air_pres_grad_x"
    ANA_air_pres_grad_y = "ANA_air_pres_grad_y"
    ANA_tide_pot_grad_x = "ANA_tide_pot_grad_x"
    ANA_tide_pot_grad_y = "ANA_tide_pot_grad_y"
    ANA_hor_viscosity_x = "ANA_hor_viscosity_x"
    ANA_hor_viscosity_y = "ANA_hor_viscosity_y"
    ANA_bclinic_force_x = "ANA_bclinic_force_x"
    ANA_bclinic_force_y = "ANA_bclinic_force_y"
    ANA_vert_viscosity_x = "ANA_vert_viscosity_x"
    ANA_vert_viscosity_y = "ANA_vert_viscosity_y"
    ANA_mom_advection_x = "ANA_mom_advection_x"
    ANA_mom_advection_y = "ANA_mom_advection_y"
    ANA_Richardson = "ANA_Richardson"
    ANA_transport_min_dt_elem = "ANA_transport_min_dt_elem"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid ANA output variable.")


class StationOutputVariables(Enum):
    """Enumeration for stations output variables"""

    ELEVATION = "elev"
    AIR_PRESSURE = "air_pressure"
    WINDX = "windx"
    WINDY = "windy"
    TEMPERATURE = "T"
    SALINITY = "S"
    U = "u"
    V = "v"
    W = "w"

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f"{name} is not a valid station output variable.")


class StationOutputIndex(Enum):
    """Indexing for stations output variables"""

    ELEVATION = 0
    AIR_PRESSURE = 1
    WINDX = 2
    WINDY = 3
    TEMPERATURE = 4
    SALINITY = 5
    U = 6
    V = 7
    W = 8


class OutputVariableUnit(Enum):
    ELEVATION = "meters"


class NWSType(Enum):
    """Atmospheric forcing type required by param.nml"""

    PARAMETRIC = -1
    TIME_HISTORY = 1
    CLIMATE_AND_FORECAST = 2
    HEAT_CONSERVATION = 3
    UNSTRUCTURED = 4

    @classmethod
    def _missing_(self, name):
        raise ValueError(f"{name} is not a valid NWS type.")


# class ForecastProduct(Enum):
#     GDAS = 'gdas_0p25'
#     GDAS_0P25 = 'gdas_0p25'
#     GFS = 'gfs_0p25_1hr'
#     GFS_0P25 = 'gfs_0p25'
#     GFS_0P25_1HR = 'gfs_0p25_1hr'
#     GFS_0P50 = 'gfs_0p50'
#     GFS_1P00 = 'gfs_1p00'

#     @classmethod
#     def _missing_(self, name):
#         ValueError(f'{name} is not a known atmospheric forecast product for '
#                    'air.')


class GFSProduct(Enum):
    GFS_0P25 = "gfs_0p25"
    GFS_0P25_1HR = "gfs_0p25_1hr"
    GFS_0P50 = "gfs_0p50"
    GFS_1P00 = "gfs_1p00"


class iof_hydro(Enum):
    elev = 0
    air_pressure = 1
    air_temperature = 2
    specific_humidity = 3
    solar_radiation = 4
    sensible_flux = 5
    latent_heat = 6
    upward_longwave = 7
    downward_longwave = 8
    total_heat_flux = 9
    evaporation = 10
    precipitation = 11
    bottom_stress = 12
    wind_speed = 13
    wind_stress = 14
    dahv = 15
    vertical_velocity = 16
    temp = 17
    salt = 18
    water_density = 19
    diffusivity = 20
    viscosity = 21
    TKE = 22
    mixing_length = 23
    hvel = 24
    hvel_side = 25
    wvel_elem = 26
    temp_elem = 27
    salt_elem = 28
    pressure_gradient = 29


class iof_dvd(Enum):
    DVD_1 = 0


class iof_wwm(Enum):
    WWM_1 = 0
    WWM_2 = 1
    WWM_3 = 2
    WWM_4 = 3
    WWM_5 = 4
    WWM_6 = 5
    WWM_9 = 6
    WWM_10 = 7
    WWM_11 = 8
    WWM_12 = 9
    WWM_13 = 10
    WWM_14 = 11
    WWM_15 = 12
    WWM_16 = 13
    WWM_17 = 14
    WWM_18 = 15
    WWM_19 = 16
    WWM_20 = 17
    WWM_21 = 18
    WWM_22 = 19
    WWM_23 = 20
    WWM_24 = 21
    WWM_25 = 22
    WWM_26 = 23
    WWM_27 = 24
    WWM_28 = 25
    WWM_energy_dir = 26
    wave_force = 27


class iof_gen(Enum):
    GEN_1 = 0
    GEN_2 = 1


class iof_age(Enum):
    AGE_1 = 0
    AGE_2 = 1


class iof_sed(Enum):
    SED_depth_change = 0
    SED_D50 = 1
    SED_bed_stress = 2
    SED_bed_roughness = 3
    SED_TSC = 4
    bed_thickness = 5
    bed_age = 6
    z0st = 7
    z0cr = 8
    z0sw = 9
    z0wr = 10
    SED3D_1 = 11
    SED_bdld_1 = 12
    SED_bedfrac_1 = 13
    SED3D_2 = 14
    SED_bdld_2 = 15
    SED_bedfrac_3 = 16


class iof_eco(Enum):
    ECO_1 = 0


class iof_icm(Enum):
    ICM_Chl = 0
    ICM_pH = 1
    ICM_PrmPrdt = 2
    ICM_DIN = 3
    ICM_PON = 4
    ICM_SED_BENDOC = 5
    ICM_SED_BENNH4 = 6
    ICM_SED_BENNO3 = 7
    ICM_SED_BENPO4 = 8
    ICM_SED_BENCOD = 9
    ICM_SED_BENDO = 10
    ICM_SED_BENSA = 11
    ICM_lfsav = 12
    ICM_stsav = 13
    ICM_rtsav = 14
    ICM_tlfsav = 15
    ICM_tstsav = 16
    ICM_trtsav = 17
    ICM_hcansav = 18
    ICM_CNH4 = 19
    ICM_CNH3 = 20
    ICM_CPIP = 21
    ICM_CPOS = 22
    ICM_CCH4 = 23
    ICM_CSO4 = 24
    ICM_CH2S = 25
    ICM_SEDPON1 = 26
    ICM_SEDPON2 = 27
    ICM_SEDPON3 = 28
    ICM_SEDPOP1 = 29
    ICM_SEDPOP2 = 30
    ICM_SEDPOP3 = 31
    ICM_SEDPOC1 = 32
    ICM_SEDPOC2 = 33
    ICM_SEDPOC3 = 34
    ICM_EROH2S = 35
    ICM_EROLPOC = 36
    ICM_ERORPOC = 37
    ICM_DO_consumption = 38
    ICM_GP1 = 39
    ICM_GP2 = 40
    ICM_GP3 = 41
    ICM_1 = 42
    ICM_2 = 43
    ICM_3 = 44
    ICM_4 = 45
    ICM_5 = 46
    ICM_6 = 47
    ICM_7 = 48
    ICM_8 = 49
    ICM_9 = 50
    ICM_10 = 51
    ICM_11 = 52
    ICM_12 = 53
    ICM_13 = 54
    ICM_14 = 55
    ICM_15 = 56
    ICM_16 = 57
    ICM_17 = 58
    ICM_18 = 59
    ICM_19 = 60
    ICM_20 = 61
    ICM_21 = 62
    ICM_22 = 63
    ICM_23 = 64
    ICM_24 = 65
    ICM_25 = 66


class iof_cos(Enum):
    COS_1 = 0
    COS_2 = 1
    COS_3 = 2
    COS_4 = 3
    COS_5 = 4
    COS_6 = 5
    COS_7 = 6
    COS_8 = 7
    COS_9 = 8
    COS_10 = 9
    COS_11 = 10
    COS_12 = 11
    COS_13 = 12


class iof_fib(Enum):
    FIB_1 = 0


class iof_sed2d(Enum):
    SED2D_depth_change = 0
    SED2D_drag_coefficient = 1
    SED2D_cflsed = 2
    SED2D_d50 = 3
    SED2D_total_transport = 4
    SED2D_susp_load = 5
    SED2D_bed_load = 6
    SED2D_average_transport = 7
    SED2D_bottom_slope = 8
    z0eq = 9
    z0cr2d = 10
    z0sw2d = 11
    z0wr2d = 12


class iof_marsh(Enum):
    marsh_flag = 0


class iof_ice(Enum):
    ICE_velocity = 0
    ICE_strain_rate = 1
    ICE_net_heat_flux = 2
    ICE_fresh_water_flux = 3
    ICE_top_T = 4
    ICE_tracer_1 = 5
    ICE_tracer_2 = 6
    ICE_tracer_3 = 7


class iof_ana(Enum):
    ANA_air_pres_grad_x = 0
    ANA_air_pres_grad_y = 1
    ANA_tide_pot_grad_x = 2
    ANA_tide_pot_grad_y = 3
    ANA_hor_viscosity_x = 4
    ANA_hor_viscosity_y = 5
    ANA_bclinic_force_x = 6
    ANA_bclinic_force_y = 7
    ANA_vert_viscosity_x = 8
    ANA_vert_viscosity_y = 9
    ANA_mom_advection_x = 10
    ANA_mom_advection_y = 11
    ANA_Richardson = 12
    ANA_transport_min_dt_elem = 13


class SchoutType(Enum):
    iof_hydro = iof_hydro
    iof_dvd = iof_dvd
    iof_wwm = iof_wwm
    iof_gen = iof_gen
    iof_age = iof_age
    iof_sed = iof_sed
    iof_eco = iof_eco
    iof_icm = iof_icm
    iof_cos = iof_cos
    iof_fib = iof_fib
    iof_sed2d = iof_sed2d
    iof_marsh = iof_marsh
    iof_ice = iof_ice
    iof_ana = iof_ana


class SchoutIofType(str, Enum):
    iof_hydro = "iof_hydro"
    iof_dvd = "iof_dvd"
    iof_wwm = "iof_wwm"
    iof_gen = "iof_gen"
    iof_age = "iof_age"
    iof_sed = "iof_sed"
    iof_eco = "iof_eco"
    iof_icm = "iof_icm"
    iof_cos = "iof_cos"
    iof_fib = "iof_fib"
    iof_sed2d = "iof_sed2d"
    iof_marsh = "iof_marsh"
    iof_ice = "iof_ice"
    iof_ana = "iof_ana"

    @classmethod
    def _missing_(self, name):
        raise ValueError(f"{name} is not a valid SCHISM output type.")


class OutputVariableShortName(Enum):

    ELEVATION = "elev"

    @classmethod
    def _missing_(self, name):
        raise ValueError(f"{name} is not a valid SCHISM output variable short name")


class NationalWaterModelDataSource(Enum):
    AWS = "AWS"
    FTP = "FTP"
    NOMADS = "NOMADS"

    @classmethod
    def _missing_(self, name):
        raise ValueError(f"{name} is not a valid National Water Model data source.")


class Sflux1Types(Enum):

    from pyschism.forcing import nws
    # GDAS = GDAS
    # GDAS_0P25 = GDAS
    GFS = nws.GFS
    GFS_0P25 = nws.GFS
    GFS_0P25_1HR = nws.GFS
    GFS_0P50 = nws.GFS
    GFS_1P00 = nws.GFS

    @classmethod
    def _missing_(cls, name):
        f = [
            f"{name} is not a valid sflux_1 type. Valid values are: ",
        ]
        for sflux_type in cls:
            f.append(sflux_type.name.lower())
        f.append(".")
        raise ValueError("".join(f))


class Sflux2Types(Enum):

    from pyschism.forcing import nws
    HRRR = nws.HRRR

    @classmethod
    def _missing_(cls, name):
        f = [
            f"{name} is not a valid sflux_2 type. Valid values are: ",
        ]
        for sflux_type in cls:
            f.append(sflux_type.name.lower())
        f.append(".")
        raise ValueError("".join(f))
