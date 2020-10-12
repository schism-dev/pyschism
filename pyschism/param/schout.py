import pathlib
from enum import Enum

import f90nml  # type: ignore[import]

PARAM_TEMPLATE = pathlib.Path(__file__).parent / 'param.nml.template'
PARAM_DEFAULTS = f90nml.read(PARAM_TEMPLATE)['schout']


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
    I6CM_25 = 66


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


class SchoutStrType(Enum):
    iof_hydro = 'iof_hydro'
    iof_dvd = 'iof_dvd'
    iof_wwm = 'iof_wwm'
    iof_gen = 'iof_gen'
    iof_age = 'iof_age'
    iof_sed = 'iof_sed'
    iof_eco = 'iof_eco'
    iof_icm = 'iof_icm'
    iof_cos = 'iof_cos'
    iof_fib = 'iof_fib'
    iof_sed2d = 'iof_sed2d'
    iof_marsh = 'iof_marsh'
    iof_ice = 'iof_ice'
    iof_ana = 'iof_ana'


allowed_vars = [
    "elev",
    "air_pressure",
    "air_temperature",
    "specific_humidity",
    "solar_radiation",
    "sensible_flux",
    "latent_heat",
    "upward_longwave",
    "downward_longwave",
    "total_heat_flux",
    "evaporation",
    "precipitation",
    "bottom_stress",
    "wind_speed",
    "wind_stress",
    "dahv",
    "vertical_velocity",
    "temp",
    "salt",
    "water_density",
    "diffusivity",
    "viscosity",
    "TKE",
    "mixing_length",
    "hvel",
    "hvel_side",
    "wvel_elem",
    "temp_elem",
    "salt_elem",
    "pressure_gradient",
    "DVD_1",
    "WWM_1",
    "WWM_2",
    "WWM_3",
    "WWM_4",
    "WWM_5",
    "WWM_6",
    "WWM_9",
    "WWM_10",
    "WWM_11",
    "WWM_12",
    "WWM_13",
    "WWM_14",
    "WWM_15",
    "WWM_16",
    "WWM_17",
    "WWM_18",
    "WWM_19",
    "WWM_20",
    "WWM_21",
    "WWM_22",
    "WWM_23",
    "WWM_24",
    "WWM_25",
    "WWM_26",
    "WWM_27",
    "WWM_28",
    "WWM_energy_dir",
    "wave_force",
    "GEN_1",
    "GEN_2",
    "AGE_1",
    "AGE_2",
    "SED_depth_change",
    "SED_D50",
    "SED_bed_stress",
    "SED_bed_roughness",
    "SED_TSC",
    "bed_thickness",
    "bed_age",
    "z0st",
    "z0cr",
    "z0sw",
    "z0wr",
    "SED3D_1",
    "SED_bdld_1",
    "SED_bedfrac_1",
    "SED3D_2",
    "SED_bdld_2",
    "SED_bedfrac_3",
    "ECO_1",
    "ICM_Chl",
    "ICM_pH",
    "ICM_PrmPrdt",
    "ICM_DIN",
    "ICM_PON",
    "ICM_SED_BENDOC",
    "ICM_SED_BENNH4",
    "ICM_SED_BENNO3",
    "ICM_SED_BENPO4",
    "ICM_SED_BENCOD",
    "ICM_SED_BENDO",
    "ICM_SED_BENSA",
    "ICM_lfsav",
    "ICM_stsav",
    "ICM_rtsav",
    "ICM_tlfsav",
    "ICM_tstsav",
    "ICM_trtsav",
    "ICM_hcansav",
    "ICM_CNH4",
    "ICM_CNH3",
    "ICM_CPIP",
    "ICM_CPOS",
    "ICM_CCH4",
    "ICM_CSO4",
    "ICM_CH2S",
    "ICM_SEDPON1",
    "ICM_SEDPON2",
    "ICM_SEDPON3",
    "ICM_SEDPOP1",
    "ICM_SEDPOP2",
    "ICM_SEDPOP3",
    "ICM_SEDPOC1",
    "ICM_SEDPOC2",
    "ICM_SEDPOC3",
    "ICM_EROH2S",
    "ICM_EROLPOC",
    "ICM_ERORPOC",
    "ICM_DO_consumption",
    "ICM_GP1",
    "ICM_GP2",
    "ICM_GP3",
    "ICM_1",
    "ICM_2",
    "ICM_3",
    "ICM_4",
    "ICM_5",
    "ICM_6",
    "ICM_7",
    "ICM_8",
    "ICM_9",
    "ICM_10",
    "ICM_11",
    "ICM_12",
    "ICM_13",
    "ICM_14",
    "ICM_15",
    "ICM_16",
    "ICM_17",
    "ICM_18",
    "ICM_19",
    "ICM_20",
    "ICM_21",
    "ICM_22",
    "ICM_23",
    "ICM_24",
    "I6CM_25",
    "COS_1",
    "COS_2",
    "COS_3",
    "COS_4",
    "COS_5",
    "COS_6",
    "COS_7",
    "COS_8",
    "COS_9",
    "COS_10",
    "COS_11",
    "COS_12",
    "COS_13",
    "FIB_1",
    "SED2D_depth_change",
    "SED2D_drag_coefficient",
    "SED2D_cflsed",
    "SED2D_d50",
    "SED2D_total_transport",
    "SED2D_susp_load",
    "SED2D_bed_load",
    "SED2D_average_transport",
    "SED2D_bottom_slope",
    "z0eq",
    "z0cr2d",
    "z0sw2d",
    "z0wr2d",
    "marsh_flag",
    "ICE_velocity",
    "ICE_strain_rate",
    "ICE_net_heat_flux",
    "ICE_fresh_water_flux",
    "ICE_top_T",
    "ICE_tracer_1",
    "ICE_tracer_2",
    "ICE_tracer_3",
    "ANA_air_pres_grad_x",
    "ANA_air_pres_grad_y",
    "ANA_tide_pot_grad_x",
    "ANA_tide_pot_grad_y",
    "ANA_hor_viscosity_x",
    "ANA_hor_viscosity_y",
    "ANA_bclinic_force_x",
    "ANA_bclinic_force_y",
    "ANA_vert_viscosity_x",
    "ANA_vert_viscosity_y",
    "ANA_mom_advection_x",
    "ANA_mom_advection_y",
    "ANA_Richardson",
    "ANA_transport_min_dt_elem",
]


class SCHOUT:
    """ Provides error checking implementation for SCHOUT group """

    def __init__(self, **outputs):
        self.__schout: dict = {}
        for key, value in PARAM_DEFAULTS.items():
            if isinstance(value, list):
                self.__schout[key] = len(value)*[0]
            else:
                self.__schout[key] = None

        self.__iof_hydro = self['iof_hydro']
        self.__iof_dvd = self['iof_dvd']
        self.__iof_wwm = self['iof_wwm']
        self.__iof_gen = self['iof_gen']
        self.__iof_age = self['iof_age']
        self.__iof_sed = self['iof_sed']
        self.__iof_eco = self['iof_eco']
        self.__iof_icm = self['iof_icm']
        self.__iof_cos = self['iof_cos']
        self.__iof_fib = self['iof_fib']
        self.__iof_sed2d = self['iof_sed2d']
        self.__iof_marsh = self['iof_marsh']
        self.__iof_ice = self['iof_ice']
        self.__iof_ana = self['iof_ana']

        for var in outputs.copy():
            if var.lower() in allowed_vars:
                exec(f"self.{var.lower()}={bool(outputs.pop(var))}")

        # check for remaining arguments.
        if len(outputs) > 0:
            raise TypeError(f'Unknown outputs: {list(outputs.keys())}')

    def __getitem__(self, key):
        return self.__schout[key]

    def __setitem__(self, key, value):
        self.__schout[key] = value

    def __iter__(self):
        for key, value in self.__schout.items():
            yield key, value

    @property
    def nhot(self):
        return self['nhot']

    @nhot.setter
    def nhot(self, nhot):
        self['nhot'] = nhot

    @property
    def nhot_write(self):
        return self['nhot_write']

    @nhot_write.setter
    def nhot_write(self, nhot_write):
        self['nhot_write'] = nhot_write

    @property
    def iout_sta(self):
        return self['iout_sta']

    @iout_sta.setter
    def iout_sta(self, iout_sta):
        assert iout_sta in [0, 1]
        self['iout_sta'] = iout_sta

    @property
    def nspool_sta(self):
        return self['nspool_sta']

    @nspool_sta.setter
    def nspool_sta(self, nspool_sta: int):
        assert isinstance(nspool_sta, int), 'Can only be set using an int'
        self['nspool_sta'] = nspool_sta

    @property
    def elev(self):
        return bool(self.__iof_hydro[iof_hydro.elev.value])

    @elev.setter
    def elev(self, elev):
        assert isinstance(elev, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.elev.value] = int(elev)

    @property
    def air_pressure(self):
        return bool(self.__iof_hydro[iof_hydro.air_pressure.value])

    @air_pressure.setter
    def air_pressure(self, air_pressure):
        assert isinstance(air_pressure, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.air_pressure.value] = int(air_pressure)

    @property
    def air_temperature(self):
        return bool(self.__iof_hydro[iof_hydro.air_temperature.value])

    @air_temperature.setter
    def air_temperature(self, air_temperature):
        assert isinstance(air_temperature, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.air_temperature.value] = int(air_temperature)

    @property
    def specific_humidity(self):
        return bool(self.__iof_hydro[iof_hydro.specific_humidity.value])

    @specific_humidity.setter
    def specific_humidity(self, specific_humidity):
        assert isinstance(specific_humidity, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.specific_humidity.value] = int(specific_humidity)

    @property
    def solar_radiation(self):
        return bool(self.__iof_hydro[iof_hydro.solar_radiation.value])

    @solar_radiation.setter
    def solar_radiation(self, solar_radiation):
        assert isinstance(solar_radiation, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.solar_radiation.value] = int(solar_radiation)

    @property
    def sensible_flux(self):
        return bool(self.__iof_hydro[iof_hydro.sensible_flux.value])

    @sensible_flux.setter
    def sensible_flux(self, sensible_flux):
        assert isinstance(sensible_flux, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.sensible_flux.value] = int(sensible_flux)

    @property
    def latent_heat(self):
        return bool(self.__iof_hydro[iof_hydro.latent_heat.value])

    @latent_heat.setter
    def latent_heat(self, latent_heat):
        assert isinstance(latent_heat, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.latent_heat.value] = int(latent_heat)

    @property
    def upward_longwave(self):
        return bool(self.__iof_hydro[iof_hydro.upward_longwave.value])

    @upward_longwave.setter
    def upward_longwave(self, upward_longwave):
        assert isinstance(upward_longwave, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.upward_longwave.value] = int(upward_longwave)

    @property
    def downward_longwave(self):
        return bool(self.__iof_hydro[iof_hydro.downward_longwave.value])

    @downward_longwave.setter
    def downward_longwave(self, downward_longwave):
        assert isinstance(downward_longwave, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.downward_longwave.value] = int(downward_longwave)

    @property
    def total_heat_flux(self):
        return bool(self.__iof_hydro[iof_hydro.total_heat_flux.value])

    @total_heat_flux.setter
    def total_heat_flux(self, total_heat_flux):
        assert isinstance(total_heat_flux, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.total_heat_flux.value] = int(total_heat_flux)

    @property
    def evaporation(self):
        return bool(self.__iof_hydro[iof_hydro.evaporation.value])

    @evaporation.setter
    def evaporation(self, evaporation):
        assert isinstance(evaporation, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.evaporation.value] = int(evaporation)

    @property
    def precipitation(self):
        return bool(self.__iof_hydro[iof_hydro.precipitation.value])

    @precipitation.setter
    def precipitation(self, precipitation):
        assert isinstance(precipitation, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.precipitation.value] = int(precipitation)

    @property
    def bottom_stress(self):
        return bool(self.__iof_hydro[iof_hydro.bottom_stress.value])

    @bottom_stress.setter
    def bottom_stress(self, bottom_stress):
        assert isinstance(bottom_stress, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.bottom_stress.value] = int(bottom_stress)

    @property
    def wind_speed(self):
        return bool(self.__iof_hydro[iof_hydro.wind_speed.value])

    @wind_speed.setter
    def wind_speed(self, wind_speed):
        assert isinstance(wind_speed, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.wind_speed.value] = int(wind_speed)

    @property
    def wind_stress(self):
        return bool(self.__iof_hydro[iof_hydro.wind_stress.value])

    @wind_stress.setter
    def wind_stress(self, wind_stress):
        assert isinstance(wind_stress, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.wind_stress.value] = int(wind_stress)

    @property
    def dahv(self):
        return bool(self.__iof_hydro[iof_hydro.dahv.value])

    @dahv.setter
    def dahv(self, dahv):
        assert isinstance(dahv, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.dahv.value] = int(dahv)

    @property
    def vertical_velocity(self):
        return bool(self.__iof_hydro[iof_hydro.vertical_velocity.value])

    @vertical_velocity.setter
    def vertical_velocity(self, vertical_velocity):
        assert isinstance(vertical_velocity, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.vertical_velocity.value] = int(vertical_velocity)

    @property
    def temp(self):
        return bool(self.__iof_hydro[iof_hydro.temp.value])

    @temp.setter
    def temp(self, temp):
        assert isinstance(temp, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.temp.value] = int(temp)

    @property
    def salt(self):
        return bool(self.__iof_hydro[iof_hydro.salt.value])

    @salt.setter
    def salt(self, salt):
        assert isinstance(salt, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.salt.value] = int(salt)

    @property
    def water_density(self):
        return bool(self.__iof_hydro[iof_hydro.water_density.value])

    @water_density.setter
    def water_density(self, water_density):
        assert isinstance(water_density, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.water_density.value] = int(water_density)

    @property
    def diffusivity(self):
        return bool(self.__iof_hydro[iof_hydro.diffusivity.value])

    @diffusivity.setter
    def diffusivity(self, diffusivity):
        assert isinstance(diffusivity, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.diffusivity.value] = int(diffusivity)

    @property
    def viscosity(self):
        return bool(self.__iof_hydro[iof_hydro.viscosity.value])

    @viscosity.setter
    def viscosity(self, viscosity):
        assert isinstance(viscosity, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.viscosity.value] = int(viscosity)

    @property
    def TKE(self):
        return bool(self.__iof_hydro[iof_hydro.TKE.value])

    @TKE.setter
    def TKE(self, TKE):
        assert isinstance(TKE, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.TKE.value] = int(TKE)

    @property
    def mixing_length(self):
        return bool(self.__iof_hydro[iof_hydro.mixing_length.value])

    @mixing_length.setter
    def mixing_length(self, mixing_length):
        assert isinstance(mixing_length, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.mixing_length.value] = int(mixing_length)

    @property
    def hvel(self):
        return bool(self.__iof_hydro[iof_hydro.hvel.value])

    @hvel.setter
    def hvel(self, hvel):
        assert isinstance(hvel, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.hvel.value] = int(hvel)

    @property
    def hvel_side(self):
        return bool(self.__iof_hydro[iof_hydro.hvel_side.value])

    @hvel_side.setter
    def hvel_side(self, hvel_side):
        assert isinstance(hvel_side, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.hvel_side.value] = int(hvel_side)

    @property
    def wvel_elem(self):
        return bool(self.__iof_hydro[iof_hydro.wvel_elem.value])

    @wvel_elem.setter
    def wvel_elem(self, wvel_elem):
        assert isinstance(wvel_elem, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.wvel_elem.value] = int(wvel_elem)

    @property
    def temp_elem(self):
        return bool(self.__iof_hydro[iof_hydro.temp_elem.value])

    @temp_elem.setter
    def temp_elem(self, temp_elem):
        assert isinstance(temp_elem, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.temp_elem.value] = int(temp_elem)

    @property
    def salt_elem(self):
        return bool(self.__iof_hydro[iof_hydro.salt_elem.value])

    @salt_elem.setter
    def salt_elem(self, salt_elem):
        assert isinstance(salt_elem, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.salt_elem.value] = int(salt_elem)

    @property
    def pressure_gradient(self):
        return bool(self.__iof_hydro[iof_hydro.pressure_gradient.value])

    @pressure_gradient.setter
    def pressure_gradient(self, pressure_gradient):
        assert isinstance(pressure_gradient, bool), 'property must be a bool'
        self.__iof_hydro[iof_hydro.pressure_gradient.value] = int(pressure_gradient)

    @property
    def DVD_1(self):
        return bool(self.__iof_dvd[iof_dvd.DVD_1.value])

    @DVD_1.setter
    def DVD_1(self, DVD_1):
        assert isinstance(DVD_1, bool), 'property must be a bool'
        self.__iof_dvd[iof_dvd.DVD_1.value] = int(DVD_1)

    @property
    def WWM_1(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_1.value])

    @WWM_1.setter
    def WWM_1(self, WWM_1):
        assert isinstance(WWM_1, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_1.value] = int(WWM_1)

    @property
    def WWM_2(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_2.value])

    @WWM_2.setter
    def WWM_2(self, WWM_2):
        assert isinstance(WWM_2, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_2.value] = int(WWM_2)

    @property
    def WWM_3(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_3.value])

    @WWM_3.setter
    def WWM_3(self, WWM_3):
        assert isinstance(WWM_3, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_3.value] = int(WWM_3)

    @property
    def WWM_4(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_4.value])

    @WWM_4.setter
    def WWM_4(self, WWM_4):
        assert isinstance(WWM_4, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_4.value] = int(WWM_4)

    @property
    def WWM_5(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_5.value])

    @WWM_5.setter
    def WWM_5(self, WWM_5):
        assert isinstance(WWM_5, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_5.value] = int(WWM_5)

    @property
    def WWM_6(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_6.value])

    @WWM_6.setter
    def WWM_6(self, WWM_6):
        assert isinstance(WWM_6, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_6.value] = int(WWM_6)

    @property
    def WWM_9(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_9.value])

    @WWM_9.setter
    def WWM_9(self, WWM_9):
        assert isinstance(WWM_9, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_9.value] = int(WWM_9)

    @property
    def WWM_10(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_10.value])

    @WWM_10.setter
    def WWM_10(self, WWM_10):
        assert isinstance(WWM_10, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_10.value] = int(WWM_10)

    @property
    def WWM_11(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_11.value])

    @WWM_11.setter
    def WWM_11(self, WWM_11):
        assert isinstance(WWM_11, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_11.value] = int(WWM_11)

    @property
    def WWM_12(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_12.value])

    @WWM_12.setter
    def WWM_12(self, WWM_12):
        assert isinstance(WWM_12, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_12.value] = int(WWM_12)

    @property
    def WWM_13(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_13.value])

    @WWM_13.setter
    def WWM_13(self, WWM_13):
        assert isinstance(WWM_13, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_13.value] = int(WWM_13)

    @property
    def WWM_14(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_14.value])

    @WWM_14.setter
    def WWM_14(self, WWM_14):
        assert isinstance(WWM_14, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_14.value] = int(WWM_14)

    @property
    def WWM_15(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_15.value])

    @WWM_15.setter
    def WWM_15(self, WWM_15):
        assert isinstance(WWM_15, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_15.value] = int(WWM_15)

    @property
    def WWM_16(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_16.value])

    @WWM_16.setter
    def WWM_16(self, WWM_16):
        assert isinstance(WWM_16, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_16.value] = int(WWM_16)

    @property
    def WWM_17(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_17.value])

    @WWM_17.setter
    def WWM_17(self, WWM_17):
        assert isinstance(WWM_17, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_17.value] = int(WWM_17)

    @property
    def WWM_18(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_18.value])

    @WWM_18.setter
    def WWM_18(self, WWM_18):
        assert isinstance(WWM_18, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_18.value] = int(WWM_18)

    @property
    def WWM_19(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_19.value])

    @WWM_19.setter
    def WWM_19(self, WWM_19):
        assert isinstance(WWM_19, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_19.value] = int(WWM_19)

    @property
    def WWM_20(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_20.value])

    @WWM_20.setter
    def WWM_20(self, WWM_20):
        assert isinstance(WWM_20, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_20.value] = int(WWM_20)

    @property
    def WWM_21(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_21.value])

    @WWM_21.setter
    def WWM_21(self, WWM_21):
        assert isinstance(WWM_21, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_21.value] = int(WWM_21)

    @property
    def WWM_22(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_22.value])

    @WWM_22.setter
    def WWM_22(self, WWM_22):
        assert isinstance(WWM_22, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_22.value] = int(WWM_22)

    @property
    def WWM_23(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_23.value])

    @WWM_23.setter
    def WWM_23(self, WWM_23):
        assert isinstance(WWM_23, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_23.value] = int(WWM_23)

    @property
    def WWM_24(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_24.value])

    @WWM_24.setter
    def WWM_24(self, WWM_24):
        assert isinstance(WWM_24, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_24.value] = int(WWM_24)

    @property
    def WWM_25(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_25.value])

    @WWM_25.setter
    def WWM_25(self, WWM_25):
        assert isinstance(WWM_25, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_25.value] = int(WWM_25)

    @property
    def WWM_26(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_26.value])

    @WWM_26.setter
    def WWM_26(self, WWM_26):
        assert isinstance(WWM_26, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_26.value] = int(WWM_26)

    @property
    def WWM_27(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_27.value])

    @WWM_27.setter
    def WWM_27(self, WWM_27):
        assert isinstance(WWM_27, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_27.value] = int(WWM_27)

    @property
    def WWM_28(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_28.value])

    @WWM_28.setter
    def WWM_28(self, WWM_28):
        assert isinstance(WWM_28, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_28.value] = int(WWM_28)

    @property
    def WWM_energy_dir(self):
        return bool(self.__iof_wwm[iof_wwm.WWM_energy_dir.value])

    @WWM_energy_dir.setter
    def WWM_energy_dir(self, WWM_energy_dir):
        assert isinstance(WWM_energy_dir, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.WWM_energy_dir.value] = int(WWM_energy_dir)

    @property
    def wave_force(self):
        return bool(self.__iof_wwm[iof_wwm.wave_force.value])

    @wave_force.setter
    def wave_force(self, wave_force):
        assert isinstance(wave_force, bool), 'property must be a bool'
        self.__iof_wwm[iof_wwm.wave_force.value] = int(wave_force)

    @property
    def GEN_1(self):
        return bool(self.__iof_gen[iof_gen.GEN_1.value])

    @GEN_1.setter
    def GEN_1(self, GEN_1):
        assert isinstance(GEN_1, bool), 'property must be a bool'
        self.__iof_gen[iof_gen.GEN_1.value] = int(GEN_1)

    @property
    def GEN_2(self):
        return bool(self.__iof_gen[iof_gen.GEN_2.value])

    @GEN_2.setter
    def GEN_2(self, GEN_2):
        assert isinstance(GEN_2, bool), 'property must be a bool'
        self.__iof_gen[iof_gen.GEN_2.value] = int(GEN_2)

    @property
    def AGE_1(self):
        return bool(self.__iof_age[iof_age.AGE_1.value])

    @AGE_1.setter
    def AGE_1(self, AGE_1):
        assert isinstance(AGE_1, bool), 'property must be a bool'
        self.__iof_age[iof_age.AGE_1.value] = int(AGE_1)

    @property
    def AGE_2(self):
        return bool(self.__iof_age[iof_age.AGE_2.value])

    @AGE_2.setter
    def AGE_2(self, AGE_2):
        assert isinstance(AGE_2, bool), 'property must be a bool'
        self.__iof_age[iof_age.AGE_2.value] = int(AGE_2)

    @property
    def SED_depth_change(self):
        return bool(self.__iof_sed[iof_sed.SED_depth_change.value])

    @SED_depth_change.setter
    def SED_depth_change(self, SED_depth_change):
        assert isinstance(SED_depth_change, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.SED_depth_change.value] = int(SED_depth_change)

    @property
    def SED_D50(self):
        return bool(self.__iof_sed[iof_sed.SED_D50.value])

    @SED_D50.setter
    def SED_D50(self, SED_D50):
        assert isinstance(SED_D50, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.SED_D50.value] = int(SED_D50)

    @property
    def SED_bed_stress(self):
        return bool(self.__iof_sed[iof_sed.SED_bed_stress.value])

    @SED_bed_stress.setter
    def SED_bed_stress(self, SED_bed_stress):
        assert isinstance(SED_bed_stress, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.SED_bed_stress.value] = int(SED_bed_stress)

    @property
    def SED_bed_roughness(self):
        return bool(self.__iof_sed[iof_sed.SED_bed_roughness.value])

    @SED_bed_roughness.setter
    def SED_bed_roughness(self, SED_bed_roughness):
        assert isinstance(SED_bed_roughness, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.SED_bed_roughness.value] = int(SED_bed_roughness)

    @property
    def SED_TSC(self):
        return bool(self.__iof_sed[iof_sed.SED_TSC.value])

    @SED_TSC.setter
    def SED_TSC(self, SED_TSC):
        assert isinstance(SED_TSC, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.SED_TSC.value] = int(SED_TSC)

    @property
    def bed_thickness(self):
        return bool(self.__iof_sed[iof_sed.bed_thickness.value])

    @bed_thickness.setter
    def bed_thickness(self, bed_thickness):
        assert isinstance(bed_thickness, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.bed_thickness.value] = int(bed_thickness)

    @property
    def bed_age(self):
        return bool(self.__iof_sed[iof_sed.bed_age.value])

    @bed_age.setter
    def bed_age(self, bed_age):
        assert isinstance(bed_age, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.bed_age.value] = int(bed_age)

    @property
    def z0st(self):
        return bool(self.__iof_sed[iof_sed.z0st.value])

    @z0st.setter
    def z0st(self, z0st):
        assert isinstance(z0st, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.z0st.value] = int(z0st)

    @property
    def z0cr(self):
        return bool(self.__iof_sed[iof_sed.z0cr.value])

    @z0cr.setter
    def z0cr(self, z0cr):
        assert isinstance(z0cr, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.z0cr.value] = int(z0cr)

    @property
    def z0sw(self):
        return bool(self.__iof_sed[iof_sed.z0sw.value])

    @z0sw.setter
    def z0sw(self, z0sw):
        assert isinstance(z0sw, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.z0sw.value] = int(z0sw)

    @property
    def z0wr(self):
        return bool(self.__iof_sed[iof_sed.z0wr.value])

    @z0wr.setter
    def z0wr(self, z0wr):
        assert isinstance(z0wr, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.z0wr.value] = int(z0wr)

    @property
    def SED3D_1(self):
        return bool(self.__iof_sed[iof_sed.SED3D_1.value])

    @SED3D_1.setter
    def SED3D_1(self, SED3D_1):
        assert isinstance(SED3D_1, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.SED3D_1.value] = int(SED3D_1)

    @property
    def SED_bdld_1(self):
        return bool(self.__iof_sed[iof_sed.SED_bdld_1.value])

    @SED_bdld_1.setter
    def SED_bdld_1(self, SED_bdld_1):
        assert isinstance(SED_bdld_1, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.SED_bdld_1.value] = int(SED_bdld_1)

    @property
    def SED_bedfrac_1(self):
        return bool(self.__iof_sed[iof_sed.SED_bedfrac_1.value])

    @SED_bedfrac_1.setter
    def SED_bedfrac_1(self, SED_bedfrac_1):
        assert isinstance(SED_bedfrac_1, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.SED_bedfrac_1.value] = int(SED_bedfrac_1)

    @property
    def SED3D_2(self):
        return bool(self.__iof_sed[iof_sed.SED3D_2.value])

    @SED3D_2.setter
    def SED3D_2(self, SED3D_2):
        assert isinstance(SED3D_2, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.SED3D_2.value] = int(SED3D_2)

    @property
    def SED_bdld_2(self):
        return bool(self.__iof_sed[iof_sed.SED_bdld_2.value])

    @SED_bdld_2.setter
    def SED_bdld_2(self, SED_bdld_2):
        assert isinstance(SED_bdld_2, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.SED_bdld_2.value] = int(SED_bdld_2)

    @property
    def SED_bedfrac_3(self):
        return bool(self.__iof_sed[iof_sed.SED_bedfrac_3.value])

    @SED_bedfrac_3.setter
    def SED_bedfrac_3(self, SED_bedfrac_3):
        assert isinstance(SED_bedfrac_3, bool), 'property must be a bool'
        self.__iof_sed[iof_sed.SED_bedfrac_3.value] = int(SED_bedfrac_3)

    @property
    def ECO_1(self):
        return bool(self.__iof_eco[iof_eco.ECO_1.value])

    @ECO_1.setter
    def ECO_1(self, ECO_1):
        assert isinstance(ECO_1, bool), 'property must be a bool'
        self.__iof_eco[iof_eco.ECO_1.value] = int(ECO_1)

    @property
    def ICM_Chl(self):
        return bool(self.__iof_icm[iof_icm.ICM_ChlICM_Chl.value])

    @ICM_Chl.setter
    def ICM_Chl(self, ICM_Chl):
        assert isinstance(ICM_Chl, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_ChloICM_Chl.value] = int(ICM_Chl)

    @property
    def ICM_pH(self):
        return bool(self.__iof_icm[iof_icm.ICM_pHICM_pH.value])

    @ICM_pH.setter
    def ICM_pH(self, ICM_pH):
        assert isinstance(ICM_pH, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_pHioICM_pH.value] = int(ICM_pH)

    @property
    def ICM_PrmPrdt(self):
        return bool(self.__iof_icm[iof_icm.ICM_PrmPrdt.value])

    @ICM_PrmPrdt.setter
    def ICM_PrmPrdt(self, ICM_PrmPrdt):
        assert isinstance(ICM_PrmPrdt, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_PrmPrdt.value] = int(ICM_PrmPrdt)

    @property
    def ICM_DIN(self):
        return bool(self.__iof_icm[iof_icm.ICM_DINICM_DIN.value])

    @ICM_DIN.setter
    def ICM_DIN(self, ICM_DIN):
        assert isinstance(ICM_DIN, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_DIN.value] = int(ICM_DIN)

    @property
    def ICM_PON(self):
        return bool(self.__iof_icm[iof_icm.ICM_PONICM_PON.value])

    @ICM_PON.setter
    def ICM_PON(self, ICM_PON):
        assert isinstance(ICM_PON, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_PON.value] = int(ICM_PON)

    @property
    def ICM_SED_BENDOC(self):
        return bool(self.__iof_icm[iof_icm.ICM_SED_BENDOCICM_SED_BENDOC.value])

    @ICM_SED_BENDOC.setter
    def ICM_SED_BENDOC(self, ICM_SED_BENDOC):
        assert isinstance(ICM_SED_BENDOC, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SED_ICM_SED_BENDOCBENDOC.value] = int(ICM_SED_BENDOC)

    @property
    def ICM_SED_BENNH4(self):
        return bool(self.__iof_icm[iof_icm.ICM_SED_BENNH4ICM_SED_BENNH4.value])

    @ICM_SED_BENNH4.setter
    def ICM_SED_BENNH4(self, ICM_SED_BENNH4):
        assert isinstance(ICM_SED_BENNH4, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SED_ICM_SED_BENNH4BENNH4.value] = int(ICM_SED_BENNH4)

    @property
    def ICM_SED_BENNO3(self):
        return bool(self.__iof_icm[iof_icm.ICM_SED_BENNO3ICM_SED_BENNO3.value])

    @ICM_SED_BENNO3.setter
    def ICM_SED_BENNO3(self, ICM_SED_BENNO3):
        assert isinstance(ICM_SED_BENNO3, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SED_ICM_SED_BENNO3BENNO3.value] = int(ICM_SED_BENNO3)

    @property
    def ICM_SED_BENPO4(self):
        return bool(self.__iof_icm[iof_icm.ICM_SED_BENPO4ICM_SED_BENPO4.value])

    @ICM_SED_BENPO4.setter
    def ICM_SED_BENPO4(self, ICM_SED_BENPO4):
        assert isinstance(ICM_SED_BENPO4, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SED_ICM_SED_BENPO4BENPO4.value] = int(ICM_SED_BENPO4)

    @property
    def ICM_SED_BENCOD(self):
        return bool(self.__iof_icm[iof_icm.ICM_SED_BENCODICM_SED_BENCOD.value])

    @ICM_SED_BENCOD.setter
    def ICM_SED_BENCOD(self, ICM_SED_BENCOD):
        assert isinstance(ICM_SED_BENCOD, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SED_ICM_SED_BENCODBENCOD.value] = int(ICM_SED_BENCOD)

    @property
    def ICM_SED_BENDO(self):
        return bool(self.__iof_icm[iof_icm.ICM_SED_BENDO.value])

    @ICM_SED_BENDO.setter
    def ICM_SED_BENDO(self, ICM_SED_BENDO):
        assert isinstance(ICM_SED_BENDO, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SED_BENDO.value] = int(ICM_SED_BENDO)

    @property
    def ICM_SED_BENSA(self):
        return bool(self.__iof_icm[iof_icm.ICM_SED_BENSA.value])

    @ICM_SED_BENSA.setter
    def ICM_SED_BENSA(self, ICM_SED_BENSA):
        assert isinstance(ICM_SED_BENSA, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SED_BENSA.value] = int(ICM_SED_BENSA)

    @property
    def ICM_lfsav(self):
        return bool(self.__iof_icm[iof_icm.ICM_lfsav.value])

    @ICM_lfsav.setter
    def ICM_lfsav(self, ICM_lfsav):
        assert isinstance(ICM_lfsav, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_lfsav.value] = int(ICM_lfsav)

    @property
    def ICM_stsav(self):
        return bool(self.__iof_icm[iof_icm.ICM_stsav.value])

    @ICM_stsav.setter
    def ICM_stsav(self, ICM_stsav):
        assert isinstance(ICM_stsav, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_stsav.value] = int(ICM_stsav)

    @property
    def ICM_rtsav(self):
        return bool(self.__iof_icm[iof_icm.ICM_rtsav.value])

    @ICM_rtsav.setter
    def ICM_rtsav(self, ICM_rtsav):
        assert isinstance(ICM_rtsav, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_rtsav.value] = int(ICM_rtsav)

    @property
    def ICM_tlfsav(self):
        return bool(self.__iof_icm[iof_icm.ICM_tlfsav.value])

    @ICM_tlfsav.setter
    def ICM_tlfsav(self, ICM_tlfsav):
        assert isinstance(ICM_tlfsav, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_tlfsav.value] = int(ICM_tlfsav)

    @property
    def ICM_tstsav(self):
        return bool(self.__iof_icm[iof_icm.ICM_tstsav.value])

    @ICM_tstsav.setter
    def ICM_tstsav(self, ICM_tstsav):
        assert isinstance(ICM_tstsav, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_tstsav.value] = int(ICM_tstsav)

    @property
    def ICM_trtsav(self):
        return bool(self.__iof_icm[iof_icm.ICM_trtsav.value])

    @ICM_trtsav.setter
    def ICM_trtsav(self, ICM_trtsav):
        assert isinstance(ICM_trtsav, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_trtsav.value] = int(ICM_trtsav)

    @property
    def ICM_hcansav(self):
        return bool(self.__iof_icm[iof_icm.ICM_hcansav.value])

    @ICM_hcansav.setter
    def ICM_hcansav(self, ICM_hcansav):
        assert isinstance(ICM_hcansav, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_hcansav.value] = int(ICM_hcansav)

    @property
    def ICM_CNH4(self):
        return bool(self.__iof_icm[iof_icm.ICM_CNH4.value])

    @ICM_CNH4.setter
    def ICM_CNH4(self, ICM_CNH4):
        assert isinstance(ICM_CNH4, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_CNH4.value] = int(ICM_CNH4)

    @property
    def ICM_CNH3(self):
        return bool(self.__iof_icm[iof_icm.ICM_CNH3.value])

    @ICM_CNH3.setter
    def ICM_CNH3(self, ICM_CNH3):
        assert isinstance(ICM_CNH3, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_CNH3.value] = int(ICM_CNH3)

    @property
    def ICM_CPIP(self):
        return bool(self.__iof_icm[iof_icm.ICM_CPIP.value])

    @ICM_CPIP.setter
    def ICM_CPIP(self, ICM_CPIP):
        assert isinstance(ICM_CPIP, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_CPIP.value] = int(ICM_CPIP)

    @property
    def ICM_CPOS(self):
        return bool(self.__iof_icm[iof_icm.ICM_CPOS.value])

    @ICM_CPOS.setter
    def ICM_CPOS(self, ICM_CPOS):
        assert isinstance(ICM_CPOS, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_CPOS.value] = int(ICM_CPOS)

    @property
    def ICM_CCH4(self):
        return bool(self.__iof_icm[iof_icm.ICM_CCH4.value])

    @ICM_CCH4.setter
    def ICM_CCH4(self, ICM_CCH4):
        assert isinstance(ICM_CCH4, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_CCH4.value] = int(ICM_CCH4)

    @property
    def ICM_CSO4(self):
        return bool(self.__iof_icm[iof_icm.ICM_CSO4.value])

    @ICM_CSO4.setter
    def ICM_CSO4(self, ICM_CSO4):
        assert isinstance(ICM_CSO4, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_CSO4.value] = int(ICM_CSO4)

    @property
    def ICM_CH2S(self):
        return bool(self.__iof_icm[iof_icm.ICM_CH2S.value])

    @ICM_CH2S.setter
    def ICM_CH2S(self, ICM_CH2S):
        assert isinstance(ICM_CH2S, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_CH2S.value] = int(ICM_CH2S)

    @property
    def ICM_SEDPON1(self):
        return bool(self.__iof_icm[iof_icm.ICM_SEDPON1.value])

    @ICM_SEDPON1.setter
    def ICM_SEDPON1(self, ICM_SEDPON1):
        assert isinstance(ICM_SEDPON1, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SEDPON1.value] = int(ICM_SEDPON1)

    @property
    def ICM_SEDPON2(self):
        return bool(self.__iof_icm[iof_icm.ICM_SEDPON2.value])

    @ICM_SEDPON2.setter
    def ICM_SEDPON2(self, ICM_SEDPON2):
        assert isinstance(ICM_SEDPON2, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SEDPON2.value] = int(ICM_SEDPON2)

    @property
    def ICM_SEDPON3(self):
        return bool(self.__iof_icm[iof_icm.ICM_SEDPON3.value])

    @ICM_SEDPON3.setter
    def ICM_SEDPON3(self, ICM_SEDPON3):
        assert isinstance(ICM_SEDPON3, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SEDPON3.value] = int(ICM_SEDPON3)

    @property
    def ICM_SEDPOP1(self):
        return bool(self.__iof_icm[iof_icm.ICM_SEDPOP1.value])

    @ICM_SEDPOP1.setter
    def ICM_SEDPOP1(self, ICM_SEDPOP1):
        assert isinstance(ICM_SEDPOP1, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SEDPOP1.value] = int(ICM_SEDPOP1)

    @property
    def ICM_SEDPOP2(self):
        return bool(self.__iof_icm[iof_icm.ICM_SEDPOP2.value])

    @ICM_SEDPOP2.setter
    def ICM_SEDPOP2(self, ICM_SEDPOP2):
        assert isinstance(ICM_SEDPOP2, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SEDPOP2.value] = int(ICM_SEDPOP2)

    @property
    def ICM_SEDPOP3(self):
        return bool(self.__iof_icm[iof_icm.ICM_SEDPOP3.value])

    @ICM_SEDPOP3.setter
    def ICM_SEDPOP3(self, ICM_SEDPOP3):
        assert isinstance(ICM_SEDPOP3, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SEDPOP3.value] = int(ICM_SEDPOP3)

    @property
    def ICM_SEDPOC1(self):
        return bool(self.__iof_icm[iof_icm.ICM_SEDPOC1.value])

    @ICM_SEDPOC1.setter
    def ICM_SEDPOC1(self, ICM_SEDPOC1):
        assert isinstance(ICM_SEDPOC1, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SEDPOC1.value] = int(ICM_SEDPOC1)

    @property
    def ICM_SEDPOC2(self):
        return bool(self.__iof_icm[iof_icm.ICM_SEDPOC2.value])

    @ICM_SEDPOC2.setter
    def ICM_SEDPOC2(self, ICM_SEDPOC2):
        assert isinstance(ICM_SEDPOC2, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SEDPOC2.value] = int(ICM_SEDPOC2)

    @property
    def ICM_SEDPOC3(self):
        return bool(self.__iof_icm[iof_icm.ICM_SEDPOC3.value])

    @ICM_SEDPOC3.setter
    def ICM_SEDPOC3(self, ICM_SEDPOC3):
        assert isinstance(ICM_SEDPOC3, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_SEDPOC3.value] = int(ICM_SEDPOC3)

    @property
    def ICM_EROH2S(self):
        return bool(self.__iof_icm[iof_icm.ICM_EROH2S.value])

    @ICM_EROH2S.setter
    def ICM_EROH2S(self, ICM_EROH2S):
        assert isinstance(ICM_EROH2S, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_EROH2S.value] = int(ICM_EROH2S)

    @property
    def ICM_EROLPOC(self):
        return bool(self.__iof_icm[iof_icm.ICM_EROLPOC.value])

    @ICM_EROLPOC.setter
    def ICM_EROLPOC(self, ICM_EROLPOC):
        assert isinstance(ICM_EROLPOC, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_EROLPOC.value] = int(ICM_EROLPOC)

    @property
    def ICM_ERORPOC(self):
        return bool(self.__iof_icm[iof_icm.ICM_ERORPOC.value])

    @ICM_ERORPOC.setter
    def ICM_ERORPOC(self, ICM_ERORPOC):
        assert isinstance(ICM_ERORPOC, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_ERORPOC.value] = int(ICM_ERORPOC)

    @property
    def ICM_DO_consumption(self):
        return bool(self.__iof_icm[iof_icm.ICM_DO_consumption.value])

    @ICM_DO_consumption.setter
    def ICM_DO_consumption(self, ICM_DO_consumption):
        assert isinstance(ICM_DO_consumption, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_DO_consumption.value] = int(ICM_DO_consumption)

    @property
    def ICM_GP1(self):
        return bool(self.__iof_icm[iof_icm.ICM_GP1.value])

    @ICM_GP1.setter
    def ICM_GP1(self, ICM_GP1):
        assert isinstance(ICM_GP1, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_GP1.value] = int(ICM_GP1)

    @property
    def ICM_GP2(self):
        return bool(self.__iof_icm[iof_icm.ICM_GP2.value])

    @ICM_GP2.setter
    def ICM_GP2(self, ICM_GP2):
        assert isinstance(ICM_GP2, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_GP2.value] = int(ICM_GP2)

    @property
    def ICM_GP3(self):
        return bool(self.__iof_icm[iof_icm.ICM_GP3.value])

    @ICM_GP3.setter
    def ICM_GP3(self, ICM_GP3):
        assert isinstance(ICM_GP3, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_GP3.value] = int(ICM_GP3)

    @property
    def ICM_1(self):
        return bool(self.__iof_icm[iof_icm.ICM_1.value])

    @ICM_1.setter
    def ICM_1(self, ICM_1):
        assert isinstance(ICM_1, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_1.value] = int(ICM_1)

    @property
    def ICM_2(self):
        return bool(self.__iof_icm[iof_icm.ICM_2.value])

    @ICM_2.setter
    def ICM_2(self, ICM_2):
        assert isinstance(ICM_2, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_2.value] = int(ICM_2)

    @property
    def ICM_3(self):
        return bool(self.__iof_icm[iof_icm.ICM_3.value])

    @ICM_3.setter
    def ICM_3(self, ICM_3):
        assert isinstance(ICM_3, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_3.value] = int(ICM_3)

    @property
    def ICM_4(self):
        return bool(self.__iof_icm[iof_icm.ICM_4.value])

    @ICM_4.setter
    def ICM_4(self, ICM_4):
        assert isinstance(ICM_4, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_4.value] = int(ICM_4)

    @property
    def ICM_5(self):
        return bool(self.__iof_icm[iof_icm.ICM_5.value])

    @ICM_5.setter
    def ICM_5(self, ICM_5):
        assert isinstance(ICM_5, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_5.value] = int(ICM_5)

    @property
    def ICM_6(self):
        return bool(self.__iof_icm[iof_icm.ICM_6.value])

    @ICM_6.setter
    def ICM_6(self, ICM_6):
        assert isinstance(ICM_6, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_6.value] = int(ICM_6)

    @property
    def ICM_7(self):
        return bool(self.__iof_icm[iof_icm.ICM_7.value])

    @ICM_7.setter
    def ICM_7(self, ICM_7):
        assert isinstance(ICM_7, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_7.value] = int(ICM_7)

    @property
    def ICM_8(self):
        return bool(self.__iof_icm[iof_icm.ICM_8.value])

    @ICM_8.setter
    def ICM_8(self, ICM_8):
        assert isinstance(ICM_8, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_8.value] = int(ICM_8)

    @property
    def ICM_9(self):
        return bool(self.__iof_icm[iof_icm.ICM_9.value])

    @ICM_9.setter
    def ICM_9(self, ICM_9):
        assert isinstance(ICM_9, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_9.value] = int(ICM_9)

    @property
    def ICM_10(self):
        return bool(self.__iof_icm[iof_icm.ICM_10.value])

    @ICM_10.setter
    def ICM_10(self, ICM_10):
        assert isinstance(ICM_10, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_10.value] = int(ICM_10)

    @property
    def ICM_11(self):
        return bool(self.__iof_icm[iof_icm.ICM_11.value])

    @ICM_11.setter
    def ICM_11(self, ICM_11):
        assert isinstance(ICM_11, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_11.value] = int(ICM_11)

    @property
    def ICM_12(self):
        return bool(self.__iof_icm[iof_icm.ICM_12.value])

    @ICM_12.setter
    def ICM_12(self, ICM_12):
        assert isinstance(ICM_12, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_12.value] = int(ICM_12)

    @property
    def ICM_13(self):
        return bool(self.__iof_icm[iof_icm.ICM_13.value])

    @ICM_13.setter
    def ICM_13(self, ICM_13):
        assert isinstance(ICM_13, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_13.value] = int(ICM_13)

    @property
    def ICM_14(self):
        return bool(self.__iof_icm[iof_icm.ICM_14.value])

    @ICM_14.setter
    def ICM_14(self, ICM_14):
        assert isinstance(ICM_14, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_14.value] = int(ICM_14)

    @property
    def ICM_15(self):
        return bool(self.__iof_icm[iof_icm.ICM_15.value])

    @ICM_15.setter
    def ICM_15(self, ICM_15):
        assert isinstance(ICM_15, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_15.value] = int(ICM_15)

    @property
    def ICM_16(self):
        return bool(self.__iof_icm[iof_icm.ICM_16.value])

    @ICM_16.setter
    def ICM_16(self, ICM_16):
        assert isinstance(ICM_16, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_16.value] = int(ICM_16)

    @property
    def ICM_17(self):
        return bool(self.__iof_icm[iof_icm.ICM_17.value])

    @ICM_17.setter
    def ICM_17(self, ICM_17):
        assert isinstance(ICM_17, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_17.value] = int(ICM_17)

    @property
    def ICM_18(self):
        return bool(self.__iof_icm[iof_icm.ICM_18.value])

    @ICM_18.setter
    def ICM_18(self, ICM_18):
        assert isinstance(ICM_18, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_18.value] = int(ICM_18)

    @property
    def ICM_19(self):
        return bool(self.__iof_icm[iof_icm.ICM_19.value])

    @ICM_19.setter
    def ICM_19(self, ICM_19):
        assert isinstance(ICM_19, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_19.value] = int(ICM_19)

    @property
    def ICM_20(self):
        return bool(self.__iof_icm[iof_icm.ICM_20.value])

    @ICM_20.setter
    def ICM_20(self, ICM_20):
        assert isinstance(ICM_20, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_20.value] = int(ICM_20)

    @property
    def ICM_21(self):
        return bool(self.__iof_icm[iof_icm.ICM_21.value])

    @ICM_21.setter
    def ICM_21(self, ICM_21):
        assert isinstance(ICM_21, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_21.value] = int(ICM_21)

    @property
    def ICM_22(self):
        return bool(self.__iof_icm[iof_icm.ICM_22.value])

    @ICM_22.setter
    def ICM_22(self, ICM_22):
        assert isinstance(ICM_22, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_22.value] = int(ICM_22)

    @property
    def ICM_23(self):
        return bool(self.__iof_icm[iof_icm.ICM_23.value])

    @ICM_23.setter
    def ICM_23(self, ICM_23):
        assert isinstance(ICM_23, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_23.value] = int(ICM_23)

    @property
    def ICM_24(self):
        return bool(self.__iof_icm[iof_icm.ICM_24.value])

    @ICM_24.setter
    def ICM_24(self, ICM_24):
        assert isinstance(ICM_24, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_24.value] = int(ICM_24)

    @property
    def ICM_25(self):
        return bool(self.__iof_icm[iof_icm.ICM_25.value])

    @ICM_25.setter
    def ICM_25(self, ICM_25):
        assert isinstance(ICM_25, bool), 'property must be a bool'
        self.__iof_icm[iof_icm.ICM_25.value] = int(ICM_25)

    @property
    def COS_1(self):
        return bool(self.__iof_cos[iof_cos.COS_1.value])

    @COS_1.setter
    def COS_1(self, COS_1):
        assert isinstance(COS_1, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_1.value] = int(COS_1)

    @property
    def COS_2(self):
        return bool(self.__iof_cos[iof_cos.COS_2.value])

    @COS_2.setter
    def COS_2(self, COS_2):
        assert isinstance(COS_2, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_2.value] = int(COS_2)

    @property
    def COS_3(self):
        return bool(self.__iof_cos[iof_cos.COS_3.value])

    @COS_3.setter
    def COS_3(self, COS_3):
        assert isinstance(COS_3, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_3.value] = int(COS_3)

    @property
    def COS_4(self):
        return bool(self.__iof_cos[iof_cos.COS_4.value])

    @COS_4.setter
    def COS_4(self, COS_4):
        assert isinstance(COS_4, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_4.value] = int(COS_4)

    @property
    def COS_5(self):
        return bool(self.__iof_cos[iof_cos.COS_5.value])

    @COS_5.setter
    def COS_5(self, COS_5):
        assert isinstance(COS_5, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_5.value] = int(COS_5)

    @property
    def COS_6(self):
        return bool(self.__iof_cos[iof_cos.COS_6.value])

    @COS_6.setter
    def COS_6(self, COS_6):
        assert isinstance(COS_6, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_6.value] = int(COS_6)

    @property
    def COS_7(self):
        return bool(self.__iof_cos[iof_cos.COS_7.value])

    @COS_7.setter
    def COS_7(self, COS_7):
        assert isinstance(COS_7, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_7.value] = int(COS_7)

    @property
    def COS_8(self):
        return bool(self.__iof_cos[iof_cos.COS_8.value])

    @COS_8.setter
    def COS_8(self, COS_8):
        assert isinstance(COS_8, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_8.value] = int(COS_8)

    @property
    def COS_9(self):
        return bool(self.__iof_cos[iof_cos.COS_9.value])

    @COS_9.setter
    def COS_9(self, COS_9):
        assert isinstance(COS_9, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_9.value] = int(COS_9)

    @property
    def COS_10(self):
        return bool(self.__iof_cos[iof_cos.COS_10.value])

    @COS_10.setter
    def COS_10(self, COS_10):
        assert isinstance(COS_10, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_10.value] = int(COS_10)

    @property
    def COS_11(self):
        return bool(self.__iof_cos[iof_cos.COS_11.value])

    @COS_11.setter
    def COS_11(self, COS_11):
        assert isinstance(COS_11, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_11.value] = int(COS_11)

    @property
    def COS_12(self):
        return bool(self.__iof_cos[iof_cos.COS_12.value])

    @COS_12.setter
    def COS_12(self, COS_12):
        assert isinstance(COS_12, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_12.value] = int(COS_12)

    @property
    def COS_13(self):
        return bool(self.__iof_cos[iof_cos.COS_13.value])

    @COS_13.setter
    def COS_13(self, COS_13):
        assert isinstance(COS_13, bool), 'property must be a bool'
        self.__iof_cos[iof_cos.COS_13.value] = int(COS_13)

    @property
    def FIB_1(self):
        return bool(self.__iof_fib[iof_fib.FIB_1.value])

    @FIB_1.setter
    def FIB_1(self, FIB_1):
        assert isinstance(FIB_1, bool), 'property must be a bool'
        self.__iof_fib[iof_fib.FIB_1.value] = int(FIB_1)

    @property
    def SED2D_depth_change(self):
        return bool(self.__iof_sed2d[iof_sed2d.SED2D_depth_change.value])

    @SED2D_depth_change.setter
    def SED2D_depth_change(self, SED2D_depth_change):
        assert isinstance(SED2D_depth_change, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.SED2D_depth_change.value] = int(SED2D_depth_change)

    @property
    def SED2D_cflsed(self):
        return bool(self.__iof_sed2d[iof_sed2d.SED2D_cflsed.value])

    @SED2D_cflsed.setter
    def SED2D_cflsed(self, SED2D_cflsed):
        assert isinstance(SED2D_cflsed, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.SED2D_cflsed.value] = int(SED2D_cflsed)

    @property
    def SED2D_d50(self):
        return bool(self.__iof_sed2d[iof_sed2d.SED2D_d50.value])

    @SED2D_d50.setter
    def SED2D_d50(self, SED2D_d50):
        assert isinstance(SED2D_d50, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.SED2D_d50.value] = int(SED2D_d50)

    @property
    def SED2D_total_transport(self):
        return bool(self.__iof_sed2d[iof_sed2d.SED2D_total_transport.value])

    @SED2D_total_transport.setter
    def SED2D_total_transport(self, SED2D_total_transport):
        assert isinstance(SED2D_total_transport, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.SED2D_total_transport.value] = int(SED2D_total_transport)

    @property
    def SED2D_susp_load(self):
        return bool(self.__iof_sed2d[iof_sed2d.SED2D_susp_load.value])

    @SED2D_susp_load.setter
    def SED2D_susp_load(self, SED2D_susp_load):
        assert isinstance(SED2D_susp_load, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.SED2D_susp_load.value] = int(SED2D_susp_load)

    @property
    def SED2D_bed_load(self):
        return bool(self.__iof_sed2d[iof_sed2d.SED2D_bed_load.value])

    @SED2D_bed_load.setter
    def SED2D_bed_load(self, SED2D_bed_load):
        assert isinstance(SED2D_bed_load, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.SED2D_bed_load.value] = int(SED2D_bed_load)

    @property
    def SED2D_average_transport(self):
        return bool(self.__iof_sed2d[iof_sed2d.SED2D_average_transport.value])

    @SED2D_average_transport.setter
    def SED2D_average_transport(self, SED2D_average_transport):
        assert isinstance(SED2D_average_transport, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.SED2D_average_transport.value] = int(SED2D_average_transport)

    @property
    def SED2D_bottom_slope(self):
        return bool(self.__iof_sed2d[iof_sed2d.SED2D_bottom_slope.value])

    @SED2D_bottom_slope.setter
    def SED2D_bottom_slope(self, SED2D_bottom_slope):
        assert isinstance(SED2D_bottom_slope, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.SED2D_bottom_slope.value] = int(SED2D_bottom_slope)

    @property
    def z0eq(self):
        return bool(self.__iof_sed2d[iof_sed2d.z0eq.value])

    @z0eq.setter
    def z0eq(self, z0eq):
        assert isinstance(z0eq, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.z0eq.value] = int(z0eq)

    @property
    def z0cr2d(self):
        return bool(self.__iof_sed2d[iof_sed2d.z0cr2d.value])

    @z0cr2d.setter
    def z0cr2d(self, z0cr):
        assert isinstance(z0cr, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.z0cr2d.value] = int(z0cr)

    @property
    def z0sw2d(self):
        return bool(self.__iof_sed2d[iof_sed2d.z0sw2d.value])

    @z0sw2d.setter
    def z0sw2d(self, z0sw):
        assert isinstance(z0sw, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.z0sw2d.value] = int(z0sw)

    @property
    def z0wr2d(self):
        return bool(self.__iof_sed2d[iof_sed2d.z0wr2d.value])

    @z0wr2d.setter
    def z0wr2d(self, z0wr):
        assert isinstance(z0wr, bool), 'property must be a bool'
        self.__iof_sed2d[iof_sed2d.z0wr2d.value] = int(z0wr)

    @property
    def marsh_flag(self):
        return bool(self.__iof_marsh[iof_marsh.marsh_flag.value])

    @marsh_flag.setter
    def marsh_flag(self, marsh_flag):
        assert isinstance(marsh_flag, bool), 'property must be a bool'
        self.__iof_marsh[iof_marsh.marsh_flag.value] = int(marsh_flag)

    @property
    def ICE_velocity(self):
        return bool(self.__iof_ice[iof_ice.ICE_velocity.value])

    @ICE_velocity.setter
    def ICE_velocity(self, ICE_velocity):
        assert isinstance(ICE_velocity, bool), 'property must be a bool'
        self.__iof_ice[iof_ice.ICE_velocity.value] = int(ICE_velocity)

    @property
    def ICE_strain_rate(self):
        return bool(self.__iof_ice[iof_ice.ICE_strain_rate.value])

    @ICE_strain_rate.setter
    def ICE_strain_rate(self, ICE_strain_rate):
        assert isinstance(ICE_strain_rate, bool), 'property must be a bool'
        self.__iof_ice[iof_ice.ICE_strain_rate.value] = int(ICE_strain_rate)

    @property
    def ICE_net_heat_flux(self):
        return bool(self.__iof_ice[iof_ice.ICE_net_heat_flux.value])

    @ICE_net_heat_flux.setter
    def ICE_net_heat_flux(self, ICE_net_heat_flux):
        assert isinstance(ICE_net_heat_flux, bool), 'property must be a bool'
        self.__iof_ice[iof_ice.ICE_net_heat_flux.value] = int(ICE_net_heat_flux)

    @property
    def ICE_fresh_water_flux(self):
        return bool(self.__iof_ice[iof_ice.ICE_fresh_water_flux.value])

    @ICE_fresh_water_flux.setter
    def ICE_fresh_water_flux(self, ICE_fresh_water_flux):
        assert isinstance(ICE_fresh_water_flux, bool), 'property must be a bool'
        self.__iof_ice[iof_ice.ICE_fresh_water_flux.value] = int(ICE_fresh_water_flux)

    @property
    def ICE_top_T(self):
        return bool(self.__iof_ice[iof_ice.ICE_top_T.value])

    @ICE_top_T.setter
    def ICE_top_T(self, ICE_top_T):
        assert isinstance(ICE_top_T, bool), 'property must be a bool'
        self.__iof_ice[iof_ice.ICE_top_T.value] = int(ICE_top_T)

    @property
    def ICE_tracer_1(self):
        return bool(self.__iof_ice[iof_ice.ICE_tracer_1.value])

    @ICE_tracer_1.setter
    def ICE_tracer_1(self, ICE_tracer_1):
        assert isinstance(ICE_tracer_1, bool), 'property must be a bool'
        self.__iof_ice[iof_ice.ICE_tracer_1.value] = int(ICE_tracer_1)

    @property
    def ICE_tracer_2(self):
        return bool(self.__iof_ice[iof_ice.ICE_tracer_2.value])

    @ICE_tracer_2.setter
    def ICE_tracer_2(self, ICE_tracer_2):
        assert isinstance(ICE_tracer_2, bool), 'property must be a bool'
        self.__iof_ice[iof_ice.ICE_tracer_2.value] = int(ICE_tracer_2)

    @property
    def ICE_tracer_3(self):
        return bool(self.__iof_ice[iof_ice.ICE_tracer_3.value])

    @ICE_tracer_3.setter
    def ICE_tracer_3(self, ICE_tracer_3):
        assert isinstance(ICE_tracer_3, bool), 'property must be a bool'
        self.__iof_ice[iof_ice.ICE_tracer_3.value] = int(ICE_tracer_3)

    @property
    def ANA_air_pres_grad_x(self):
        return bool(self.__iof_ana[iof_ana.ANA_air_pres_grad_x.value])

    @ANA_air_pres_grad_x.setter
    def ANA_air_pres_grad_x(self, ANA_air_pres_grad_x):
        assert isinstance(ANA_air_pres_grad_x, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_air_pres_grad_x.value] = int(ANA_air_pres_grad_x)

    @property
    def ANA_air_pres_grad_y(self):
        return bool(self.__iof_ana[iof_ana.ANA_air_pres_grad_y.value])

    @ANA_air_pres_grad_y.setter
    def ANA_air_pres_grad_y(self, ANA_air_pres_grad_y):
        assert isinstance(ANA_air_pres_grad_y, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_air_pres_grad_y.value] = int(ANA_air_pres_grad_y)

    @property
    def ANA_tide_pot_grad_x(self):
        return bool(self.__iof_ana[iof_ana.ANA_tide_pot_grad_x.value])

    @ANA_tide_pot_grad_x.setter
    def ANA_tide_pot_grad_x(self, ANA_tide_pot_grad_x):
        assert isinstance(ANA_tide_pot_grad_x, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_tide_pot_grad_x.value] = int(ANA_tide_pot_grad_x)

    @property
    def ANA_tide_pot_grad_y(self):
        return bool(self.__iof_ana[iof_ana.ANA_tide_pot_grad_y.value])

    @ANA_tide_pot_grad_y.setter
    def ANA_tide_pot_grad_y(self, ANA_tide_pot_grad_y):
        assert isinstance(ANA_tide_pot_grad_y, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_tide_pot_grad_y.value] = int(ANA_tide_pot_grad_y)

    @property
    def ANA_hor_viscosity_x(self):
        return bool(self.__iof_ana[iof_ana.ANA_hor_viscosity_x.value])

    @ANA_hor_viscosity_x.setter
    def ANA_hor_viscosity_x(self, ANA_hor_viscosity_x):
        assert isinstance(ANA_hor_viscosity_x, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_hor_viscosity_x.value] = int(ANA_hor_viscosity_x)

    @property
    def ANA_hor_viscosity_y(self):
        return bool(self.__iof_ana[iof_ana.ANA_hor_viscosity_y.value])

    @ANA_hor_viscosity_y.setter
    def ANA_hor_viscosity_y(self, ANA_hor_viscosity_y):
        assert isinstance(ANA_hor_viscosity_y, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_hor_viscosity_y.value] = int(ANA_hor_viscosity_y)

    @property
    def ANA_bclinic_force_x(self):
        return bool(self.__iof_ana[iof_ana.ANA_bclinic_force_x.value])

    @ANA_bclinic_force_x.setter
    def ANA_bclinic_force_x(self, ANA_bclinic_force_x):
        assert isinstance(ANA_bclinic_force_x, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_bclinic_force_x.value] = int(ANA_bclinic_force_x)

    @property
    def ANA_bclinic_force_y(self):
        return bool(self.__iof_ana[iof_ana.ANA_bclinic_force_y.value])

    @ANA_bclinic_force_y.setter
    def ANA_bclinic_force_y(self, ANA_bclinic_force_y):
        assert isinstance(ANA_bclinic_force_y, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_bclinic_force_y.value] = int(ANA_bclinic_force_y)

    @property
    def ANA_vert_viscosity_x(self):
        return bool(self.__iof_ana[iof_ana.ANA_vert_viscosity_x.value])

    @ANA_vert_viscosity_x.setter
    def ANA_vert_viscosity_x(self, ANA_vert_viscosity_x):
        assert isinstance(ANA_vert_viscosity_x, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_vert_viscosity_x.value] = int(ANA_vert_viscosity_x)

    @property
    def ANA_vert_viscosity_y(self):
        return bool(self.__iof_ana[iof_ana.ANA_vert_viscosity_y.value])

    @ANA_vert_viscosity_y.setter
    def ANA_vert_viscosity_y(self, ANA_vert_viscosity_y):
        assert isinstance(ANA_vert_viscosity_y, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_vert_viscosity_y.value] = int(ANA_vert_viscosity_y)

    @property
    def ANA_mom_advection_x(self):
        return bool(self.__iof_ana[iof_ana.ANA_mom_advection_x.value])

    @ANA_mom_advection_x.setter
    def ANA_mom_advection_x(self, ANA_mom_advection_x):
        assert isinstance(ANA_mom_advection_x, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_mom_advection_x.value] = int(ANA_mom_advection_x)

    @property
    def ANA_mom_advection_y(self):
        return bool(self.__iof_ana[iof_ana.ANA_mom_advection_y.value])

    @ANA_mom_advection_y.setter
    def ANA_mom_advection_y(self, ANA_mom_advection_y):
        assert isinstance(ANA_mom_advection_y, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_mom_advection_y.value] = int(ANA_mom_advection_y)

    @property
    def ANA_Richardson(self):
        return bool(self.__iof_ana[iof_ana.ANA_Richardson.value])

    @ANA_Richardson.setter
    def ANA_Richardson(self, ANA_Richardson):
        assert isinstance(ANA_Richardson, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_Richardson.value] = int(ANA_Richardson)

    @property
    def ANA_transport_min_dt_elem(self):
        return bool(self.__iof_ana[iof_ana.ANA_transport_min_dt_elem.value])

    @ANA_transport_min_dt_elem.setter
    def ANA_transport_min_dt_elem(self, ANA_transport_min_dt_elem):
        assert isinstance(ANA_transport_min_dt_elem, bool), 'property must be a bool'
        self.__iof_ana[iof_ana.ANA_transport_min_dt_elem.value] = int(ANA_transport_min_dt_elem)
