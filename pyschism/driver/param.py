import pathlib
from functools import lru_cache
import f90nml
from datetime import datetime, timedelta, timezone
from math import modf
from pyschism.driver.core import CORE
from pyschism.driver.opt import OPT
from pyschism.driver.schout import SCHOUT


class Param:

    @staticmethod
    def open(path):
        nml = f90nml.read(path)
        param = Param()
        for group, attrs in nml.items():
            for attr, value in attrs.items():
                setattr(param, f"{group}.{attr}", value)
                # exec(f"param.{group}.{attr}={value}")

        # year = nml['opt']['year'] if 'year' in nml['opt'] else 1
        # month = nml['opt']['month'] if 'month' in nml['opt'] else 1
        # day = nml['opt']['day'] if 'day' in nml['opt'] else 1
        # hour = nml['opt']['hour'] if 'hour' in nml['opt'] else 0
        # param.start_date = datetime(year, month, day, hour)
        # param.
        # param.end_date = param.start_date + \
            # timedelta(days=nml['core']['rnday'])
        # adjust time for utc
        # dt = nml['opt']['utc_start'] if 'utc_start' in nml['opt'] else 0
        # param.start_date = 
        # param.run_time = 
        # param.spinup_time = 
        return param

    def write(self, path, overwrite=False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            msg = f"File {path} exists and overwrite=False"
            raise IOError(msg)
        f90nml.patch(self.src, self._nml, path)

    @property
    @lru_cache(maxsize=None)
    def core(self):
        return CORE(self._nml)

    @property
    @lru_cache(maxsize=None)
    def opt(self):
        return OPT(self._nml)

    @property
    @lru_cache(maxsize=None)
    def schout(self):
        return SCHOUT(self._nml)

    @property
    @lru_cache(maxsize=None)
    def src(self):
        return pathlib.Path(__file__).parent / 'param.nml'

    @property
    def start_date(self):
        return self._forcing_start_date + self._spinup_time

    @property
    def end_date(self):
        return self.start_date + self._run_time

    @property
    def run_time(self):
        return timedelta(days=self.core.rnday)

    @property
    def spinup_time(self):
        return self._spinup_time

    @property
    def forcing_start_date(self):
        minutes, hour = modf(self.opt.start_hour)
        return datetime(
            self.opt.start_year,
            self.opt.start_month,
            self.opt.start_day,
            int(hour % 24),
            int(minutes % 60),
            tzinfo=timezone(timedelta(hours=self.opt.utc_start))
            )

    @property
    def timestep(self):
        return self.core.dt

    @property
    def output_sampling_rate(self):
        return timedelta(seconds=self.opt.nspool / self.timestep)

    @property
    def nspool(self):
        return self.core.nspool

    @property
    def ihfskip(self):
        return self.core.ihfskip

    @start_date.setter
    def start_date(self, start_date):
        self._start_date = start_date

    @run_time.setter
    def run_time(self, run_time):
        assert isinstance(run_time, timedelta)
        msg = "run_time must be positive."
        assert run_time.total_seconds() >= 0, msg
        self._run_time = run_time

    @spinup_time.setter
    def spinup_time(self, spinup_time):
        assert isinstance(spinup_time, timedelta)
        total_seconds = spinup_time.total_seconds()
        msg = "spinup_time must be larger or equal than zero."
        assert total_seconds >= 0, msg
        if total_seconds == 0:
            self.opt.nramp = 0
            self.opt.dramp = 0.
        else:
            self.opt.nramp = 1
            self.opt.dramp = float(total_seconds / (60. * 60. * 24.))
        self._forcing_start_date = self.start_date - spinup_time
        self.__spinup_time = spinup_time

    @nspool.setter
    def nspool(self, nspool):
        if isinstance(nspool, timedelta):
            nspool = int(round(nspool.total_seconds() / self.timestep))
        msg = "nspool must be a timedelta object, or an integer greater or "
        msg += "equal to zero."
        assert isinstance(nspool, (int, timedelta))
        self.core.nspool = nspool

    @ihfskip.setter
    def ihfskip(self, ihfskip):
        self.core.ihfskip = ihfskip

    @timestep.setter
    def timestep(self, timestep):
        self.core.dt = timestep

    @property
    @lru_cache(maxsize=None)
    def _nml(self):
        return f90nml.read(self.src.resolve())

    @property
    def _start_date(self):
        return self.__start_date

    @property
    def _end_date(self):
        return self._start_date + self._run_time

    @property
    def _forcing_start_date(self):
        return self.__forcing_start_date

    @property
    def _run_time(self):
        try:
            return self.__run_time
        except AttributeError:
            return timedelta(0.)

    @property
    def _spinup_time(self):
        try:
            return self.__spinup_time
        except AttributeError:
            return timedelta(0.)

    @_start_date.setter
    def _start_date(self, start_date):
        # https://stackoverflow.com/questions/5802108/how-to-check-if-a-datetime-object-is-localized-with-pytz
        assert isinstance(start_date, datetime)
        self._forcing_start_date = start_date - self.spinup_time
        self.__start_date = start_date

    @_forcing_start_date.setter
    def _forcing_start_date(self, forcing_start_date):
        assert isinstance(forcing_start_date, datetime)
        if forcing_start_date.tzinfo is not None \
                and forcing_start_date.tzinfo.utcoffset(forcing_start_date) is not None:
            self.opt.utc_start = forcing_start_date.utcoffset().total_seconds()/3600
        self.opt.start_year = forcing_start_date.year
        self.opt.start_month = forcing_start_date.month
        self.opt.start_day = forcing_start_date.day
        self.opt.start_hour = forcing_start_date.hour
        self.opt.start_hour += forcing_start_date.minute / 60.
        self.__forcing_start_date = forcing_start_date

    @_run_time.setter
    def _run_time(self, run_time):
        assert isinstance(run_time, timedelta)
        self.__run_time = run_time

    # --------------------------------------------------
    # The following are pointers to self.schout.iof_*[*]
    # --------------------------------------------------

    @property
    def elev(self):
        return self.schout.elev

    @property
    def air_pressure(self):
        return self.schout.air_pressure

    @property
    def air_temperature(self):
        return self.schout.air_temperature

    @property
    def specific_humidity(self):
        return self.schout.specific_humidity

    @property
    def solar_radiation(self):
        return self.schout.solar_radiation

    @property
    def sensible_flux(self):
        return self.schout.sensible_flux

    @property
    def latent_heat(self):
        return self.schout.latent_heat

    @property
    def upward_longwave(self):
        return self.schout.upward_longwave

    @property
    def downward_longwave(self):
        return self.schout.downward_longwave

    @property
    def total_heat_flux(self):
        return self.schout.total_heat_flux

    @property
    def evaporation(self):
        return self.schout.evaporation

    @property
    def precipitation(self):
        return self.schout.precipitation

    @property
    def bottom_stress(self):
        return self.schout.bottom_stress

    @property
    def wind_speed(self):
        return self.schout.wind_speed

    @property
    def wind_stress(self):
        return self.schout.wind_stress

    @property
    def dahv(self):
        return self.schout.dahv

    @property
    def vertical_velocity(self):
        return self.schout.vertical_velocity

    @property
    def temp(self):
        return self.schout.temp

    @property
    def salt(self):
        return self.schout.salt

    @property
    def water_density(self):
        return self.schout.water_density

    @property
    def diffusivity(self):
        return self.schout.diffusivity

    @property
    def viscosity(self):
        return self.schout.viscosity

    @property
    def TKE(self):
        return self.schout.TKE

    @property
    def mixing_length(self):
        return self.schout.mixing_length

    @property
    def hvel(self):
        return self.schout.hvel

    @property
    def hvel_side(self):
        return self.schout.hvel_side

    @property
    def wvel_elem(self):
        return self.schout.wvel_elem

    @property
    def temp_elem(self):
        return self.schout.temp_elem

    @property
    def salt_elem(self):
        return self.schout.salt_elem

    @property
    def pressure_gradient(self):
        return self.schout.pressure_gradient

    @property
    def WWM_1(self):
        return self.schout.WWM_1

    @property
    def WWM_2(self):
        return self.schout.WWM_2

    @property
    def WWM_3(self):
        return self.schout.WWM_3

    @property
    def WWM_4(self):
        return self.schout.WWM_4

    @property
    def WWM_5(self):
        return self.schout.WWM_5

    @property
    def WWM_6(self):
        return self.schout.WWM_6

    @property
    def WWM_9(self):
        return self.schout.WWM_9

    @property
    def WWM_10(self):
        return self.schout.WWM_10

    @property
    def WWM_11(self):
        return self.schout.WWM_11

    @property
    def WWM_12(self):
        return self.schout.WWM_12

    @property
    def WWM_13(self):
        return self.schout.WWM_13

    @property
    def WWM_14(self):
        return self.schout.WWM_14

    @property
    def WWM_15(self):
        return self.schout.WWM_15

    @property
    def WWM_16(self):
        return self.schout.WWM_16

    @property
    def WWM_17(self):
        return self.schout.WWM_17

    @property
    def WWM_18(self):
        return self.schout.WWM_18

    @property
    def WWM_19(self):
        return self.schout.WWM_19

    @property
    def WWM_20(self):
        return self.schout.WWM_20

    @property
    def WWM_21(self):
        return self.schout.WWM_21

    @property
    def WWM_22(self):
        return self.schout.WWM_22

    @property
    def WWM_23(self):
        return self.schout.WWM_23

    @property
    def WWM_24(self):
        return self.schout.WWM_24

    @property
    def WWM_25(self):
        return self.schout.WWM_25

    @property
    def WWM_26(self):
        return self.schout.WWM_26

    @property
    def WWM_27(self):
        return self.schout.WWM_27

    @property
    def WWM_28(self):
        return self.schout.WWM_28

    @property
    def WWM_energy_dir(self):
        return self.schout.WWM_energy_dir

    @property
    def wave_force(self):
        return self.schout.wave_force

    @property
    def GEN_1(self):
        return self.schout.GEN_1

    @property
    def GEN_2(self):
        return self.schout.GEN_2

    @property
    def AGE_1(self):
        return self.schout.AGE_1

    @property
    def AGE_2(self):
        return self.schout.AGE_2

    @property
    def SED_depth_change(self):
        return self.schout.SED_depth_change

    @property
    def SED_D50(self):
        return self.schout.SED_D50

    @property
    def SED_bed_stress(self):
        return self.schout.SED_bed_stress

    @property
    def SED_bed_roughness(self):
        return self.schout.SED_bed_roughness

    @property
    def SED_TSC(self):
        return self.schout.SED_TSC

    @property
    def bed_thickness(self):
        return self.schout.bed_thickness

    @property
    def bed_age(self):
        return self.schout.bed_age

    @property
    def z0st(self):
        return self.schout.z0st

    @property
    def z0cr(self):
        return self.schout.z0cr

    @property
    def z0sw(self):
        return self.schout.z0sw

    @property
    def z0wr(self):
        return self.schout.z0wr

    @property
    def SED3D_1(self):
        return self.schout.SED3D_1

    @property
    def SED_bdld_1(self):
        return self.schout.SED_bdld_1

    @property
    def SED_bedfrac_1(self):
        return self.schout.SED_bedfrac_1

    @property
    def SED3D_2(self):
        return self.schout.SED3D_2

    @property
    def SED_bdld_2(self):
        return self.schout.SED_bdld_2

    @property
    def SED_bedfrac_3(self):
        return self.schout.SED_bedfrac_3

    @property
    def ECO_1(self):
        return self.schout.ECO_1

    @property
    def ICM_Chl(self):
        return self.schout.ICM_Chl

    @property
    def ICM_pH(self):
        return self.schout.ICM_pH

    @property
    def ICM_PrmPrdt(self):
        return self.schout.ICM_PrmPrdt

    @property
    def ICM_DIN(self):
        return self.schout.ICM_DIN

    @property
    def ICM_PON(self):
        return self.schout.ICM_PON

    @property
    def ICM_SED_BENDOC(self):
        return self.schout.ICM_SED_BENDOC

    @property
    def ICM_SED_BENNH4(self):
        return self.schout.ICM_SED_BENNH4

    @property
    def ICM_SED_BENNO3(self):
        return self.schout.ICM_SED_BENNO3

    @property
    def ICM_SED_BENPO4(self):
        return self.schout.ICM_SED_BENPO4

    @property
    def ICM_SED_BENCOD(self):
        return self.schout.ICM_SED_BENCOD

    @property
    def ICM_SED_BENDO(self):
        return self.schout.ICM_SED_BENDO

    @property
    def ICM_SED_BENSA(self):
        return self.schout.ICM_SED_BENSA

    @property
    def ICM_lfsav(self):
        return self.schout.ICM_lfsav

    @property
    def ICM_stsav(self):
        return self.schout.ICM_stsav

    @property
    def ICM_rtsav(self):
        return self.schout.ICM_rtsav

    @property
    def ICM_tlfsav(self):
        return self.schout.ICM_tlfsav

    @property
    def ICM_tstsav(self):
        return self.schout.ICM_tstsav

    @property
    def ICM_trtsav(self):
        return self.schout.ICM_trtsav

    @property
    def ICM_hcansav(self):
        return self.schout.ICM_hcansav

    @property
    def ICM_CNH4(self):
        return self.schout.ICM_CNH4

    @property
    def ICM_CNH3(self):
        return self.schout.ICM_CNH3

    @property
    def ICM_CPIP(self):
        return self.schout.ICM_CPIP

    @property
    def ICM_CPOS(self):
        return self.schout.ICM_CPOS

    @property
    def ICM_CCH4(self):
        return self.schout.ICM_CCH4

    @property
    def ICM_CSO4(self):
        return self.schout.ICM_CSO4

    @property
    def ICM_CH2S(self):
        return self.schout.ICM_CH2S

    @property
    def ICM_SEDPON1(self):
        return self.schout.ICM_SEDPON1

    @property
    def ICM_SEDPON2(self):
        return self.schout.ICM_SEDPON2

    @property
    def ICM_SEDPON3(self):
        return self.schout.ICM_SEDPON3

    @property
    def ICM_SEDPOP1(self):
        return self.schout.ICM_SEDPOP1

    @property
    def ICM_SEDPOP2(self):
        return self.schout.ICM_SEDPOP2

    @property
    def ICM_SEDPOP3(self):
        return self.schout.ICM_SEDPOP3

    @property
    def ICM_SEDPOC1(self):
        return self.schout.ICM_SEDPOC1

    @property
    def ICM_SEDPOC2(self):
        return self.schout.ICM_SEDPOC2

    @property
    def ICM_SEDPOC3(self):
        return self.schout.ICM_SEDPOC3

    @property
    def ICM_EROH2S(self):
        return self.schout.ICM_EROH2S

    @property
    def ICM_EROLPOC(self):
        return self.schout.ICM_EROLPOC

    @property
    def ICM_ERORPOC(self):
        return self.schout.ICM_ERORPOC

    @property
    def ICM_DO_consumption(self):
        return self.schout.ICM_DO_consumption

    @property
    def ICM_GP1(self):
        return self.schout.ICM_GP1

    @property
    def ICM_GP2(self):
        return self.schout.ICM_GP2

    @property
    def ICM_GP3(self):
        return self.schout.ICM_GP3

    @property
    def ICM_1(self):
        return self.schout.ICM_1

    @property
    def ICM_2(self):
        return self.schout.ICM_2

    @property
    def ICM_3(self):
        return self.schout.ICM_3

    @property
    def ICM_4(self):
        return self.schout.ICM_4

    @property
    def ICM_5(self):
        return self.schout.ICM_5

    @property
    def ICM_6(self):
        return self.schout.ICM_6

    @property
    def ICM_7(self):
        return self.schout.ICM_7

    @property
    def ICM_8(self):
        return self.schout.ICM_8

    @property
    def ICM_9(self):
        return self.schout.ICM_9

    @property
    def ICM_10(self):
        return self.schout.ICM_10

    @property
    def ICM_11(self):
        return self.schout.ICM_11

    @property
    def ICM_12(self):
        return self.schout.ICM_12

    @property
    def ICM_13(self):
        return self.schout.ICM_13

    @property
    def ICM_14(self):
        return self.schout.ICM_14

    @property
    def ICM_15(self):
        return self.schout.ICM_15

    @property
    def ICM_16(self):
        return self.schout.ICM_16

    @property
    def ICM_17(self):
        return self.schout.ICM_17

    @property
    def ICM_18(self):
        return self.schout.ICM_18

    @property
    def ICM_19(self):
        return self.schout.ICM_19

    @property
    def ICM_20(self):
        return self.schout.ICM_20

    @property
    def ICM_21(self):
        return self.schout.ICM_21

    @property
    def ICM_22(self):
        return self.schout.ICM_22

    @property
    def ICM_23(self):
        return self.schout.ICM_23

    @property
    def ICM_24(self):
        return self.schout.ICM_24

    @property
    def ICM_25(self):
        return self.schout.ICM_25

    @property
    def COS_1(self):
        return self.schout.COS_1

    @property
    def COS_2(self):
        return self.schout.COS_2

    @property
    def COS_3(self):
        return self.schout.COS_3

    @property
    def COS_4(self):
        return self.schout.COS_4

    @property
    def COS_5(self):
        return self.schout.COS_5

    @property
    def COS_6(self):
        return self.schout.COS_6

    @property
    def COS_7(self):
        return self.schout.COS_7

    @property
    def COS_8(self):
        return self.schout.COS_8

    @property
    def COS_9(self):
        return self.schout.COS_9

    @property
    def COS_10(self):
        return self.schout.COS_10

    @property
    def COS_11(self):
        return self.schout.COS_11

    @property
    def COS_12(self):
        return self.schout.COS_12

    @property
    def COS_13(self):
        return self.schout.COS_13

    @property
    def FIB_1(self):
        return self.schout.FIB_1

    @property
    def SED2D_depth_change(self):
        return self.schout.SED2D_depth_change

    @property
    def SED2D_cflsed(self):
        return self.schout.SED2D_cflsed

    @property
    def SED2D_d50(self):
        return self.schout.SED2D_d50

    @property
    def SED2D_total_transport(self):
        return self.schout.SED2D_total_transport

    @property
    def SED2D_susp_load(self):
        return self.schout.SED2D_susp_load

    @property
    def SED2D_bed_load(self):
        return self.schout.SED2D_bed_load

    @property
    def SED2D_average_transport(self):
        return self.schout.SED2D_average_transport

    @property
    def SED2D_bottom_slope(self):
        return self.schout.SED2D_bottom_slope

    @property
    def z0eq(self):
        return self.schout.z0eq

    @property
    def z0cr2d(self):
        return self.schout.z0cr2d

    @property
    def z0sw2d(self):
        return self.schout.z0sw2d

    @property
    def z0wr2d(self):
        return self.schout.z0wr2d

    @property
    def marsh_flag(self):
        return self.schout.marsh_flag

    @property
    def ICE_velocity(self):
        return self.schout.ICE_velocity

    @property
    def ICE_strain_rate(self):
        return self.schout.ICE_strain_rate

    @property
    def ICE_net_heat_flux(self):
        return self.schout.ICE_net_heat_flux

    @property
    def ICE_fresh_water_flux(self):
        return self.schout.ICE_fresh_water_flux

    @property
    def ICE_top_T(self):
        return self.schout.ICE_top_T

    @property
    def ICE_tracer_1(self):
        return self.schout.ICE_tracer_1

    @property
    def ICE_tracer_2(self):
        return self.schout.ICE_tracer_2

    @property
    def ICE_tracer_3(self):
        return self.schout.ICE_tracer_3

    @property
    def ANA_air_pres_grad_x(self):
        return self.schout.ANA_air_pres_grad_x

    @property
    def ANA_air_pres_grad_y(self):
        return self.schout.ANA_air_pres_grad_y

    @property
    def ANA_tide_pot_grad_x(self):
        return self.schout.ANA_tide_pot_grad_x

    @property
    def ANA_tide_pot_grad_y(self):
        return self.schout.ANA_tide_pot_grad_y

    @property
    def ANA_hor_viscosity_x(self):
        return self.schout.ANA_hor_viscosity_x

    @property
    def ANA_hor_viscosity_y(self):
        return self.schout.ANA_hor_viscosity_y

    @property
    def ANA_bclinic_force_x(self):
        return self.schout.ANA_bclinic_force_x

    @property
    def ANA_bclinic_force_y(self):
        return self.schout.ANA_bclinic_force_y

    @property
    def ANA_vert_viscosity_x(self):
        return self.schout.ANA_vert_viscosity_x

    @property
    def ANA_vert_viscosity_y(self):
        return self.schout.ANA_vert_viscosity_y

    @property
    def ANA_mom_advection_x(self):
        return self.schout.ANA_mom_advection_x

    @property
    def ANA_mom_advection_y(self):
        return self.schout.ANA_mom_advection_y

    @property
    def ANA_Richardson(self):
        return self.schout.ANA_Richardson

    @property
    def ANA_transport_min_dt_elem(self):
        return self.schout.ANA_transport_min_dt_elem

    @elev.setter
    def elev(self, elev):
        self.schout.elev = elev

    @air_pressure.setter
    def air_pressure(self, air_pressure):
        self.schout.air_pressure = air_pressure

    @air_temperature.setter
    def air_temperature(self, air_temperature):
        self.schout.air_temperature = air_temperature

    @specific_humidity.setter
    def specific_humidity(self, specific_humidity):
        self.schout.specific_humidity = specific_humidity

    @solar_radiation.setter
    def solar_radiation(self, solar_radiation):
        self.schout.solar_radiation = solar_radiation

    @sensible_flux.setter
    def sensible_flux(self, sensible_flux):
        self.schout.sensible_flux = sensible_flux

    @latent_heat.setter
    def latent_heat(self, latent_heat):
        self.schout.latent_heat = latent_heat

    @upward_longwave.setter
    def upward_longwave(self, upward_longwave):
        self.schout.upward_longwave = upward_longwave

    @downward_longwave.setter
    def downward_longwave(self, downward_longwave):
        self.schout.downward_longwave = downward_longwave

    @total_heat_flux.setter
    def total_heat_flux(self, total_heat_flux):
        self.schout.total_heat_flux = total_heat_flux

    @evaporation.setter
    def evaporation(self, evaporation):
        self.schout.evaporation = evaporation

    @precipitation.setter
    def precipitation(self, precipitation):
        self.schout.precipitation = precipitation

    @bottom_stress.setter
    def bottom_stress(self, bottom_stress):
        self.schout.bottom_stress = bottom_stress

    @wind_speed.setter
    def wind_speed(self, wind_speed):
        self.schout.wind_speed = wind_speed

    @wind_stress.setter
    def wind_stress(self, wind_stress):
        self.schout.wind_stress = wind_stress

    @dahv.setter
    def dahv(self, dahv):
        self.schout.dahv = dahv

    @vertical_velocity.setter
    def vertical_velocity(self, vertical_velocity):
        self.schout.vertical_velocity = vertical_velocity

    @temp.setter
    def temp(self, temp):
        self.schout.temp = temp

    @salt.setter
    def salt(self, salt):
        self.schout.salt = salt

    @water_density.setter
    def water_density(self, water_density):
        self.schout.water_density = water_density

    @diffusivity.setter
    def diffusivity(self, diffusivity):
        self.schout.diffusivity = diffusivity

    @viscosity.setter
    def viscosity(self, viscosity):
        self.schout.viscosity = viscosity

    @TKE.setter
    def TKE(self, TKE):
        self.schout.TKE = TKE

    @mixing_length.setter
    def mixing_length(self, mixing_length):
        self.schout.mixing_length = mixing_length

    @hvel.setter
    def hvel(self, hvel):
        self.schout.hvel = hvel

    @hvel_side.setter
    def hvel_side(self, hvel_side):
        self.schout.hvel_side = hvel_side

    @wvel_elem.setter
    def wvel_elem(self, wvel_elem):
        self.schout.wvel_elem = wvel_elem

    @temp_elem.setter
    def temp_elem(self, temp_elem):
        self.schout.temp_elem = temp_elem

    @salt_elem.setter
    def salt_elem(self, salt_elem):
        self.schout.salt_elem = salt_elem

    @pressure_gradient.setter
    def pressure_gradient(self, pressure_gradient):
        self.schout.pressure_gradient = pressure_gradient

    @WWM_1.setter
    def WWM_1(self, WWM_1):
        self.schout.WWM_1 = WWM_1

    @WWM_2.setter
    def WWM_2(self, WWM_2):
        self.schout.WWM_2 = WWM_2

    @WWM_3.setter
    def WWM_3(self, WWM_3):
        self.schout.WWM_3 = WWM_3

    @WWM_4.setter
    def WWM_4(self, WWM_4):
        self.schout.WWM_4 = WWM_4

    @WWM_5.setter
    def WWM_5(self, WWM_5):
        self.schout.WWM_5 = WWM_5

    @WWM_6.setter
    def WWM_6(self, WWM_6):
        self.schout.WWM_6 = WWM_6

    @WWM_9.setter
    def WWM_9(self, WWM_9):
        self.schout.WWM_9 = WWM_9

    @WWM_10.setter
    def WWM_10(self, WWM_10):
        self.schout.WWM_10 = WWM_10

    @WWM_11.setter
    def WWM_11(self, WWM_11):
        self.schout.WWM_11 = WWM_11

    @WWM_12.setter
    def WWM_12(self, WWM_12):
        self.schout.WWM_12 = WWM_12

    @WWM_13.setter
    def WWM_13(self, WWM_13):
        self.schout.WWM_13 = WWM_13

    @WWM_14.setter
    def WWM_14(self, WWM_14):
        self.schout.WWM_14 = WWM_14

    @WWM_15.setter
    def WWM_15(self, WWM_15):
        self.schout.WWM_15 = WWM_15

    @WWM_16.setter
    def WWM_16(self, WWM_16):
        self.schout.WWM_16 = WWM_16

    @WWM_17.setter
    def WWM_17(self, WWM_17):
        self.schout.WWM_17 = WWM_17

    @WWM_18.setter
    def WWM_18(self, WWM_18):
        self.schout.WWM_18 = WWM_18

    @WWM_19.setter
    def WWM_19(self, WWM_19):
        self.schout.WWM_19 = WWM_19

    @WWM_20.setter
    def WWM_20(self, WWM_20):
        self.schout.WWM_20 = WWM_20

    @WWM_21.setter
    def WWM_21(self, WWM_21):
        self.schout.WWM_21 = WWM_21

    @WWM_22.setter
    def WWM_22(self, WWM_22):
        self.schout.WWM_22 = WWM_22

    @WWM_23.setter
    def WWM_23(self, WWM_23):
        self.schout.WWM_23 = WWM_23

    @WWM_24.setter
    def WWM_24(self, WWM_24):
        self.schout.WWM_24 = WWM_24

    @WWM_25.setter
    def WWM_25(self, WWM_25):
        self.schout.WWM_25 = WWM_25

    @WWM_26.setter
    def WWM_26(self, WWM_26):
        self.schout.WWM_26 = WWM_26

    @WWM_27.setter
    def WWM_27(self, WWM_27):
        self.schout.WWM_27 = WWM_27

    @WWM_28.setter
    def WWM_28(self, WWM_28):
        self.schout.WWM_28 = WWM_28

    @WWM_energy_dir.setter
    def WWM_energy_dir(self, WWM_energy_dir):
        self.schout.WWM_energy_dir = WWM_energy_dir

    @wave_force.setter
    def wave_force(self, wave_force):
        self.schout.wave_force = wave_force

    @GEN_1.setter
    def GEN_1(self, GEN_1):
        self.schout.GEN_1 = GEN_1

    @GEN_2.setter
    def GEN_2(self, GEN_2):
        self.schout.GEN_2 = GEN_2

    @AGE_1.setter
    def AGE_1(self, AGE_1):
        self.schout.AGE_1 = AGE_1

    @AGE_2.setter
    def AGE_2(self, AGE_2):
        self.schout.AGE_2 = AGE_2

    @SED_depth_change.setter
    def SED_depth_change(self, SED_depth_change):
        self.schout.SED_depth_change = SED_depth_change

    @SED_D50.setter
    def SED_D50(self, SED_D50):
        self.schout.SED_D50 = SED_D50

    @SED_bed_stress.setter
    def SED_bed_stress(self, SED_bed_stress):
        self.schout.SED_bed_stress = SED_bed_stress

    @SED_bed_roughness.setter
    def SED_bed_roughness(self, SED_bed_roughness):
        self.schout.SED_bed_roughness = SED_bed_roughness

    @SED_TSC.setter
    def SED_TSC(self, SED_TSC):
        self.schout.SED_TSC = SED_TSC

    @bed_thickness.setter
    def bed_thickness(self, bed_thickness):
        self.schout.bed_thickness = bed_thickness

    @bed_age.setter
    def bed_age(self, bed_age):
        self.schout.bed_age = bed_age

    @z0st.setter
    def z0st(self, z0st):
        self.schout.z0st = z0st

    @z0cr.setter
    def z0cr(self, z0cr):
        self.schout.z0cr = z0cr

    @z0sw.setter
    def z0sw(self, z0sw):
        self.schout.z0sw = z0sw

    @z0wr.setter
    def z0wr(self, z0wr):
        self.schout.z0wr = z0wr

    @SED3D_1.setter
    def SED3D_1(self, SED3D_1):
        self.schout.SED3D_1 = SED3D_1

    @SED_bdld_1.setter
    def SED_bdld_1(self, SED_bdld_1):
        self.schout.SED_bdld_1 = SED_bdld_1

    @SED_bedfrac_1.setter
    def SED_bedfrac_1(self, SED_bedfrac_1):
        self.schout.SED_bedfrac_1 = SED_bedfrac_1

    @SED3D_2.setter
    def SED3D_2(self, SED3D_2):
        self.schout.SED3D_2 = SED3D_2

    @SED_bdld_2.setter
    def SED_bdld_2(self, SED_bdld_2):
        self.schout.SED_bdld_2 = SED_bdld_2

    @SED_bedfrac_3.setter
    def SED_bedfrac_3(self, SED_bedfrac_3):
        self.schout.SED_bedfrac_3 = SED_bedfrac_3

    @ECO_1.setter
    def ECO_1(self, ECO_1):
        self.schout.ECO_1 = ECO_1

    @ICM_Chl.setter
    def ICM_Chl(self, ICM_Chl):
        self.schout.ICM_Chl = ICM_Chl

    @ICM_pH.setter
    def ICM_pH(self, ICM_pH):
        self.schout.ICM_pH = ICM_pH

    @ICM_PrmPrdt.setter
    def ICM_PrmPrdt(self, ICM_PrmPrdt):
        self.schout.ICM_PrmPrdt = ICM_PrmPrdt

    @ICM_DIN.setter
    def ICM_DIN(self, ICM_DIN):
        self.schout.ICM_DIN = ICM_DIN

    @ICM_PON.setter
    def ICM_PON(self, ICM_PON):
        self.schout.ICM_PON = ICM_PON

    @ICM_SED_BENDOC.setter
    def ICM_SED_BENDOC(self, ICM_SED_BENDOC):
        self.schout.ICM_SED_BENDOC = ICM_SED_BENDOC

    @ICM_SED_BENNH4.setter
    def ICM_SED_BENNH4(self, ICM_SED_BENNH4):
        self.schout.ICM_SED_BENNH4 = ICM_SED_BENNH4

    @ICM_SED_BENNO3.setter
    def ICM_SED_BENNO3(self, ICM_SED_BENNO3):
        self.schout.ICM_SED_BENNO3 = ICM_SED_BENNO3

    @ICM_SED_BENPO4.setter
    def ICM_SED_BENPO4(self, ICM_SED_BENPO4):
        self.schout.ICM_SED_BENPO4 = ICM_SED_BENPO4

    @ICM_SED_BENCOD.setter
    def ICM_SED_BENCOD(self, ICM_SED_BENCOD):
        self.schout.ICM_SED_BENCOD = ICM_SED_BENCOD

    @ICM_SED_BENDO.setter
    def ICM_SED_BENDO(self, ICM_SED_BENDO):
        self.schout.ICM_SED_BENDO = ICM_SED_BENDO

    @ICM_SED_BENSA.setter
    def ICM_SED_BENSA(self, ICM_SED_BENSA):
        self.schout.ICM_SED_BENSA = ICM_SED_BENSA

    @ICM_lfsav.setter
    def ICM_lfsav(self, ICM_lfsav):
        self.schout.ICM_lfsav = ICM_lfsav

    @ICM_stsav.setter
    def ICM_stsav(self, ICM_stsav):
        self.schout.ICM_stsav = ICM_stsav

    @ICM_rtsav.setter
    def ICM_rtsav(self, ICM_rtsav):
        self.schout.ICM_rtsav = ICM_rtsav

    @ICM_tlfsav.setter
    def ICM_tlfsav(self, ICM_tlfsav):
        self.schout.ICM_tlfsav = ICM_tlfsav

    @ICM_tstsav.setter
    def ICM_tstsav(self, ICM_tstsav):
        self.schout.ICM_tstsav = ICM_tstsav

    @ICM_trtsav.setter
    def ICM_trtsav(self, ICM_trtsav):
        self.schout.ICM_trtsav = ICM_trtsav

    @ICM_hcansav.setter
    def ICM_hcansav(self, ICM_hcansav):
        self.schout.ICM_hcansav = ICM_hcansav

    @ICM_CNH4.setter
    def ICM_CNH4(self, ICM_CNH4):
        self.schout.ICM_CNH4 = ICM_CNH4

    @ICM_CNH3.setter
    def ICM_CNH3(self, ICM_CNH3):
        self.schout.ICM_CNH3 = ICM_CNH3

    @ICM_CPIP.setter
    def ICM_CPIP(self, ICM_CPIP):
        self.schout.ICM_CPIP = ICM_CPIP

    @ICM_CPOS.setter
    def ICM_CPOS(self, ICM_CPOS):
        self.schout.ICM_CPOS = ICM_CPOS

    @ICM_CCH4.setter
    def ICM_CCH4(self, ICM_CCH4):
        self.schout.ICM_CCH4 = ICM_CCH4

    @ICM_CSO4.setter
    def ICM_CSO4(self, ICM_CSO4):
        self.schout.ICM_CSO4 = ICM_CSO4

    @ICM_CH2S.setter
    def ICM_CH2S(self, ICM_CH2S):
        self.schout.ICM_CH2S = ICM_CH2S

    @ICM_SEDPON1.setter
    def ICM_SEDPON1(self, ICM_SEDPON1):
        self.schout.ICM_SEDPON1 = ICM_SEDPON1

    @ICM_SEDPON2.setter
    def ICM_SEDPON2(self, ICM_SEDPON2):
        self.schout.ICM_SEDPON2 = ICM_SEDPON2

    @ICM_SEDPON3.setter
    def ICM_SEDPON3(self, ICM_SEDPON3):
        self.schout.ICM_SEDPON3 = ICM_SEDPON3

    @ICM_SEDPOP1.setter
    def ICM_SEDPOP1(self, ICM_SEDPOP1):
        self.schout.ICM_SEDPOP1 = ICM_SEDPOP1

    @ICM_SEDPOP2.setter
    def ICM_SEDPOP2(self, ICM_SEDPOP2):
        self.schout.ICM_SEDPOP2 = ICM_SEDPOP2

    @ICM_SEDPOP3.setter
    def ICM_SEDPOP3(self, ICM_SEDPOP3):
        self.schout.ICM_SEDPOP3 = ICM_SEDPOP3

    @ICM_SEDPOC1.setter
    def ICM_SEDPOC1(self, ICM_SEDPOC1):
        self.schout.ICM_SEDPOC1 = ICM_SEDPOC1

    @ICM_SEDPOC2.setter
    def ICM_SEDPOC2(self, ICM_SEDPOC2):
        self.schout.ICM_SEDPOC2 = ICM_SEDPOC2

    @ICM_SEDPOC3.setter
    def ICM_SEDPOC3(self, ICM_SEDPOC3):
        self.schout.ICM_SEDPOC3 = ICM_SEDPOC3

    @ICM_EROH2S.setter
    def ICM_EROH2S(self, ICM_EROH2S):
        self.schout.ICM_EROH2S = ICM_EROH2S

    @ICM_EROLPOC.setter
    def ICM_EROLPOC(self, ICM_EROLPOC):
        self.schout.ICM_EROLPOC = ICM_EROLPOC

    @ICM_ERORPOC.setter
    def ICM_ERORPOC(self, ICM_ERORPOC):
        self.schout.ICM_ERORPOC = ICM_ERORPOC

    @ICM_DO_consumption.setter
    def ICM_DO_consumption(self, ICM_DO_consumption):
        self.schout.ICM_DO_consumption = ICM_DO_consumption

    @ICM_GP1.setter
    def ICM_GP1(self, ICM_GP1):
        self.schout.ICM_GP1 = ICM_GP1

    @ICM_GP2.setter
    def ICM_GP2(self, ICM_GP2):
        self.schout.ICM_GP2 = ICM_GP2

    @ICM_GP3.setter
    def ICM_GP3(self, ICM_GP3):
        self.schout.ICM_GP3 = ICM_GP3

    @ICM_1.setter
    def ICM_1(self, ICM_1):
        self.schout.ICM_1 = ICM_1

    @ICM_2.setter
    def ICM_2(self, ICM_2):
        self.schout.ICM_2 = ICM_2

    @ICM_3.setter
    def ICM_3(self, ICM_3):
        self.schout.ICM_3 = ICM_3

    @ICM_4.setter
    def ICM_4(self, ICM_4):
        self.schout.ICM_4 = ICM_4

    @ICM_5.setter
    def ICM_5(self, ICM_5):
        self.schout.ICM_5 = ICM_5

    @ICM_6.setter
    def ICM_6(self, ICM_6):
        self.schout.ICM_6 = ICM_6

    @ICM_7.setter
    def ICM_7(self, ICM_7):
        self.schout.ICM_7 = ICM_7

    @ICM_8.setter
    def ICM_8(self, ICM_8):
        self.schout.ICM_8 = ICM_8

    @ICM_9.setter
    def ICM_9(self, ICM_9):
        self.schout.ICM_9 = ICM_9

    @ICM_10.setter
    def ICM_10(self, ICM_10):
        self.schout.ICM_10 = ICM_10

    @ICM_11.setter
    def ICM_11(self, ICM_11):
        self.schout.ICM_11 = ICM_11

    @ICM_12.setter
    def ICM_12(self, ICM_12):
        self.schout.ICM_12 = ICM_12

    @ICM_13.setter
    def ICM_13(self, ICM_13):
        self.schout.ICM_13 = ICM_13

    @ICM_14.setter
    def ICM_14(self, ICM_14):
        self.schout.ICM_14 = ICM_14

    @ICM_15.setter
    def ICM_15(self, ICM_15):
        self.schout.ICM_15 = ICM_15

    @ICM_16.setter
    def ICM_16(self, ICM_16):
        self.schout.ICM_16 = ICM_16

    @ICM_17.setter
    def ICM_17(self, ICM_17):
        self.schout.ICM_17 = ICM_17

    @ICM_18.setter
    def ICM_18(self, ICM_18):
        self.schout.ICM_18 = ICM_18

    @ICM_19.setter
    def ICM_19(self, ICM_19):
        self.schout.ICM_19 = ICM_19

    @ICM_20.setter
    def ICM_20(self, ICM_20):
        self.schout.ICM_20 = ICM_20

    @ICM_21.setter
    def ICM_21(self, ICM_21):
        self.schout.ICM_21 = ICM_21

    @ICM_22.setter
    def ICM_22(self, ICM_22):
        self.schout.ICM_22 = ICM_22

    @ICM_23.setter
    def ICM_23(self, ICM_23):
        self.schout.ICM_23 = ICM_23

    @ICM_24.setter
    def ICM_24(self, ICM_24):
        self.schout.ICM_24 = ICM_24

    @ICM_25.setter
    def ICM_25(self, ICM_25):
        self.schout.ICM_25 = ICM_25

    @COS_1.setter
    def COS_1(self, COS_1):
        self.schout.COS_1 = COS_1

    @COS_2.setter
    def COS_2(self, COS_2):
        self.schout.COS_2 = COS_2

    @COS_3.setter
    def COS_3(self, COS_3):
        self.schout.COS_3 = COS_3

    @COS_4.setter
    def COS_4(self, COS_4):
        self.schout.COS_4 = COS_4

    @COS_5.setter
    def COS_5(self, COS_5):
        self.schout.COS_5 = COS_5

    @COS_6.setter
    def COS_6(self, COS_6):
        self.schout.COS_6 = COS_6

    @COS_7.setter
    def COS_7(self, COS_7):
        self.schout.COS_7 = COS_7

    @COS_8.setter
    def COS_8(self, COS_8):
        self.schout.COS_8 = COS_8

    @COS_9.setter
    def COS_9(self, COS_9):
        self.schout.COS_9 = COS_9

    @COS_10.setter
    def COS_10(self, COS_10):
        self.schout.COS_10 = COS_10

    @COS_11.setter
    def COS_11(self, COS_11):
        self.schout.COS_11 = COS_11

    @COS_12.setter
    def COS_12(self, COS_12):
        self.schout.COS_12 = COS_12

    @COS_13.setter
    def COS_13(self, COS_13):
        self.schout.COS_13 = COS_13

    @FIB_1.setter
    def FIB_1(self, FIB_1):
        self.schout.FIB_1 = FIB_1

    @SED2D_depth_change.setter
    def SED2D_depth_change(self, SED2D_depth_change):
        self.schout.SED2D_depth_change = SED2D_depth_change

    @SED2D_cflsed.setter
    def SED2D_cflsed(self, SED2D_cflsed):
        self.schout.SED2D_cflsed = SED2D_cflsed

    @SED2D_d50.setter
    def SED2D_d50(self, SED2D_d50):
        self.schout.SED2D_d50 = SED2D_d50

    @SED2D_total_transport.setter
    def SED2D_total_transport(self, SED2D_total_transport):
        self.schout.SED2D_total_transport = SED2D_total_transport

    @SED2D_susp_load.setter
    def SED2D_susp_load(self, SED2D_susp_load):
        self.schout.SED2D_susp_load = SED2D_susp_load

    @SED2D_bed_load.setter
    def SED2D_bed_load(self, SED2D_bed_load):
        self.schout.SED2D_bed_load = SED2D_bed_load

    @SED2D_average_transport.setter
    def SED2D_average_transport(self, SED2D_average_transport):
        self.schout.SED2D_average_transport = SED2D_average_transport

    @SED2D_bottom_slope.setter
    def SED2D_bottom_slope(self, SED2D_bottom_slope):
        self.schout.SED2D_bottom_slope = SED2D_bottom_slope

    @z0eq.setter
    def z0eq(self, z0eq):
        self.schout.z0eq = z0eq

    @z0cr2d.setter
    def z0cr2d(self, z0cr):
        self.schout.z0cr2d = z0cr

    @z0sw2d.setter
    def z0sw2d(self, z0sw):
        self.schout.z0sw2d = z0sw

    @z0wr2d.setter
    def z0wr2d(self, z0wr):
        self.schout.z0wr2d = z0wr

    @marsh_flag.setter
    def marsh_flag(self, marsh_flag):
        self.schout.marsh_flag = marsh_flag

    @ICE_velocity.setter
    def ICE_velocity(self, ICE_velocity):
        self.schout.ICE_velocity = ICE_velocity

    @ICE_strain_rate.setter
    def ICE_strain_rate(self, ICE_strain_rate):
        self.schout.ICE_strain_rate = ICE_strain_rate

    @ICE_net_heat_flux.setter
    def ICE_net_heat_flux(self, ICE_net_heat_flux):
        self.schout.ICE_net_heat_flux = ICE_net_heat_flux

    @ICE_fresh_water_flux.setter
    def ICE_fresh_water_flux(self, ICE_fresh_water_flux):
        self.schout.ICE_fresh_water_flux = ICE_fresh_water_flux

    @ICE_top_T.setter
    def ICE_top_T(self, ICE_top_T):
        self.schout.ICE_top_T = ICE_top_T

    @ICE_tracer_1.setter
    def ICE_tracer_1(self, ICE_tracer_1):
        self.schout.ICE_tracer_1 = ICE_tracer_1

    @ICE_tracer_2.setter
    def ICE_tracer_2(self, ICE_tracer_2):
        self.schout.ICE_tracer_2 = ICE_tracer_2

    @ICE_tracer_3.setter
    def ICE_tracer_3(self, ICE_tracer_3):
        self.schout.ICE_tracer_3 = ICE_tracer_3

    @ANA_air_pres_grad_x.setter
    def ANA_air_pres_grad_x(self, ANA_air_pres_grad_x):
        self.schout.ANA_air_pres_grad_x = ANA_air_pres_grad_x

    @ANA_air_pres_grad_y.setter
    def ANA_air_pres_grad_y(self, ANA_air_pres_grad_y):
        self.schout.ANA_air_pres_grad_y = ANA_air_pres_grad_y

    @ANA_tide_pot_grad_x.setter
    def ANA_tide_pot_grad_x(self, ANA_tide_pot_grad_x):
        self.schout.ANA_tide_pot_grad_x = ANA_tide_pot_grad_x

    @ANA_tide_pot_grad_y.setter
    def ANA_tide_pot_grad_y(self, ANA_tide_pot_grad_y):
        self.schout.ANA_tide_pot_grad_y = ANA_tide_pot_grad_y

    @ANA_hor_viscosity_x.setter
    def ANA_hor_viscosity_x(self, ANA_hor_viscosity_x):
        self.schout.ANA_hor_viscosity_x = ANA_hor_viscosity_x

    @ANA_hor_viscosity_y.setter
    def ANA_hor_viscosity_y(self, ANA_hor_viscosity_y):
        self.schout.ANA_hor_viscosity_y = ANA_hor_viscosity_y

    @ANA_bclinic_force_x.setter
    def ANA_bclinic_force_x(self, ANA_bclinic_force_x):
        self.schout.ANA_bclinic_force_x = ANA_bclinic_force_x

    @ANA_bclinic_force_y.setter
    def ANA_bclinic_force_y(self, ANA_bclinic_force_y):
        self.schout.ANA_bclinic_force_y = ANA_bclinic_force_y

    @ANA_vert_viscosity_x.setter
    def ANA_vert_viscosity_x(self, ANA_vert_viscosity_x):
        self.schout.ANA_vert_viscosity_x = ANA_vert_viscosity_x

    @ANA_vert_viscosity_y.setter
    def ANA_vert_viscosity_y(self, ANA_vert_viscosity_y):
        self.schout.ANA_vert_viscosity_y = ANA_vert_viscosity_y

    @ANA_mom_advection_x.setter
    def ANA_mom_advection_x(self, ANA_mom_advection_x):
        self.schout.ANA_mom_advection_x = ANA_mom_advection_x

    @ANA_mom_advection_y.setter
    def ANA_mom_advection_y(self, ANA_mom_advection_y):
        self.schout.ANA_mom_advection_y = ANA_mom_advection_y

    @ANA_Richardson.setter
    def ANA_Richardson(self, ANA_Richardson):
        self.schout.ANA_Richardson = ANA_Richardson

    @ANA_transport_min_dt_elem.setter
    def ANA_transport_min_dt_elem(self, ANA_transport_min_dt_elem):
        self.schout.ANA_transport_min_dt_elem = ANA_transport_min_dt_elem
