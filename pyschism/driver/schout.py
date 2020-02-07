

class SCHOUT:
    """ SCHOUT group """

    def __init__(self, nml):
        self._nml = nml

    @property
    def nhot(self):
        return self._nml['schout']['nhot']

    @property
    def nhot_write(self):
        return self._nml['schout']['nhot_write']

    @property
    def iout_sta(self):
        return self._nml['schout']['iout_sta']

    @property
    def nspool_sta(self):
        return self._nml['schout']['nspool_sta']

    @property
    def iof_hydro(self):
        return self._nml['schout']['iof_hydro']

    @property
    def iof_wwm(self):
        return self._nml['schout']['iof_wwm']

    @property
    def iof_gen(self):
        return self._nml['schout']['iof_gen']

    @property
    def iof_age(self):
        return self._nml['schout']['iof_age']

    @property
    def iof_sed(self):
        return self._nml['schout']['iof_sed']

    @property
    def iof_eco(self):
        return self._nml['schout']['iof_eco']

    @property
    def iof_icm(self):
        return self._nml['schout']['iof_icm']

    @property
    def iof_cos(self):
        return self._nml['schout']['iof_cos']

    @property
    def iof_fib(self):
        return self._nml['schout']['iof_fib']

    @property
    def iof_sed2d(self):
        return self._nml['schout']['iof_sed2d']

    @property
    def iof_marsh(self):
        return self._nml['schout']['iof_marsh']

    @property
    def iof_ice(self):
        return self._nml['schout']['iof_ice']

    @property
    def iof_ana(self):
        return self._nml['schout']['iof_ana']

    @property
    def elev(self):
        return bool(self._nml['schout']['iof_hydro'][0])

    @property
    def air_pressure(self):
        return bool(self._nml['schout']['iof_hydro'][1])

    @property
    def air_temperature(self):
        return bool(self._nml['schout']['iof_hydro'][2])

    @property
    def specific_humidity(self):
        return bool(self._nml['schout']['iof_hydro'][3])

    @property
    def solar_radiation(self):
        return bool(self._nml['schout']['iof_hydro'][4])

    @property
    def sensible_flux(self):
        return bool(self._nml['schout']['iof_hydro'][5])

    @property
    def latent_heat(self):
        return bool(self._nml['schout']['iof_hydro'][6])

    @property
    def upward_longwave(self):
        return bool(self._nml['schout']['iof_hydro'][7])

    @property
    def downward_longwave(self):
        return bool(self._nml['schout']['iof_hydro'][8])

    @property
    def total_heat_flux(self):
        return bool(self._nml['schout']['iof_hydro'][9])

    @property
    def evaporation(self):
        return bool(self._nml['schout']['iof_hydro'][10])

    @property
    def precipitation(self):
        return bool(self._nml['schout']['iof_hydro'][11])

    @property
    def bottom_stress(self):
        return bool(self._nml['schout']['iof_hydro'][12])

    @property
    def wind_speed(self):
        return bool(self._nml['schout']['iof_hydro'][13])

    @property
    def wind_stress(self):
        return bool(self._nml['schout']['iof_hydro'][14])

    @property
    def dahv(self):
        return bool(self._nml['schout']['iof_hydro'][15])

    @property
    def vertical_velocity(self):
        return bool(self._nml['schout']['iof_hydro'][16])

    @property
    def temp(self):
        return bool(self._nml['schout']['iof_hydro'][17])

    @property
    def salt(self):
        return bool(self._nml['schout']['iof_hydro'][18])

    @property
    def water_density(self):
        return bool(self._nml['schout']['iof_hydro'][19])

    @property
    def diffusivity(self):
        return bool(self._nml['schout']['iof_hydro'][20])

    @property
    def viscosity(self):
        return bool(self._nml['schout']['iof_hydro'][21])

    @property
    def TKE(self):
        return bool(self._nml['schout']['iof_hydro'][22])

    @property
    def mixing_length(self):
        return bool(self._nml['schout']['iof_hydro'][23])

    @property
    def hvel(self):
        return bool(self._nml['schout']['iof_hydro'][24])

    @property
    def hvel_side(self):
        return bool(self._nml['schout']['iof_hydro'][25])

    @property
    def wvel_elem(self):
        return bool(self._nml['schout']['iof_hydro'][26])

    @property
    def temp_elem(self):
        return bool(self._nml['schout']['iof_hydro'][27])

    @property
    def salt_elem(self):
        return bool(self._nml['schout']['iof_hydro'][28])

    @property
    def pressure_gradient(self):
        return bool(self._nml['schout']['iof_hydro'][29])

    @property
    def WWM_1(self):
        return bool(self._nml['schout']['iof_wwm'][0])

    @property
    def WWM_2(self):
        return bool(self._nml['schout']['iof_wwm'][1])

    @property
    def WWM_3(self):
        return bool(self._nml['schout']['iof_wwm'][2])

    @property
    def WWM_4(self):
        return bool(self._nml['schout']['iof_wwm'][3])

    @property
    def WWM_5(self):
        return bool(self._nml['schout']['iof_wwm'][4])

    @property
    def WWM_6(self):
        return bool(self._nml['schout']['iof_wwm'][5])

    @property
    def WWM_9(self):
        return bool(self._nml['schout']['iof_wwm'][6])

    @property
    def WWM_10(self):
        return bool(self._nml['schout']['iof_wwm'][7])

    @property
    def WWM_11(self):
        return bool(self._nml['schout']['iof_wwm'][8])

    @property
    def WWM_12(self):
        return bool(self._nml['schout']['iof_wwm'][9])

    @property
    def WWM_13(self):
        return bool(self._nml['schout']['iof_wwm'][10])

    @property
    def WWM_14(self):
        return bool(self._nml['schout']['iof_wwm'][11])

    @property
    def WWM_15(self):
        return bool(self._nml['schout']['iof_wwm'][12])

    @property
    def WWM_16(self):
        return bool(self._nml['schout']['iof_wwm'][13])

    @property
    def WWM_17(self):
        return bool(self._nml['schout']['iof_wwm'][14])

    @property
    def WWM_18(self):
        return bool(self._nml['schout']['iof_wwm'][15])

    @property
    def WWM_19(self):
        return bool(self._nml['schout']['iof_wwm'][16])

    @property
    def WWM_20(self):
        return bool(self._nml['schout']['iof_wwm'][17])

    @property
    def WWM_21(self):
        return bool(self._nml['schout']['iof_wwm'][18])

    @property
    def WWM_22(self):
        return bool(self._nml['schout']['iof_wwm'][19])

    @property
    def WWM_23(self):
        return bool(self._nml['schout']['iof_wwm'][20])

    @property
    def WWM_24(self):
        return bool(self._nml['schout']['iof_wwm'][21])

    @property
    def WWM_25(self):
        return bool(self._nml['schout']['iof_wwm'][22])

    @property
    def WWM_26(self):
        return bool(self._nml['schout']['iof_wwm'][23])

    @property
    def WWM_27(self):
        return bool(self._nml['schout']['iof_wwm'][24])

    @property
    def WWM_28(self):
        return bool(self._nml['schout']['iof_wwm'][25])

    @property
    def WWM_energy_dir(self):
        return bool(self._nml['schout']['iof_wwm'][26])

    @property
    def wave_force(self):
        return bool(self._nml['schout']['iof_wwm'][27])

    @property
    def GEN_1(self):
        return bool(self._nml['schout']['iof_gen'][0])

    @property
    def GEN_2(self):
        return bool(self._nml['schout']['iof_gen'][1])

    @property
    def AGE_1(self):
        return bool(self._nml['schout']['iof_age'][0])

    @property
    def AGE_2(self):
        return bool(self._nml['schout']['iof_age'][1])

    @property
    def SED_depth_change(self):
        return bool(self._nml['schout']['iof_sed'][0])

    @property
    def SED_D50(self):
        return bool(self._nml['schout']['iof_sed'][1])

    @property
    def SED_bed_stress(self):
        return bool(self._nml['schout']['iof_sed'][2])

    @property
    def SED_bed_roughness(self):
        return bool(self._nml['schout']['iof_sed'][3])

    @property
    def SED_TSC(self):
        return bool(self._nml['schout']['iof_sed'][4])

    @property
    def bed_thickness(self):
        return bool(self._nml['schout']['iof_sed'][5])

    @property
    def bed_age(self):
        return bool(self._nml['schout']['iof_sed'][6])

    @property
    def z0st(self):
        return bool(self._nml['schout']['iof_sed'][7])

    @property
    def z0cr(self):
        return bool(self._nml['schout']['iof_sed'][8])

    @property
    def z0sw(self):
        return bool(self._nml['schout']['iof_sed'][9])

    @property
    def z0wr(self):
        return bool(self._nml['schout']['iof_sed'][10])

    @property
    def SED3D_1(self):
        return bool(self._nml['schout']['iof_sed'][11])

    @property
    def SED_bdld_1(self):
        return bool(self._nml['schout']['iof_sed'][12])

    @property
    def SED_bedfrac_1(self):
        return bool(self._nml['schout']['iof_sed'][13])

    @property
    def SED3D_2(self):
        return bool(self._nml['schout']['iof_sed'][14])

    @property
    def SED_bdld_2(self):
        return bool(self._nml['schout']['iof_sed'][15])

    @property
    def SED_bedfrac_3(self):
        return bool(self._nml['schout']['iof_sed'][16])

    @property
    def ECO_1(self):
        return bool(self._nml['schout']['iof_eco'][0])

    @property
    def ICM_Chl(self):
        return bool(self._nml['schout']['iof_icm'][0])

    @property
    def ICM_pH(self):
        return bool(self._nml['schout']['iof_icm'][1])

    @property
    def ICM_PrmPrdt(self):
        return bool(self._nml['schout']['iof_icm'][2])

    @property
    def ICM_DIN(self):
        return bool(self._nml['schout']['iof_icm'][3])

    @property
    def ICM_PON(self):
        return bool(self._nml['schout']['iof_icm'][4])

    @property
    def ICM_SED_BENDOC(self):
        return bool(self._nml['schout']['iof_icm'][5])

    @property
    def ICM_SED_BENNH4(self):
        return bool(self._nml['schout']['iof_icm'][6])

    @property
    def ICM_SED_BENNO3(self):
        return bool(self._nml['schout']['iof_icm'][7])

    @property
    def ICM_SED_BENPO4(self):
        return bool(self._nml['schout']['iof_icm'][8])

    @property
    def ICM_SED_BENCOD(self):
        return bool(self._nml['schout']['iof_icm'][9])

    @property
    def ICM_SED_BENDO(self):
        return bool(self._nml['schout']['iof_icm'][10])

    @property
    def ICM_SED_BENSA(self):
        return bool(self._nml['schout']['iof_icm'][11])

    @property
    def ICM_lfsav(self):
        return bool(self._nml['schout']['iof_icm'][12])

    @property
    def ICM_stsav(self):
        return bool(self._nml['schout']['iof_icm'][13])

    @property
    def ICM_rtsav(self):
        return bool(self._nml['schout']['iof_icm'][14])

    @property
    def ICM_tlfsav(self):
        return bool(self._nml['schout']['iof_icm'][15])

    @property
    def ICM_tstsav(self):
        return bool(self._nml['schout']['iof_icm'][16])

    @property
    def ICM_trtsav(self):
        return bool(self._nml['schout']['iof_icm'][17])

    @property
    def ICM_hcansav(self):
        return bool(self._nml['schout']['iof_icm'][18])

    @property
    def ICM_CNH4(self):
        return bool(self._nml['schout']['iof_icm'][19])

    @property
    def ICM_CNH3(self):
        return bool(self._nml['schout']['iof_icm'][20])

    @property
    def ICM_CPIP(self):
        return bool(self._nml['schout']['iof_icm'][21])

    @property
    def ICM_CPOS(self):
        return bool(self._nml['schout']['iof_icm'][22])

    @property
    def ICM_CCH4(self):
        return bool(self._nml['schout']['iof_icm'][23])

    @property
    def ICM_CSO4(self):
        return bool(self._nml['schout']['iof_icm'][24])

    @property
    def ICM_CH2S(self):
        return bool(self._nml['schout']['iof_icm'][25])

    @property
    def ICM_SEDPON1(self):
        return bool(self._nml['schout']['iof_icm'][26])

    @property
    def ICM_SEDPON2(self):
        return bool(self._nml['schout']['iof_icm'][27])

    @property
    def ICM_SEDPON3(self):
        return bool(self._nml['schout']['iof_icm'][28])

    @property
    def ICM_SEDPOP1(self):
        return bool(self._nml['schout']['iof_icm'][29])

    @property
    def ICM_SEDPOP2(self):
        return bool(self._nml['schout']['iof_icm'][30])

    @property
    def ICM_SEDPOP3(self):
        return bool(self._nml['schout']['iof_icm'][31])

    @property
    def ICM_SEDPOC1(self):
        return bool(self._nml['schout']['iof_icm'][32])

    @property
    def ICM_SEDPOC2(self):
        return bool(self._nml['schout']['iof_icm'][33])

    @property
    def ICM_SEDPOC3(self):
        return bool(self._nml['schout']['iof_icm'][34])

    @property
    def ICM_EROH2S(self):
        return bool(self._nml['schout']['iof_icm'][35])

    @property
    def ICM_EROLPOC(self):
        return bool(self._nml['schout']['iof_icm'][36])

    @property
    def ICM_ERORPOC(self):
        return bool(self._nml['schout']['iof_icm'][37])

    @property
    def ICM_DO_consumption(self):
        return bool(self._nml['schout']['iof_icm'][38])

    @property
    def ICM_GP1(self):
        return bool(self._nml['schout']['iof_icm'][39])

    @property
    def ICM_GP2(self):
        return bool(self._nml['schout']['iof_icm'][40])

    @property
    def ICM_GP3(self):
        return bool(self._nml['schout']['iof_icm'][41])

    @property
    def ICM_1(self):
        return bool(self._nml['schout']['iof_icm'][42])

    @property
    def ICM_2(self):
        return bool(self._nml['schout']['iof_icm'][43])

    @property
    def ICM_3(self):
        return bool(self._nml['schout']['iof_icm'][44])

    @property
    def ICM_4(self):
        return bool(self._nml['schout']['iof_icm'][45])

    @property
    def ICM_5(self):
        return bool(self._nml['schout']['iof_icm'][46])

    @property
    def ICM_6(self):
        return bool(self._nml['schout']['iof_icm'][47])

    @property
    def ICM_7(self):
        return bool(self._nml['schout']['iof_icm'][48])

    @property
    def ICM_8(self):
        return bool(self._nml['schout']['iof_icm'][49])

    @property
    def ICM_9(self):
        return bool(self._nml['schout']['iof_icm'][50])

    @property
    def ICM_10(self):
        return bool(self._nml['schout']['iof_icm'][51])

    @property
    def ICM_11(self):
        return bool(self._nml['schout']['iof_icm'][52])

    @property
    def ICM_12(self):
        return bool(self._nml['schout']['iof_icm'][53])

    @property
    def ICM_13(self):
        return bool(self._nml['schout']['iof_icm'][54])

    @property
    def ICM_14(self):
        return bool(self._nml['schout']['iof_icm'][55])

    @property
    def ICM_15(self):
        return bool(self._nml['schout']['iof_icm'][56])

    @property
    def ICM_16(self):
        return bool(self._nml['schout']['iof_icm'][57])

    @property
    def ICM_17(self):
        return bool(self._nml['schout']['iof_icm'][58])

    @property
    def ICM_18(self):
        return bool(self._nml['schout']['iof_icm'][59])

    @property
    def ICM_19(self):
        return bool(self._nml['schout']['iof_icm'][60])

    @property
    def ICM_20(self):
        return bool(self._nml['schout']['iof_icm'][61])

    @property
    def ICM_21(self):
        return bool(self._nml['schout']['iof_icm'][62])

    @property
    def ICM_22(self):
        return bool(self._nml['schout']['iof_icm'][63])

    @property
    def ICM_23(self):
        return bool(self._nml['schout']['iof_icm'][64])

    @property
    def ICM_24(self):
        return bool(self._nml['schout']['iof_icm'][65])

    @property
    def ICM_25(self):
        return bool(self._nml['schout']['iof_icm'][66])

    @property
    def COS_1(self):
        return bool(self._nml['schout']['iof_cos'][0])

    @property
    def COS_2(self):
        return bool(self._nml['schout']['iof_cos'][1])

    @property
    def COS_3(self):
        return bool(self._nml['schout']['iof_cos'][2])

    @property
    def COS_4(self):
        return bool(self._nml['schout']['iof_cos'][3])

    @property
    def COS_5(self):
        return bool(self._nml['schout']['iof_cos'][4])

    @property
    def COS_6(self):
        return bool(self._nml['schout']['iof_cos'][5])

    @property
    def COS_7(self):
        return bool(self._nml['schout']['iof_cos'][6])

    @property
    def COS_8(self):
        return bool(self._nml['schout']['iof_cos'][7])

    @property
    def COS_9(self):
        return bool(self._nml['schout']['iof_cos'][8])

    @property
    def COS_10(self):
        return bool(self._nml['schout']['iof_cos'][9])

    @property
    def COS_11(self):
        return bool(self._nml['schout']['iof_cos'][10])

    @property
    def COS_12(self):
        return bool(self._nml['schout']['iof_cos'][11])

    @property
    def COS_13(self):
        return bool(self._nml['schout']['iof_cos'][12])

    @property
    def FIB_1(self):
        return bool(self._nml['schout']['iof_fib'][0])

    @property
    def SED2D_depth_change(self):
        return bool(self._nml['schout']['iof_sed2d'][1])

    @property
    def SED2D_cflsed(self):
        return bool(self._nml['schout']['iof_sed2d'][2])

    @property
    def SED2D_d50(self):
        return bool(self._nml['schout']['iof_sed2d'][3])

    @property
    def SED2D_total_transport(self):
        return bool(self._nml['schout']['iof_sed2d'][4])

    @property
    def SED2D_susp_load(self):
        return bool(self._nml['schout']['iof_sed2d'][5])

    @property
    def SED2D_bed_load(self):
        return bool(self._nml['schout']['iof_sed2d'][6])

    @property
    def SED2D_average_transport(self):
        return bool(self._nml['schout']['iof_sed2d'][7])

    @property
    def SED2D_bottom_slope(self):
        return bool(self._nml['schout']['iof_sed2d'][8])

    @property
    def z0eq(self):
        return bool(self._nml['schout']['iof_sed2d'][9])

    @property
    def z0cr2d(self):
        return bool(self._nml['schout']['iof_sed2d'][10])

    @property
    def z0sw2d(self):
        return bool(self._nml['schout']['iof_sed2d'][11])

    @property
    def z0wr2d(self):
        return bool(self._nml['schout']['iof_sed2d'][12])

    @property
    def marsh_flag(self):
        return bool(self._nml['schout']['iof_marsh'][0])

    @property
    def ICE_velocity(self):
        return bool(self._nml['schout']['iof_ice'][0])

    @property
    def ICE_strain_rate(self):
        return bool(self._nml['schout']['iof_ice'][1])

    @property
    def ICE_net_heat_flux(self):
        return bool(self._nml['schout']['iof_ice'][2])

    @property
    def ICE_fresh_water_flux(self):
        return bool(self._nml['schout']['iof_ice'][3])

    @property
    def ICE_top_T(self):
        return bool(self._nml['schout']['iof_ice'][4])

    @property
    def ICE_tracer_1(self):
        return bool(self._nml['schout']['iof_ice'][5])

    @property
    def ICE_tracer_2(self):
        return bool(self._nml['schout']['iof_ice'][6])

    @property
    def ICE_tracer_3(self):
        return bool(self._nml['schout']['iof_ice'][7])

    @property
    def ANA_air_pres_grad_x(self):
        return bool(self._nml['schout']['iof_ana'][0])

    @property
    def ANA_air_pres_grad_y(self):
        return bool(self._nml['schout']['iof_ana'][1])

    @property
    def ANA_tide_pot_grad_x(self):
        return bool(self._nml['schout']['iof_ana'][2])

    @property
    def ANA_tide_pot_grad_y(self):
        return bool(self._nml['schout']['iof_ana'][3])

    @property
    def ANA_hor_viscosity_x(self):
        return bool(self._nml['schout']['iof_ana'][4])

    @property
    def ANA_hor_viscosity_y(self):
        return bool(self._nml['schout']['iof_ana'][5])

    @property
    def ANA_bclinic_force_x(self):
        return bool(self._nml['schout']['iof_ana'][6])

    @property
    def ANA_bclinic_force_y(self):
        return bool(self._nml['schout']['iof_ana'][7])

    @property
    def ANA_vert_viscosity_x(self):
        return bool(self._nml['schout']['iof_ana'][8])

    @property
    def ANA_vert_viscosity_y(self):
        return bool(self._nml['schout']['iof_ana'][9])

    @property
    def ANA_mom_advection_x(self):
        return bool(self._nml['schout']['iof_ana'][10])

    @property
    def ANA_mom_advection_y(self):
        return bool(self._nml['schout']['iof_ana'][11])

    @property
    def ANA_Richardson(self):
        return bool(self._nml['schout']['iof_ana'][12])

    @property
    def ANA_transport_min_dt_elem(self):
        return bool(self._nml['schout']['iof_ana'][13])

    @nhot.setter
    def nhot(self, nhot):
        self._nml['schout']['nhot'] = nhot

    @nhot_write.setter
    def nhot_write(self, nhot_write):
        self._nml['schout']['nhot_write'] = nhot_write

    @iout_sta.setter
    def iout_sta(self, iout_sta):
        self._nml['schout']['iout_sta'] = iout_sta

    @nspool_sta.setter
    def nspool_sta(self, nspool_sta):
        self._nml['schout']['nspool_sta'] = nspool_sta

    @iof_hydro.setter
    def iof_hydro(self, iof_hydro):
        self._nml['schout']['iof_hydro'] = iof_hydro

    @iof_wwm.setter
    def iof_wwm(self, iof_wwm):
        self._nml['schout']['iof_wwm'] = iof_wwm

    @iof_gen.setter
    def iof_gen(self, iof_gen):
        self._nml['schout']['iof_gen'] = iof_gen

    @iof_age.setter
    def iof_age(self, iof_age):
        self._nml['schout']['iof_age'] = iof_age

    @iof_sed.setter
    def iof_sed(self, iof_sed):
        self._nml['schout']['iof_sed'] = iof_sed

    @iof_eco.setter
    def iof_eco(self, iof_eco):
        self._nml['schout']['iof_eco'] = iof_eco

    @iof_icm.setter
    def iof_icm(self, iof_icm):
        self._nml['schout']['iof_icm'] = iof_icm

    @iof_cos.setter
    def iof_cos(self, iof_cos):
        self._nml['schout']['iof_cos'] = iof_cos

    @iof_fib.setter
    def iof_fib(self, iof_fib):
        self._nml['schout']['iof_fib'] = iof_fib

    @iof_sed2d.setter
    def iof_sed2d(self, iof_sed2d):
        self._nml['schout']['iof_sed2d'] = iof_sed2d

    @iof_marsh.setter
    def iof_marsh(self, iof_marsh):
        self._nml['schout']['iof_marsh'] = iof_marsh

    @iof_ice.setter
    def iof_ice(self, iof_ice):
        self._nml['schout']['iof_ice'] = iof_ice

    @iof_ana.setter
    def iof_ana(self, iof_ana):
        self._nml['schout']['iof_ana'] = iof_ana

    @elev.setter
    def elev(self, elev):
        self._nml['schout']['iof_hydro'][0] = int(bool(elev))

    @air_pressure.setter
    def air_pressure(self, air_pressure):
        self._nml['schout']['iof_hydro'][1] = int(bool(air_pressure))

    @air_temperature.setter
    def air_temperature(self, air_temperature):
        self._nml['schout']['iof_hydro'][2] = int(bool(air_temperature))

    @specific_humidity.setter
    def specific_humidity(self, specific_humidity):
        self._nml['schout']['iof_hydro'][3] = int(bool(specific_humidity))

    @solar_radiation.setter
    def solar_radiation(self, solar_radiation):
        self._nml['schout']['iof_hydro'][4] = int(bool(solar_radiation))

    @sensible_flux.setter
    def sensible_flux(self, sensible_flux):
        self._nml['schout']['iof_hydro'][5] = int(bool(sensible_flux))

    @latent_heat.setter
    def latent_heat(self, latent_heat):
        self._nml['schout']['iof_hydro'][6] = int(bool(latent_heat))

    @upward_longwave.setter
    def upward_longwave(self, upward_longwave):
        self._nml['schout']['iof_hydro'][7] = int(bool(upward_longwave))

    @downward_longwave.setter
    def downward_longwave(self, downward_longwave):
        self._nml['schout']['iof_hydro'][8] = int(bool(downward_longwave))

    @total_heat_flux.setter
    def total_heat_flux(self, total_heat_flux):
        self._nml['schout']['iof_hydro'][9] = int(bool(total_heat_flux))

    @evaporation.setter
    def evaporation(self, evaporation):
        self._nml['schout']['iof_hydro'][10] = int(bool(evaporation))

    @precipitation.setter
    def precipitation(self, precipitation):
        self._nml['schout']['iof_hydro'][11] = int(bool(precipitation))

    @bottom_stress.setter
    def bottom_stress(self, bottom_stress):
        self._nml['schout']['iof_hydro'][12] = int(bool(bottom_stress))

    @wind_speed.setter
    def wind_speed(self, wind_speed):
        self._nml['schout']['iof_hydro'][13] = int(bool(wind_speed))

    @wind_stress.setter
    def wind_stress(self, wind_stress):
        self._nml['schout']['iof_hydro'][14] = int(bool(wind_stress))

    @dahv.setter
    def dahv(self, dahv):
        self._nml['schout']['iof_hydro'][15] = int(bool(dahv))

    @vertical_velocity.setter
    def vertical_velocity(self, vertical_velocity):
        self._nml['schout']['iof_hydro'][16] = int(bool(vertical_velocity))

    @temp.setter
    def temp(self, temp):
        self._nml['schout']['iof_hydro'][17] = int(bool(temp))

    @salt.setter
    def salt(self, salt):
        self._nml['schout']['iof_hydro'][18] = int(bool(salt))

    @water_density.setter
    def water_density(self, water_density):
        self._nml['schout']['iof_hydro'][19] = int(bool(water_density))

    @diffusivity.setter
    def diffusivity(self, diffusivity):
        self._nml['schout']['iof_hydro'][20] = int(bool(diffusivity))

    @viscosity.setter
    def viscosity(self, viscosity):
        self._nml['schout']['iof_hydro'][21] = int(bool(viscosity))

    @TKE.setter
    def TKE(self, TKE):
        self._nml['schout']['iof_hydro'][22] = int(bool(TKE))

    @mixing_length.setter
    def mixing_length(self, mixing_length):
        self._nml['schout']['iof_hydro'][23] = int(bool(mixing_length))

    @hvel.setter
    def hvel(self, hvel):
        self._nml['schout']['iof_hydro'][24] = int(bool(hvel))

    @hvel_side.setter
    def hvel_side(self, hvel_side):
        self._nml['schout']['iof_hydro'][25] = int(bool(hvel_side))

    @wvel_elem.setter
    def wvel_elem(self, wvel_elem):
        self._nml['schout']['iof_hydro'][26] = int(bool(wvel_elem))

    @temp_elem.setter
    def temp_elem(self, temp_elem):
        self._nml['schout']['iof_hydro'][27] = int(bool(temp_elem))

    @salt_elem.setter
    def salt_elem(self, salt_elem):
        self._nml['schout']['iof_hydro'][28] = int(bool(salt_elem))

    @pressure_gradient.setter
    def pressure_gradient(self, pressure_gradient):
        self._nml['schout']['iof_hydro'][29] = int(bool(pressure_gradient))

    @WWM_1.setter
    def WWM_1(self, WWM_1):
        self._nml['schout']['iof_wwm'][0] = int(bool(WWM_1))

    @WWM_2.setter
    def WWM_2(self, WWM_2):
        self._nml['schout']['iof_wwm'][1] = int(bool(WWM_2))

    @WWM_3.setter
    def WWM_3(self, WWM_3):
        self._nml['schout']['iof_wwm'][2] = int(bool(WWM_3))

    @WWM_4.setter
    def WWM_4(self, WWM_4):
        self._nml['schout']['iof_wwm'][3] = int(bool(WWM_4))

    @WWM_5.setter
    def WWM_5(self, WWM_5):
        self._nml['schout']['iof_wwm'][4] = int(bool(WWM_5))

    @WWM_6.setter
    def WWM_6(self, WWM_6):
        self._nml['schout']['iof_wwm'][5] = int(bool(WWM_6))

    @WWM_9.setter
    def WWM_9(self, WWM_9):
        self._nml['schout']['iof_wwm'][6] = int(bool(WWM_9))

    @WWM_10.setter
    def WWM_10(self, WWM_10):
        self._nml['schout']['iof_wwm'][7] = int(bool(WWM_10))

    @WWM_11.setter
    def WWM_11(self, WWM_11):
        self._nml['schout']['iof_wwm'][8] = int(bool(WWM_11))

    @WWM_12.setter
    def WWM_12(self, WWM_12):
        self._nml['schout']['iof_wwm'][9] = int(bool(WWM_12))

    @WWM_13.setter
    def WWM_13(self, WWM_13):
        self._nml['schout']['iof_wwm'][10] = int(bool(WWM_13))

    @WWM_14.setter
    def WWM_14(self, WWM_14):
        self._nml['schout']['iof_wwm'][11] = int(bool(WWM_14))

    @WWM_15.setter
    def WWM_15(self, WWM_15):
        self._nml['schout']['iof_wwm'][12] = int(bool(WWM_15))

    @WWM_16.setter
    def WWM_16(self, WWM_16):
        self._nml['schout']['iof_wwm'][13] = int(bool(WWM_16))

    @WWM_17.setter
    def WWM_17(self, WWM_17):
        self._nml['schout']['iof_wwm'][14] = int(bool(WWM_17))

    @WWM_18.setter
    def WWM_18(self, WWM_18):
        self._nml['schout']['iof_wwm'][15] = int(bool(WWM_18))

    @WWM_19.setter
    def WWM_19(self, WWM_19):
        self._nml['schout']['iof_wwm'][16] = int(bool(WWM_19))

    @WWM_20.setter
    def WWM_20(self, WWM_20):
        self._nml['schout']['iof_wwm'][17] = int(bool(WWM_20))

    @WWM_21.setter
    def WWM_21(self, WWM_21):
        self._nml['schout']['iof_wwm'][18] = int(bool(WWM_21))

    @WWM_22.setter
    def WWM_22(self, WWM_22):
        self._nml['schout']['iof_wwm'][19] = int(bool(WWM_22))

    @WWM_23.setter
    def WWM_23(self, WWM_23):
        self._nml['schout']['iof_wwm'][20] = int(bool(WWM_23))

    @WWM_24.setter
    def WWM_24(self, WWM_24):
        self._nml['schout']['iof_wwm'][21] = int(bool(WWM_24))

    @WWM_25.setter
    def WWM_25(self, WWM_25):
        self._nml['schout']['iof_wwm'][22] = int(bool(WWM_25))

    @WWM_26.setter
    def WWM_26(self, WWM_26):
        self._nml['schout']['iof_wwm'][23] = int(bool(WWM_26))

    @WWM_27.setter
    def WWM_27(self, WWM_27):
        self._nml['schout']['iof_wwm'][24] = int(bool(WWM_27))

    @WWM_28.setter
    def WWM_28(self, WWM_28):
        self._nml['schout']['iof_wwm'][25] = int(bool(WWM_28))

    @WWM_energy_dir.setter
    def WWM_energy_dir(self, WWM_energy_dir):
        self._nml['schout']['iof_wwm'][26] = int(bool(WWM_energy_dir))

    @wave_force.setter
    def wave_force(self, wave_force):
        self._nml['schout']['iof_wwm'][27] = int(bool(wave_force))

    @GEN_1.setter
    def GEN_1(self, GEN_1):
        self._nml['schout']['iof_gen'][0] = int(bool(GEN_1))

    @GEN_2.setter
    def GEN_2(self, GEN_2):
        self._nml['schout']['iof_gen'][1] = int(bool(GEN_2))

    @AGE_1.setter
    def AGE_1(self, AGE_1):
        self._nml['schout']['iof_age'][0] = int(bool(AGE_1))

    @AGE_2.setter
    def AGE_2(self, AGE_2):
        self._nml['schout']['iof_age'][1] = int(bool(AGE_2))

    @SED_depth_change.setter
    def SED_depth_change(self, SED_depth_change):
        self._nml['schout']['iof_sed'][0] = int(bool(SED_depth_change))

    @SED_D50.setter
    def SED_D50(self, SED_D50):
        self._nml['schout']['iof_sed'][1] = int(bool(SED_D50))

    @SED_bed_stress.setter
    def SED_bed_stress(self, SED_bed_stress):
        self._nml['schout']['iof_sed'][2] = int(bool(SED_bed_stress))

    @SED_bed_roughness.setter
    def SED_bed_roughness(self, SED_bed_roughness):
        self._nml['schout']['iof_sed'][3] = int(bool(SED_bed_roughness))

    @SED_TSC.setter
    def SED_TSC(self, SED_TSC):
        self._nml['schout']['iof_sed'][4] = int(bool(SED_TSC))

    @bed_thickness.setter
    def bed_thickness(self, bed_thickness):
        self._nml['schout']['iof_sed'][5] = int(bool(bed_thickness))

    @bed_age.setter
    def bed_age(self, bed_age):
        self._nml['schout']['iof_sed'][6] = int(bool(bed_age))

    @z0st.setter
    def z0st(self, z0st):
        self._nml['schout']['iof_sed'][7] = int(bool(z0st))

    @z0cr.setter
    def z0cr(self, z0cr):
        self._nml['schout']['iof_sed'][8] = int(bool(z0cr))

    @z0sw.setter
    def z0sw(self, z0sw):
        self._nml['schout']['iof_sed'][9] = int(bool(z0sw))

    @z0wr.setter
    def z0wr(self, z0wr):
        self._nml['schout']['iof_sed'][10] = int(bool(z0wr))

    @SED3D_1.setter
    def SED3D_1(self, SED3D_1):
        self._nml['schout']['iof_sed'][11] = int(bool(SED3D_1))

    @SED_bdld_1.setter
    def SED_bdld_1(self, SED_bdld_1):
        self._nml['schout']['iof_sed'][12] = int(bool(SED_bdld_1))

    @SED_bedfrac_1.setter
    def SED_bedfrac_1(self, SED_bedfrac_1):
        self._nml['schout']['iof_sed'][13] = int(bool(SED_bedfrac_1))

    @SED3D_2.setter
    def SED3D_2(self, SED3D_2):
        self._nml['schout']['iof_sed'][14] = int(bool(SED3D_2))

    @SED_bdld_2.setter
    def SED_bdld_2(self, SED_bdld_2):
        self._nml['schout']['iof_sed'][15] = int(bool(SED_bdld_2))

    @SED_bedfrac_3.setter
    def SED_bedfrac_3(self, SED_bedfrac_3):
        self._nml['schout']['iof_sed'][16] = int(bool(SED_bedfrac_3))

    @ECO_1.setter
    def ECO_1(self, ECO_1):
        self._nml['schout']['iof_eco'][0] = int(bool(ECO_1))

    @ICM_Chl.setter
    def ICM_Chl(self, ICM_Chl):
        self._nml['schout']['iof_icm'][0] = int(bool(ICM_Chl))

    @ICM_pH.setter
    def ICM_pH(self, ICM_pH):
        self._nml['schout']['iof_icm'][1] = int(bool(ICM_pH))

    @ICM_PrmPrdt.setter
    def ICM_PrmPrdt(self, ICM_PrmPrdt):
        self._nml['schout']['iof_icm'][2] = int(bool(ICM_PrmPrdt))

    @ICM_DIN.setter
    def ICM_DIN(self, ICM_DIN):
        self._nml['schout']['iof_icm'][3] = int(bool(ICM_DIN))

    @ICM_PON.setter
    def ICM_PON(self, ICM_PON):
        self._nml['schout']['iof_icm'][4] = int(bool(ICM_PON))

    @ICM_SED_BENDOC.setter
    def ICM_SED_BENDOC(self, ICM_SED_BENDOC):
        self._nml['schout']['iof_icm'][5] = int(bool(ICM_SED_BENDOC))

    @ICM_SED_BENNH4.setter
    def ICM_SED_BENNH4(self, ICM_SED_BENNH4):
        self._nml['schout']['iof_icm'][6] = int(bool(ICM_SED_BENNH4))

    @ICM_SED_BENNO3.setter
    def ICM_SED_BENNO3(self, ICM_SED_BENNO3):
        self._nml['schout']['iof_icm'][7] = int(bool(ICM_SED_BENNO3))

    @ICM_SED_BENPO4.setter
    def ICM_SED_BENPO4(self, ICM_SED_BENPO4):
        self._nml['schout']['iof_icm'][8] = int(bool(ICM_SED_BENPO4))

    @ICM_SED_BENCOD.setter
    def ICM_SED_BENCOD(self, ICM_SED_BENCOD):
        self._nml['schout']['iof_icm'][9] = int(bool(ICM_SED_BENCOD))

    @ICM_SED_BENDO.setter
    def ICM_SED_BENDO(self, ICM_SED_BENDO):
        self._nml['schout']['iof_icm'][10] = int(bool(ICM_SED_BENDO))

    @ICM_SED_BENSA.setter
    def ICM_SED_BENSA(self, ICM_SED_BENSA):
        self._nml['schout']['iof_icm'][11] = int(bool(ICM_SED_BENSA))

    @ICM_lfsav.setter
    def ICM_lfsav(self, ICM_lfsav):
        self._nml['schout']['iof_icm'][12] = int(bool(ICM_lfsav))

    @ICM_stsav.setter
    def ICM_stsav(self, ICM_stsav):
        self._nml['schout']['iof_icm'][13] = int(bool(ICM_stsav))

    @ICM_rtsav.setter
    def ICM_rtsav(self, ICM_rtsav):
        self._nml['schout']['iof_icm'][14] = int(bool(ICM_rtsav))

    @ICM_tlfsav.setter
    def ICM_tlfsav(self, ICM_tlfsav):
        self._nml['schout']['iof_icm'][15] = int(bool(ICM_tlfsav))

    @ICM_tstsav.setter
    def ICM_tstsav(self, ICM_tstsav):
        self._nml['schout']['iof_icm'][16] = int(bool(ICM_tstsav))

    @ICM_trtsav.setter
    def ICM_trtsav(self, ICM_trtsav):
        self._nml['schout']['iof_icm'][17] = int(bool(ICM_trtsav))

    @ICM_hcansav.setter
    def ICM_hcansav(self, ICM_hcansav):
        self._nml['schout']['iof_icm'][18] = int(bool(ICM_hcansav))

    @ICM_CNH4.setter
    def ICM_CNH4(self, ICM_CNH4):
        self._nml['schout']['iof_icm'][19] = int(bool(ICM_CNH4))

    @ICM_CNH3.setter
    def ICM_CNH3(self, ICM_CNH3):
        self._nml['schout']['iof_icm'][20] = int(bool(ICM_CNH3))

    @ICM_CPIP.setter
    def ICM_CPIP(self, ICM_CPIP):
        self._nml['schout']['iof_icm'][21] = int(bool(ICM_CPIP))

    @ICM_CPOS.setter
    def ICM_CPOS(self, ICM_CPOS):
        self._nml['schout']['iof_icm'][22] = int(bool(ICM_CPOS))

    @ICM_CCH4.setter
    def ICM_CCH4(self, ICM_CCH4):
        self._nml['schout']['iof_icm'][23] = int(bool(ICM_CCH4))

    @ICM_CSO4.setter
    def ICM_CSO4(self, ICM_CSO4):
        self._nml['schout']['iof_icm'][24] = int(bool(ICM_CSO4))

    @ICM_CH2S.setter
    def ICM_CH2S(self, ICM_CH2S):
        self._nml['schout']['iof_icm'][25] = int(bool(ICM_CH2S))

    @ICM_SEDPON1.setter
    def ICM_SEDPON1(self, ICM_SEDPON1):
        self._nml['schout']['iof_icm'][26] = int(bool(ICM_SEDPON1))

    @ICM_SEDPON2.setter
    def ICM_SEDPON2(self, ICM_SEDPON2):
        self._nml['schout']['iof_icm'][27] = int(bool(ICM_SEDPON2))

    @ICM_SEDPON3.setter
    def ICM_SEDPON3(self, ICM_SEDPON3):
        self._nml['schout']['iof_icm'][28] = int(bool(ICM_SEDPON3))

    @ICM_SEDPOP1.setter
    def ICM_SEDPOP1(self, ICM_SEDPOP1):
        self._nml['schout']['iof_icm'][29] = int(bool(ICM_SEDPOP1))

    @ICM_SEDPOP2.setter
    def ICM_SEDPOP2(self, ICM_SEDPOP2):
        self._nml['schout']['iof_icm'][30] = int(bool(ICM_SEDPOP2))

    @ICM_SEDPOP3.setter
    def ICM_SEDPOP3(self, ICM_SEDPOP3):
        self._nml['schout']['iof_icm'][31] = int(bool(ICM_SEDPOP3))

    @ICM_SEDPOC1.setter
    def ICM_SEDPOC1(self, ICM_SEDPOC1):
        self._nml['schout']['iof_icm'][32] = int(bool(ICM_SEDPOC1))

    @ICM_SEDPOC2.setter
    def ICM_SEDPOC2(self, ICM_SEDPOC2):
        self._nml['schout']['iof_icm'][33] = int(bool(ICM_SEDPOC2))

    @ICM_SEDPOC3.setter
    def ICM_SEDPOC3(self, ICM_SEDPOC3):
        self._nml['schout']['iof_icm'][34] = int(bool(ICM_SEDPOC3))

    @ICM_EROH2S.setter
    def ICM_EROH2S(self, ICM_EROH2S):
        self._nml['schout']['iof_icm'][35] = int(bool(ICM_EROH2S))

    @ICM_EROLPOC.setter
    def ICM_EROLPOC(self, ICM_EROLPOC):
        self._nml['schout']['iof_icm'][36] = int(bool(ICM_EROLPOC))

    @ICM_ERORPOC.setter
    def ICM_ERORPOC(self, ICM_ERORPOC):
        self._nml['schout']['iof_icm'][37] = int(bool(ICM_ERORPOC))

    @ICM_DO_consumption.setter
    def ICM_DO_consumption(self, ICM_DO_consumption):
        self._nml['schout']['iof_icm'][38] = int(bool(ICM_DO_consumption))

    @ICM_GP1.setter
    def ICM_GP1(self, ICM_GP1):
        self._nml['schout']['iof_icm'][39] = int(bool(ICM_GP1))

    @ICM_GP2.setter
    def ICM_GP2(self, ICM_GP2):
        self._nml['schout']['iof_icm'][40] = int(bool(ICM_GP2))

    @ICM_GP3.setter
    def ICM_GP3(self, ICM_GP3):
        self._nml['schout']['iof_icm'][41] = int(bool(ICM_GP3))

    @ICM_1.setter
    def ICM_1(self, ICM_1):
        self._nml['schout']['iof_icm'][42] = int(bool(ICM_1))

    @ICM_2.setter
    def ICM_2(self, ICM_2):
        self._nml['schout']['iof_icm'][43] = int(bool(ICM_2))

    @ICM_3.setter
    def ICM_3(self, ICM_3):
        self._nml['schout']['iof_icm'][44] = int(bool(ICM_3))

    @ICM_4.setter
    def ICM_4(self, ICM_4):
        self._nml['schout']['iof_icm'][45] = int(bool(ICM_4))

    @ICM_5.setter
    def ICM_5(self, ICM_5):
        self._nml['schout']['iof_icm'][46] = int(bool(ICM_5))

    @ICM_6.setter
    def ICM_6(self, ICM_6):
        self._nml['schout']['iof_icm'][47] = int(bool(ICM_6))

    @ICM_7.setter
    def ICM_7(self, ICM_7):
        self._nml['schout']['iof_icm'][48] = int(bool(ICM_7))

    @ICM_8.setter
    def ICM_8(self, ICM_8):
        self._nml['schout']['iof_icm'][49] = int(bool(ICM_8))

    @ICM_9.setter
    def ICM_9(self, ICM_9):
        self._nml['schout']['iof_icm'][50] = int(bool(ICM_9))

    @ICM_10.setter
    def ICM_10(self, ICM_10):
        self._nml['schout']['iof_icm'][51] = int(bool(ICM_10))

    @ICM_11.setter
    def ICM_11(self, ICM_11):
        self._nml['schout']['iof_icm'][52] = int(bool(ICM_11))

    @ICM_12.setter
    def ICM_12(self, ICM_12):
        self._nml['schout']['iof_icm'][53] = int(bool(ICM_12))

    @ICM_13.setter
    def ICM_13(self, ICM_13):
        self._nml['schout']['iof_icm'][54] = int(bool(ICM_13))

    @ICM_14.setter
    def ICM_14(self, ICM_14):
        self._nml['schout']['iof_icm'][55] = int(bool(ICM_14))

    @ICM_15.setter
    def ICM_15(self, ICM_15):
        self._nml['schout']['iof_icm'][56] = int(bool(ICM_15))

    @ICM_16.setter
    def ICM_16(self, ICM_16):
        self._nml['schout']['iof_icm'][57] = int(bool(ICM_16))

    @ICM_17.setter
    def ICM_17(self, ICM_17):
        self._nml['schout']['iof_icm'][58] = int(bool(ICM_17))

    @ICM_18.setter
    def ICM_18(self, ICM_18):
        self._nml['schout']['iof_icm'][59] = int(bool(ICM_18))

    @ICM_19.setter
    def ICM_19(self, ICM_19):
        self._nml['schout']['iof_icm'][60] = int(bool(ICM_19))

    @ICM_20.setter
    def ICM_20(self, ICM_20):
        self._nml['schout']['iof_icm'][61] = int(bool(ICM_20))

    @ICM_21.setter
    def ICM_21(self, ICM_21):
        self._nml['schout']['iof_icm'][62] = int(bool(ICM_21))

    @ICM_22.setter
    def ICM_22(self, ICM_22):
        self._nml['schout']['iof_icm'][63] = int(bool(ICM_22))

    @ICM_23.setter
    def ICM_23(self, ICM_23):
        self._nml['schout']['iof_icm'][64] = int(bool(ICM_23))

    @ICM_24.setter
    def ICM_24(self, ICM_24):
        self._nml['schout']['iof_icm'][65] = int(bool(ICM_24))

    @ICM_25.setter
    def ICM_25(self, ICM_25):
        self._nml['schout']['iof_icm'][66] = int(bool(ICM_25))

    @COS_1.setter
    def COS_1(self, COS_1):
        self._nml['schout']['iof_cos'][0] = int(bool(COS_1))

    @COS_2.setter
    def COS_2(self, COS_2):
        self._nml['schout']['iof_cos'][1] = int(bool(COS_2))

    @COS_3.setter
    def COS_3(self, COS_3):
        self._nml['schout']['iof_cos'][2] = int(bool(COS_3))

    @COS_4.setter
    def COS_4(self, COS_4):
        self._nml['schout']['iof_cos'][3] = int(bool(COS_4))

    @COS_5.setter
    def COS_5(self, COS_5):
        self._nml['schout']['iof_cos'][4] = int(bool(COS_5))

    @COS_6.setter
    def COS_6(self, COS_6):
        self._nml['schout']['iof_cos'][5] = int(bool(COS_6))

    @COS_7.setter
    def COS_7(self, COS_7):
        self._nml['schout']['iof_cos'][6] = int(bool(COS_7))

    @COS_8.setter
    def COS_8(self, COS_8):
        self._nml['schout']['iof_cos'][7] = int(bool(COS_8))

    @COS_9.setter
    def COS_9(self, COS_9):
        self._nml['schout']['iof_cos'][8] = int(bool(COS_9))

    @COS_10.setter
    def COS_10(self, COS_10):
        self._nml['schout']['iof_cos'][9] = int(bool(COS_10))

    @COS_11.setter
    def COS_11(self, COS_11):
        self._nml['schout']['iof_cos'][10] = int(bool(COS_11))

    @COS_12.setter
    def COS_12(self, COS_12):
        self._nml['schout']['iof_cos'][11] = int(bool(COS_12))

    @COS_13.setter
    def COS_13(self, COS_13):
        self._nml['schout']['iof_cos'][12] = int(bool(COS_13))

    @FIB_1.setter
    def FIB_1(self, FIB_1):
        self._nml['schout']['iof_fib'][0] = int(bool(FIB_1))

    @SED2D_depth_change.setter
    def SED2D_depth_change(self, SED2D_depth_change):
        self._nml['schout']['iof_sed2d'][1] = int(bool(SED2D_depth_change))

    @SED2D_cflsed.setter
    def SED2D_cflsed(self, SED2D_cflsed):
        self._nml['schout']['iof_sed2d'][2] = int(bool(SED2D_cflsed))

    @SED2D_d50.setter
    def SED2D_d50(self, SED2D_d50):
        self._nml['schout']['iof_sed2d'][3] = int(bool(SED2D_d50))

    @SED2D_total_transport.setter
    def SED2D_total_transport(self, SED2D_total_transport):
        self._nml['schout']['iof_sed2d'][4] = int(bool(SED2D_total_transport))

    @SED2D_susp_load.setter
    def SED2D_susp_load(self, SED2D_susp_load):
        self._nml['schout']['iof_sed2d'][5] = int(bool(SED2D_susp_load))

    @SED2D_bed_load.setter
    def SED2D_bed_load(self, SED2D_bed_load):
        self._nml['schout']['iof_sed2d'][6] = int(bool(SED2D_bed_load))

    @SED2D_average_transport.setter
    def SED2D_average_transport(self, SED2D_average_transport):
        self._nml['schout']['iof_sed2d'][7] = int(bool(SED2D_average_transport))

    @SED2D_bottom_slope.setter
    def SED2D_bottom_slope(self, SED2D_bottom_slope):
        self._nml['schout']['iof_sed2d'][8] = int(bool(SED2D_bottom_slope))

    @z0eq.setter
    def z0eq(self, z0eq):
        self._nml['schout']['iof_sed2d'][9] = int(bool(z0eq))

    @z0cr2d.setter
    def z0cr2d(self, z0cr):
        self._nml['schout']['iof_sed2d'][10] = int(bool(z0cr))

    @z0sw2d.setter
    def z0sw2d(self, z0sw):
        self._nml['schout']['iof_sed2d'][11] = int(bool(z0sw))

    @z0wr2d.setter
    def z0wr2d(self, z0wr):
        self._nml['schout']['iof_sed2d'][12] = int(bool(z0wr))

    @marsh_flag.setter
    def marsh_flag(self, marsh_flag):
        self._nml['schout']['iof_marsh'][0] = int(bool(marsh_flag))

    @ICE_velocity.setter
    def ICE_velocity(self, ICE_velocity):
        self._nml['schout']['iof_ice'][0] = int(bool(ICE_velocity))

    @ICE_strain_rate.setter
    def ICE_strain_rate(self, ICE_strain_rate):
        self._nml['schout']['iof_ice'][1] = int(bool(ICE_strain_rate))

    @ICE_net_heat_flux.setter
    def ICE_net_heat_flux(self, ICE_net_heat_flux):
        self._nml['schout']['iof_ice'][2] = int(bool(ICE_net_heat_flux))

    @ICE_fresh_water_flux.setter
    def ICE_fresh_water_flux(self, ICE_fresh_water_flux):
        self._nml['schout']['iof_ice'][3] = int(bool(ICE_fresh_water_flux))

    @ICE_top_T.setter
    def ICE_top_T(self, ICE_top_T):
        self._nml['schout']['iof_ice'][4] = int(bool(ICE_top_T))

    @ICE_tracer_1.setter
    def ICE_tracer_1(self, ICE_tracer_1):
        self._nml['schout']['iof_ice'][5] = int(bool(ICE_tracer_1))

    @ICE_tracer_2.setter
    def ICE_tracer_2(self, ICE_tracer_2):
        self._nml['schout']['iof_ice'][6] = int(bool(ICE_tracer_2))

    @ICE_tracer_3.setter
    def ICE_tracer_3(self, ICE_tracer_3):
        self._nml['schout']['iof_ice'][7] = int(bool(ICE_tracer_3))

    @ANA_air_pres_grad_x.setter
    def ANA_air_pres_grad_x(self, ANA_air_pres_grad_x):
        self._nml['schout']['iof_ana'][0] = int(bool(ANA_air_pres_grad_x))

    @ANA_air_pres_grad_y.setter
    def ANA_air_pres_grad_y(self, ANA_air_pres_grad_y):
        self._nml['schout']['iof_ana'][1] = int(bool(ANA_air_pres_grad_y))

    @ANA_tide_pot_grad_x.setter
    def ANA_tide_pot_grad_x(self, ANA_tide_pot_grad_x):
        self._nml['schout']['iof_ana'][2] = int(bool(ANA_tide_pot_grad_x))

    @ANA_tide_pot_grad_y.setter
    def ANA_tide_pot_grad_y(self, ANA_tide_pot_grad_y):
        self._nml['schout']['iof_ana'][3] = int(bool(ANA_tide_pot_grad_y))

    @ANA_hor_viscosity_x.setter
    def ANA_hor_viscosity_x(self, ANA_hor_viscosity_x):
        self._nml['schout']['iof_ana'][4] = int(bool(ANA_hor_viscosity_x))

    @ANA_hor_viscosity_y.setter
    def ANA_hor_viscosity_y(self, ANA_hor_viscosity_y):
        self._nml['schout']['iof_ana'][5] = int(bool(ANA_hor_viscosity_y))

    @ANA_bclinic_force_x.setter
    def ANA_bclinic_force_x(self, ANA_bclinic_force_x):
        self._nml['schout']['iof_ana'][6] = int(bool(ANA_bclinic_force_x))

    @ANA_bclinic_force_y.setter
    def ANA_bclinic_force_y(self, ANA_bclinic_force_y):
        self._nml['schout']['iof_ana'][7] = int(bool(ANA_bclinic_force_y))

    @ANA_vert_viscosity_x.setter
    def ANA_vert_viscosity_x(self, ANA_vert_viscosity_x):
        self._nml['schout']['iof_ana'][8] = int(bool(ANA_vert_viscosity_x))

    @ANA_vert_viscosity_y.setter
    def ANA_vert_viscosity_y(self, ANA_vert_viscosity_y):
        self._nml['schout']['iof_ana'][9] = int(bool(ANA_vert_viscosity_y))

    @ANA_mom_advection_x.setter
    def ANA_mom_advection_x(self, ANA_mom_advection_x):
        self._nml['schout']['iof_ana'][10] = int(bool(ANA_mom_advection_x))

    @ANA_mom_advection_y.setter
    def ANA_mom_advection_y(self, ANA_mom_advection_y):
        self._nml['schout']['iof_ana'][11] = int(bool(ANA_mom_advection_y))

    @ANA_Richardson.setter
    def ANA_Richardson(self, ANA_Richardson):
        self._nml['schout']['iof_ana'][12] = int(bool(ANA_Richardson))

    @ANA_transport_min_dt_elem.setter
    def ANA_transport_min_dt_elem(self, ANA_transport_min_dt_elem):
        self._nml['schout']['iof_ana'][13] = int(bool(ANA_transport_min_dt_elem))
