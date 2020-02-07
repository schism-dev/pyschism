#! /usr/bin/env python
import pathlib
import tempfile
from pyschism.driver.param import Param
import unittest


PARAM = (pathlib.Path(__file__)
         / "../../../pyschism/driver/param.nml").resolve()


class ParamTestCase(unittest.TestCase):

    def setUp(self):
        self.param = Param()

    def test_read(self):
        self.assertIsInstance(Param.open(PARAM), Param)

    def test_write(self):
        p = Param()
        tmpdir = tempfile.TemporaryDirectory()
        self.assertIsNone(p.write(pathlib.Path(tmpdir.name) / 'param.nml'))

    def test_write_overwrite_False_raises(self):
        file = tempfile.NamedTemporaryFile()
        p = Param()
        self.assertRaises(OSError, p.write, file.name)

    def test_core_ipre(self):
        self.param.ipre = 0
        self.assertEqual(self.param.core.ipre, 0)

    def test_core_ibc(self):
        self.param.ibc = 0
        self.assertEqual(self.param.core.ibc, 0)

    def test_core_ibtp(self):
        self.param.ibtp = 1
        self.assertEqual(self.param.core.ibtp, 1)

    def test_core_rnday(self):
        self.param.core.rnday = 30
        self.assertEqual(self.param.core.rnday, 30)

    def test_core_dt(self):
        self.param.timestep = 100.0
        self.assertEqual(self.param.core.dt, 100.0)

    def test_core_msc2(self):
        self.param.msc2 = 24
        self.assertEqual(self.param.core.msc2, 24)

    def test_core_mdc2(self):
        self.param.mdc2 = 30
        self.assertEqual(self.param.core.mdc2, 30)

    def test_core_ntracer_gen(self):
        self.param.ntracer_gen = 2
        self.assertEqual(self.param.core.ntracer_gen, 2)

    def test_core_ntracer_age(self):
        self.param.ntracer_age = 4
        self.assertEqual(self.param.core.ntracer_age, 4)

    def test_core_sed_class(self):
        self.param.sed_class = 5
        self.assertEqual(self.param.core.sed_class, 5)

    def test_core_eco_class(self):
        self.param.eco_class = 27
        self.assertEqual(self.param.core.eco_class, 27)

    def test_core_nspool(self):
        self.param.nspool = 36
        self.assertEqual(self.param.core.nspool, 36)

    def test_core_ihfskip(self):
        self.param.ihfskip = 864
        self.assertEqual(self.param.core.ihfskip, 864)

    def test_opt_ipre2(self):
        self.param.opt.ipre2 = 0
        self.assertEqual(self.param.opt.ipre2, 0)

    def test_opt_start_year(self):
        self.param.opt.start_year = 2000
        self.assertEqual(self.param.opt.start_year, 2000)

    def test_opt_start_month(self):
        self.param.opt.start_month = 1
        self.assertEqual(self.param.opt.start_month, 1)

    def test_opt_start_day(self):
        self.param.opt.start_day = 1
        self.assertEqual(self.param.opt.start_day, 1)

    def test_opt_start_hour(self):
        self.param.opt.start_hour = 0
        self.assertEqual(self.param.opt.start_hour, 0)

    def test_opt_utc_start(self):
        self.param.opt.utc_start = 8
        self.assertEqual(self.param.opt.utc_start, 8)

    def test_opt_ics(self):
        self.param.opt.ics = 1
        self.assertEqual(self.param.opt.ics, 1)

    def test_opt_ihot(self):
        self.param.opt.ihot = 0
        self.assertEqual(self.param.opt.ihot, 0)

    def test_opt_ieos_type(self):
        self.param.opt.ieos_type = 0
        self.assertEqual(self.param.opt.ieos_type, 0)

    def test_opt_ieos_pres(self):
        self.param.opt.ieos_pres = 0
        self.assertEqual(self.param.opt.ieos_pres, 0)

    def test_opt_eos_a(self):
        self.param.opt.eos_a = -0.1
        self.assertEqual(self.param.opt.eos_a, -0.1)

    def test_opt_eos_b(self):
        self.param.opt.eos_b = 1001.0
        self.assertEqual(self.param.opt.eos_b, 1001.0)

    def test_opt_nramp(self):
        self.param.opt.nramp = 1
        self.assertEqual(self.param.opt.nramp, 1)

    def test_opt_dramp(self):
        self.param.opt.dramp = 1.0
        self.assertEqual(self.param.opt.dramp, 1.0)

    def test_opt_nrampbc(self):
        self.param.opt.nrampbc = 0
        self.assertEqual(self.param.opt.nrampbc, 0)

    def test_opt_drampbc(self):
        self.param.opt.drampbc = 1.0
        self.assertEqual(self.param.opt.drampbc, 1.0)

    def test_opt_iupwind_mom(self):
        self.param.opt.iupwind_mom = 0
        self.assertEqual(self.param.opt.iupwind_mom, 0)

    def test_opt_indvel(self):
        self.param.opt.indvel = 0
        self.assertEqual(self.param.opt.indvel, 0)

    def test_opt_ihorcon(self):
        self.param.opt.ihorcon = 0
        self.assertEqual(self.param.opt.ihorcon, 0)

    def test_opt_hvis_coef0(self):
        self.param.opt.hvis_coef0 = 0.025
        self.assertEqual(self.param.opt.hvis_coef0, 0.025)

    def test_opt_ishapiro(self):
        self.param.opt.ishapiro = 1
        self.assertEqual(self.param.opt.ishapiro, 1)

    def test_opt_shapiro0(self):
        self.param.opt.shapiro0 = 0.5
        self.assertEqual(self.param.opt.shapiro0, 0.5)

    def test_opt_niter_shap(self):
        self.param.opt.niter_shap = 1
        self.assertEqual(self.param.opt.niter_shap, 1)

    def test_opt_thetai(self):
        self.param.opt.thetai = 0.6
        self.assertEqual(self.param.opt.thetai, 0.6)

    def test_opt_icou_elfe_wwm(self):
        self.param.opt.icou_elfe_wwm = 0
        self.assertEqual(self.param.opt.icou_elfe_wwm, 0)

    def test_opt_nstep_wwm(self):
        self.param.opt.nstep_wwm = 1
        self.assertEqual(self.param.opt.nstep_wwm, 1)

    def test_opt_iwbl(self):
        self.param.opt.iwbl = 0
        self.assertEqual(self.param.opt.iwbl, 0)

    def test_opt_hmin_radstress(self):
        self.param.opt.hmin_radstress = 1.0
        self.assertEqual(self.param.opt.hmin_radstress, 1.0)

    def test_opt_nrampwafo(self):
        self.param.opt.nrampwafo = 0
        self.assertEqual(self.param.opt.nrampwafo, 0)

    def test_opt_drampwafo(self):
        self.param.opt.drampwafo = 1.0
        self.assertEqual(self.param.opt.drampwafo, 1.0)

    def test_opt_turbinj(self):
        self.param.opt.turbinj = 0.15
        self.assertEqual(self.param.opt.turbinj, 0.15)

    def test_opt_imm(self):
        self.param.opt.imm = 0
        self.assertEqual(self.param.opt.imm, 0)

    def test_opt_ibdef(self):
        self.param.opt.ibdef = 10
        self.assertEqual(self.param.opt.ibdef, 10)

    def test_opt_slam0(self):
        self.param.opt.slam0 = -124
        self.assertEqual(self.param.opt.slam0, -124)

    def test_opt_sfea0(self):
        self.param.opt.sfea0 = 45
        self.assertEqual(self.param.opt.sfea0, 45)

    def test_opt_iunder_deep(self):
        self.param.opt.iunder_deep = 0
        self.assertEqual(self.param.opt.iunder_deep, 0)

    def test_opt_h1_bcc(self):
        self.param.opt.h1_bcc = 50.0
        self.assertEqual(self.param.opt.h1_bcc, 50.0)

    def test_opt_h2_bcc(self):
        self.param.opt.h2_bcc = 100.0
        self.assertEqual(self.param.opt.h2_bcc, 100.0)

    def test_opt_hw_depth(self):
        self.param.opt.hw_depth = 1000000.0
        self.assertEqual(self.param.opt.hw_depth, 1000000.0)

    def test_opt_hw_ratio(self):
        self.param.opt.hw_ratio = 0.5
        self.assertEqual(self.param.opt.hw_ratio, 0.5)

    def test_opt_ihydraulics(self):
        self.param.opt.ihydraulics = 0
        self.assertEqual(self.param.opt.ihydraulics, 0)

    def test_opt_if_source(self):
        self.param.opt.if_source = 0
        self.assertEqual(self.param.opt.if_source, 0)

    def test_opt_nramp_ss(self):
        self.param.opt.nramp_ss = 1
        self.assertEqual(self.param.opt.nramp_ss, 1)

    def test_opt_dramp_ss(self):
        self.param.opt.dramp_ss = 2
        self.assertEqual(self.param.opt.dramp_ss, 2)

    def test_opt_ihdif(self):
        self.param.opt.ihdif = 0
        self.assertEqual(self.param.opt.ihdif, 0)

    def test_opt_nchi(self):
        self.param.opt.nchi = 0
        self.assertEqual(self.param.opt.nchi, 0)

    def test_opt_dzb_min(self):
        self.param.opt.dzb_min = 0.5
        self.assertEqual(self.param.opt.dzb_min, 0.5)

    def test_opt_dzb_decay(self):
        self.param.opt.dzb_decay = 0.0
        self.assertEqual(self.param.opt.dzb_decay, 0.0)

    def test_opt_hmin_man(self):
        self.param.opt.hmin_man = 1.0
        self.assertEqual(self.param.opt.hmin_man, 1.0)

    def test_opt_ncor(self):
        self.param.opt.ncor = 0
        self.assertEqual(self.param.opt.ncor, 0)

    def test_opt_rlatitude(self):
        self.param.opt.rlatitude = 46
        self.assertEqual(self.param.opt.rlatitude, 46)

    def test_opt_coricoef(self):
        self.param.opt.coricoef = 0
        self.assertEqual(self.param.opt.coricoef, 0)

    def test_opt_ic_elev(self):
        self.param.opt.ic_elev = 0
        self.assertEqual(self.param.opt.ic_elev, 0)

    def test_opt_nramp_elev(self):
        self.param.opt.nramp_elev = 0
        self.assertEqual(self.param.opt.nramp_elev, 0)

    def test_opt_inv_atm_bnd(self):
        self.param.opt.inv_atm_bnd = 0
        self.assertEqual(self.param.opt.inv_atm_bnd, 0)

    def test_opt_prmsl_ref(self):
        self.param.opt.prmsl_ref = 101325.0
        self.assertEqual(self.param.opt.prmsl_ref, 101325.0)

    def test_opt_(self):
        self.param.opt.flag_ic = [1, 1, 0, None, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(self.param.opt.flag_ic, [1, 1, 0, None, 0, 0, 0, 0, 0, 0, 0])

    def test_opt_gen_wsett(self):
        self.param.opt.gen_wsett = 0.0001
        self.assertEqual(self.param.opt.gen_wsett, 0.0001)

    def test_opt_ibcc_mean(self):
        self.param.opt.ibcc_mean = 0
        self.assertEqual(self.param.opt.ibcc_mean, 0)

    def test_opt_rmaxvel(self):
        self.param.opt.rmaxvel = 10.0
        self.assertEqual(self.param.opt.rmaxvel, 10.0)

    def test_opt_velmin_btrack(self):
        self.param.opt.velmin_btrack = 0.0001
        self.assertEqual(self.param.opt.velmin_btrack, 0.0001)

    def test_opt_btrack_nudge(self):
        self.param.opt.btrack_nudge = 0.009013
        self.assertEqual(self.param.opt.btrack_nudge, 0.009013)

    def test_opt_ibtrack_openbnd(self):
        self.param.opt.ibtrack_openbnd = 1
        self.assertEqual(self.param.opt.ibtrack_openbnd, 1)

    def test_opt_ihhat(self):
        self.param.opt.ihhat = 1
        self.assertEqual(self.param.opt.ihhat, 1)

    def test_opt_inunfl(self):
        self.param.opt.inunfl = 0
        self.assertEqual(self.param.opt.inunfl, 0)

    def test_opt_h0(self):
        self.param.opt.h0 = 0.01
        self.assertEqual(self.param.opt.h0, 0.01)

    def test_opt_shorewafo(self):
        self.param.opt.shorewafo = 0
        self.assertEqual(self.param.opt.shorewafo, 0)

    def test_opt_moitn0(self):
        self.param.opt.moitn0 = 50
        self.assertEqual(self.param.opt.moitn0, 50)

    def test_opt_mxitn0(self):
        self.param.opt.mxitn0 = 1500
        self.assertEqual(self.param.opt.mxitn0, 1500)

    def test_opt_rtol0(self):
        self.param.opt.rtol0 = 1e-12
        self.assertEqual(self.param.opt.rtol0, 1e-12)

    def test_opt_nadv(self):
        self.param.opt.nadv = 1
        self.assertEqual(self.param.opt.nadv, 1)

    def test_opt_dtb_max(self):
        self.param.opt.dtb_max = 30.0
        self.assertEqual(self.param.opt.dtb_max, 30.0)

    def test_opt_dtb_min(self):
        self.param.opt.dtb_min = 10.0
        self.assertEqual(self.param.opt.dtb_min, 10.0)

    def test_opt_inter_mom(self):
        self.param.opt.inter_mom = 0
        self.assertEqual(self.param.opt.inter_mom, 0)

    def test_opt_kr_co(self):
        self.param.opt.kr_co = 1
        self.assertEqual(self.param.opt.kr_co, 1)

    def test_opt_itr_met(self):
        self.param.opt.itr_met = 3
        self.assertEqual(self.param.opt.itr_met, 3)

    def test_opt_h_tvd(self):
        self.param.opt.h_tvd = 5.0
        self.assertEqual(self.param.opt.h_tvd, 5.0)

    def test_opt_eps1_tvd_imp(self):
        self.param.opt.eps1_tvd_imp = 0.0001
        self.assertEqual(self.param.opt.eps1_tvd_imp, 0.0001)

    def test_opt_eps2_tvd_imp(self):
        self.param.opt.eps2_tvd_imp = 1e-14
        self.assertEqual(self.param.opt.eps2_tvd_imp, 1e-14)

    def test_opt_ip_weno(self):
        self.param.opt.ip_weno = 2
        self.assertEqual(self.param.opt.ip_weno, 2)

    def test_opt_courant_weno(self):
        self.param.opt.courant_weno = 0.5
        self.assertEqual(self.param.opt.courant_weno, 0.5)

    def test_opt_nquad(self):
        self.param.opt.nquad = 2
        self.assertEqual(self.param.opt.nquad, 2)

    def test_opt_ntd_weno(self):
        self.param.opt.ntd_weno = 1
        self.assertEqual(self.param.opt.ntd_weno, 1)

    def test_opt_epsilon1(self):
        self.param.opt.epsilon1 = 0.001
        self.assertEqual(self.param.opt.epsilon1, 0.001)

    def test_opt_epsilon2(self):
        self.param.opt.epsilon2 = 1e-10
        self.assertEqual(self.param.opt.epsilon2, 1e-10)

    def test_opt_i_prtnftl_weno(self):
        self.param.opt.i_prtnftl_weno = 0
        self.assertEqual(self.param.opt.i_prtnftl_weno, 0)

    def test_opt_epsilon3(self):
        self.param.opt.epsilon3 = 1e-25
        self.assertEqual(self.param.opt.epsilon3, 1e-25)

    def test_opt_ielad_weno(self):
        self.param.opt.ielad_weno = 0
        self.assertEqual(self.param.opt.ielad_weno, 0)

    def test_opt_small_elad(self):
        self.param.opt.small_elad = 0.0001
        self.assertEqual(self.param.opt.small_elad, 0.0001)

    def test_opt_nws(self):
        self.param.opt.nws = 0
        self.assertEqual(self.param.opt.nws, 0)

    def test_opt_wtiminc(self):
        self.param.opt.wtiminc = 150.0
        self.assertEqual(self.param.opt.wtiminc, 150.0)

    def test_opt_nrampwind(self):
        self.param.opt.nrampwind = 1
        self.assertEqual(self.param.opt.nrampwind, 1)

    def test_opt_drampwind(self):
        self.param.opt.drampwind = 1.0
        self.assertEqual(self.param.opt.drampwind, 1.0)

    def test_opt_iwindoff(self):
        self.param.opt.iwindoff = 0
        self.assertEqual(self.param.opt.iwindoff, 0)

    def test_opt_iwind_form(self):
        self.param.opt.iwind_form = -1
        self.assertEqual(self.param.opt.iwind_form, -1)

    def test_opt_impose_net_flux(self):
        self.param.opt.impose_net_flux = 0
        self.assertEqual(self.param.opt.impose_net_flux, 0)

    def test_opt_ihconsv(self):
        self.param.opt.ihconsv = 0
        self.assertEqual(self.param.opt.ihconsv, 0)

    def test_opt_isconsv(self):
        self.param.opt.isconsv = 0
        self.assertEqual(self.param.opt.isconsv, 0)

    def test_opt_itur(self):
        self.param.opt.itur = 3
        self.assertEqual(self.param.opt.itur, 3)

    def test_opt_dfv0(self):
        self.param.opt.dfv0 = 0.01
        self.assertEqual(self.param.opt.dfv0, 0.01)

    def test_opt_dfh0(self):
        self.param.opt.dfh0 = 0.0001
        self.assertEqual(self.param.opt.dfh0, 0.0001)

    def test_opt_mid(self):
        self.param.opt.mid = 'KL'
        self.assertEqual(self.param.opt.mid, 'KL')

    def test_opt_stab(self):
        self.param.opt.stab = 'KC'
        self.assertEqual(self.param.opt.stab, 'KC')

    def test_opt_xlsc0(self):
        self.param.opt.xlsc0 = 0.1
        self.assertEqual(self.param.opt.xlsc0, 0.1)

    def test_opt_inu_elev(self):
        self.param.opt.inu_elev = 0
        self.assertEqual(self.param.opt.inu_elev, 0)

    def test_opt_inu_uv(self):
        self.param.opt.inu_uv = 0
        self.assertEqual(self.param.opt.inu_uv, 0)

    def test_opt_inu_tr(self):
        self.param.opt.inu_tr = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(self.param.opt.inu_tr, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_opt_vnh1(self):
        self.param.opt.vnh1 = 400
        self.assertEqual(self.param.opt.vnh1, 400)

    def test_opt_vnf1(self):
        self.param.opt.vnf1 = 0.0
        self.assertEqual(self.param.opt.vnf1, 0.0)

    def test_opt_vnh2(self):
        self.param.opt.vnh2 = 500
        self.assertEqual(self.param.opt.vnh2, 500)

    def test_opt_vnf2(self):
        self.param.opt.vnf2 = 0.0
        self.assertEqual(self.param.opt.vnf2, 0.0)

    def test_opt_step_nu_tr(self):
        self.param.opt.step_nu_tr = 86400.0
        self.assertEqual(self.param.opt.step_nu_tr, 86400.0)

    def test_opt_h_bcc1(self):
        self.param.opt.h_bcc1 = 100.0
        self.assertEqual(self.param.opt.h_bcc1, 100.0)

    def test_opt_s1_mxnbt(self):
        self.param.opt.s1_mxnbt = 0.5
        self.assertEqual(self.param.opt.s1_mxnbt, 0.5)

    def test_opt_s2_mxnbt(self):
        self.param.opt.s2_mxnbt = 3.5
        self.assertEqual(self.param.opt.s2_mxnbt, 3.5)

    def test_opt_iharind(self):
        self.param.opt.iharind = 0
        self.assertEqual(self.param.opt.iharind, 0)

    def test_opt_iflux(self):
        self.param.opt.iflux = 0
        self.assertEqual(self.param.opt.iflux, 0)

    def test_opt_izonal5(self):
        self.param.opt.izonal5 = 0
        self.assertEqual(self.param.opt.izonal5, 0)

    def test_opt_ibtrack_test(self):
        self.param.opt.ibtrack_test = 0
        self.assertEqual(self.param.opt.ibtrack_test, 0)

    def test_opt_irouse_test(self):
        self.param.opt.irouse_test = 0
        self.assertEqual(self.param.opt.irouse_test, 0)

    def test_opt_flag_fib(self):
        self.param.opt.flag_fib = 1
        self.assertEqual(self.param.opt.flag_fib, 1)

    def test_opt_slr_rate(self):
        self.param.opt.slr_rate = 120.0
        self.assertEqual(self.param.opt.slr_rate, 120.0)

    def test_opt_isav(self):
        self.param.opt.isav = 0
        self.assertEqual(self.param.opt.isav, 0)

    def test_opt_sav_cd(self):
        self.param.opt.sav_cd = 1.13
        self.assertEqual(self.param.opt.sav_cd, 1.13)

    def test_opt_nstep_ice(self):
        self.param.opt.nstep_ice = 1
        self.assertEqual(self.param.opt.nstep_ice, 1)

    def test_opt_level_age(self):
        self.param.opt.level_age = [9, -999]
        self.assertEqual(self.param.opt.level_age, [9, -999])

    def test_opt_rearth_pole(self):
        self.param.opt.rearth_pole = 6378206.4
        self.assertEqual(self.param.opt.rearth_pole, 6378206.4)

    def test_opt_rearth_eq(self):
        self.param.opt.rearth_eq = 6378206.4
        self.assertEqual(self.param.opt.rearth_eq, 6378206.4)

    def test_opt_shw(self):
        self.param.opt.shw = 4184.0
        self.assertEqual(self.param.opt.shw, 4184.0)

    def test_opt_rho0(self):
        self.param.opt.rho0 = 1000.0
        self.assertEqual(self.param.opt.rho0, 1000.0)

    def test_opt_vclose_surf_frac(self):
        self.param.opt.vclose_surf_frac = 1.0
        self.assertEqual(self.param.opt.vclose_surf_frac, 1.0)

    def test_opt_iadjust_mass_consv0(self):
        self.param.opt.iadjust_mass_consv0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(self.param.opt.iadjust_mass_consv0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_schout_nhot(self):
        self.param.schout.nhot = 0
        self.assertEqual(self.param.schout.nhot, 0)

    def test_schout_nhot_write(self):
        self.param.schout.nhot_write = 8640
        self.assertEqual(self.param.schout.nhot_write, 8640)

    def test_schout_iout_sta(self):
        self.param.schout.iout_sta = 0
        self.assertEqual(self.param.schout.iout_sta, 0)

    def test_schout_nspool_sta(self):
        self.param.schout.nspool_sta = 10
        self.assertEqual(self.param.schout.nspool_sta, 10)

    def test_schout_iof_hydro(self):
        self.param.schout.iof_hydro = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(self.param.schout.iof_hydro, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_schout_iof_wwm(self):
        self.param.schout.iof_wwm = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(self.param.schout.iof_wwm, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_schout_iof_gen(self):
        self.param.schout.iof_gen = [0, 0]
        self.assertEqual(self.param.schout.iof_gen, [0, 0])

    def test_schout_iof_age(self):
        self.param.schout.iof_age = [0, 0]
        self.assertEqual(self.param.schout.iof_age, [0, 0])

    def test_schout_iof_sed(self):
        self.param.schout.iof_sed = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(self.param.schout.iof_sed, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_schout_iof_eco(self):
        self.param.schout.iof_eco = [0]
        self.assertEqual(self.param.schout.iof_eco, [0])

    def test_schout_iof_icm(self):
        self.param.schout.iof_icm = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(self.param.schout.iof_icm, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_schout_iof_cos(self):
        self.param.schout.iof_cos = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(self.param.schout.iof_cos, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_schout_iof_fib(self):
        self.param.schout.iof_fib = [0]
        self.assertEqual(self.param.schout.iof_fib, [0])

    def test_schout_iof_sed2d(self):
        self.param.schout.iof_sed2d = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(self.param.schout.iof_sed2d, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_schout_iof_marsh(self):
        self.param.schout.iof_marsh = [0]
        self.assertEqual(self.param.schout.iof_marsh, [0])

    def test_schout_iof_ice(self):
        self.param.schout.iof_ice = [0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(self.param.schout.iof_ice, [0, 0, 0, 0, 0, 0, 0, 0])

    def test_schout_iof_ana(self):
        self.param.schout.iof_ana = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(self.param.schout.iof_ana, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_elev(self):
        self.param.elev = 1
        self.assertEqual(self.param.elev, 1)

    def test_air_pressure(self):
        self.param.air_pressure = 1
        self.assertEqual(self.param.air_pressure, 1)

    def test_air_temperature(self):
        self.param.air_temperature = 1
        self.assertEqual(self.param.air_temperature, 1)

    def test_specific_humidity(self):
        self.param.specific_humidity = 1
        self.assertEqual(self.param.specific_humidity, 1)

    def test_solar_radiation(self):
        self.param.solar_radiation = 1
        self.assertEqual(self.param.solar_radiation, 1)

    def test_sensible_flux(self):
        self.param.sensible_flux = 1
        self.assertEqual(self.param.sensible_flux, 1)

    def test_latent_heat(self):
        self.param.latent_heat = 1
        self.assertEqual(self.param.latent_heat, 1)

    def test_upward_longwave(self):
        self.param.upward_longwave = 1
        self.assertEqual(self.param.upward_longwave, 1)

    def test_downward_longwave(self):
        self.param.downward_longwave = 1
        self.assertEqual(self.param.downward_longwave, 1)

    def test_total_heat_flux(self):
        self.param.total_heat_flux = 1
        self.assertEqual(self.param.total_heat_flux, 1)

    def test_evaporation(self):
        self.param.evaporation = 1
        self.assertEqual(self.param.evaporation, 1)

    def test_precipitation(self):
        self.param.precipitation = 1
        self.assertEqual(self.param.precipitation, 1)

    def test_bottom_stress(self):
        self.param.bottom_stress = 1
        self.assertEqual(self.param.bottom_stress, 1)

    def test_wind_speed(self):
        self.param.wind_speed = 1
        self.assertEqual(self.param.wind_speed, 1)

    def test_wind_stress(self):
        self.param.wind_stress = 1
        self.assertEqual(self.param.wind_stress, 1)

    def test_dahv(self):
        self.param.dahv = 1
        self.assertEqual(self.param.dahv, 1)

    def test_vertical_velocity(self):
        self.param.vertical_velocity = 1
        self.assertEqual(self.param.vertical_velocity, 1)

    def test_temp(self):
        self.param.temp = 1
        self.assertEqual(self.param.temp, 1)

    def test_salt(self):
        self.param.salt = 1
        self.assertEqual(self.param.salt, 1)

    def test_water_density(self):
        self.param.water_density = 1
        self.assertEqual(self.param.water_density, 1)

    def test_diffusivity(self):
        self.param.diffusivity = 1
        self.assertEqual(self.param.diffusivity, 1)

    def test_viscosity(self):
        self.param.viscosity = 1
        self.assertEqual(self.param.viscosity, 1)

    def test_TKE(self):
        self.param.TKE = 1
        self.assertEqual(self.param.TKE, 1)

    def test_mixing_length(self):
        self.param.mixing_length = 1
        self.assertEqual(self.param.mixing_length, 1)

    def test_hvel(self):
        self.param.hvel = 1
        self.assertEqual(self.param.hvel, 1)

    def test_hvel_side(self):
        self.param.hvel_side = 1
        self.assertEqual(self.param.hvel_side, 1)

    def test_wvel_elem(self):
        self.param.wvel_elem = 1
        self.assertEqual(self.param.wvel_elem, 1)

    def test_temp_elem(self):
        self.param.temp_elem = 1
        self.assertEqual(self.param.temp_elem, 1)

    def test_salt_elem(self):
        self.param.salt_elem = 1
        self.assertEqual(self.param.salt_elem, 1)

    def test_pressure_gradient(self):
        self.param.pressure_gradient = 1
        self.assertEqual(self.param.pressure_gradient, 1)

    def test_WWM_1(self):
        self.param.WWM_1 = 1
        self.assertEqual(self.param.WWM_1, 1)

    def test_WWM_2(self):
        self.param.WWM_2 = 1
        self.assertEqual(self.param.WWM_2, 1)

    def test_WWM_3(self):
        self.param.WWM_3 = 1
        self.assertEqual(self.param.WWM_3, 1)

    def test_WWM_4(self):
        self.param.WWM_4 = 1
        self.assertEqual(self.param.WWM_4, 1)

    def test_WWM_5(self):
        self.param.WWM_5 = 1
        self.assertEqual(self.param.WWM_5, 1)

    def test_WWM_6(self):
        self.param.WWM_6 = 1
        self.assertEqual(self.param.WWM_6, 1)

    def test_WWM_9(self):
        self.param.WWM_9 = 1
        self.assertEqual(self.param.WWM_9, 1)

    def test_WWM_10(self):
        self.param.WWM_10 = 1
        self.assertEqual(self.param.WWM_10, 1)

    def test_WWM_11(self):
        self.param.WWM_11 = 1
        self.assertEqual(self.param.WWM_11, 1)

    def test_WWM_12(self):
        self.param.WWM_12 = 1
        self.assertEqual(self.param.WWM_12, 1)

    def test_WWM_13(self):
        self.param.WWM_13 = 1
        self.assertEqual(self.param.WWM_13, 1)

    def test_WWM_14(self):
        self.param.WWM_14 = 1
        self.assertEqual(self.param.WWM_14, 1)

    def test_WWM_15(self):
        self.param.WWM_15 = 1
        self.assertEqual(self.param.WWM_15, 1)

    def test_WWM_16(self):
        self.param.WWM_16 = 1
        self.assertEqual(self.param.WWM_16, 1)

    def test_WWM_17(self):
        self.param.WWM_17 = 1
        self.assertEqual(self.param.WWM_17, 1)

    def test_WWM_18(self):
        self.param.WWM_18 = 1
        self.assertEqual(self.param.WWM_18, 1)

    def test_WWM_19(self):
        self.param.WWM_19 = 1
        self.assertEqual(self.param.WWM_19, 1)

    def test_WWM_20(self):
        self.param.WWM_20 = 1
        self.assertEqual(self.param.WWM_20, 1)

    def test_WWM_21(self):
        self.param.WWM_21 = 1
        self.assertEqual(self.param.WWM_21, 1)

    def test_WWM_22(self):
        self.param.WWM_22 = 1
        self.assertEqual(self.param.WWM_22, 1)

    def test_WWM_23(self):
        self.param.WWM_23 = 1
        self.assertEqual(self.param.WWM_23, 1)

    def test_WWM_24(self):
        self.param.WWM_24 = 1
        self.assertEqual(self.param.WWM_24, 1)

    def test_WWM_25(self):
        self.param.WWM_25 = 1
        self.assertEqual(self.param.WWM_25, 1)

    def test_WWM_26(self):
        self.param.WWM_26 = 1
        self.assertEqual(self.param.WWM_26, 1)

    def test_WWM_27(self):
        self.param.WWM_27 = 1
        self.assertEqual(self.param.WWM_27, 1)

    def test_WWM_28(self):
        self.param.WWM_28 = 1
        self.assertEqual(self.param.WWM_28, 1)

    def test_WWM_energy_dir(self):
        self.param.WWM_energy_dir = 1
        self.assertEqual(self.param.WWM_energy_dir, 1)

    def test_wave_force(self):
        self.param.wave_force = 1
        self.assertEqual(self.param.wave_force, 1)

    def test_GEN_1(self):
        self.param.GEN_1 = 1
        self.assertEqual(self.param.GEN_1, 1)

    def test_GEN_2(self):
        self.param.GEN_2 = 1
        self.assertEqual(self.param.GEN_2, 1)

    def test_AGE_1(self):
        self.param.AGE_1 = 1
        self.assertEqual(self.param.AGE_1, 1)

    def test_AGE_2(self):
        self.param.AGE_2 = 1
        self.assertEqual(self.param.AGE_2, 1)

    def test_SED_depth_change(self):
        self.param.SED_depth_change = 1
        self.assertEqual(self.param.SED_depth_change, 1)

    def test_SED_D50(self):
        self.param.SED_D50 = 1
        self.assertEqual(self.param.SED_D50, 1)

    def test_SED_bed_stress(self):
        self.param.SED_bed_stress = 1
        self.assertEqual(self.param.SED_bed_stress, 1)

    def test_SED_bed_roughness(self):
        self.param.SED_bed_roughness = 1
        self.assertEqual(self.param.SED_bed_roughness, 1)

    def test_SED_TSC(self):
        self.param.SED_TSC = 1
        self.assertEqual(self.param.SED_TSC, 1)

    def test_bed_thickness(self):
        self.param.bed_thickness = 1
        self.assertEqual(self.param.bed_thickness, 1)

    def test_bed_age(self):
        self.param.bed_age = 1
        self.assertEqual(self.param.bed_age, 1)

    def test_z0st(self):
        self.param.z0st = 1
        self.assertEqual(self.param.z0st, 1)

    def test_z0cr(self):
        self.param.z0cr = 1
        self.assertEqual(self.param.z0cr, 1)

    def test_z0sw(self):
        self.param.z0sw = 1
        self.assertEqual(self.param.z0sw, 1)

    def test_z0wr(self):
        self.param.z0wr = 1
        self.assertEqual(self.param.z0wr, 1)

    def test_SED3D_1(self):
        self.param.SED3D_1 = 1
        self.assertEqual(self.param.SED3D_1, 1)

    def test_SED_bdld_1(self):
        self.param.SED_bdld_1 = 1
        self.assertEqual(self.param.SED_bdld_1, 1)

    def test_SED_bedfrac_1(self):
        self.param.SED_bedfrac_1 = 1
        self.assertEqual(self.param.SED_bedfrac_1, 1)

    def test_SED3D_2(self):
        self.param.SED3D_2 = 1
        self.assertEqual(self.param.SED3D_2, 1)

    def test_SED_bdld_2(self):
        self.param.SED_bdld_2 = 1
        self.assertEqual(self.param.SED_bdld_2, 1)

    def test_SED_bedfrac_3(self):
        self.param.SED_bedfrac_3 = 1
        self.assertEqual(self.param.SED_bedfrac_3, 1)

    def test_ECO_1(self):
        self.param.ECO_1 = 1
        self.assertEqual(self.param.ECO_1, 1)

    def test_ICM_Chl(self):
        self.param.ICM_Chl = 1
        self.assertEqual(self.param.ICM_Chl, 1)

    def test_ICM_pH(self):
        self.param.ICM_pH = 1
        self.assertEqual(self.param.ICM_pH, 1)

    def test_ICM_PrmPrdt(self):
        self.param.ICM_PrmPrdt = 1
        self.assertEqual(self.param.ICM_PrmPrdt, 1)

    def test_ICM_DIN(self):
        self.param.ICM_DIN = 1
        self.assertEqual(self.param.ICM_DIN, 1)

    def test_ICM_PON(self):
        self.param.ICM_PON = 1
        self.assertEqual(self.param.ICM_PON, 1)

    def test_ICM_SED_BENDOC(self):
        self.param.ICM_SED_BENDOC = 1
        self.assertEqual(self.param.ICM_SED_BENDOC, 1)

    def test_ICM_SED_BENNH4(self):
        self.param.ICM_SED_BENNH4 = 1
        self.assertEqual(self.param.ICM_SED_BENNH4, 1)

    def test_ICM_SED_BENNO3(self):
        self.param.ICM_SED_BENNO3 = 1
        self.assertEqual(self.param.ICM_SED_BENNO3, 1)

    def test_ICM_SED_BENPO4(self):
        self.param.ICM_SED_BENPO4 = 1
        self.assertEqual(self.param.ICM_SED_BENPO4, 1)

    def test_ICM_SED_BENCOD(self):
        self.param.ICM_SED_BENCOD = 1
        self.assertEqual(self.param.ICM_SED_BENCOD, 1)

    def test_ICM_SED_BENDO(self):
        self.param.ICM_SED_BENDO = 1
        self.assertEqual(self.param.ICM_SED_BENDO, 1)

    def test_ICM_SED_BENSA(self):
        self.param.ICM_SED_BENSA = 1
        self.assertEqual(self.param.ICM_SED_BENSA, 1)

    def test_ICM_lfsav(self):
        self.param.ICM_lfsav = 1
        self.assertEqual(self.param.ICM_lfsav, 1)

    def test_ICM_stsav(self):
        self.param.ICM_stsav = 1
        self.assertEqual(self.param.ICM_stsav, 1)

    def test_ICM_rtsav(self):
        self.param.ICM_rtsav = 1
        self.assertEqual(self.param.ICM_rtsav, 1)

    def test_ICM_tlfsav(self):
        self.param.ICM_tlfsav = 1
        self.assertEqual(self.param.ICM_tlfsav, 1)

    def test_ICM_tstsav(self):
        self.param.ICM_tstsav = 1
        self.assertEqual(self.param.ICM_tstsav, 1)

    def test_ICM_trtsav(self):
        self.param.ICM_trtsav = 1
        self.assertEqual(self.param.ICM_trtsav, 1)

    def test_ICM_hcansav(self):
        self.param.ICM_hcansav = 1
        self.assertEqual(self.param.ICM_hcansav, 1)

    def test_ICM_CNH4(self):
        self.param.ICM_CNH4 = 1
        self.assertEqual(self.param.ICM_CNH4, 1)

    def test_ICM_CNH3(self):
        self.param.ICM_CNH3 = 1
        self.assertEqual(self.param.ICM_CNH3, 1)

    def test_ICM_CPIP(self):
        self.param.ICM_CPIP = 1
        self.assertEqual(self.param.ICM_CPIP, 1)

    def test_ICM_CPOS(self):
        self.param.ICM_CPOS = 1
        self.assertEqual(self.param.ICM_CPOS, 1)

    def test_ICM_CCH4(self):
        self.param.ICM_CCH4 = 1
        self.assertEqual(self.param.ICM_CCH4, 1)

    def test_ICM_CSO4(self):
        self.param.ICM_CSO4 = 1
        self.assertEqual(self.param.ICM_CSO4, 1)

    def test_ICM_CH2S(self):
        self.param.ICM_CH2S = 1
        self.assertEqual(self.param.ICM_CH2S, 1)

    def test_ICM_SEDPON1(self):
        self.param.ICM_SEDPON1 = 1
        self.assertEqual(self.param.ICM_SEDPON1, 1)

    def test_ICM_SEDPON2(self):
        self.param.ICM_SEDPON2 = 1
        self.assertEqual(self.param.ICM_SEDPON2, 1)

    def test_ICM_SEDPON3(self):
        self.param.ICM_SEDPON3 = 1
        self.assertEqual(self.param.ICM_SEDPON3, 1)

    def test_ICM_SEDPOP1(self):
        self.param.ICM_SEDPOP1 = 1
        self.assertEqual(self.param.ICM_SEDPOP1, 1)

    def test_ICM_SEDPOP2(self):
        self.param.ICM_SEDPOP2 = 1
        self.assertEqual(self.param.ICM_SEDPOP2, 1)

    def test_ICM_SEDPOP3(self):
        self.param.ICM_SEDPOP3 = 1
        self.assertEqual(self.param.ICM_SEDPOP3, 1)

    def test_ICM_SEDPOC1(self):
        self.param.ICM_SEDPOC1 = 1
        self.assertEqual(self.param.ICM_SEDPOC1, 1)

    def test_ICM_SEDPOC2(self):
        self.param.ICM_SEDPOC2 = 1
        self.assertEqual(self.param.ICM_SEDPOC2, 1)

    def test_ICM_SEDPOC3(self):
        self.param.ICM_SEDPOC3 = 1
        self.assertEqual(self.param.ICM_SEDPOC3, 1)

    def test_ICM_EROH2S(self):
        self.param.ICM_EROH2S = 1
        self.assertEqual(self.param.ICM_EROH2S, 1)

    def test_ICM_EROLPOC(self):
        self.param.ICM_EROLPOC = 1
        self.assertEqual(self.param.ICM_EROLPOC, 1)

    def test_ICM_ERORPOC(self):
        self.param.ICM_ERORPOC = 1
        self.assertEqual(self.param.ICM_ERORPOC, 1)

    def test_ICM_DO_consumption(self):
        self.param.ICM_DO_consumption = 1
        self.assertEqual(self.param.ICM_DO_consumption, 1)

    def test_ICM_GP1(self):
        self.param.ICM_GP1 = 1
        self.assertEqual(self.param.ICM_GP1, 1)

    def test_ICM_GP2(self):
        self.param.ICM_GP2 = 1
        self.assertEqual(self.param.ICM_GP2, 1)

    def test_ICM_GP3(self):
        self.param.ICM_GP3 = 1
        self.assertEqual(self.param.ICM_GP3, 1)

    def test_ICM_1(self):
        self.param.ICM_1 = 1
        self.assertEqual(self.param.ICM_1, 1)

    def test_ICM_2(self):
        self.param.ICM_2 = 1
        self.assertEqual(self.param.ICM_2, 1)

    def test_ICM_3(self):
        self.param.ICM_3 = 1
        self.assertEqual(self.param.ICM_3, 1)

    def test_ICM_4(self):
        self.param.ICM_4 = 1
        self.assertEqual(self.param.ICM_4, 1)

    def test_ICM_5(self):
        self.param.ICM_5 = 1
        self.assertEqual(self.param.ICM_5, 1)

    def test_ICM_6(self):
        self.param.ICM_6 = 1
        self.assertEqual(self.param.ICM_6, 1)

    def test_ICM_7(self):
        self.param.ICM_7 = 1
        self.assertEqual(self.param.ICM_7, 1)

    def test_ICM_8(self):
        self.param.ICM_8 = 1
        self.assertEqual(self.param.ICM_8, 1)

    def test_ICM_9(self):
        self.param.ICM_9 = 1
        self.assertEqual(self.param.ICM_9, 1)

    def test_ICM_10(self):
        self.param.ICM_10 = 1
        self.assertEqual(self.param.ICM_10, 1)

    def test_ICM_11(self):
        self.param.ICM_11 = 1
        self.assertEqual(self.param.ICM_11, 1)

    def test_ICM_12(self):
        self.param.ICM_12 = 1
        self.assertEqual(self.param.ICM_12, 1)

    def test_ICM_13(self):
        self.param.ICM_13 = 1
        self.assertEqual(self.param.ICM_13, 1)

    def test_ICM_14(self):
        self.param.ICM_14 = 1
        self.assertEqual(self.param.ICM_14, 1)

    def test_ICM_15(self):
        self.param.ICM_15 = 1
        self.assertEqual(self.param.ICM_15, 1)

    def test_ICM_16(self):
        self.param.ICM_16 = 1
        self.assertEqual(self.param.ICM_16, 1)

    def test_ICM_17(self):
        self.param.ICM_17 = 1
        self.assertEqual(self.param.ICM_17, 1)

    def test_ICM_18(self):
        self.param.ICM_18 = 1
        self.assertEqual(self.param.ICM_18, 1)

    def test_ICM_19(self):
        self.param.ICM_19 = 1
        self.assertEqual(self.param.ICM_19, 1)

    def test_ICM_20(self):
        self.param.ICM_20 = 1
        self.assertEqual(self.param.ICM_20, 1)

    def test_ICM_21(self):
        self.param.ICM_21 = 1
        self.assertEqual(self.param.ICM_21, 1)

    def test_ICM_22(self):
        self.param.ICM_22 = 1
        self.assertEqual(self.param.ICM_22, 1)

    def test_ICM_23(self):
        self.param.ICM_23 = 1
        self.assertEqual(self.param.ICM_23, 1)

    def test_ICM_24(self):
        self.param.ICM_24 = 1
        self.assertEqual(self.param.ICM_24, 1)

    def test_ICM_25(self):
        self.param.ICM_25 = 1
        self.assertEqual(self.param.ICM_25, 1)

    def test_COS_1(self):
        self.param.COS_1 = 1
        self.assertEqual(self.param.COS_1, 1)

    def test_COS_2(self):
        self.param.COS_2 = 1
        self.assertEqual(self.param.COS_2, 1)

    def test_COS_3(self):
        self.param.COS_3 = 1
        self.assertEqual(self.param.COS_3, 1)

    def test_COS_4(self):
        self.param.COS_4 = 1
        self.assertEqual(self.param.COS_4, 1)

    def test_COS_5(self):
        self.param.COS_5 = 1
        self.assertEqual(self.param.COS_5, 1)

    def test_COS_6(self):
        self.param.COS_6 = 1
        self.assertEqual(self.param.COS_6, 1)

    def test_COS_7(self):
        self.param.COS_7 = 1
        self.assertEqual(self.param.COS_7, 1)

    def test_COS_8(self):
        self.param.COS_8 = 1
        self.assertEqual(self.param.COS_8, 1)

    def test_COS_9(self):
        self.param.COS_9 = 1
        self.assertEqual(self.param.COS_9, 1)

    def test_COS_10(self):
        self.param.COS_10 = 1
        self.assertEqual(self.param.COS_10, 1)

    def test_COS_11(self):
        self.param.COS_11 = 1
        self.assertEqual(self.param.COS_11, 1)

    def test_COS_12(self):
        self.param.COS_12 = 1
        self.assertEqual(self.param.COS_12, 1)

    def test_COS_13(self):
        self.param.COS_13 = 1
        self.assertEqual(self.param.COS_13, 1)

    def test_FIB_1(self):
        self.param.FIB_1 = 1
        self.assertEqual(self.param.FIB_1, 1)

    def test_SED2D_depth_change(self):
        self.param.SED2D_depth_change = 1
        self.assertEqual(self.param.SED2D_depth_change, 1)

    def test_SED2D_cflsed(self):
        self.param.SED2D_cflsed = 1
        self.assertEqual(self.param.SED2D_cflsed, 1)

    def test_SED2D_d50(self):
        self.param.SED2D_d50 = 1
        self.assertEqual(self.param.SED2D_d50, 1)

    def test_SED2D_total_transport(self):
        self.param.SED2D_total_transport = 1
        self.assertEqual(self.param.SED2D_total_transport, 1)

    def test_SED2D_susp_load(self):
        self.param.SED2D_susp_load = 1
        self.assertEqual(self.param.SED2D_susp_load, 1)

    def test_SED2D_bed_load(self):
        self.param.SED2D_bed_load = 1
        self.assertEqual(self.param.SED2D_bed_load, 1)

    def test_SED2D_average_transport(self):
        self.param.SED2D_average_transport = 1
        self.assertEqual(self.param.SED2D_average_transport, 1)

    def test_SED2D_bottom_slope(self):
        self.param.SED2D_bottom_slope = 1
        self.assertEqual(self.param.SED2D_bottom_slope, 1)

    def test_z0eq(self):
        self.param.z0eq = 1
        self.assertEqual(self.param.z0eq, 1)

    def test_z0cr2d(self):
        self.param.z0cr2d = 1
        self.assertEqual(self.param.z0cr2d, 1)

    def test_z0sw2d(self):
        self.param.z0sw2d = 1
        self.assertEqual(self.param.z0sw2d, 1)

    def test_z0wr2d(self):
        self.param.z0wr2d = 1
        self.assertEqual(self.param.z0wr2d, 1)

    def test_marsh_flag(self):
        self.param.marsh_flag = 1
        self.assertEqual(self.param.marsh_flag, 1)

    def test_ICE_velocity(self):
        self.param.ICE_velocity = 1
        self.assertEqual(self.param.ICE_velocity, 1)

    def test_ICE_strain_rate(self):
        self.param.ICE_strain_rate = 1
        self.assertEqual(self.param.ICE_strain_rate, 1)

    def test_ICE_net_heat_flux(self):
        self.param.ICE_net_heat_flux = 1
        self.assertEqual(self.param.ICE_net_heat_flux, 1)

    def test_ICE_fresh_water_flux(self):
        self.param.ICE_fresh_water_flux = 1
        self.assertEqual(self.param.ICE_fresh_water_flux, 1)

    def test_ICE_top_T(self):
        self.param.ICE_top_T = 1
        self.assertEqual(self.param.ICE_top_T, 1)

    def test_ICE_tracer_1(self):
        self.param.ICE_tracer_1 = 1
        self.assertEqual(self.param.ICE_tracer_1, 1)

    def test_ICE_tracer_2(self):
        self.param.ICE_tracer_2 = 1
        self.assertEqual(self.param.ICE_tracer_2, 1)

    def test_ICE_tracer_3(self):
        self.param.ICE_tracer_3 = 1
        self.assertEqual(self.param.ICE_tracer_3, 1)

    def test_ANA_air_pres_grad_x(self):
        self.param.ANA_air_pres_grad_x = 1
        self.assertEqual(self.param.ANA_air_pres_grad_x, 1)

    def test_ANA_air_pres_grad_y(self):
        self.param.ANA_air_pres_grad_y = 1
        self.assertEqual(self.param.ANA_air_pres_grad_y, 1)

    def test_ANA_tide_pot_grad_x(self):
        self.param.ANA_tide_pot_grad_x = 1
        self.assertEqual(self.param.ANA_tide_pot_grad_x, 1)

    def test_ANA_tide_pot_grad_y(self):
        self.param.ANA_tide_pot_grad_y = 1
        self.assertEqual(self.param.ANA_tide_pot_grad_y, 1)

    def test_ANA_hor_viscosity_x(self):
        self.param.ANA_hor_viscosity_x = 1
        self.assertEqual(self.param.ANA_hor_viscosity_x, 1)

    def test_ANA_hor_viscosity_y(self):
        self.param.ANA_hor_viscosity_y = 1
        self.assertEqual(self.param.ANA_hor_viscosity_y, 1)

    def test_ANA_bclinic_force_x(self):
        self.param.ANA_bclinic_force_x = 1
        self.assertEqual(self.param.ANA_bclinic_force_x, 1)

    def test_ANA_bclinic_force_y(self):
        self.param.ANA_bclinic_force_y = 1
        self.assertEqual(self.param.ANA_bclinic_force_y, 1)

    def test_ANA_vert_viscosity_x(self):
        self.param.ANA_vert_viscosity_x = 1
        self.assertEqual(self.param.ANA_vert_viscosity_x, 1)

    def test_ANA_vert_viscosity_y(self):
        self.param.ANA_vert_viscosity_y = 1
        self.assertEqual(self.param.ANA_vert_viscosity_y, 1)

    def test_ANA_mom_advection_x(self):
        self.param.ANA_mom_advection_x = 1
        self.assertEqual(self.param.ANA_mom_advection_x, 1)

    def test_ANA_mom_advection_y(self):
        self.param.ANA_mom_advection_y = 1
        self.assertEqual(self.param.ANA_mom_advection_y, 1)

    def test_ANA_Richardson(self):
        self.param.ANA_Richardson = 1
        self.assertEqual(self.param.ANA_Richardson, 1)

    def test_ANA_transport_min_dt_elem(self):
        self.param.ANA_transport_min_dt_elem = 1
        self.assertEqual(self.param.ANA_transport_min_dt_elem, 1)


if __name__ == '__main__':
    unittest.main()
