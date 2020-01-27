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
        self.param.core.ipre = 0
        self.assertEqual(self.param.core.ipre, 0)

    def test_core_ibc(self):
        self.param.core.ibc = 0
        self.assertEqual(self.param.core.ibc, 0)

    def test_core_ibtp(self):
        self.param.core.ibtp = 1
        self.assertEqual(self.param.core.ibtp, 1)

    def test_core_rnday(self):
        self.param.core.rnday = 30
        self.assertEqual(self.param.core.rnday, 30)

    def test_core_dt(self):
        self.param.core.dt = 100.0
        self.assertEqual(self.param.core.dt, 100.0)

    def test_core_msc2(self):
        self.param.core.msc2 = 24
        self.assertEqual(self.param.core.msc2, 24)

    def test_core_mdc2(self):
        self.param.core.mdc2 = 30
        self.assertEqual(self.param.core.mdc2, 30)

    def test_core_ntracer_gen(self):
        self.param.core.ntracer_gen = 2
        self.assertEqual(self.param.core.ntracer_gen, 2)

    def test_core_ntracer_age(self):
        self.param.core.ntracer_age = 4
        self.assertEqual(self.param.core.ntracer_age, 4)

    def test_core_sed_class(self):
        self.param.core.sed_class = 5
        self.assertEqual(self.param.core.sed_class, 5)

    def test_core_eco_class(self):
        self.param.core.eco_class = 27
        self.assertEqual(self.param.core.eco_class, 27)

    def test_core_nspool(self):
        self.param.core.nspool = 36
        self.assertEqual(self.param.core.nspool, 36)

    def test_core_ihfskip(self):
        self.param.core.ihfskip = 864
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


if __name__ == '__main__':
    unittest.main()
