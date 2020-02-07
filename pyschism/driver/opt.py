class OPT:
    """ OPT group for Param namelist """

    def __init__(self, nml):
        self._nml = nml

    @property
    def ipre2(self):
        return self._nml['opt']['ipre2']

    @property
    def start_year(self):
        return int(self._nml['opt']['start_year'])

    @property
    def start_month(self):
        return int(self._nml['opt']['start_month'])

    @property
    def start_day(self):
        return int(self._nml['opt']['start_day'])

    @property
    def start_hour(self):
        return float(self._nml['opt']['start_hour'])

    @property
    def utc_start(self):
        return self._nml['opt']['utc_start']

    @property
    def ics(self):
        return self._nml['opt']['ics']

    @property
    def ihot(self):
        return self._nml['opt']['ihot']

    @property
    def ieos_type(self):
        return self._nml['opt']['ieos_type']

    @property
    def ieos_pres(self):
        return self._nml['opt']['ieos_pres']

    @property
    def eos_a(self):
        return self._nml['opt']['eos_a']

    @property
    def eos_b(self):
        return self._nml['opt']['eos_b']

    @property
    def nramp(self):
        return self._nml['opt']['nramp']

    @property
    def dramp(self):
        return self._nml['opt']['dramp']

    @property
    def nrampbc(self):
        return self._nml['opt']['nrampbc']

    @property
    def drampbc(self):
        return self._nml['opt']['drampbc']

    @property
    def iupwind_mom(self):
        return self._nml['opt']['iupwind_mom']

    @property
    def indvel(self):
        return self._nml['opt']['indvel']

    @property
    def ihorcon(self):
        return self._nml['opt']['ihorcon']

    @property
    def hvis_coef0(self):
        return self._nml['opt']['hvis_coef0']

    @property
    def ishapiro(self):
        return self._nml['opt']['ishapiro']

    @property
    def shapiro0(self):
        return self._nml['opt']['shapiro0']

    @property
    def niter_shap(self):
        return self._nml['opt']['niter_shap']

    @property
    def thetai(self):
        return self._nml['opt']['thetai']

    @property
    def icou_elfe_wwm(self):
        return self._nml['opt']['icou_elfe_wwm']

    @property
    def nstep_wwm(self):
        return self._nml['opt']['nstep_wwm']

    @property
    def iwbl(self):
        return self._nml['opt']['iwbl']

    @property
    def hmin_radstress(self):
        return self._nml['opt']['hmin_radstress']

    @property
    def nrampwafo(self):
        return self._nml['opt']['nrampwafo']

    @property
    def drampwafo(self):
        return self._nml['opt']['drampwafo']

    @property
    def turbinj(self):
        return self._nml['opt']['turbinj']

    @property
    def imm(self):
        return self._nml['opt']['imm']

    @property
    def ibdef(self):
        return self._nml['opt']['ibdef']

    @property
    def slam0(self):
        return self._nml['opt']['slam0']

    @property
    def sfea0(self):
        return self._nml['opt']['sfea0']

    @property
    def iunder_deep(self):
        return self._nml['opt']['iunder_deep']

    @property
    def h1_bcc(self):
        return self._nml['opt']['h1_bcc']

    @property
    def h2_bcc(self):
        return self._nml['opt']['h2_bcc']

    @property
    def hw_depth(self):
        return self._nml['opt']['hw_depth']

    @property
    def hw_ratio(self):
        return self._nml['opt']['hw_ratio']

    @property
    def ihydraulics(self):
        return self._nml['opt']['ihydraulics']

    @property
    def if_source(self):
        return self._nml['opt']['if_source']

    @property
    def nramp_ss(self):
        return self._nml['opt']['nramp_ss']

    @property
    def dramp_ss(self):
        return self._nml['opt']['dramp_ss']

    @property
    def ihdif(self):
        return self._nml['opt']['ihdif']

    @property
    def nchi(self):
        return self._nml['opt']['nchi']

    @property
    def dzb_min(self):
        return self._nml['opt']['dzb_min']

    @property
    def dzb_decay(self):
        return self._nml['opt']['dzb_decay']

    @property
    def hmin_man(self):
        return self._nml['opt']['hmin_man']

    @property
    def ncor(self):
        return self._nml['opt']['ncor']

    @property
    def rlatitude(self):
        return self._nml['opt']['rlatitude']

    @property
    def coricoef(self):
        return self._nml['opt']['coricoef']

    @property
    def ic_elev(self):
        return self._nml['opt']['ic_elev']

    @property
    def nramp_elev(self):
        return self._nml['opt']['nramp_elev']

    @property
    def inv_atm_bnd(self):
        return self._nml['opt']['inv_atm_bnd']

    @property
    def prmsl_ref(self):
        return self._nml['opt']['prmsl_ref']

    @property
    def flag_ic(self):
        return self._nml['opt']['flag_ic']

    @property
    def gen_wsett(self):
        return self._nml['opt']['gen_wsett']

    @property
    def ibcc_mean(self):
        return self._nml['opt']['ibcc_mean']

    @property
    def rmaxvel(self):
        return self._nml['opt']['rmaxvel']

    @property
    def velmin_btrack(self):
        return self._nml['opt']['velmin_btrack']

    @property
    def btrack_nudge(self):
        return self._nml['opt']['btrack_nudge']

    @property
    def ibtrack_openbnd(self):
        return self._nml['opt']['ibtrack_openbnd']

    @property
    def ihhat(self):
        return self._nml['opt']['ihhat']

    @property
    def inunfl(self):
        return self._nml['opt']['inunfl']

    @property
    def h0(self):
        return self._nml['opt']['h0']

    @property
    def shorewafo(self):
        return self._nml['opt']['shorewafo']

    @property
    def moitn0(self):
        return self._nml['opt']['moitn0']

    @property
    def mxitn0(self):
        return self._nml['opt']['mxitn0']

    @property
    def rtol0(self):
        return self._nml['opt']['rtol0']

    @property
    def nadv(self):
        return self._nml['opt']['nadv']

    @property
    def dtb_max(self):
        return self._nml['opt']['dtb_max']

    @property
    def dtb_min(self):
        return self._nml['opt']['dtb_min']

    @property
    def inter_mom(self):
        return self._nml['opt']['inter_mom']

    @property
    def kr_co(self):
        return self._nml['opt']['kr_co']

    @property
    def itr_met(self):
        return self._nml['opt']['itr_met']

    @property
    def h_tvd(self):
        return self._nml['opt']['h_tvd']

    @property
    def eps1_tvd_imp(self):
        return self._nml['opt']['eps1_tvd_imp']

    @property
    def eps2_tvd_imp(self):
        return self._nml['opt']['eps2_tvd_imp']

    @property
    def ip_weno(self):
        return self._nml['opt']['ip_weno']

    @property
    def courant_weno(self):
        return self._nml['opt']['courant_weno']

    @property
    def nquad(self):
        return self._nml['opt']['nquad']

    @property
    def ntd_weno(self):
        return self._nml['opt']['ntd_weno']

    @property
    def epsilon1(self):
        return self._nml['opt']['epsilon1']

    @property
    def epsilon2(self):
        return self._nml['opt']['epsilon2']

    @property
    def i_prtnftl_weno(self):
        return self._nml['opt']['i_prtnftl_weno']

    @property
    def epsilon3(self):
        return self._nml['opt']['epsilon3']

    @property
    def ielad_weno(self):
        return self._nml['opt']['ielad_weno']

    @property
    def small_elad(self):
        return self._nml['opt']['small_elad']

    @property
    def nws(self):
        return self._nml['opt']['nws']

    @property
    def wtiminc(self):
        return self._nml['opt']['wtiminc']

    @property
    def nrampwind(self):
        return self._nml['opt']['nrampwind']

    @property
    def drampwind(self):
        return self._nml['opt']['drampwind']

    @property
    def iwindoff(self):
        return self._nml['opt']['iwindoff']

    @property
    def iwind_form(self):
        return self._nml['opt']['iwind_form']

    @property
    def impose_net_flux(self):
        return self._nml['opt']['impose_net_flux']

    @property
    def ihconsv(self):
        return self._nml['opt']['ihconsv']

    @property
    def isconsv(self):
        return self._nml['opt']['isconsv']

    @property
    def itur(self):
        return self._nml['opt']['itur']

    @property
    def dfv0(self):
        return self._nml['opt']['dfv0']

    @property
    def dfh0(self):
        return self._nml['opt']['dfh0']

    @property
    def mid(self):
        return self._nml['opt']['mid']

    @property
    def stab(self):
        return self._nml['opt']['stab']

    @property
    def xlsc0(self):
        return self._nml['opt']['xlsc0']

    @property
    def inu_elev(self):
        return self._nml['opt']['inu_elev']

    @property
    def inu_uv(self):
        return self._nml['opt']['inu_uv']

    @property
    def inu_tr(self):
        return self._nml['opt']['inu_tr']

    @property
    def vnh1(self):
        return self._nml['opt']['vnh1']

    @property
    def vnf1(self):
        return self._nml['opt']['vnf1']

    @property
    def vnh2(self):
        return self._nml['opt']['vnh2']

    @property
    def vnf2(self):
        return self._nml['opt']['vnf2']

    @property
    def step_nu_tr(self):
        return self._nml['opt']['step_nu_tr']

    @property
    def h_bcc1(self):
        return self._nml['opt']['h_bcc1']

    @property
    def s1_mxnbt(self):
        return self._nml['opt']['s1_mxnbt']

    @property
    def s2_mxnbt(self):
        return self._nml['opt']['s2_mxnbt']

    @property
    def iharind(self):
        return self._nml['opt']['iharind']

    @property
    def iflux(self):
        return self._nml['opt']['iflux']

    @property
    def izonal5(self):
        return self._nml['opt']['izonal5']

    @property
    def ibtrack_test(self):
        return self._nml['opt']['ibtrack_test']

    @property
    def irouse_test(self):
        return self._nml['opt']['irouse_test']

    @property
    def flag_fib(self):
        return self._nml['opt']['flag_fib']

    @property
    def slr_rate(self):
        return self._nml['opt']['slr_rate']

    @property
    def isav(self):
        return self._nml['opt']['isav']

    @property
    def sav_cd(self):
        return self._nml['opt']['sav_cd']

    @property
    def nstep_ice(self):
        return self._nml['opt']['nstep_ice']

    @property
    def level_age(self):
        return self._nml['opt']['level_age']

    @property
    def rearth_pole(self):
        return self._nml['opt']['rearth_pole']

    @property
    def rearth_eq(self):
        return self._nml['opt']['rearth_eq']

    @property
    def shw(self):
        return self._nml['opt']['shw']

    @property
    def rho0(self):
        return self._nml['opt']['rho0']

    @property
    def vclose_surf_frac(self):
        return self._nml['opt']['vclose_surf_frac']

    @property
    def iadjust_mass_consv0(self):
        return self._nml['opt']['iadjust_mass_consv0']

    @ipre2.setter
    def ipre2(self, ipre2):
        self._nml['opt']['ipre2'] = ipre2

    @start_year.setter
    def start_year(self, start_year):
        self._nml['opt']['start_year'] = start_year

    @start_month.setter
    def start_month(self, start_month):
        assert isinstance(start_month, int)
        self._nml['opt']['start_month'] = start_month

    @start_day.setter
    def start_day(self, start_day):
        self._nml['opt']['start_day'] = start_day

    @start_hour.setter
    def start_hour(self, start_hour):
        self._nml['opt']['start_hour'] = start_hour

    @utc_start.setter
    def utc_start(self, utc_start):
        self._nml['opt']['utc_start'] = utc_start

    @ics.setter
    def ics(self, ics):
        self._nml['opt']['ics'] = ics

    @ihot.setter
    def ihot(self, ihot):
        self._nml['opt']['ihot'] = ihot

    @ieos_type.setter
    def ieos_type(self, ieos_type):
        self._nml['opt']['ieos_type'] = ieos_type

    @ieos_pres.setter
    def ieos_pres(self, ieos_pres):
        self._nml['opt']['ieos_pres'] = ieos_pres

    @eos_a.setter
    def eos_a(self, eos_a):
        self._nml['opt']['eos_a'] = eos_a

    @eos_b.setter
    def eos_b(self, eos_b):
        self._nml['opt']['eos_b'] = eos_b

    @nramp.setter
    def nramp(self, nramp):
        self._nml['opt']['nramp'] = nramp

    @dramp.setter
    def dramp(self, dramp):
        self._nml['opt']['dramp'] = dramp

    @nrampbc.setter
    def nrampbc(self, nrampbc):
        self._nml['opt']['nrampbc'] = nrampbc

    @drampbc.setter
    def drampbc(self, drampbc):
        self._nml['opt']['drampbc'] = drampbc

    @iupwind_mom.setter
    def iupwind_mom(self, iupwind_mom):
        self._nml['opt']['iupwind_mom'] = iupwind_mom

    @indvel.setter
    def indvel(self, indvel):
        self._nml['opt']['indvel'] = indvel

    @ihorcon.setter
    def ihorcon(self, ihorcon):
        self._nml['opt']['ihorcon'] = ihorcon

    @hvis_coef0.setter
    def hvis_coef0(self, hvis_coef0):
        self._nml['opt']['hvis_coef0'] = hvis_coef0

    @ishapiro.setter
    def ishapiro(self, ishapiro):
        self._nml['opt']['ishapiro'] = ishapiro

    @shapiro0.setter
    def shapiro0(self, shapiro0):
        self._nml['opt']['shapiro0'] = shapiro0

    @niter_shap.setter
    def niter_shap(self, niter_shap):
        self._nml['opt']['niter_shap'] = niter_shap

    @thetai.setter
    def thetai(self, thetai):
        self._nml['opt']['thetai'] = thetai

    @icou_elfe_wwm.setter
    def icou_elfe_wwm(self, icou_elfe_wwm):
        self._nml['opt']['icou_elfe_wwm'] = icou_elfe_wwm

    @nstep_wwm.setter
    def nstep_wwm(self, nstep_wwm):
        self._nml['opt']['nstep_wwm'] = nstep_wwm

    @iwbl.setter
    def iwbl(self, iwbl):
        self._nml['opt']['iwbl'] = iwbl

    @hmin_radstress.setter
    def hmin_radstress(self, hmin_radstress):
        self._nml['opt']['hmin_radstress'] = hmin_radstress

    @nrampwafo.setter
    def nrampwafo(self, nrampwafo):
        self._nml['opt']['nrampwafo'] = nrampwafo

    @drampwafo.setter
    def drampwafo(self, drampwafo):
        self._nml['opt']['drampwafo'] = drampwafo

    @turbinj.setter
    def turbinj(self, turbinj):
        self._nml['opt']['turbinj'] = turbinj

    @imm.setter
    def imm(self, imm):
        self._nml['opt']['imm'] = imm

    @ibdef.setter
    def ibdef(self, ibdef):
        self._nml['opt']['ibdef'] = ibdef

    @slam0.setter
    def slam0(self, slam0):
        self._nml['opt']['slam0'] = slam0

    @sfea0.setter
    def sfea0(self, sfea0):
        self._nml['opt']['sfea0'] = sfea0

    @iunder_deep.setter
    def iunder_deep(self, iunder_deep):
        self._nml['opt']['iunder_deep'] = iunder_deep

    @h1_bcc.setter
    def h1_bcc(self, h1_bcc):
        self._nml['opt']['h1_bcc'] = h1_bcc

    @h2_bcc.setter
    def h2_bcc(self, h2_bcc):
        self._nml['opt']['h2_bcc'] = h2_bcc

    @hw_depth.setter
    def hw_depth(self, hw_depth):
        self._nml['opt']['hw_depth'] = hw_depth

    @hw_ratio.setter
    def hw_ratio(self, hw_ratio):
        self._nml['opt']['hw_ratio'] = hw_ratio

    @ihydraulics.setter
    def ihydraulics(self, ihydraulics):
        self._nml['opt']['ihydraulics'] = ihydraulics

    @if_source.setter
    def if_source(self, if_source):
        self._nml['opt']['if_source'] = if_source

    @nramp_ss.setter
    def nramp_ss(self, nramp_ss):
        self._nml['opt']['nramp_ss'] = nramp_ss

    @dramp_ss.setter
    def dramp_ss(self, dramp_ss):
        self._nml['opt']['dramp_ss'] = dramp_ss

    @ihdif.setter
    def ihdif(self, ihdif):
        self._nml['opt']['ihdif'] = ihdif

    @nchi.setter
    def nchi(self, nchi):
        self._nml['opt']['nchi'] = nchi

    @dzb_min.setter
    def dzb_min(self, dzb_min):
        self._nml['opt']['dzb_min'] = dzb_min

    @dzb_decay.setter
    def dzb_decay(self, dzb_decay):
        self._nml['opt']['dzb_decay'] = dzb_decay

    @hmin_man.setter
    def hmin_man(self, hmin_man):
        self._nml['opt']['hmin_man'] = hmin_man

    @ncor.setter
    def ncor(self, ncor):
        self._nml['opt']['ncor'] = ncor

    @rlatitude.setter
    def rlatitude(self, rlatitude):
        self._nml['opt']['rlatitude'] = rlatitude

    @coricoef.setter
    def coricoef(self, coricoef):
        self._nml['opt']['coricoef'] = coricoef

    @ic_elev.setter
    def ic_elev(self, ic_elev):
        self._nml['opt']['ic_elev'] = ic_elev

    @nramp_elev.setter
    def nramp_elev(self, nramp_elev):
        self._nml['opt']['nramp_elev'] = nramp_elev

    @inv_atm_bnd.setter
    def inv_atm_bnd(self, inv_atm_bnd):
        self._nml['opt']['inv_atm_bnd'] = inv_atm_bnd

    @prmsl_ref.setter
    def prmsl_ref(self, prmsl_ref):
        self._nml['opt']['prmsl_ref'] = prmsl_ref

    @flag_ic.setter
    def flag_ic(self, flag_ic):
        self._nml['opt']['flag_ic'] = flag_ic

    @gen_wsett.setter
    def gen_wsett(self, gen_wsett):
        self._nml['opt']['gen_wsett'] = gen_wsett

    @ibcc_mean.setter
    def ibcc_mean(self, ibcc_mean):
        self._nml['opt']['ibcc_mean'] = ibcc_mean

    @rmaxvel.setter
    def rmaxvel(self, rmaxvel):
        self._nml['opt']['rmaxvel'] = rmaxvel

    @velmin_btrack.setter
    def velmin_btrack(self, velmin_btrack):
        self._nml['opt']['velmin_btrack'] = velmin_btrack

    @btrack_nudge.setter
    def btrack_nudge(self, btrack_nudge):
        self._nml['opt']['btrack_nudge'] = btrack_nudge

    @ibtrack_openbnd.setter
    def ibtrack_openbnd(self, ibtrack_openbnd):
        self._nml['opt']['ibtrack_openbnd'] = ibtrack_openbnd

    @ihhat.setter
    def ihhat(self, ihhat):
        self._nml['opt']['ihhat'] = ihhat

    @inunfl.setter
    def inunfl(self, inunfl):
        self._nml['opt']['inunfl'] = inunfl

    @h0.setter
    def h0(self, h0):
        self._nml['opt']['h0'] = h0

    @shorewafo.setter
    def shorewafo(self, shorewafo):
        self._nml['opt']['shorewafo'] = shorewafo

    @moitn0.setter
    def moitn0(self, moitn0):
        self._nml['opt']['moitn0'] = moitn0

    @mxitn0.setter
    def mxitn0(self, mxitn0):
        self._nml['opt']['mxitn0'] = mxitn0

    @rtol0.setter
    def rtol0(self, rtol0):
        self._nml['opt']['rtol0'] = rtol0

    @nadv.setter
    def nadv(self, nadv):
        self._nml['opt']['nadv'] = nadv

    @dtb_max.setter
    def dtb_max(self, dtb_max):
        self._nml['opt']['dtb_max'] = dtb_max

    @dtb_min.setter
    def dtb_min(self, dtb_min):
        self._nml['opt']['dtb_min'] = dtb_min

    @inter_mom.setter
    def inter_mom(self, inter_mom):
        self._nml['opt']['inter_mom'] = inter_mom

    @kr_co.setter
    def kr_co(self, kr_co):
        self._nml['opt']['kr_co'] = kr_co

    @itr_met.setter
    def itr_met(self, itr_met):
        self._nml['opt']['itr_met'] = itr_met

    @h_tvd.setter
    def h_tvd(self, h_tvd):
        self._nml['opt']['h_tvd'] = h_tvd

    @eps1_tvd_imp.setter
    def eps1_tvd_imp(self, eps1_tvd_imp):
        self._nml['opt']['eps1_tvd_imp'] = eps1_tvd_imp

    @eps2_tvd_imp.setter
    def eps2_tvd_imp(self, eps2_tvd_imp):
        self._nml['opt']['eps2_tvd_imp'] = eps2_tvd_imp

    @ip_weno.setter
    def ip_weno(self, ip_weno):
        self._nml['opt']['ip_weno'] = ip_weno

    @courant_weno.setter
    def courant_weno(self, courant_weno):
        self._nml['opt']['courant_weno'] = courant_weno

    @nquad.setter
    def nquad(self, nquad):
        self._nml['opt']['nquad'] = nquad

    @ntd_weno.setter
    def ntd_weno(self, ntd_weno):
        self._nml['opt']['ntd_weno'] = ntd_weno

    @epsilon1.setter
    def epsilon1(self, epsilon1):
        self._nml['opt']['epsilon1'] = epsilon1

    @epsilon2.setter
    def epsilon2(self, epsilon2):
        self._nml['opt']['epsilon2'] = epsilon2

    @i_prtnftl_weno.setter
    def i_prtnftl_weno(self, i_prtnftl_weno):
        self._nml['opt']['i_prtnftl_weno'] = i_prtnftl_weno

    @epsilon3.setter
    def epsilon3(self, epsilon3):
        self._nml['opt']['epsilon3'] = epsilon3

    @ielad_weno.setter
    def ielad_weno(self, ielad_weno):
        self._nml['opt']['ielad_weno'] = ielad_weno

    @small_elad.setter
    def small_elad(self, small_elad):
        self._nml['opt']['small_elad'] = small_elad

    @nws.setter
    def nws(self, nws):
        self._nml['opt']['nws'] = nws

    @wtiminc.setter
    def wtiminc(self, wtiminc):
        self._nml['opt']['wtiminc'] = wtiminc

    @nrampwind.setter
    def nrampwind(self, nrampwind):
        self._nml['opt']['nrampwind'] = nrampwind

    @drampwind.setter
    def drampwind(self, drampwind):
        self._nml['opt']['drampwind'] = drampwind

    @iwindoff.setter
    def iwindoff(self, iwindoff):
        self._nml['opt']['iwindoff'] = iwindoff

    @iwind_form.setter
    def iwind_form(self, iwind_form):
        self._nml['opt']['iwind_form'] = iwind_form

    @impose_net_flux.setter
    def impose_net_flux(self, impose_net_flux):
        self._nml['opt']['impose_net_flux'] = impose_net_flux

    @ihconsv.setter
    def ihconsv(self, ihconsv):
        self._nml['opt']['ihconsv'] = ihconsv

    @isconsv.setter
    def isconsv(self, isconsv):
        self._nml['opt']['isconsv'] = isconsv

    @itur.setter
    def itur(self, itur):
        self._nml['opt']['itur'] = itur

    @dfv0.setter
    def dfv0(self, dfv0):
        self._nml['opt']['dfv0'] = dfv0

    @dfh0.setter
    def dfh0(self, dfh0):
        self._nml['opt']['dfh0'] = dfh0

    @mid.setter
    def mid(self, mid):
        self._nml['opt']['mid'] = mid

    @stab.setter
    def stab(self, stab):
        self._nml['opt']['stab'] = stab

    @xlsc0.setter
    def xlsc0(self, xlsc0):
        self._nml['opt']['xlsc0'] = xlsc0

    @inu_elev.setter
    def inu_elev(self, inu_elev):
        self._nml['opt']['inu_elev'] = inu_elev

    @inu_uv.setter
    def inu_uv(self, inu_uv):
        self._nml['opt']['inu_uv'] = inu_uv

    @inu_tr.setter
    def inu_tr(self, inu_tr):
        self._nml['opt']['inu_tr'] = inu_tr

    @vnh1.setter
    def vnh1(self, vnh1):
        self._nml['opt']['vnh1'] = vnh1

    @vnf1.setter
    def vnf1(self, vnf1):
        self._nml['opt']['vnf1'] = vnf1

    @vnh2.setter
    def vnh2(self, vnh2):
        self._nml['opt']['vnh2'] = vnh2

    @vnf2.setter
    def vnf2(self, vnf2):
        self._nml['opt']['vnf2'] = vnf2

    @step_nu_tr.setter
    def step_nu_tr(self, step_nu_tr):
        self._nml['opt']['step_nu_tr'] = step_nu_tr

    @h_bcc1.setter
    def h_bcc1(self, h_bcc1):
        self._nml['opt']['h_bcc1'] = h_bcc1

    @s1_mxnbt.setter
    def s1_mxnbt(self, s1_mxnbt):
        self._nml['opt']['s1_mxnbt'] = s1_mxnbt

    @s2_mxnbt.setter
    def s2_mxnbt(self, s2_mxnbt):
        self._nml['opt']['s2_mxnbt'] = s2_mxnbt

    @iharind.setter
    def iharind(self, iharind):
        self._nml['opt']['iharind'] = iharind

    @iflux.setter
    def iflux(self, iflux):
        self._nml['opt']['iflux'] = iflux

    @izonal5.setter
    def izonal5(self, izonal5):
        self._nml['opt']['izonal5'] = izonal5

    @ibtrack_test.setter
    def ibtrack_test(self, ibtrack_test):
        self._nml['opt']['ibtrack_test'] = ibtrack_test

    @irouse_test.setter
    def irouse_test(self, irouse_test):
        self._nml['opt']['irouse_test'] = irouse_test

    @flag_fib.setter
    def flag_fib(self, flag_fib):
        self._nml['opt']['flag_fib'] = flag_fib

    @slr_rate.setter
    def slr_rate(self, slr_rate):
        self._nml['opt']['slr_rate'] = slr_rate

    @isav.setter
    def isav(self, isav):
        self._nml['opt']['isav'] = isav

    @sav_cd.setter
    def sav_cd(self, sav_cd):
        self._nml['opt']['sav_cd'] = sav_cd

    @nstep_ice.setter
    def nstep_ice(self, nstep_ice):
        self._nml['opt']['nstep_ice'] = nstep_ice

    @level_age.setter
    def level_age(self, level_age):
        self._nml['opt']['level_age'] = level_age

    @rearth_pole.setter
    def rearth_pole(self, rearth_pole):
        self._nml['opt']['rearth_pole'] = rearth_pole

    @rearth_eq.setter
    def rearth_eq(self, rearth_eq):
        self._nml['opt']['rearth_eq'] = rearth_eq

    @shw.setter
    def shw(self, shw):
        self._nml['opt']['shw'] = shw

    @rho0.setter
    def rho0(self, rho0):
        self._nml['opt']['rho0'] = rho0

    @vclose_surf_frac.setter
    def vclose_surf_frac(self, vclose_surf_frac):
        self._nml['opt']['vclose_surf_frac'] = vclose_surf_frac

    @iadjust_mass_consv0.setter
    def iadjust_mass_consv0(self, iadjust_mass_consv0):
        self._nml['opt']['iadjust_mass_consv0'] = iadjust_mass_consv0
