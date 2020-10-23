from datetime import datetime, timedelta, timezone
import pathlib
from typing import Union

import f90nml  # type: ignore[import]

PARAM_TEMPLATE = pathlib.Path(__file__).parent / 'param.nml.template'
PARAM_DEFAULTS = f90nml.read(PARAM_TEMPLATE)['opt']


class OPT:
    """ Provides error checking implementation for OPT group """

    def __init__(self):
        self.__opt: dict = {}
        for key, value in PARAM_DEFAULTS.items():
            if isinstance(value, list):
                self.__opt[key] = len(value)*[0]
            else:
                self.__opt[key] = None
        self.__start_date = None

    def __getitem__(self, key):
        return self.__opt[key]

    def __setitem__(self, key, value):
        self.__opt[key] = value

    def __iter__(self):
        for key, value in self.__opt.items():
            yield key, value

    @property
    def start_date(self):
        return self.__start_date

    @start_date.setter
    def start_date(self, start_date: datetime):
        if start_date.tzinfo is None or \
                start_date.tzinfo.utcoffset(start_date) is None:
            start_date = start_date.replace(tzinfo=timezone(timedelta(0)))
        self['start_year'] = start_date.year
        self['start_month'] = start_date.month
        self['start_day'] = start_date.day
        self['start_hour'] = start_date.hour
        self['start_hour'] += start_date.minute / 60.
        self['utc_start'] = -start_date.utcoffset().total_seconds() / 3600  # type: ignore[union-attr]  # noqa: E501
        # get rid of "negative" zero
        self['utc_start'] = +0. if self['utc_start'] == -0. \
            else self['utc_start']
        self.__start_date = start_date

    @property
    def ipre2(self):
        return self['ipre2']

    @ipre2.setter
    def ipre2(self, ipre2):
        self['ipre2'] = ipre2

    @property
    def start_year(self):
        return int(self['start_year'])

    @start_year.setter
    def start_year(self, start_year):
        raise AttributeError(
            'start_year must be set using start_date attribute.')

    @property
    def start_month(self):
        return int(self['start_month'])

    @start_month.setter
    def start_month(self, start_month):
        raise AttributeError(
            'start_month must be set using start_date attribute.')

    @property
    def start_day(self):
        return int(self['start_day'])

    @start_day.setter
    def start_day(self, start_day):
        raise AttributeError(
            'start_day must be set using start_date attribute.')

    @property
    def start_hour(self):
        return float(self['start_hour'])

    @start_hour.setter
    def start_hour(self, start_hour):
        raise AttributeError(
            'start_hour must be set using start_date attribute.')

    @property
    def utc_start(self):
        return self['utc_start']

    @utc_start.setter
    def utc_start(self, utc_start):
        raise AttributeError(
            'utc_start must be set using start_date attribute.')

    @property
    def ics(self):
        return self['ics']

    @ics.setter
    def ics(self, ics):
        self['ics'] = ics

    @property
    def ihot(self):
        return self['ihot']

    @ihot.setter
    def ihot(self, ihot):
        self['ihot'] = ihot

    @property
    def ieos_type(self):
        return self['ieos_type']

    @ieos_type.setter
    def ieos_type(self, ieos_type):
        self['ieos_type'] = ieos_type

    @property
    def ieos_pres(self):
        return self['ieos_pres']

    @ieos_pres.setter
    def ieos_pres(self, ieos_pres):
        self['ieos_pres'] = ieos_pres

    @property
    def eos_a(self):
        return self['eos_a']

    @eos_a.setter
    def eos_a(self, eos_a):
        self['eos_a'] = eos_a

    @property
    def eos_b(self):
        return self['eos_b']

    @eos_b.setter
    def eos_b(self, eos_b):
        self['eos_b'] = eos_b

    @property
    def nramp(self):
        return self['nramp']

    # No setter for nramp

    @property
    def dramp(self):
        return timedelta(days=self._dramp)

    @dramp.setter
    def dramp(self, dramp: Union[None, float, int, timedelta]):

        if dramp is None:
            self['dramp'] = None
            self['nramp'] = None
            return

        elif isinstance(dramp, timedelta):
            dramp = dramp.total_seconds() / 86400.

        self['dramp'] = float(dramp)
        self['nramp'] = 1

    @property
    def nrampbc(self):
        return self['nrampbc']

    # No setter for nrampbc

    @property
    def drampbc(self):
        return self['drampbc']

    @drampbc.setter
    def drampbc(self, drampbc: Union[None, float, int, timedelta]):

        if drampbc is None:
            self['drampbc'] = None
            self['nrampbc'] = None
            return

        elif isinstance(drampbc, timedelta):
            dramp = drampbc.total_seconds() / 86400.

        self['drampbc'] = float(dramp)
        self['nrampbc'] = 1

    @property
    def iupwind_mom(self):
        return self['iupwind_mom']

    @iupwind_mom.setter
    def iupwind_mom(self, iupwind_mom):
        self['iupwind_mom'] = iupwind_mom

    @property
    def indvel(self):
        return self['indvel']

    @indvel.setter
    def indvel(self, indvel):
        self['indvel'] = indvel

    @property
    def ihorcon(self):
        return self['ihorcon']

    @ihorcon.setter
    def ihorcon(self, ihorcon):
        self['ihorcon'] = ihorcon

    @property
    def hvis_coef0(self):
        return self['hvis_coef0']

    @hvis_coef0.setter
    def hvis_coef0(self, hvis_coef0):
        self['hvis_coef0'] = hvis_coef0

    @property
    def ishapiro(self):
        return self['ishapiro']

    @ishapiro.setter
    def ishapiro(self, ishapiro):
        self['ishapiro'] = ishapiro

    @property
    def shapiro0(self):
        return self['shapiro0']

    @shapiro0.setter
    def shapiro0(self, shapiro0):
        self['shapiro0'] = shapiro0

    @property
    def niter_shap(self):
        return self['niter_shap']

    @niter_shap.setter
    def niter_shap(self, niter_shap):
        self['niter_shap'] = niter_shap

    @property
    def thetai(self):
        return self['thetai']

    @thetai.setter
    def thetai(self, thetai):
        self['thetai'] = thetai

    @property
    def icou_elfe_wwm(self):
        return self['icou_elfe_wwm']

    @icou_elfe_wwm.setter
    def icou_elfe_wwm(self, icou_elfe_wwm):
        self['icou_elfe_wwm'] = icou_elfe_wwm

    @property
    def nstep_wwm(self):
        return self['nstep_wwm']

    @nstep_wwm.setter
    def nstep_wwm(self, nstep_wwm):
        self['nstep_wwm'] = nstep_wwm

    @property
    def iwbl(self):
        return self['iwbl']

    @iwbl.setter
    def iwbl(self, iwbl):
        self['iwbl'] = iwbl

    @property
    def hmin_radstress(self):
        return self['hmin_radstress']

    @hmin_radstress.setter
    def hmin_radstress(self, hmin_radstress):
        self['hmin_radstress'] = hmin_radstress

    @property
    def nrampwafo(self):
        return self['nrampwafo']

    @nrampwafo.setter
    def nrampwafo(self, nrampwafo):
        self['nrampwafo'] = nrampwafo

    @property
    def drampwafo(self):
        return self['drampwafo']

    @drampwafo.setter
    def drampwafo(self, drampwafo):
        self['drampwafo'] = drampwafo

    @property
    def turbinj(self):
        return self['turbinj']

    @turbinj.setter
    def turbinj(self, turbinj):
        self['turbinj'] = turbinj

    @property
    def imm(self):
        return self['imm']

    @imm.setter
    def imm(self, imm):
        self['imm'] = imm

    @property
    def ibdef(self):
        return self['ibdef']

    @ibdef.setter
    def ibdef(self, ibdef):
        self['ibdef'] = ibdef

    @property
    def slam0(self):
        return self['slam0']

    @slam0.setter
    def slam0(self, slam0):
        self['slam0'] = slam0

    @property
    def sfea0(self):
        return self['sfea0']

    @sfea0.setter
    def sfea0(self, sfea0):
        self['sfea0'] = sfea0

    @property
    def iunder_deep(self):
        return self['iunder_deep']

    @iunder_deep.setter
    def iunder_deep(self, iunder_deep):
        self['iunder_deep'] = iunder_deep

    @property
    def h1_bcc(self):
        return self['h1_bcc']

    @h1_bcc.setter
    def h1_bcc(self, h1_bcc):
        self['h1_bcc'] = h1_bcc

    @property
    def h2_bcc(self):
        return self['h2_bcc']

    @h2_bcc.setter
    def h2_bcc(self, h2_bcc):
        self['h2_bcc'] = h2_bcc

    @property
    def hw_depth(self):
        return self['hw_depth']

    @hw_depth.setter
    def hw_depth(self, hw_depth):
        self['hw_depth'] = hw_depth

    @property
    def hw_ratio(self):
        return self['hw_ratio']

    @hw_ratio.setter
    def hw_ratio(self, hw_ratio):
        self['hw_ratio'] = hw_ratio

    @property
    def ihydraulics(self):
        return self['ihydraulics']

    @ihydraulics.setter
    def ihydraulics(self, ihydraulics):
        self['ihydraulics'] = ihydraulics

    @property
    def if_source(self):
        return self['if_source']

    @if_source.setter
    def if_source(self, if_source):
        self['if_source'] = if_source

    @property
    def nramp_ss(self):
        return self['nramp_ss']

    @nramp_ss.setter
    def nramp_ss(self, nramp_ss):
        self['nramp_ss'] = nramp_ss

    @property
    def dramp_ss(self):
        return self['dramp_ss']

    @dramp_ss.setter
    def dramp_ss(self, dramp_ss):
        self['dramp_ss'] = dramp_ss

    @property
    def ihdif(self):
        return self['ihdif']

    @ihdif.setter
    def ihdif(self, ihdif):
        self['ihdif'] = ihdif

    @property
    def nchi(self):
        return self['nchi']

    @nchi.setter
    def nchi(self, nchi):
        self['nchi'] = nchi

    @property
    def dzb_min(self):
        return self['dzb_min']

    @dzb_min.setter
    def dzb_min(self, dzb_min):
        self['dzb_min'] = dzb_min

    @property
    def dzb_decay(self):
        return self['dzb_decay']

    @dzb_decay.setter
    def dzb_decay(self, dzb_decay):
        self['dzb_decay'] = dzb_decay

    @property
    def hmin_man(self):
        return self['hmin_man']

    @hmin_man.setter
    def hmin_man(self, hmin_man):
        self['hmin_man'] = hmin_man

    @property
    def ncor(self):
        return self['ncor']

    @ncor.setter
    def ncor(self, ncor):
        raise AttributeError('ncor attribute must be set using set_ncor()'
                             'method.')

    # ncor needs arguments
    def set_ncor(self, ncor, coricoef=PARAM_DEFAULTS['coricoef'],
                 rlatitude=PARAM_DEFAULTS['rlatitude'],
                 sfea0=None):

        # clear any possible values
        self['coricoef'] = None
        self['rlatitude'] = None
        self['sfea0'] = None

        assert ncor in [0, -1, 1, 'coricoef', 'rlatitude', 'mesh']
        self['ncor'] = ncor

        if ncor in [0, 'coricoef']:
            if coricoef is None:
                raise AttributeError('coricoef is required when ncor=0')
            self['coricoef'] = coricoef

        elif ncor in [-1, 'rlatitude']:
            if rlatitude is None:
                raise AttributeError('rlatitude is required when ncor=-1')
            self['coricoef'] = rlatitude

        elif ncor in [1, 'mesh']:
            if self.ics == 1:
                if sfea0 is None:
                    raise AttributeError('sfea0 is required when ncor=-1 and '
                                         'ics=1')
                self['sfea0'] = sfea0

    @property
    def rlatitude(self):
        return self['rlatitude']

    @rlatitude.setter
    def rlatitude(self, rlatitude):
        raise AttributeError('rlatitude attribute must be set using set_ncor()'
                             'method.')

    @property
    def coricoef(self):
        return self['coricoef']

    @coricoef.setter
    def coricoef(self, coricoef):
        raise AttributeError('coricoef attribute must be set using set_ncor()'
                             'method.')

    @property
    def ic_elev(self):
        return self['ic_elev']

    @ic_elev.setter
    def ic_elev(self, ic_elev):
        self['ic_elev'] = ic_elev

    @property
    def nramp_elev(self):
        return self['nramp_elev']

    @nramp_elev.setter
    def nramp_elev(self, nramp_elev):
        self['nramp_elev'] = nramp_elev

    @property
    def inv_atm_bnd(self):
        return self['inv_atm_bnd']

    @inv_atm_bnd.setter
    def inv_atm_bnd(self, inv_atm_bnd):
        self['inv_atm_bnd'] = inv_atm_bnd

    @property
    def prmsl_ref(self):
        return self['prmsl_ref']

    @prmsl_ref.setter
    def prmsl_ref(self, prmsl_ref):
        self['prmsl_ref'] = prmsl_ref

    @property
    def flag_ic(self):
        return self['flag_ic']

    @flag_ic.setter
    def flag_ic(self, flag_ic):
        self['flag_ic'] = flag_ic

    @property
    def gen_wsett(self):
        return self['gen_wsett']

    @gen_wsett.setter
    def gen_wsett(self, gen_wsett):
        self['gen_wsett'] = gen_wsett

    @property
    def ibcc_mean(self):
        return self['ibcc_mean']

    @ibcc_mean.setter
    def ibcc_mean(self, ibcc_mean):
        self['ibcc_mean'] = ibcc_mean

    @property
    def rmaxvel(self):
        return self['rmaxvel']

    @rmaxvel.setter
    def rmaxvel(self, rmaxvel):
        self['rmaxvel'] = rmaxvel

    @property
    def velmin_btrack(self):
        return self['velmin_btrack']

    @velmin_btrack.setter
    def velmin_btrack(self, velmin_btrack):
        self['velmin_btrack'] = velmin_btrack

    @property
    def btrack_nudge(self):
        return self['btrack_nudge']

    @btrack_nudge.setter
    def btrack_nudge(self, btrack_nudge):
        self['btrack_nudge'] = btrack_nudge

    @property
    def ibtrack_openbnd(self):
        return self['ibtrack_openbnd']

    @ibtrack_openbnd.setter
    def ibtrack_openbnd(self, ibtrack_openbnd):
        self['ibtrack_openbnd'] = ibtrack_openbnd

    @property
    def ihhat(self):
        return self['ihhat']

    @ihhat.setter
    def ihhat(self, ihhat):
        self['ihhat'] = ihhat

    @property
    def inunfl(self):
        return self['inunfl']

    @inunfl.setter
    def inunfl(self, inunfl):
        self['inunfl'] = inunfl

    @property
    def h0(self):
        return self['h0']

    @h0.setter
    def h0(self, h0):
        self['h0'] = h0

    @property
    def shorewafo(self):
        return self['shorewafo']

    @shorewafo.setter
    def shorewafo(self, shorewafo):
        self['shorewafo'] = shorewafo

    @property
    def moitn0(self):
        return self['moitn0']

    @moitn0.setter
    def moitn0(self, moitn0):
        self['moitn0'] = moitn0

    @property
    def mxitn0(self):
        return self['mxitn0']

    @mxitn0.setter
    def mxitn0(self, mxitn0):
        self['mxitn0'] = mxitn0

    @property
    def rtol0(self):
        return self['rtol0']

    @rtol0.setter
    def rtol0(self, rtol0):
        self['rtol0'] = rtol0

    @property
    def nadv(self):
        return self['nadv']

    @nadv.setter
    def nadv(self, nadv):
        self['nadv'] = nadv

    @property
    def dtb_max(self):
        return self['dtb_max']

    @dtb_max.setter
    def dtb_max(self, dtb_max):
        self['dtb_max'] = dtb_max

    @property
    def dtb_min(self):
        return self['dtb_min']

    @dtb_min.setter
    def dtb_min(self, dtb_min):
        self['dtb_min'] = dtb_min

    @property
    def inter_mom(self):
        return self['inter_mom']

    @inter_mom.setter
    def inter_mom(self, inter_mom):
        self['inter_mom'] = inter_mom

    @property
    def kr_co(self):
        return self['kr_co']

    @kr_co.setter
    def kr_co(self, kr_co):
        self['kr_co'] = kr_co

    @property
    def itr_met(self):
        return self['itr_met']

    @itr_met.setter
    def itr_met(self, itr_met):
        self['itr_met'] = itr_met

    @property
    def h_tvd(self):
        return self['h_tvd']

    @h_tvd.setter
    def h_tvd(self, h_tvd):
        self['h_tvd'] = h_tvd

    @property
    def eps1_tvd_imp(self):
        return self['eps1_tvd_imp']

    @eps1_tvd_imp.setter
    def eps1_tvd_imp(self, eps1_tvd_imp):
        self['eps1_tvd_imp'] = eps1_tvd_imp

    @property
    def eps2_tvd_imp(self):
        return self['eps2_tvd_imp']

    @eps2_tvd_imp.setter
    def eps2_tvd_imp(self, eps2_tvd_imp):
        self['eps2_tvd_imp'] = eps2_tvd_imp

    @property
    def ip_weno(self):
        return self['ip_weno']

    @ip_weno.setter
    def ip_weno(self, ip_weno):
        self['ip_weno'] = ip_weno

    @property
    def courant_weno(self):
        return self['courant_weno']

    @courant_weno.setter
    def courant_weno(self, courant_weno):
        self['courant_weno'] = courant_weno

    @property
    def nquad(self):
        return self['nquad']

    @nquad.setter
    def nquad(self, nquad):
        self['nquad'] = nquad

    @property
    def ntd_weno(self):
        return self['ntd_weno']

    @ntd_weno.setter
    def ntd_weno(self, ntd_weno):
        self['ntd_weno'] = ntd_weno

    @property
    def epsilon1(self):
        return self['epsilon1']

    @epsilon1.setter
    def epsilon1(self, epsilon1):
        self['epsilon1'] = epsilon1

    @property
    def epsilon2(self):
        return self['epsilon2']

    @epsilon2.setter
    def epsilon2(self, epsilon2):
        self['epsilon2'] = epsilon2

    @property
    def i_prtnftl_weno(self):
        return self['i_prtnftl_weno']

    @i_prtnftl_weno.setter
    def i_prtnftl_weno(self, i_prtnftl_weno):
        self['i_prtnftl_weno'] = i_prtnftl_weno

    @property
    def epsilon3(self):
        return self['epsilon3']

    @epsilon3.setter
    def epsilon3(self, epsilon3):
        self['epsilon3'] = epsilon3

    @property
    def ielad_weno(self):
        return self['ielad_weno']

    @ielad_weno.setter
    def ielad_weno(self, ielad_weno):
        self['ielad_weno'] = ielad_weno

    @property
    def small_elad(self):
        return self['small_elad']

    @small_elad.setter
    def small_elad(self, small_elad):
        self['small_elad'] = small_elad

    @property
    def nws(self):
        return self['nws']

    @nws.setter
    def nws(self, nws):
        self['nws'] = nws

    @property
    def wtiminc(self):
        return self['wtiminc']

    @wtiminc.setter
    def wtiminc(self, wtiminc):
        self['wtiminc'] = wtiminc

    @property
    def nrampwind(self):
        return self['nrampwind']

    @nrampwind.setter
    def nrampwind(self, nrampwind):
        self['nrampwind'] = nrampwind

    @property
    def drampwind(self):
        return self['drampwind']

    @drampwind.setter
    def drampwind(self, drampwind):
        self['drampwind'] = drampwind

    @property
    def iwindoff(self):
        return self['iwindoff']

    @iwindoff.setter
    def iwindoff(self, iwindoff):
        self['iwindoff'] = iwindoff

    @property
    def iwind_form(self):
        return self['iwind_form']

    @iwind_form.setter
    def iwind_form(self, iwind_form):
        self['iwind_form'] = iwind_form

    @property
    def impose_net_flux(self):
        return self['impose_net_flux']

    @impose_net_flux.setter
    def impose_net_flux(self, impose_net_flux):
        self['impose_net_flux'] = impose_net_flux

    @property
    def ihconsv(self):
        return self['ihconsv']

    @ihconsv.setter
    def ihconsv(self, ihconsv):
        self['ihconsv'] = ihconsv

    @property
    def isconsv(self):
        return self['isconsv']

    @isconsv.setter
    def isconsv(self, isconsv):
        self['isconsv'] = isconsv

    @property
    def itur(self):
        return self['itur']

    @itur.setter
    def itur(self, itur):
        self['itur'] = itur

    @property
    def dfv0(self):
        return self['dfv0']

    @dfv0.setter
    def dfv0(self, dfv0):
        self['dfv0'] = dfv0

    @property
    def dfh0(self):
        return self['dfh0']

    @dfh0.setter
    def dfh0(self, dfh0):
        self['dfh0'] = dfh0

    @property
    def mid(self):
        return self['mid']

    @mid.setter
    def mid(self, mid):
        self['mid'] = mid

    @property
    def stab(self):
        return self['stab']

    @stab.setter
    def stab(self, stab):
        self['stab'] = stab

    @property
    def xlsc0(self):
        return self['xlsc0']

    @xlsc0.setter
    def xlsc0(self, xlsc0):
        self['xlsc0'] = xlsc0

    @property
    def inu_elev(self):
        return self['inu_elev']

    @inu_elev.setter
    def inu_elev(self, inu_elev):
        self['inu_elev'] = inu_elev

    @property
    def inu_uv(self):
        return self['inu_uv']

    @inu_uv.setter
    def inu_uv(self, inu_uv):
        self['inu_uv'] = inu_uv

    @property
    def inu_tr(self):
        return self['inu_tr']

    @inu_tr.setter
    def inu_tr(self, inu_tr):
        self['inu_tr'] = inu_tr

    @property
    def vnh1(self):
        return self['vnh1']

    @vnh1.setter
    def vnh1(self, vnh1):
        self['vnh1'] = vnh1

    @property
    def vnf1(self):
        return self['vnf1']

    @vnf1.setter
    def vnf1(self, vnf1):
        self['vnf1'] = vnf1

    @property
    def vnh2(self):
        return self['vnh2']

    @vnh2.setter
    def vnh2(self, vnh2):
        self['vnh2'] = vnh2

    @property
    def vnf2(self):
        return self['vnf2']

    @vnf2.setter
    def vnf2(self, vnf2):
        self['vnf2'] = vnf2

    @property
    def step_nu_tr(self):
        return self['step_nu_tr']

    @step_nu_tr.setter
    def step_nu_tr(self, step_nu_tr):
        self['step_nu_tr'] = step_nu_tr

    @property
    def h_bcc1(self):
        return self['h_bcc1']

    @h_bcc1.setter
    def h_bcc1(self, h_bcc1):
        self['h_bcc1'] = h_bcc1

    @property
    def s1_mxnbt(self):
        return self['s1_mxnbt']

    @s1_mxnbt.setter
    def s1_mxnbt(self, s1_mxnbt):
        self['s1_mxnbt'] = s1_mxnbt

    @property
    def s2_mxnbt(self):
        return self['s2_mxnbt']

    @s2_mxnbt.setter
    def s2_mxnbt(self, s2_mxnbt):
        self['s2_mxnbt'] = s2_mxnbt

    @property
    def iharind(self):
        return self['iharind']

    @iharind.setter
    def iharind(self, iharind):
        self['iharind'] = iharind

    @property
    def iflux(self):
        return self['iflux']

    @iflux.setter
    def iflux(self, iflux):
        self['iflux'] = iflux

    @property
    def izonal5(self):
        return self['izonal5']

    @izonal5.setter
    def izonal5(self, izonal5):
        self['izonal5'] = izonal5

    @property
    def ibtrack_test(self):
        return self['ibtrack_test']

    @ibtrack_test.setter
    def ibtrack_test(self, ibtrack_test):
        self['ibtrack_test'] = ibtrack_test

    @property
    def irouse_test(self):
        return self['irouse_test']

    @irouse_test.setter
    def irouse_test(self, irouse_test):
        self['irouse_test'] = irouse_test

    @property
    def flag_fib(self):
        return self['flag_fib']

    @flag_fib.setter
    def flag_fib(self, flag_fib):
        self['flag_fib'] = flag_fib

    @property
    def slr_rate(self):
        return self['slr_rate']

    @slr_rate.setter
    def slr_rate(self, slr_rate):
        self['slr_rate'] = slr_rate

    @property
    def isav(self):
        return self['isav']

    @isav.setter
    def isav(self, isav):
        self['isav'] = isav

    @property
    def sav_cd(self):
        return self['sav_cd']

    @sav_cd.setter
    def sav_cd(self, sav_cd):
        self['sav_cd'] = sav_cd

    @property
    def nstep_ice(self):
        return self['nstep_ice']

    @nstep_ice.setter
    def nstep_ice(self, nstep_ice):
        self['nstep_ice'] = nstep_ice

    @property
    def level_age(self):
        return self['level_age']

    @level_age.setter
    def level_age(self, level_age):
        self['level_age'] = level_age

    @property
    def rearth_pole(self):
        return self['rearth_pole']

    @rearth_pole.setter
    def rearth_pole(self, rearth_pole):
        self['rearth_pole'] = rearth_pole

    @property
    def rearth_eq(self):
        return self['rearth_eq']

    @rearth_eq.setter
    def rearth_eq(self, rearth_eq):
        self['rearth_eq'] = rearth_eq

    @property
    def shw(self):
        return self['shw']

    @shw.setter
    def shw(self, shw):
        self['shw'] = shw

    @property
    def rho0(self):
        return self['rho0']

    @rho0.setter
    def rho0(self, rho0):
        self['rho0'] = rho0

    @property
    def vclose_surf_frac(self):
        return self['vclose_surf_frac']

    @vclose_surf_frac.setter
    def vclose_surf_frac(self, vclose_surf_frac):
        self['vclose_surf_frac'] = vclose_surf_frac

    @property
    def iadjust_mass_consv0(self):
        return self['iadjust_mass_consv0']

    @iadjust_mass_consv0.setter
    def iadjust_mass_consv0(self, iadjust_mass_consv0):
        self['iadjust_mass_consv0'] = iadjust_mass_consv0
