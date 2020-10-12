from datetime import timedelta
from enum import Enum
import pathlib
from typing import Union

import f90nml  # type: ignore[import]

PARAM_TEMPLATE = pathlib.Path(__file__).parent / 'param.nml.template'
PARAM_DEFAULTS = f90nml.read(PARAM_TEMPLATE)['core']


class IntIbcType(Enum):
    BAROCLINIC = 0
    BAROTROPIC = 1


class StrIbcType(Enum):
    BAROCLINIC = 'baroclinic'
    BAROTROPIC = 'barotropic'


class Stratification(Enum):
    BAROCLINIC = 0
    BAROTROPIC = 1


class CORE:
    """ Provides error checking implementation for CORE group """

    def __init__(self):
        self.__core: dict = {}
        for key, value in PARAM_DEFAULTS.items():
            if isinstance(value, list):
                self.__core[key] = len(value)*[0]
            else:
                self.__core[key] = None

    def __getitem__(self, key):
        return self.__core[key]

    def __setitem__(self, key, value):
        self.__core[key] = value

    def __iter__(self):
        for key, value in self.__core.items():
            yield key, value

    @property
    def ipre(self):
        return self.__core['ipre']

    @ipre.setter
    def ipre(self, ipre):
        assert ipre in [0, 1], 'ipre must be 0 or 1'
        self.__core['ipre'] = ipre

    @ipre.deleter
    def ipre(self):
        self.__core['ipre'] = None

    @property
    def ibc(self):
        return Stratification(self.__core['ibc'])

    @ibc.setter
    def ibc(self, ibc: Union[Stratification, int, str]):
        if isinstance(ibc, int):
            assert ibc in [0, 1]
            ibc = Stratification[IntIbcType(ibc).name].value

        elif isinstance(ibc, str):
            ibc = Stratification[StrIbcType(ibc.lower()).name].value

        else:
            assert isinstance(ibc, Stratification)
            ibc = ibc.value
        self.__core['ibtp'] = 0
        self.__core['ibc'] = ibc

    @property
    def ibtp(self):
        return self.__core['ibtp']

    @ibtp.setter
    def ibtp(self, ibtp):
        assert ibtp in [0, 1], 'ibtp must be 0 or 1'
        if self.ibc.value == 1 and ibtp == 1:
            raise AttributeError(
                'ERROR: Cannot set ibtp=1: ibc must be equal to 0 but ibc is '
                'currently equal to 1')
        self.__core['ibtp'] = ibtp

    @property
    def rnday(self):
        return timedelta(days=self._rnday)

    @rnday.setter
    def rnday(self, rnday: Union[int, float, timedelta]):
        assert isinstance(rnday, (int, float, timedelta))
        if isinstance(rnday, int):
            rnday = float(rnday)
        elif isinstance(rnday, timedelta):
            rnday = rnday.days
        assert rnday > 0.
        self.__core['rnday'] = float(rnday)

    @property
    def _rnday(self):
        return self.__core['rnday']

    @property
    def dt(self):
        return timedelta(seconds=self._dt)

    @dt.setter
    def dt(self, dt: Union[int, float, timedelta]):
        if isinstance(dt, timedelta):
            dt = dt.total_seconds()
        assert dt > 0., "dt (timestep) must a positive float or integer."
        self.__core['dt'] = float(dt)

    @property
    def _dt(self):
        return self.__core['dt']

    @property
    def msc2(self):
        return self.__core['msc2']

    @msc2.setter
    def msc2(self, msc2):
        self.__core['msc2'] = msc2

    @property
    def mdc2(self):
        return self.__core['mdc2']

    @mdc2.setter
    def mdc2(self, mdc2):
        self.__core['mdc2'] = mdc2

    @property
    def ntracer_gen(self):
        return self.__core['ntracer_gen']

    @ntracer_gen.setter
    def ntracer_gen(self, ntracer_gen):
        self.__core['ntracer_gen'] = ntracer_gen

    @property
    def ntracer_age(self):
        return self.__core['ntracer_age']

    @ntracer_age.setter
    def ntracer_age(self, ntracer_age):
        self.__core['ntracer_age'] = ntracer_age

    @property
    def sed_class(self):
        return self.__core['sed_class']

    @sed_class.setter
    def sed_class(self, sed_class):
        self.__core['sed_class'] = sed_class

    @property
    def eco_class(self):
        return self.__core['eco_class']

    @eco_class.setter
    def eco_class(self, eco_class):
        self.__core['eco_class'] = eco_class

    @property
    def nspool(self):
        return self.__core['nspool'] * self.dt

    @nspool.setter
    def nspool(self, nspool):
        if nspool is None:
            if self.rnday is None:
                raise AttributeError('Cannot auto-set nspool: rnday is not '
                                     'set.')
            nspool = timedelta(days=self.rnday)
        if isinstance(nspool, timedelta):
            nspool = int(round(nspool / self.dt))
        assert nspool > 0, "nspool a positive integer."
        self.__core['nspool'] = int(nspool)

    @property
    def _nspool(self):
        return self.__core['nspool']

    @property
    def ihfskip(self):
        return self._ihfskip * self.dt

    @ihfskip.setter
    def ihfskip(self, ihfskip):
        if ihfskip is None:
            if self.rnday is None:
                raise AttributeError('Cannot auto-set ihfskip: rnday is not '
                                     'set.')
            ihfskip = timedelta(days=self.__core['rnday']) / self.dt

        else:
            assert (ihfskip / self.dt).is_integer(), \
                "ihfskip / dt must be an integer but:\n" \
                f"{ihfskip} / {self.dt} = {ihfskip / self.dt}"
        self.__core['ihfskip'] = int(ihfskip)

    @property
    def _ihfskip(self):
        return self.__core['ihfskip']
