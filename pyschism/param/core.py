from datetime import timedelta
from enum import Enum
from typing import Union


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

    def __init__(self, nml):
        self.__nml = nml

    @property
    def ipre(self):
        return self.__nml['core']['ipre']

    @ipre.setter
    def ipre(self, ipre):
        assert ipre in [0, 1]
        self.__nml['core']['ipre'] = ipre

    @property
    def ibc(self):
        return self.__nml['core']['ibc']

    @ibc.setter
    def ibc(self, ibc):
        if isinstance(ibc, int):
            assert ibc in [0, 1]
            ibc = Stratification[IntIbcType(ibc).name].value

        elif isinstance(ibc, str):
            ibc = Stratification[StrIbcType(ibc.lower()).name].value

        else:
            assert isinstance(ibc, Stratification)
            ibc = ibc.value
        self.ibtp = 0
        self.__nml['core']['ibc'] = ibc

    @property
    def ibtp(self):
        return self.__nml['core']['ibtp']

    @ibtp.setter
    def ibtp(self, ibtp):
        assert ibtp in [0, 1], 'ibtp must be 0 or 1'
        if self.ibc == 1 and ibtp == 1:
            raise AttributeError(
                'ERROR: Cannot set ibtp=1: ibc must be equal to 0 but ibc is '
                'currently equal to 1')
        self.__nml['core']['ibtp'] = ibtp

    @property
    def rnday(self):
        return self.__nml['core']['rnday']

    @rnday.setter
    def rnday(self, rnday: Union[int, float, timedelta]):
        assert isinstance(rnday, (int, float, timedelta))
        if isinstance(rnday, int):
            rnday = float(rnday)
        elif isinstance(rnday, timedelta):
            rnday = rnday.days
        assert rnday > 0.
        self.__nml['core']['rnday'] = float(rnday)

    @property
    def dt(self):
        return self.__nml['core']['dt']

    @dt.setter
    def dt(self, dt: Union[int, float, timedelta]):
        if isinstance(dt, timedelta):
            dt = dt.total_seconds()
        assert dt > 0., "dt (timestep) must a positive float or integer."
        self.__nml['core']['dt'] = dt

    @property
    def msc2(self):
        return self.__nml['core']['msc2']

    @msc2.setter
    def msc2(self, msc2):
        self.__nml['core']['msc2'] = msc2

    @property
    def mdc2(self):
        return self.__nml['core']['mdc2']

    @mdc2.setter
    def mdc2(self, mdc2):
        self.__nml['core']['mdc2'] = mdc2

    @property
    def ntracer_gen(self):
        return self.__nml['core']['ntracer_gen']

    @ntracer_gen.setter
    def ntracer_gen(self, ntracer_gen):
        self.__nml['core']['ntracer_gen'] = ntracer_gen

    @property
    def ntracer_age(self):
        return self.__nml['core']['ntracer_age']

    @ntracer_age.setter
    def ntracer_age(self, ntracer_age):
        self.__nml['core']['ntracer_age'] = ntracer_age

    @property
    def sed_class(self):
        return self.__nml['core']['sed_class']

    @sed_class.setter
    def sed_class(self, sed_class):
        self.__nml['core']['sed_class'] = sed_class

    @property
    def eco_class(self):
        return self.__nml['core']['eco_class']

    @eco_class.setter
    def eco_class(self, eco_class):
        self.__nml['core']['eco_class'] = eco_class

    @property
    def nspool(self):
        return self.__nml['core']['nspool']

    @nspool.setter
    def nspool(self, nspool):
        if isinstance(nspool, timedelta):
            nspool = int(round(nspool.total_seconds() / self.dt))
        assert nspool >= 0, "nspool must be zero or positive integer."
        self.__nml['core']['nspool'] = nspool

    @property
    def ihfskip(self):
        return self.__nml['core']['ihfskip']

    @ihfskip.setter
    def ihfskip(self, ihfskip):
        if ihfskip is None:
            ihfskip = timedelta(days=self.rnday).total_seconds() / self.dt
        else:
            assert (ihfskip / self.dt).is_integer(), \
                "ihfskip / dt must be an integer but:\n" \
                f"{ihfskip} / {self.dt} = {ihfskip / self.dt}"
        self.__nml['core']['ihfskip'] = int(ihfskip)
