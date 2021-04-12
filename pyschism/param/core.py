from datetime import timedelta
from enum import Enum
import pathlib
from typing import Union

import f90nml

from pyschism.enums import Stratification


PARAM_TEMPLATE = pathlib.Path(__file__).parent / 'param.nml.template'
PARAM_DEFAULTS = f90nml.read(PARAM_TEMPLATE)['core']


class IbcType(Enum):
    BAROCLINIC = 0
    BAROTROPIC = 1

    @classmethod
    def _missing_(self, name):
        raise ValueError(f'{name} is not a valid integer for ibc. '
                         'Valid integers are 0 or 1.')


class CoreMeta(type):

    def __new__(meta, name, bases, attrs):
        attrs.update(**PARAM_DEFAULTS)
        return type(name, bases, attrs)


class CORE(metaclass=CoreMeta):
    """ Provides error checking implementation for CORE group """

    mandatory = ['ipre', 'ibc', 'ibtp', 'rnday', 'dt', 'nspool', 'ihfskip']

    def __init__(
            self,
            ipre: int = 0,
            ibc: int = 0,
            ibtp: int = 0,
            rnday: Union[float, timedelta] = 0.,
            dt: Union[float, timedelta] = 150.,
            nspool: Union[int, float, timedelta] = None,
            ihfskip: Union[int, timedelta] = None
    ):
        self.ipre = ipre
        self.ibc = ibc
        self.ibtp = ibtp
        self.rnday = rnday
        self.dt = dt
        self.nspool = nspool
        self.ihfskip = ihfskip

    def __str__(self):
        data = '\n'.join([
            f"  {key}={str(getattr(self, key))}" for key in self.to_dict()])
        return f"&CORE\n{data}\n/"

    def to_dict(self):
        output = {}
        for key, default in PARAM_DEFAULTS.items():
            current = getattr(self, key)
            if key == 'rnday':
                current = current.total_seconds() / (60.*60.*24)
            if key == 'dt':
                current = current.total_seconds()
            if key in self.mandatory:
                output[key] = current
            elif default != current:
                output[key] = current
        return output

    @property
    def ipre(self) -> int:
        return self._ipre

    @ipre.setter
    def ipre(self, ipre: int):
        if ipre not in [0, 1]:
            raise ValueError("Argument to ipre attribute must be 0 or 1")
        self._ipre = ipre

    @property
    def ibc(self):
        return self._ibc

    @ibc.setter
    def ibc(self, ibc: Union[Stratification, int, str]):

        if isinstance(ibc, str):
            ibc = IbcType[ibc.upper()].value

        if isinstance(ibc, Stratification):
            ibc = ibc.value

        if isinstance(ibc, int):
            if ibc not in [0, 1]:
                raise ValueError(
                    'Argument to attribute ibc must be of type '
                    f"{Stratification} or an 0, 1 integer or a string "
                    "'barotropic', 'baroclinic', not type "
                    f"{type(ibc)}.")

        self._ibc = ibc

    @property
    def ibtp(self):
        return self._ibtp

    @ibtp.setter
    def ibtp(self, ibtp: int):

        if ibtp not in [0, 1]:
            raise TypeError('Argument to attribute ibtp must be 0 or 1, not '
                            f"{ibtp}.")
        if ibtp == 1 and self.ibc == 0:
            raise ValueError("ibtp cannot be set to 1 because ibc is equal to "
                             'zero.')
        self._ibtp = ibtp

    @property
    def rnday(self) -> float:
        return self._rnday

    @rnday.setter
    def rnday(self, rnday: Union[float, timedelta]):
        if not isinstance(rnday, timedelta):
            rnday = timedelta(days=float(rnday))
        self._rnday = rnday.total_seconds() / timedelta(days=1).total_seconds()

    @property
    def dt(self):
        return self._dt

    @dt.setter
    def dt(self, dt: Union[float, timedelta]):
        if not isinstance(dt, timedelta):
            dt = timedelta(seconds=float(dt))
        self._dt = dt.total_seconds()

    @property
    def nspool(self):
        return self._nspool

    @nspool.setter
    def nspool(self, nspool: Union[int, float, timedelta, None]):
        if nspool is None:
            nspool = int(round(self.rnday / self.dt))
        if isinstance(nspool, timedelta):
            nspool = int(round(nspool / self.dt))
        if isinstance(nspool, float):
            nspool = int(round(timedelta(hours=nspool) / self.dt))
        if isinstance(nspool, (int, float)):
            if nspool <= 0:
                raise ValueError("nspool must be positive.")
        self._nspool = int(nspool)

    @property
    def ihfskip(self):
        return self._ihfskip

    @ihfskip.setter
    def ihfskip(self, ihfskip: Union[int, timedelta, None]):

        if not isinstance(ihfskip, (int, timedelta, type(None))):
            raise TypeError('Argument ihfskip must be int, timedelta or None.')

        if ihfskip is None:
            ihfskip = int(round(self.rnday / self.dt))

        if isinstance(ihfskip, timedelta):
            ihfskip = int(round(ihfskip / self.dt))

        if not (ihfskip / self.nspool).is_integer():
            raise ValueError("ihfskip/nspool must be an integer but got  "
                             f"{ihfskip}/{self.nspool}={ihfskip/self.nspool}")

        self._ihfskip = ihfskip
