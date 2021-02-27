from datetime import timedelta
from enum import Enum
import logging
import pathlib
from typing import Union

import f90nml  # type: ignore[import]

from pyschism.enums import Stratification

_logger = logging.getLogger(__name__)

PARAM_TEMPLATE = pathlib.Path(__file__).parent / 'param.nml.template'
PARAM_DEFAULTS = f90nml.read(PARAM_TEMPLATE)['core']


class IntIbcType(Enum):
    BAROCLINIC = 0
    BAROTROPIC = 1

    @classmethod
    def _missing_(self, name):
        raise ValueError(f'{name} is not a valid integer for ibc. '
                         'Valid integers are 0 or 1.')


class StrIbcType(Enum):
    BAROCLINIC = 'baroclinic'
    BAROTROPIC = 'barotropic'

    @classmethod
    def _missing_(self, name):
        raise ValueError(f'{name} is not a valid string for ibc. '
                         'Valid strings are "barotropic" or "baroclinic".')


class Ipre:

    def __set__(self, obj, ipre: int):
        if ipre not in [0, 1]:
            raise ValueError("Argument to ipre attribute must be 0 or 1")
        obj.__dict__['ipre'] = ipre

    def __get__(self, obj, val):
        return obj.__dict__.get('ipre', 0)


class Ibc:

    def __set__(self, obj, ibc: Union[Stratification, int, str]):
        if isinstance(ibc, int):
            ibc = Stratification(IntIbcType(ibc).value)

        if not isinstance(ibc, Stratification):
            raise TypeError('Argument to attribute ibc must be of type '
                            f"{Stratification} or an 0, 1 integer or a string "
                            "'barotropic', 'baroclinic', not type "
                            f"{type(ibc)}.")

        obj.__dict__['ibc'] = ibc

    def __get__(self, obj, val):
        return obj.__dict__['ibc'].value


class Ibtp:

    def __set__(self, obj, ibtp: int):

        if ibtp not in [0, 1]:
            raise TypeError('Argument to attribute ibtp must be 0 or 1, not '
                            f"{ibtp}.")
        if ibtp == 1 and obj.ibc == 0:
            raise ValueError("ibtp cannot be set to 1 because ibc is equal to "
                             'zero.')
        obj.__dict__['ibtp'] = ibtp

    def __get__(self, obj, val):
        return obj.__dict__.get('ibtp', 0)


class Rnday:

    def __set__(self, obj, rnday: Union[float, timedelta]):

        if isinstance(rnday, float):
            rnday = timedelta(days=rnday)

        if not isinstance(rnday, timedelta):
            raise TypeError('Argument to attribute rnday must be of type float'
                            f' (to represent days) or type {timedelta}.')

        obj.__dict__['rnday'] = rnday

    def __get__(self, obj, val):
        return obj.__dict__['rnday']


class Dt:

    def __set__(self, obj, dt: Union[int, float, timedelta]):

        if isinstance(dt, (int, float)):
            dt = timedelta(seconds=dt)

        if not isinstance(dt, timedelta):
            raise TypeError('Argument to attribute dt must be of type float'
                            f' (to represent days) or type {timedelta}.')

        obj.__dict__['dt'] = dt

    def __get__(self, obj, val):
        return obj.__dict__['dt']


class Nspool:

    def __set__(self, obj, nspool: Union[int, float, timedelta, None]):
        if nspool is None:
            nspool = int(round(obj.rnday / obj.dt))
        if isinstance(nspool, timedelta):
            nspool = int(round(nspool / obj.dt))
        if isinstance(nspool, float):
            nspool = int(round(timedelta(hours=nspool) / obj.dt))
        if isinstance(nspool, (int, float)):
            if nspool <= 0:
                raise ValueError("nspool must be positive.")
        obj.__dict__['nspool'] = int(nspool)

    def __get__(self, obj, val):
        return obj.__dict__['nspool']


class Ihfskip:

    def __set__(self, obj, ihfskip: Union[int, timedelta, None]):

        if not isinstance(ihfskip, (int, timedelta, type(None))):
            raise TypeError('Argument ihfskip must be int, timedelta or None.')

        if ihfskip is None:
            ihfskip = int(round(obj.rnday / obj.dt))

        if isinstance(ihfskip, timedelta):
            ihfskip = int(round(ihfskip / obj.dt))

        if not (ihfskip / obj.nspool).is_integer():
            raise ValueError("ihfskip/nspool must be an integer but got  "
                             f"{ihfskip}/{obj.nspool}={ihfskip/obj.nspool}")

        obj.__dict__['ihfskip'] = ihfskip

    def __get__(self, obj, val):
        return obj.__dict__['ihfskip']


class CoreMeta(type):

    def __new__(meta, name, bases, attrs):
        for key, value in PARAM_DEFAULTS.items():
            if key in attrs:
                continue
            else:
                attrs[key] = value
        return type(name, bases, attrs)


class CORE(metaclass=CoreMeta):
    """ Provides error checking implementation for CORE group """

    ipre = Ipre()
    ibc = Ibc()
    ibtp = Ibtp()
    rnday = Rnday()
    dt = Dt()
    nspool = Nspool()
    ihfskip = Ihfskip()

    def __init__(
            self,
            ibc: Union[int, str, Stratification],
            rnday: Union[int, float, timedelta],
            dt: Union[int, float, timedelta],
            nspool: Union[int, float, timedelta] = None,
            ihfskip: Union[int, timedelta] = None,
    ):
        _logger.info('Initializing CORE')
        self.ibc = ibc
        self.rnday = rnday
        self.dt = dt
        self.nspool = nspool
        self.ihfskip = ihfskip
        self._mandatory = ['ipre', 'ibc', 'ibtp', 'rnday', 'dt', 'nspool',
                           'ihfskip']

    def __str__(self):
        data = []
        for key, default in PARAM_DEFAULTS.items():
            current = getattr(self, key)
            if key == 'rnday':
                current = current.total_seconds() / (60.*60.*24)
            if key == 'dt':
                current = current.total_seconds()
            if key in self._mandatory:
                data.append(f"  {key}={str(current)}")
            elif default != current:
                data.append(f"  {key}={str(current)}")
        data = '\n'.join(data)
        return f"&CORE\n{data}\n/"
