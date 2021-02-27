from datetime import datetime, timedelta
import logging
import pathlib
from typing import Union

import f90nml
import pytz

from pyschism.enums import Coriolis
from pyschism.forcing.atmosphere.nws import NWS
from pyschism.forcing.atmosphere.nws.nws2 import NWS2

PARAM_TEMPLATE = pathlib.Path(__file__).parent / 'param.nml.template'
PARAM_DEFAULTS = f90nml.read(PARAM_TEMPLATE)['opt']

_logger = logging.getLogger(__name__)

class OptMeta(type):

    def __new__(meta, name, bases, attrs):
        for key, value in PARAM_DEFAULTS.items():
            if key not in attrs:
                if isinstance(value, list):
                    attrs[key] = len(value)*[0]
                else:
                    attrs[key] = None
        return type(name, bases, attrs)


class Dramp:

    def __set__(self, obj, dramp: Union[int, float, timedelta, None]):

        if not isinstance(dramp, (int, float, timedelta, type(None))):
            raise TypeError("Argument drampbc must be an int, float, "
                            "timedelta, or None.")

        if isinstance(dramp, (int, float)):
            dramp = timedelta(days=dramp)

        if dramp is not None:
            obj.nramp = 1

        obj.__dict__['dramp'] = dramp

    def __get__(self, obj, val):
        return obj.__dict__.get('dramp')


class Drampbc:

    def __set__(self, obj, drampbc: Union[int, float, timedelta, None]):
        if not isinstance(drampbc, (int, float, timedelta, type(None))):
            raise TypeError("Argument drampbc must be an int, float, "
                            "timedelta, or None.")

        if isinstance(drampbc, (int, float)):
            drampbc = timedelta(days=drampbc)

        if drampbc is not None:
            obj.nrampbc = 1

        obj.__dict__['drampbc'] = drampbc

    def __get__(self, obj, val):
        return obj.__dict__.get('drampbc')


class StartDate:

    def __set__(self, obj, start_date: Union[datetime, None]):
        if not isinstance(start_date, (datetime, type(None))):
            raise TypeError("Argument start_date must be of type datetime or "
                            "None. The datetime object will be assumed to be "
                            "in UTC if ")
        if start_date is not None:
            if start_date.tzinfo is None \
                    or start_date.tzinfo.utcoffset(start_date) is None:
                start_date = start_date.replace(tzinfo=pytz.utc)
            obj.start_year = start_date.year
            obj.start_month = start_date.month
            obj.start_day = start_date.day
            obj.start_hour = start_date.hour
            obj.start_hour += start_date.minute / 60.
            obj.utc_start = -start_date.utcoffset().total_seconds() / 3600  # type: ignore[union-attr]  # noqa: E501
            # get rid of "negative" zero
            obj.utc_start = +0. if obj.utc_start == -0. \
                else obj.utc_start
        obj.__dict__['start_date'] = start_date

    def __get__(self, obj, val):
        return obj.__dict__.get('start_date')


class Nchi:

    def __set__(self, obj, fgrid):
        nchi = fgrid.nchi
        if nchi == -1:
            obj.hmin_man = fgrid.hmin_man
        if obj.nchi == 1:
            obj.dbz_min = fgrid.dbz_min
            obj.dbz_decay = fgrid.dbz_decay
        obj.__dict__['nchi'] = nchi

    def __get__(self, obj, val):
        return obj.__dict__.get('nchi')


class Ics:

    def __set__(self, obj, ics: int):
        obj.__dict__['ics'] = ics

    def __get__(self, obj, val):
        return obj.__dict__.get('ics')


class Sfea0:

    def __set__(self, obj, sfea0: float):
        obj.__dict__['sfea0'] = sfea0

    def __get__(self, obj, val):
        return obj.__dict__.get('sfea0')


class Ncor:

    def __set__(self, obj, ncor: Coriolis):
        if not isinstance(ncor, Coriolis):
            raise TypeError(f"ncor must be of type {Coriolis}, not type "
                            f"{type(ncor)}.")
        obj.__dict__['ncor'] = ncor.value

    def __get__(self, obj, val):
        return obj.__dict__.get('ncor')


class Ihot:

    def __set__(self, obj, ihot: int):
        if not isinstance(ihot, int):
            raise TypeError(f"ihot must be of type {int}, not type "
                            f"{type(ihot)}.")
        if ihot not in [0, 1]:
            raise ValueError("ihot must be 0 or 1")

        obj.__dict__['ihot'] = ihot

    def __get__(self, obj, val):
        return obj.__dict__.get('ihot')


class Nws:
    def __set__(self, obj, nws: NWS):
        if not isinstance(nws, NWS):
            raise TypeError(
                f"nws must be of type {NWS}, not type {type(nws)}.")

        if obj.start_date is None:
            raise ValueError(
                "Can't initialize atmospheric data without start_date")
            nws._start_date = obj.start_date

        if isinstance(nws, NWS2):
            obj.__dict__['nws'] = 2

        else:
            raise NotImplementedError(
                f'NWS type {type(nws)} is not implemented.')

    def __get__(self, obj, val):
        return obj.__dict__.get('nws')


class OPT(metaclass=OptMeta):
    """ Provides error checking implementation for OPT group """

    dramp = Dramp()
    drampbc = Drampbc()
    start_date = StartDate()
    nchi = Nchi()
    ics = Ics()
    ncor = Ncor()
    sfea0 = Sfea0()
    ihot = Ihot()
    nws = Nws()

    def __init__(
            self,
            dramp: Union[int, float, timedelta] = None,
            drampbc: Union[int, float, timedelta] = None,
            # drampwind:,
            start_date: datetime = None):
        _logger.info('Iniatializing OPT.')
        self.dramp = dramp
        self.drampbc = drampbc
        self.start_date = start_date
        # TODO: Generate descriptors for the following properties:
        # wtiminc = 150. !time step for atmos. forcing. Default: same as dt
        # nrampwind = 1 !ramp-up option for atmos. forcing
        # drampwind = 1. !needed if nrampwind/=0; ramp-up period in days
        # iwindoff = 0 !needed only if nws/=0; '1': needs windfactor.gr3
        # iwind_form = -1 !needed if nws/=0
        # impose_net_flux = 0

    def __str__(self):
        data = []
        for key, default in PARAM_DEFAULTS.items():
            current = getattr(self, key)
            if key in ['dramp', 'drampbc']:
                if current is not None:
                    current = current.total_seconds() / (60*60*24)
            if current is not None:
                if not isinstance(current, list):
                    data.append(f'  {key}={current}')
                else:
                    for i, state in enumerate(current):
                        if state:
                            current.append(f'  {current}({i+1}) = 1')
        data = '\n'.join(data)
        return f"&OPT\n{data}\n/"
