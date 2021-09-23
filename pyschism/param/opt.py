from datetime import datetime, timedelta
import logging

# import pathlib
from typing import Union

import f90nml

# import pytz

from pyschism import dates
from pyschism.mesh.fgrid import NchiType
from pyschism.param.schism_init import GitParamTemplate

PARAM_DEFAULTS = f90nml.read(GitParamTemplate().path)["opt"]

logger = logging.getLogger(__name__)


# class FlagIcDescriptor:

#     ic_types = [
#         # (name, type, index)
#         ("ic_temp", gridgr3.TempIc, 0),
# ]

#     def __init__(self, name, ic_type, index):
#         self.type = ic_type
#         self.name = name
#         self.index = index

#     def __get__(self, obj, val):
#         return bool(obj.flag_ic[self.index])

#     def __set__(self, obj, val: bool):
#         if not isinstance(val, bool):
#             raise TypeError(
#                 f"Argument to {self.name} must be boolean, not " f"type {type(val)}."
#             )
#         obj.flag_ic[self.index] = int(val)


class OptMeta(type):
    def __new__(meta, name, bases, attrs):
        for key, value in PARAM_DEFAULTS.items():
            if key not in attrs:
                if isinstance(value, list):
                    attrs[key] = len(value) * [0]
                else:
                    attrs[key] = None
        # for ic_type in FlagIcDescriptor.ic_types:
        #     attrs[ic_type] = FlagIcDescriptor()
        return type(name, bases, attrs)


class OPT(metaclass=OptMeta):
    """Provides error checking implementation for OPT group"""

    def __init__(
        self,
        start_date: datetime = None,
        dramp: Union[int, float, timedelta] = None,
        drampbc: Union[int, float, timedelta] = None,
        dramp_ss: Union[float, timedelta] = None,
        drampwafo: Union[float, timedelta] = None,
        drampwind: Union[float, timedelta] = None,
        ics=None,
        nws=None,
        nchi: Union[int, NchiType] = None,
        hmin_man=None,
        dbz_min=None,
        dbz_decay=None,
        ic_elev=None,
    ):
        self.start_date = start_date
        self.dramp = dramp
        self.drampbc = drampbc
        self.dramp_ss = dramp_ss
        self.drampwafo = drampwafo
        self.drampwind = drampwind
        self.ics = ics
        self.nchi = nchi
        self.hmin_man = hmin_man
        self.dbz_min = dbz_min
        self.dbz_decay = dbz_decay
        self.ic_elev = ic_elev
        self.nws = nws

    def __str__(self):
        data = []
        for key, _ in PARAM_DEFAULTS.items():
            current = getattr(self, key)
            if current is not None:
                if not isinstance(current, list):
                    data.append(f"  {key}={current}")
                else:
                    for i, state in enumerate(current):
                        if state:
                            current.append(f"  {current}({i+1}) = 1")
        data = "\n".join(data)
        return f"&OPT\n{data}\n/"

    def to_dict(self):
        data = {}
        for key, _ in PARAM_DEFAULTS.items():
            current = getattr(self, key)
            if current is not None:
                if not isinstance(current, list):
                    data[key] = current
                else:
                    data[key] = len(current) * [0]
                    for i, state in enumerate(current):
                        if state:
                            data[key][i] = 1
        return data

    @property
    def start_date(self) -> Union[datetime, None]:
        return self._start_date

    @start_date.setter
    def start_date(self, start_date: Union[datetime, None]):
        if start_date is not None:
            start_date = dates.localize_datetime(start_date)
            self.start_year = start_date.year
            self.start_month = start_date.month
            self.start_day = start_date.day
            self.start_hour = start_date.hour
            self.start_hour += start_date.minute / 60.0
            self.utc_start = -start_date.utcoffset().total_seconds() / 3600.0  # type: ignore[union-attr]  # noqa: E501
        else:
            self.start_year = None
            self.start_month = None
            self.start_day = None
            self.start_hour = None
            self.utc_start = None
        self._start_date = start_date

    @property
    def start_year(self) -> Union[int, None]:
        return self._start_year

    @start_year.setter
    def start_year(self, start_year: Union[int, None]):
        if start_year is not None:
            start_year = int(start_year)
        self._start_year = start_year

    @property
    def start_month(self) -> Union[int, None]:
        return self._start_month

    @start_month.setter
    def start_month(self, start_month: Union[int, None]):
        if start_month is not None:
            start_month = int(start_month)
        self._start_month = start_month

    @property
    def start_day(self) -> Union[int, None]:
        return self._start_day

    @start_day.setter
    def start_day(self, start_day: Union[int, None]):
        if start_day is not None:
            start_day = int(start_day)
        self._start_day = start_day

    @property
    def start_hour(self) -> Union[float, None]:
        return self._start_hour

    @start_hour.setter
    def start_hour(self, start_hour: Union[float, None]):
        if start_hour is not None:
            start_hour = float(start_hour)
        self._start_hour = start_hour

    @property
    def utc_start(self) -> Union[float, None]:
        return self._utc_start

    @utc_start.setter
    def utc_start(self, utc_start: Union[float, None]):
        if utc_start is not None:
            utc_start = float(utc_start)
        self._utc_start = utc_start

    @property
    def dramp(self) -> Union[float, None]:
        return self._dramp

    @dramp.setter
    def dramp(self, dramp: Union[float, timedelta, None]):
        if dramp is not None:
            if not isinstance(dramp, timedelta):
                dramp = timedelta(days=float(dramp))
            self._dramp = float(dramp / timedelta(days=1))
        else:
            self._dramp = None

    @property
    def nramp(self) -> Union[int, None]:
        if not hasattr(self, "_nramp") and self.dramp is not None:
            return 1
        elif hasattr(self, "_nramp"):
            return self._nramp

    @nramp.setter
    def nramp(self, nramp: Union[int, None]):
        if nramp not in [0, 1]:
            raise ValueError("Argument nramp must be 0, 1.")
        self._nramp = nramp

    @nramp.deleter
    def nramp(self):
        if hasattr(self, "_nramp"):
            del self._nramp

    @property
    def drampbc(self) -> Union[float, None]:
        return self._drampbc

    @drampbc.setter
    def drampbc(self, drampbc: Union[float, timedelta, None]):
        if drampbc is not None:
            if not isinstance(drampbc, timedelta):
                drampbc = timedelta(days=float(drampbc))
            self._drampbc = float(drampbc / timedelta(days=1))
        else:
            self._drampbc = None

    @property
    def nrampbc(self) -> Union[int, None]:
        if not hasattr(self, "_nrampbc") and self.drampbc is not None:
            return 1
        elif hasattr(self, "_nrampbc"):
            return self._nrampbc

    @nrampbc.setter
    def nrampbc(self, nrampbc: Union[int, None]):
        if nrampbc not in [0, 1]:
            raise ValueError("Argument nrampbc must be 0, 1.")
        self._nrampbc = nrampbc

    @nrampbc.deleter
    def nrampbc(self):
        if hasattr(self, "_nrampbc"):
            del self._nrampbc

    @property
    def dramp_ss(self) -> Union[float, None]:
        return self._dramp_ss

    @dramp_ss.setter
    def dramp_ss(self, dramp_ss: Union[float, timedelta, None]):
        if dramp_ss is not None:
            if not isinstance(dramp_ss, timedelta):
                dramp_ss = timedelta(days=float(dramp_ss))
            self._dramp_ss = float(dramp_ss / timedelta(days=1))
        else:
            self._dramp_ss = None

    @property
    def nramp_ss(self) -> Union[int, None]:
        if not hasattr(self, "_nramp_ss") and self.dramp_ss is not None:
            return 1
        elif hasattr(self, "_nramp_ss"):
            return self._nramp_ss

    @nramp_ss.setter
    def nramp_ss(self, nramp_ss: Union[int, None]):
        if nramp_ss not in [0, 1]:
            raise ValueError("Argument nramp_ss must be 0, 1.")
        self._nramp_ss = nramp_ss

    @nramp_ss.deleter
    def nramp_ss(self):
        if hasattr(self, "_nramp_ss"):
            del self._nramp_ss

    @property
    def drampwafo(self) -> Union[float, None]:
        return self._drampwafo

    @drampwafo.setter
    def drampwafo(self, drampwafo: Union[float, timedelta, None]):
        if drampwafo is not None:
            if not isinstance(drampwafo, timedelta):
                drampwafo = timedelta(days=float(drampwafo))
            self._drampwafo = float(drampwafo / timedelta(days=1))
        else:
            self._drampwafo = None

    @property
    def nrampwafo(self) -> Union[int, None]:
        if not hasattr(self, "_nrampwafo") and self.drampwafo is not None:
            return 1
        elif hasattr(self, "_nrampwafo"):
            return self._nrampwafo

    @nrampwafo.setter
    def nrampwafo(self, nrampwafo: Union[int, None]):
        if nrampwafo not in [0, 1]:
            raise ValueError("Argument nrampwafo must be 0, 1.")
        self._nrampwafo = nrampwafo

    @nrampwafo.deleter
    def nrampwafo(self):
        if hasattr(self, "_nrampwafo"):
            del self._nrampwafo

    @property
    def drampwind(self) -> Union[float, None]:
        return self._drampwind

    @drampwind.setter
    def drampwind(self, drampwind: Union[float, timedelta, None]):
        if drampwind is not None:
            if not isinstance(drampwind, timedelta):
                drampwind = timedelta(days=float(drampwind))
            self._drampwind = float(drampwind / timedelta(days=1))
        else:
            self._drampwind = None

    @property
    def nrampwind(self) -> Union[int, None]:
        if not hasattr(self, "_nrampwind") and self.drampwind is not None:
            return 1
        elif hasattr(self, "_nrampwind"):
            return self._nrampwind

    @nrampwind.setter
    def nrampwind(self, nrampwind: Union[int, None]):
        if nrampwind not in [0, 1]:
            raise ValueError("Argument nrampwind must be 0, 1.")
        self._nrampwind = nrampwind

    @nrampwind.deleter
    def nrampwind(self):
        if hasattr(self, "_nrampwind"):
            del self._nrampwind

    @property
    def nchi(self) -> Union[int, None]:
        return self._nchi

    @nchi.setter
    def nchi(self, nchi: Union[int, NchiType, None]):
        if not isinstance(nchi, NchiType) and nchi is not None:
            nchi = NchiType(nchi).value
        self._nchi = nchi

    @property
    def ic_elev(self) -> None:
        return self._ic_elev

    @ic_elev.setter
    def ic_elev(self, ic_elev: Union[int, None]):
        assert ic_elev in [0, 1, None]
        self._ic_elev = ic_elev

    # @property
    # def flag_ic(self) -> List[int]:
    #     if not hasattr(self, "_flag_ic"):
    #         self._flag_ic = len(FlagIcDescriptor.ic_types) * [0]
    #     return self._flag_ic


# --- I think none of this is needed anymore but kept for here just in case.

# class Dramp:

#     def __set__(self, obj, dramp: Union[int, float, timedelta, None]):

#         if not isinstance(dramp, (int, float, timedelta, type(None))):
#             raise TypeError("Argument drampbc must be an int, float, "
#                             "timedelta, or None.")

#         if isinstance(dramp, (int, float)):
#             dramp = timedelta(days=dramp)

#         if dramp is not None:
#             obj.nramp = 1

#         obj.__dict__['dramp'] = dramp

#     def __get__(self, obj, val):
#         return obj.__dict__.get('dramp')


# class Drampbc:

#     def __set__(self, obj, drampbc: Union[int, float, timedelta, None]):
#         if not isinstance(drampbc, (int, float, timedelta, type(None))):
#             raise TypeError("Argument drampbc must be an int, float, "
#                             "timedelta, or None.")

#         if isinstance(drampbc, (int, float)):
#             drampbc = timedelta(days=drampbc)

#         if drampbc is not None:
#             obj.nrampbc = 1

#         obj.__dict__['drampbc'] = drampbc

#     def __get__(self, obj, val):
#         return obj.__dict__.get('drampbc')


# class StartDate:

#     def __set__(self, obj, start_date: Union[datetime, None]):
#         if not isinstance(start_date, (datetime, type(None))):
#             raise TypeError("Argument start_date must be of type datetime or "
#                             "None. The datetime object will be assumed to be "
#                             "in UTC if ")
#         if start_date is not None:
#             if start_date.tzinfo is None \
#                     or start_date.tzinfo.utcoffset(start_date) is None:
#                 start_date = start_date.replace(tzinfo=pytz.utc)
#             obj.start_year = start_date.year
#             obj.start_month = start_date.month
#             obj.start_day = start_date.day
#             obj.start_hour = start_date.hour
#             obj.start_hour += start_date.minute / 60.
#             obj.utc_start = -start_date.utcoffset().total_seconds() / 3600  # type: ignore[union-attr]  # noqa: E501
#             # get rid of "negative" zero
#             obj.utc_start = +0. if obj.utc_start == -0. \
#                 else obj.utc_start
#         obj.__dict__['start_date'] = start_date

#     def __get__(self, obj, val):
#         return obj.__dict__.get('start_date')


# class Nchi:

#     def __set__(self, obj, fgrid):
#         nchi = fgrid.nchi
#         if nchi == -1:
#             obj.hmin_man = fgrid.hmin_man
#         if obj.nchi == 1:
#             obj.dbz_min = fgrid.dbz_min
#             obj.dbz_decay = fgrid.dbz_decay
#         obj.__dict__['nchi'] = nchi

#     def __get__(self, obj, val):
#         return obj.__dict__.get('nchi')


# class Ics:

#     def __set__(self, obj, ics: int):
#         obj.__dict__['ics'] = ics

#     def __get__(self, obj, val):
#         return obj.__dict__.get('ics')


# class Sfea0:

#     def __set__(self, obj, sfea0: float):
#         obj.__dict__['sfea0'] = sfea0

#     def __get__(self, obj, val):
#         return obj.__dict__.get('sfea0')


# class Ncor:

#     def __set__(self, obj, ncor: Coriolis):
#         if not isinstance(ncor, Coriolis):
#             raise TypeError(f"ncor must be of type {Coriolis}, not type "
#                             f"{type(ncor)}.")
#         obj.__dict__['ncor'] = ncor.value

#     def __get__(self, obj, val):
#         return obj.__dict__.get('ncor')


# class Ihot:

#     def __set__(self, obj, ihot: int):
#         if not isinstance(ihot, int):
#             raise TypeError(f"ihot must be of type {int}, not type "
#                             f"{type(ihot)}.")
#         if ihot not in [0, 1]:
#             raise ValueError("ihot must be 0 or 1")

#         obj.__dict__['ihot'] = ihot

#     def __get__(self, obj, val):
#         return obj.__dict__.get('ihot')


# class Nws:
#     def __set__(self, obj, nws: NWS):
#         if not isinstance(nws, NWS):
#             raise TypeError(
#                 f"nws must be of type {NWS}, not type {type(nws)}.")

#         if obj.start_date is None:
#             raise ValueError(
#                 "Can't initialize atmospheric data without start_date")
#             nws._start_date = obj.start_date

#         if isinstance(nws, NWS2):
#             obj.__dict__['nws'] = 2

#         else:
#             raise NotImplementedError(
#                 f'NWS type {type(nws)} is not implemented.')

#     def __get__(self, obj, val):
#         return obj.__dict__.get('nws')
