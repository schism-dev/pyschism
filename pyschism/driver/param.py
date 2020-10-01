from datetime import datetime, timedelta
from enum import Enum
import pathlib
from typing import Union

import f90nml  # type: ignore[import]

from .core import CORE
from .opt import OPT
from .schout import SCHOUT
from .output_vars import OUTPUT_VARS


class IntIbcType(Enum):
    BAROCLINIC = 0
    BAROTROPIC = 1


class StrIbcType(Enum):
    BAROCLINIC = 'baroclinic'
    BAROTROPIC = 'barotropic'


class Stratification(Enum):
    BAROCLINIC = 0
    BAROTROPIC = 1


class Param:

    def __init__(
            self,
            rnday: Union[float, timedelta],
            nspool: Union[int, timedelta],
            dt: Union[float, int, timedelta] = timedelta(seconds=150.),
            start_date: Union[None, datetime] = None,
            dramp: Union[float, int, timedelta] = timedelta(seconds=0.),
            ibc: Union[Stratification, int, str] = Stratification.BAROTROPIC,
            **output_vars
    ):
        self._src = pathlib.Path(__file__).parent / 'param.nml'
        self._nml = f90nml.read(self._src)
        self.__core = CORE(self._nml)
        self.__opt = OPT(self._nml)
        self.__schout = SCHOUT(self._nml)
        self._rnday = rnday
        self._dt = dt
        self._nspool = nspool
        self._dramp = dramp
        self._start_date = start_date
        self._ibc = ibc

        # set SCHISM output requests
        allowed_vars = [var.lower() for var in OUTPUT_VARS]
        for var in output_vars.copy():
            if var.lower() in allowed_vars:
                setattr(self.schout, var, bool(output_vars.pop(var)))

        # check for remaining arguments.
        if len(output_vars) > 0:
            raise TypeError(f'Unknown output_vars: {list(output_vars.keys())}')

    def write(self, path, overwrite=False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise IOError(f"File {path} exists and overwrite=False")
        f90nml.patch(self._src, self._nml, path)

    @property
    def core(self):
        return self.__core

    @property
    def opt(self):
        return self.__opt

    @property
    def schout(self):
        return self.__schout

    @property
    def start_date(self):
        return self.__start_date

    @property
    def end_date(self):
        return self.start_date + timedelta(days=self.core.rnday)

    @property
    def _rnday(self):
        return self.__rnday

    @_rnday.setter
    def _rnday(self, rnday: Union[float, timedelta]):
        self.__rnday = rnday
        if isinstance(rnday, timedelta):
            rnday = rnday.days
        self.core.rnday = rnday

    @property
    def _dt(self):
        return self.__dt

    @_dt.setter
    def _dt(self, dt: Union[float, int, timedelta]):
        self.__dt = dt
        if isinstance(dt, int):
            dt = float(dt)
        if isinstance(dt, timedelta):
            dt = dt.total_seconds()
        self.core.dt = dt

    @property
    def _nspool(self):
        return self.__nspool

    @_nspool.setter
    def _nspool(self, nspool: Union[int, timedelta]):
        self.__nspool = nspool
        if isinstance(nspool, timedelta):
            nspool = int(round(nspool.total_seconds() / self.core.dt))
        self.core.nspool = nspool

    @property
    def _dramp(self):
        return self.__dramp

    @_dramp.setter
    def _dramp(self, dramp: Union[float, int, timedelta]):
        self.__dramp = dramp
        if isinstance(dramp, int):
            dramp = float(dramp)
        if isinstance(dramp, timedelta):
            dramp = dramp.days
        self.opt.dramp = dramp

    @property
    def _start_date(self):
        return self.__start_date

    @_start_date.setter
    def _start_date(self, start_date: Union[None, datetime]):
        self.__start_date = start_date
        if isinstance(start_date, datetime):
            self.__forcing_start_date = start_date - timedelta(
                days=self.opt.dramp)
            if start_date.tzinfo is not None \
                    and start_date.tzinfo.utcoffset(start_date) is not None:
                self.opt.utc_start = -start_date.utcoffset().total_seconds() / 3600  # type: ignore[union-attr]
            self.opt.start_year = self._forcing_start_date.year
            self.opt.start_month = self._forcing_start_date.month
            self.opt.start_day = self._forcing_start_date.day
            self.opt.start_hour = self._forcing_start_date.hour
            self.opt.start_hour += self._forcing_start_date.minute / 60.

    @property
    def _ibc(self):
        return self.__ibc

    @_ibc.setter
    def _ibc(self, ibc: Union[Stratification, int, str]):
        self.__ibc = ibc
        if isinstance(ibc, int):
            self.core.ibc = Stratification[IntIbcType(ibc).name].value

        if isinstance(ibc, str):
            self.core.ibc = Stratification[StrIbcType(ibc.lower()).name].value

        # set ibtp which depends on ibc
        if self.core.ibc == 1:
            self.core.ibtp = 0
        elif self.core.ibc == 0:
            self.core.ibtp = 1

    @property
    def _forcing_start_date(self):
        try:
            return self.__forcing_start_date
        except AttributeError:
            raise Exception('ERROR: Cannot compute start date for forcing, '
                            'no start date provided.')
