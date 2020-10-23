from datetime import datetime, timedelta
import pathlib
from typing import Union

import f90nml  # type: ignore[import]

from .core import CORE, Stratification
from .opt import OPT
from .schout import SCHOUT
from ..stations import Stations

PARAM_TEMPLATE = pathlib.Path(__file__).parent / 'param.nml.template'


class Param:

    def __init__(
            self,
            dt: Union[float, timedelta],
            rnday: Union[float, timedelta],
            ihfskip: int = None,
            dramp: Union[float, timedelta] = None,
            start_date: datetime = None,
            ibc: Union[Stratification, int, str] = Stratification.BAROTROPIC,
            drampbc: Union[float, timedelta] = None,
            stations: Stations = None,
            nspool: Union[int, timedelta] = None,
            **outputs
    ):
        """Main interface for Param class



        """
        # -----------------
        # initialize core |
        # -----------------
        self.__core = CORE()
        self.core.rnday = rnday
        self.core.dt = dt
        self.core.nspool = nspool
        self.core.ibc = ibc
        self.core.ihfskip = ihfskip

        # ---------------
        # intialize opt |
        # ---------------
        self.__opt = OPT()
        self.opt.dramp = dramp
        self.opt.drampbc = drampbc
        self.opt.start_date = start_date

        # -------------------
        # initialize schout |
        # -------------------
        self.__schout = SCHOUT(**outputs)

        # --------------------
        # initialize station |
        # --------------------
        if stations is not None:
            if not isinstance(stations, Stations):
                raise TypeError('stations argument must be of type '
                                f'{Stations}')
        self.__stations = stations

    def __str__(self):

        def append_items(f, group):
            for var, value in getattr(self, group):
                if value is not None:
                    if not isinstance(value, list):
                        f.append(f'  {var}={value}')
                    else:
                        for i, state in enumerate(value):
                            if state:
                                f.append(f'  {var}({i+1}) = 1')
            return f

        f = ['&CORE']
        append_items(f, 'core')
        f.extend(['/\n', '&OPT'])
        append_items(f, 'opt')
        f.extend(['/\n', '&SCHOUT'])
        append_items(f, 'schout')
        f.append('/\n')
        return '\n'.join(f)

    def write(self, path, overwrite=False, use_template=False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise IOError(f"File {path} exists and overwrite=False")
        if use_template:
            f90nml.patch(PARAM_TEMPLATE, self.__nml, path)
            return
        with open(path, 'w') as f:
            f.write(str(self))

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
    def stations(self):
        return self.__stations
