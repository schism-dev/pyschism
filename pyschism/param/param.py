from datetime import datetime, timedelta
import pathlib
from typing import Union

import f90nml  # type: ignore[import]

from .core import CORE, Stratification
from .opt import OPT
from .schout import SCHOUT


PARAM_TEMPLATE = pathlib.Path(__file__).parent / 'param.nml.template'
PARAM_DEFAULTS = f90nml.read(PARAM_TEMPLATE)


class Param:

    def __init__(
            self,
            dt: Union[float, int, timedelta],
            rnday: Union[float, timedelta],
            dramp: Union[int, float, timedelta],
            nspool: Union[int, timedelta],
            start_date: Union[None, datetime] = None,
            ibc: Union[Stratification, int, str] = Stratification.BAROTROPIC,
            ihfskip: Union[int, None] = None,
            **outputs
    ):
        # ---------------------------
        # initialize main container |
        # ---------------------------
        # Blank out all values in template, this program assumes no defaults
        self.__nml: dict = {'core': {}, 'opt': {}, 'schout': {}}
        for group, data in PARAM_DEFAULTS.items():
            for key, value in data.items():
                if isinstance(value, list):
                    self.__nml[group][key] = len(value)*[0]
                else:
                    self.__nml[group][key] = None

        # -----------------
        # initialize core |
        # -----------------
        self.__core = CORE(self.__nml)
        self.core.rnday = rnday
        self.core.dt = dt
        self.core.nspool = nspool
        self.core.dramp = dramp
        self.core.ibc = ibc
        self.core.ihfskip = ihfskip
        # TODO: Must set (when applicable)
        # msc2
        # mdc2
        # ntracer_gen
        # ntracer_age
        # sed_class
        # eco_class

        # ---------------
        # intialize opt |
        # ---------------
        self.__opt = OPT(self.__nml)
        self.opt.start_date = start_date

        # -------------------
        # initialize schout |
        # -------------------
        self.__schout = SCHOUT(self.__nml, **outputs)

    def __str__(self):

        def append_items(f, group):
            for var, value in self.__nml[group].items():
                if value is not None:
                    if not isinstance(value, list):
                        f.append(f'  {var}={value}')
                    else:
                        for i, state in enumerate(value):
                            if state == 1:
                                f.append(f'  {var}({i+1}) = {state}')
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
