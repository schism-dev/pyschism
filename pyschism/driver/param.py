import pathlib
from functools import lru_cache
import f90nml
from pyschism.driver.core import CORE
from pyschism.driver.opt import OPT
from pyschism.driver.schout import SCHOUT


class Param:

    @staticmethod
    def open(path):
        nml = f90nml.read(path)
        param = Param()
        for group, attrs in nml.items():
            for attr, value in attrs.items():
                setattr(param, f"{group}.{attr}", value)
        return param

    def write(self, path, overwrite=False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            msg = f"File {path} exists and overwrite=False"
            raise IOError(msg)
        f90nml.patch(self.src, self.nml, path)

    @property
    def core(self):
        return CORE(self.nml)

    @property
    def opt(self):
        return OPT(self.nml)

    @property
    def schout(self):
        return SCHOUT(self.nml)

    @property
    @lru_cache
    def src(self):
        return pathlib.Path(__file__).parent / 'param.nml'

    @property
    @lru_cache
    def nml(self):
        return f90nml.read(self.src.resolve())
