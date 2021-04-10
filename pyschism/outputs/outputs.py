import os
import pathlib
from typing import Union

import cf

from pyschism.param.schout import SurfaceOutputVars


class OutputVariable:

    def __init__(self, name, fields) -> None:
        self.name = name
        self.nx_grids = fields.select_by_ncvar('lon')
        self.ny_grids = fields.select_by_ncvar('lat')
        self.fields = fields.select_by_ncvar(self.name)


class OutputCollection:

    surface_output_vars = SurfaceOutputVars()

    def __init__(self, path: Union[str, os.PathLike]):
        self.path = pathlib.Path(path)
        self.fields = cf.read(self.path.glob('*.nc'))

        for varname in self.surface_output_vars:
            setattr(self, varname, OutputVariable(varname, self.fields))
