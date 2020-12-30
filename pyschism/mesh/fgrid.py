from enum import Enum
import os
import pathlib
from typing import Union

from pyproj import CRS  # type: ignore[import]

from pyschism.mesh.base import Gr3
from pyschism.mesh.parsers import grd


class FrictionFilename(Enum):
    MANNINGS_N = 'manning.gr3'
    DRAG_COEFFICIENT = 'drag.gr3'
    ROUGHNESS_LENGTH = 'rough.gr3'

    @classmethod
    def _missing_(self, name):
        raise ValueError(f'{name} is not a valid filename for a friction '
                         'file.')


class NchiType(Enum):
    MANNINGS_N = -1
    ROUGHNESS_LENGTH = 1
    DRAG_COEFFICIENT = 0


class Fgrid(Gr3):
    """
    Base class for all friction types (e.g. manning.grd, drag.grd, etc...)
    """

    def __init__(self, nchi: NchiType, *argv, **kwargs):
        self._nchi = nchi
        self._fname = FrictionFilename[NchiType(nchi).name]
        super().__init__(*argv, **kwargs)

    @property
    def nchi(self):
        return self._nchi.value

    @property
    def fname(self):
        return self._fname.value

    @staticmethod
    def open(file: Union[str, os.PathLike],
             crs: Union[str, CRS] = None):
        filename = pathlib.Path(file).name
        return FrictionDispatch[
            FrictionFilename(filename).name].value(
                **grd.read(pathlib.Path(file), boundaries=False, crs=crs))


class ManningsN(Fgrid):
    """  Class for representing Manning's n values.  """

    def __init__(self, *argv, **kwargs):
        self.hmin_man = 1.
        super().__init__(NchiType.MANNINGS_N, *argv, **kwargs)

    @classmethod
    def open(cls, file: Union[str, os.PathLike],
             crs: Union[str, CRS] = None):
        return super(Fgrid, cls).open(file, crs)


class RoughnessLength(Fgrid):

    def __init__(self, *argv, **kwargs):
        super().__init__(NchiType.ROUGHNESS_LENGTH, *argv, **kwargs)

    @classmethod
    def open(cls, file: Union[str, os.PathLike],
             crs: Union[str, CRS] = None):
        return super(Fgrid, cls).open(file, crs)


class DragCoefficient(Fgrid):

    def __init__(self, *argv, **kwargs):
        self.dzb_min = 0.5
        self.dzb_decay = 0.
        super().__init__(NchiType.DRAG_COEFFICIENT, *argv, **kwargs)

    @classmethod
    def open(cls, file: Union[str, os.PathLike],
             crs: Union[str, CRS] = None):
        return super(Fgrid, cls).open(file, crs)


class FrictionDispatch(Enum):
    MANNINGS_N = ManningsN
    DRAG_COEFFICIENT = DragCoefficient
    ROUGHNESS_LENGTH = RoughnessLength
