from enum import Enum
import os
import pathlib
from typing import Union
from copy import deepcopy
import numpy as np

from pyproj import CRS  # type: ignore[import]
from shapely.geometry import Polygon, MultiPolygon, Point
import geopandas as gpd

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

    @classmethod
    def open(cls, file: Union[str, os.PathLike],
             crs: Union[str, CRS] = None):
        filename = pathlib.Path(file).name
        if cls.__name__ == "Fgrid":
            return FrictionDispatch[
                FrictionFilename(filename).name].value(
                    **grd.read(pathlib.Path(file), boundaries=False, crs=crs))
        else:
            return super().open(file, crs)

    @classmethod
    def constant(cls, hgrid, value):
        obj = cls(**{k: v for k, v in hgrid.to_dict().items() if k
                     in ['nodes', 'elements', 'description', 'crs']})
        obj.values[:] = value
        return obj

    def add_region(
            self,
            region: Union[Polygon, MultiPolygon],
            value
    ):
        # Assuming input polygons are in EPSG:4326
        if isinstance(region, Polygon):
            region = [region]
        gdf1 = gpd.GeoDataFrame(
                {'geometry': region}, crs=self.crs)

        points = [Point(*coord) for coord in self.coords]
        gdf2 = gpd.GeoDataFrame(
                {'geometry': points, 'index': list(range(len(points)))},
                crs=self.crs)
        gdf_in = gpd.sjoin(gdf2, gdf1, op="within")
        picks = ([i.index for i in gdf_in.itertuples()])
        self.values[picks] = value


class ManningsN(Fgrid):
    """  Class for representing Manning's n values.  """

    def __init__(self, *argv, **kwargs):
        self.hmin_man = 1.
        super().__init__(NchiType.MANNINGS_N, *argv, **kwargs)

    @classmethod
    def linear_with_depth(
            cls,
            hgrid: Union[str, os.PathLike, Gr3],
            min_value: float = 0.02,
            max_value: float = 0.05,
            min_depth: float = None,
            max_depth: float = None):

        # Inspired by https://github.com/schism-dev/schism/blob/master/src/Utility/Pre-Processing/NWM/Manning/write_manning.py
        obj = cls.constant(hgrid, np.nan)
        hgrid_depths = hgrid.values.copy()
        if min_depth is None:
            min_depth = np.min(hgrid_depths)
        if max_depth is None:
            max_depth = np.max(hgrid_depths)

        values = (
                min_value + (hgrid_depths - min_depth)
                * (max_value - min_value) / (max_depth - min_depth))

        if min_value is not None:
            values[values < min_value] = min_value

        if max_value is not None:
            values[values > max_value] = max_value

        obj.values[:] = values

        return obj


class RoughnessLength(Fgrid):

    def __init__(self, *argv, **kwargs):
        super().__init__(NchiType.ROUGHNESS_LENGTH, *argv, **kwargs)


class DragCoefficient(Fgrid):

    def __init__(self, *argv, **kwargs):
        self.dzb_min = 0.5
        self.dzb_decay = 0.
        super().__init__(NchiType.DRAG_COEFFICIENT, *argv, **kwargs)


class FrictionDispatch(Enum):
    MANNINGS_N = ManningsN
    DRAG_COEFFICIENT = DragCoefficient
    ROUGHNESS_LENGTH = RoughnessLength
