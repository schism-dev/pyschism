from enum import Enum
import os
import pathlib
from typing import Union, Tuple
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

    @staticmethod
    def open(file: Union[str, os.PathLike],
             crs: Union[str, CRS] = None):
        filename = pathlib.Path(file).name
        return FrictionDispatch[
            FrictionFilename(filename).name].value(
                **grd.read(pathlib.Path(file), boundaries=False, crs=crs))

    @classmethod
    def from_hgrid(cls, hgrid):
        # NOTE: nchi is set by subclass calling this method
        if isinstance(hgrid, (str, os.PathLike)):
            obj = cls.open(hgrid)
        elif isinstance(hgrid, Gr3):
            hgrid_dict = hgrid.to_dict()
            create_dict = deepcopy({
                k: hgrid_dict[k] for k in [
                    "description", "nodes", "elements", "crs"]})
            create_dict["nodes"] = create_dict["nodes"].copy()
            create_dict["elements"] = create_dict["elements"].copy()
            obj = cls(**create_dict)
        else:
            raise TypeError(
                f"Invalid hgrid type passed to create constant value"
                f" Fgrid: {type(hgrid)}")

        return obj

    @classmethod
    def constant(cls, hgrid, value):
        obj = cls.from_hgrid(hgrid)
        obj.nodes.values[:] = value
        return obj

    @classmethod
    def linear_with_depth(
            cls,
            hgrid: Union[str, os.PathLike, Gr3],
            values: Tuple[float, float],
            depths: Union[Tuple[float, float], None] = None):

        # Inspired by https://github.com/schism-dev/schism/blob/master/src/Utility/Pre-Processing/NWM/Manning/write_manning.py
        obj = cls.from_hgrid(hgrid)
        hgrid_depths = obj.values.copy()
        if not depths:
            depths = (np.min(hgrid_depths.ravel()),
                      np.max(hgrid_depths.ravel()))
        obj.nodes.values[:] = (
                values[0] + (hgrid_depths - depths[0])
                * (values[1] - values[0]) / (depths[1] - depths[0]))
        obj.nodes.values[obj.nodes.values < values[0]] = values[0]
        obj.nodes.values[obj.nodes.values > values[1]] = values[1]

        return obj


    def add_region(
            self,
            region: Union[Polygon, MultiPolygon],
            value: float = 0.02):
        # Assuming input polygons are in EPSG:4326
        if isinstance(region, Polygon):
            region = [region]
        gdf1 = gpd.GeoDataFrame(
                {'geometry': region}, crs='EPSG:4326')

        points = [Point(*coord) for coord in self.coords]
        gdf2 = gpd.GeoDataFrame(
                {'geometry': points, 'index': list(range(len(points)))},
                crs='EPSG:4326')
        gdf_in = gpd.sjoin(gdf2, gdf1, op="within")
        picks = ([i.index for i in gdf_in.itertuples()])
        self.nodes.values[picks] = value


class ManningsN(Fgrid):
    """  Class for representing Manning's n values.  """

    def __init__(self, *argv, **kwargs):
        self.hmin_man = 1.
        super().__init__(NchiType.MANNINGS_N, *argv, **kwargs)

    @classmethod
    def open(cls, file: Union[str, os.PathLike],
             crs: Union[str, CRS] = None):
        return super(Fgrid, cls).open(file, crs)

    @classmethod
    def constant(cls, hgrid, value):
        return super(ManningsN, cls).constant(hgrid, value)

    @classmethod
    def linear_with_depth(
            cls,
            hgrid: Union[str, os.PathLike, Gr3],
            values: Tuple[float, float],
            depths: Union[Tuple[float, float], None] = None):
        return super(ManningsN, cls).linear_with_depth(
                hgrid, values, depths)

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
