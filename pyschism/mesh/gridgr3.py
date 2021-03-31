import os
import pathlib
from typing import Union

from pyproj import CRS
from shapely.geometry import Polygon, MultiPolygon, Point
import geopandas as gpd

from pyschism.mesh.base import Gr3
from pyschism.mesh.parsers import grd

class Gr3Filename(Enum):
    ALBEDO = 'albedo.gr3'
    DIFFMAX = 'diffmax.gr3'
    DIFFMIN = 'diffmin.gr3'
    WATERTYPE = 'watertype.gr3'
    SHAPIRO = 'shapiro.gr3'

class NchiType(Enum):
    AB = 
    

class Gr3ALL(Gr3):

    def __init__(self, *argv, **kwargs):
     
        pass

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

    def write(self, path: Union[str, os.PathLike], overwrite: bool = False):
        path = pathlib.Path(path)
        self.write(path.path.parent / , overwrite)

class Albedo(Gr3ALL):
    """ Class for writing albedo.gr3 file with constant value"""

    def __init__(self, *argv, **kwargs):
     
