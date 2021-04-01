import pathlib
import subprocess
import tempfile
from typing import Union

import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon, Point

from pyschism.mesh.base import Gr3


class Gr3Field(Gr3):

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


class Albedo(Gr3Field):
    """ Class for writing albedo.gr3 file with constant value"""


class Diffmax(Gr3Field):
    pass


class Diffmin(Gr3Field):
    pass


class Watertype(Gr3Field):
    pass


class Shapiro(Gr3Field):

    @classmethod
    def from_binary(cls, hgrid, dst_crs, ref_slope):
        _tmpdir = tempfile.TemporaryDirectory()
        tmpdir = pathlib.Path(_tmpdir.name)
        hgrid = hgrid.copy()
        hgrid.transform_to(dst_crs)
        hgrid.write(tmpdir / 'hgrid.gr3')
        subprocess.check_call(['gen_slope_filter', ref_slope], cwd=tmpdir)
        obj = cls.open(tmpdir / 'slope_filter.gr3')
        obj.description = 'shapiro'
        return obj


class Windrot(Gr3Field):
    pass
