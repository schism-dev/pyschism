import os
import pathlib
import subprocess
import tempfile
from typing import Union

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point

from pyschism.forcing.hycom import Hycom
from pyschism.mesh.base import Gr3


class Gr3Field(Gr3):
    @classmethod
    def constant(cls, hgrid, value):
        obj = cls(
            **{
                k: v
                for k, v in hgrid.to_dict().items()
                if k in ["nodes", "elements", "description", "crs"]
            }
        )
        obj.values[:] = value
        obj.description = f"{cls.__name__.lower()} {obj.crs}"
        return obj

    @classmethod
    def default(cls, hgrid):
        raise NotImplementedError(f"No default defined for {cls.__name__}.")

    def add_region(self, region: Union[Polygon, MultiPolygon], value):
        if isinstance(region, Polygon):
            region = [region]
        gdf1 = gpd.GeoDataFrame({"geometry": region}, crs=self.crs)

        points = [Point(*coord) for coord in self.coords]
        gdf2 = gpd.GeoDataFrame(
            {"geometry": points, "index": list(range(len(points)))}, crs=self.crs
        )
        gdf_in = gpd.sjoin(gdf2, gdf1, op="within")
        picks = [i.index for i in gdf_in.itertuples()]
        self.values[picks] = value

    def modify_by_region(self, hgrid, fname, value, depth1, flag):
        '''
        reset (flag==0) or add (flag==1) value to a region
        '''
        lines=[line.strip().split() for line in open(fname, 'r').readlines()]
        data=np.squeeze(np.array([lines[3:]])).astype('float')
        x=data[:,0]
        y=data[:,1]
        coords = list( zip(x, y))
        poly = Polygon(coords)

        #region is in cpp projection 
        gdf1 = gpd.GeoDataFrame(
                {'geometry': [poly]})

        points = [Point(*coord) for coord in self.coords]
        gdf2 = gpd.GeoDataFrame(
                 {'geometry': points, 'index': list(range(len(points)))})
        gdf_in = gpd.sjoin(gdf2, gdf1, op="within")
        picks = [i.index for i in gdf_in.itertuples()]
        if flag == 0:
            self.values[picks] = value
        else:
            picks2 = np.where(-hgrid.values > depth1)
            picks3 = np.intersect1d(picks, picks2)
            self.values[picks3] = self.values[picks3] + value


class Albedo(Gr3Field):
    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 0.15)


class Diffmax(Gr3Field):
    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 1.0)


class Diffmin(Gr3Field):
    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 1e-6)


class Watertype(Gr3Field):
    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 1.0)


class Shapiro(Gr3Field):
    @classmethod
    #def from_binary(cls, outdir: Union[str, os.PathLike], hgrid):
    def from_binary(cls, hgrid):
        _tmpdir = tempfile.TemporaryDirectory()
        tmpdir = pathlib.Path(_tmpdir.name)
        hgrid = hgrid.copy()
        hgrid.nodes.transform_to_cpp()
        hgrid.write(tmpdir / "hgrid.gr3")
        subprocess.check_call(["gen_slope_filter"], cwd=tmpdir)
        #outdir = pathlib.Path(outdir)
        #shutil.copy2(tmpdir / 'slope_filter.gr3', outdir / 'shapiro.gr3')
        obj = cls.open(tmpdir / "slope_filter.gr3", crs='epsg:4326')
        obj.description = "shapiro"
        return obj


class Windrot(Gr3Field):
    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 0.0)


class IcField(Gr3Field):
    pass


class ElevIc(IcField):
    @classmethod
    def default(cls, hgrid, offset=-0.1):
        obj = cls.constant(hgrid, 0.0)
        mask = np.logical_and(hgrid.values > 0.0, obj.values < hgrid.values)
        idxs = np.where(mask)
        obj.values[idxs] = hgrid.values[idxs] + offset
        return obj


class TempIc(IcField):
    @classmethod
    def from_hycom(cls, gr3: Gr3, hycom: Hycom, date):
        obj = cls.constant(gr3, np.nan)
        obj.values[:] = hycom.temperature.interpolate(obj, date)
        return obj


class SaltIc(IcField):
    @classmethod
    def from_hycom(cls, gr3: Gr3, hycom: Hycom, date):
        obj = cls.constant(gr3, np.nan)
        obj.values[:] = hycom.salinity.interpolate(obj, date)
        return obj


class Estuary(Gr3Field):
    pass
