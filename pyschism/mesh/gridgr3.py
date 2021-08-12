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
    def from_binary(cls, outdir: Union[str, os.PathLike], hgrid, dst_crs):
        _tmpdir = tempfile.TemporaryDirectory()
        tmpdir = pathlib.Path(_tmpdir.name)
        hgrid = hgrid.copy()
        hgrid.transform_to(dst_crs)
        hgrid.write(tmpdir / "hgrid.gr3")
        subprocess.check_call(["gen_slope_filter"], cwd=tmpdir)
        outdir = pathlib.Path(outdir)
        obj = cls.open(tmpdir / "slope_filter.gr3")
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
