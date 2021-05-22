import os
import pathlib
import subprocess
import tempfile
from typing import Union

import geopandas as gpd
from numba import jit, prange
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point

from pyschism.forcing.baroclinic import BaroclinicForcing
from pyschism.mesh.base import Gr3


class Gr3Field(Gr3):

    @classmethod
    def constant(cls, hgrid, value):
        obj = cls(**{k: v for k, v in hgrid.to_dict().items() if k
                     in ['nodes', 'elements', 'description', 'crs']})
        obj.values[:] = value
        obj.description = f'{cls.__name__.lower()} {obj.crs}'
        return obj

    @classmethod
    def default(cls, hgrid):
        raise NotImplementedError(f'No default defined for {cls.__name__}.')

    def add_region(
            self,
            region: Union[Polygon, MultiPolygon],
            value
    ):
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

    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 0.15)


class Diffmax(Gr3Field):

    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 1.)


class Diffmin(Gr3Field):

    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 1e-6)


class Watertype(Gr3Field):

    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 1.)


class Shapiro(Gr3Field):

    @classmethod
    def from_binary(cls, outdir: Union[str, os.PathLike], hgrid, dst_crs):
        _tmpdir = tempfile.TemporaryDirectory()
        tmpdir = pathlib.Path(_tmpdir.name)
        hgrid = hgrid.copy()
        hgrid.transform_to(dst_crs)
        hgrid.write(tmpdir / 'hgrid.gr3')
        subprocess.check_call(['gen_slope_filter'], cwd=tmpdir)
        outdir = pathlib.Path(outdir)
        obj = cls.open(tmpdir / 'slope_filter.gr3')
        obj.description = 'shapiro'
        return obj


class Windrot(Gr3Field):

    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 0.)


class ElevIc(Gr3Field):

    @classmethod
    def default(cls, hgrid, offset=-0.1):
        obj = cls.constant(hgrid, 0.)
        mask = np.logical_and(
                hgrid.values > 0.,
                obj.values < hgrid.values
            )
        idxs = np.where(mask)
        obj.values[idxs] = hgrid.values[idxs] + offset
        return obj


class TempIc(Gr3Field):

    @classmethod
    def from_forcing(cls, gr3: Gr3, forcing: BaroclinicForcing, date):
        obj = cls.constant(gr3, np.nan)
        obj.values[:] = forcing.temperature.interpolate(obj, date)
        return obj


class SaltIc(Gr3Field):

    @classmethod
    def from_forcing(cls, gr3: Gr3, forcing: BaroclinicForcing, date):
        obj = cls.constant(gr3, np.nan)
        obj.values[:] = forcing.salinity.interpolate(obj, date)
        return obj


class Estuary(Gr3Field):
    pass


class Nudge(Gr3Field):

    def __init__(self, hgrid, rlmax=1.5, rnu_day=0.25):

        @jit(nopython=True, parallel=True)
        def compute_nudge(lon, lat, opbd, out):

            nnode = lon.shape[0]

            rnu_max = 1./rnu_day/86400.

            for idn in prange(nnode):
                if idn in opbd:
                    rnu = rnu_max
                    distmin = 0.
                else:
                    distmin = np.finfo(np.float64).max
                    for j in opbd:
                        rl2 = np.sqrt(
                            np.square(lon[idn] - lon[j-1])
                            + np.square(lat[idn] - lat[j-1])
                            )
                        if rl2 < distmin:
                            distmin = rl2
                rnu = 0.
                if distmin <= rlmax:
                    rnu = (1-distmin/rlmax)*rnu_max
                out[idn] = rnu

        opbd = []
        for row in hgrid.boundaries.ocean.itertuples():
            opbd.extend(row.indexes.tolist())

        out = np.zeros(hgrid.values.shape)
        lon, lat = hgrid.get_xy(crs='epsg:4326')
        compute_nudge(lon, lat, opbd, out)
        self.values[:] = out
        self.description = f"{rlmax}, {rnu_day}"

    @classmethod
    def default(cls, hgrid):
        return cls(hgrid)
