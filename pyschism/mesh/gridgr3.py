import os
import pathlib
import subprocess
import tempfile
from typing import Union

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point

from pyschism.mesh.base import Gr3


class Gr3Field(Gr3):

    @classmethod
    def constant(cls, hgrid, value):
        obj = cls(**{k: v for k, v in hgrid.to_dict().items() if k
                     in ['nodes', 'elements', 'description', 'crs']})
        obj.values[:] = value
        obj.description = f'{cls.__name__.lower()} {obj.crs}'
        return obj

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
    """ Class for writing albedo.gr3 file with constant value"""
    pass


class Diffmax(Gr3Field):
    pass


class Diffmin(Gr3Field):
    pass


class Watertype(Gr3Field):
    pass


class Fluxflag(Gr3Field):
    pass


class Tvdflag(Gr3Field):
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


class Nudge(Gr3Field):

    """
    This class is to generate nudge.gr3 file. The time complexity is O(n^2), which
    is bad for large mesh. 
    """

    def __init__(self):

        pass

    def gen_nudge(self, outdir: Union[str, os.PathLike], hgrid):

        #self.hgrid = hgrid

        outdir = pathlib.Path(outdir)

        hgrid = hgrid.to_dict()
        nodes = hgrid['nodes']
        elements = hgrid['elements']
        NE, NP = len(elements), len(nodes)
        lon = []
        lat = []
        for id, (coords, values) in nodes.items():
            lon.append(coords[0])
            lat.append(coords[1])

        bnd = hgrid['boundaries']
        opbd = bnd[None][0]['indexes']

        # Max relax distance in degr
        rlmax = 1.5
        # Max relax strength in days
        rnu_day = 0.25

        rnu = 0
        rnu_max = 1./rnu_day/86400.
        out = [f"{rlmax}, {rnu_day}"]
        out.extend("\n")
        out.append(f"{NE} {NP}")
        out.extend("\n")
        print(f'Max relax distnce is {rlmax}')
        for idn, (coords, values) in nodes.items():
            if idn in opbd:
                rnu = rnu_max
                distmin = 0.
            else:
                distmin = np.finfo(np.float64).max
                for j in opbd:
                    tmp = np.square(lon[int(idn)-1]-lon[int(j)-1]) +  \
                        np.square(lat[int(idn)-1]-lat[int(j)-1])
                    rl2 = np.sqrt(tmp)
                    if rl2 < distmin:
                        distmin = rl2
            rnu = 0.
            if distmin <= rlmax:
                rnu = (1-distmin/rlmax)*rnu_max
            line = [f"{idn}"]
            line.extend([f"{x:<.7e}" for x in coords])
            line.extend([f"{rnu:<.7e}"])
            line.extend("\n")
            out.append(" ".join(line))
            with open(outdir / 'TEM_nudge.gr3', 'w+') as fid:
                fid.writelines(out)
