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


def transform_to_cpp(coords, lonc, latc):
    longitude = list(coords[:, 0] / 180 * np.pi)
    latitude = list(coords[:, 1] / 180 * np.pi)
    radius = 6378206.4
    loncc = lonc / 180 * np.pi
    latcc = latc / 180 * np.pi
    x = [
        radius * (longitude[i] - loncc) * np.cos(latcc)
        for i in np.arange(len(longitude))
    ]
    y = [radius * latitude[i] for i in np.arange(len(latitude))]
    return np.vstack([x, y]).T


class Shapiro(Gr3Field):

    @classmethod
    def slope_filter(cls, hgrid, shapiro_vals1, depths, shapiro_max, threshold_slope, 
        regions=None, shapiro_vals2=None, flags=None, lonc=None, latc=None):
        """
        https://github.com/wzhengui/pylibs/blob/1ee35efaa2d52fa682113126d84846ba33318f99/Utility/schism_file.py#L265
        """
        hgrid = hgrid.copy()
        hgrid.nodes.transform_to_cpp(lonc, latc)
        xy = hgrid.nodes.coords
        x = xy[:, 0]
        y = xy[:, 1]
        dp = -hgrid.nodes.values
        elnode = hgrid.elements.array

        #breakpoint()
        #fp = np.any(elnode.mask, axis=1)
        fp = np.any(elnode.mask)
        if fp:
            fpn = ~fp
        x1 = x[elnode[:, 0]]
        y1 = y[elnode[:, 0]]
        v1 = dp[elnode[:, 0]]
        x2 = x[elnode[:, 1]]
        y2 = y[elnode[:, 1]]
        v2 = dp[elnode[:, 1]]
        x3 = x[elnode[:, 2]]
        y3 = y[elnode[:, 2]]
        v3 = dp[elnode[:, 2]]
        a1 = ((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2
        if fp:
            mask = np.any(elnode.mask, axis=1)
            x4 = x[elnode[:, 3]]
            y4 = y[elnode[:, 3]]
            v4 = dp[elnode[:, 3]]
            x4[mask] = x1[mask]
            y4[mask] = y1[mask]
            v4[mask] = v1[mask]
            a2 = ((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2

        # compute gradients
        dpedx = (v1*(y2-y3)+v2*(y3-y1)+v3*(y1-y2))/(2*a1)
        dpedy = ((x3-x2)*v1+(x1-x3)*v2+(x2-x1)*v3)/(2*a1)
        dpedxy = np.sqrt(dpedx**2+dpedy**2)

        # modify quads
        if fp:
            mask = np.any(elnode.mask, axis=1)
            fpn = ~mask
            dpedx2 = (v1[fpn]*(y3[fpn]-y4[fpn])+v3[fpn]*(y4[fpn]-y1[fpn])+v4[fpn]*(y1[fpn]-y3[fpn]))/(2*a2[fpn])
            dpedy2 = ((x4[fpn]-x3[fpn])*v1[fpn]+(x1[fpn]-x4[fpn])*v3[fpn]+(x3[fpn]-x1[fpn])*v4[fpn])/(2*a2[fpn])
            dpedxy2 = np.sqrt(dpedx2**2+dpedy2**2)
            
            dpedx[fpn] = (dpedx[fpn]+dpedx2)/2
            dpedy[fpn] = (dpedy[fpn]+dpedy2)/2
            dpedxy[fpn] = (dpedxy[fpn]+dpedxy2)/2

        # get node ball information
        nne, ine = hgrid.elements.get_node_ball()
         
        # interpolate into nodes
        dpdxy = []
        for inode in np.arange(len(dp)):
            ind = ine[inode]
            dpdxy.append(np.max(dpedxy[ind]))
        slope = np.array(dpdxy)

        shapiro = shapiro_max*np.tanh(2*slope/threshold_slope)

        # further tweaks on shallow waters
        if len(depths) != len(shapiro_vals1):
            raise Exception(f'lengths of depths {len(depths)} and shapiro_vals1 {len(shapiro_vals1)} inconsistent')
        fp = dp < depths[-1]
        shapiro[fp] = np.maximum(shapiro[fp], np.interp(dp[fp], depths, shapiro_vals1))

        shapiro = cls(
            nodes={id: (xy[i, :], shapiro[i]) for i, id in enumerate(hgrid.nodes.id)},
            elements=hgrid.elements.to_dict(),
            crs=None,
            description=f"shapiro threshold_slope={threshold_slope} lonc={lonc} latc={latc}",
        )

        if regions is not None:
            for reg, value, flag in zip(regions, shapiro_vals2, flags):
                shapiro.modify_by_region(hgrid, reg, value, depths[0], flag)

        return shapiro


class Windrot(Gr3Field):
    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 0.0)


class IcField(Gr3Field):
    pass


class ElevIc(IcField):
    @classmethod
    def default(cls, hgrid, offset=0.1):
        obj = cls.constant(hgrid, 0.0)
        obj.values[:] = np.maximum(0.0, hgrid.values - offset)
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
