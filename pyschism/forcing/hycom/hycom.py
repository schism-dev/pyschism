import os
import numpy as np
from datetime import datetime,timedelta
import logging
import pathlib
from typing import Union

from netCDF4 import Dataset
from matplotlib.transforms import Bbox

from pyschism.mesh import Hgrid
from pyschism.dates import localize_datetime, nearest_cycle_date, pivot_time

logger = logging.getLogger(__name__)

class Nudge:

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

        #Max relax distance in degr
        rlmax = 1.5
        #Max relax strength in days
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
                        distmin=rl2
            rnu = 0.
            if distmin <= rlmax:
                rnu = (1-distmin/rlmax)*rnu_max
            line = [f"{idn}"]
            line.extend([f"{x:<.7e}" for x in coords])
            line.extend([f"{rnu:<.7e}"])
            line.extend("\n")
            out.append(" ".join(line))

            with open(outdir / 'TEM_nudge.gr3','w+') as fid:
                fid.writelines(out)

class HotStartInventory():

    def __init__(self):

        pass

    logger.info('Fetching RTOFS data')

    def fetch_data(self, outdir: Union[str, os.PathLike], hgrid, start_date):

        self.hgrid = hgrid
        self.start_date = start_date

        outdir = pathlib.Path(outdir)
        if outdir.name != '':
            outdir /= 'hotstart'
        outdir.mkdir(exist_ok=True)

        bbox = self.hgrid.get_bbox(crs='epsg:4326', output_type='bbox')
 
        #Go to yesterday's directory to get tomorrow's data
        dt = self.start_date - timedelta(days=1)
        nc_ssh = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_2ds_forecast_3hrly_diag')
        nc_salt = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_salt')
        nc_temp = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_temp')
        nc_uvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_uvel')
        nc_vvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_vvel')

        lon = nc_ssh['lon'][:]
        lat = nc_ssh['lat'][:]
        lev = nc_salt['lev'][:]

        xmin = bbox.x0 + 360. if bbox.x0 < 0 else bbox.x0
        xmax = bbox.x1 + 360. if bbox.x1 < 0 else bbox.x1
        bbox = Bbox.from_extents(xmin, bbox.y0, xmax, bbox.y1)

        lat_idxs = np.where((lat >= bbox.ymin) & (lat <= bbox.ymax))[0]
        lon_idxs = np.where((lon >= bbox.xmin) & (lon <= bbox.xmax))[0]

        xlon = lon[lon_idxs]
        for ilon in range(lon_idxs.size):
            if xlon[ilon] > 180.:
                xlon[ilon] = xlon[ilon] - 360.
        ylat = lat[lat_idxs]

        with Dataset(outdir / 'SSH_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat 
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:2]
            #ssh
            dst.createVariable('surf_el', 'f4', ('time', 'ylat', 'xlon',), fill_value=-30000.0)
            dst['surf_el'].long_name = "sea_surface_elevation (m)"
            dst['surf_el'].add_offset = 0.
            dst['surf_el'].scale_factor = 0.001
            dst['surf_el'][:,:,:] = nc_ssh['ssh'][8:9,0,lat_idxs,lon_idxs]

        with Dataset(outdir / 'TS_1.nc', 'w', format='NETCDF3_CLASSIC') as dst: 
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat 
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            #lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:2]
            #salt 
            dst.createVariable('salinity', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.0)
            dst['salinity'].long_name = "sea_water_salinity (psu)" 
            dst['salinity'].add_offset = 20.
            dst['salinity'].scale_factor = 0.001
            #temp
            dst.createVariable('temperature', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.0)
            dst['temperature'].long_name = "sea_water_potential_temperature (degc)"
            dst['temperature'].add_offset = 20.
            dst['temperature'].scale_factor = 0.001
        
            for k in np.arange(len(lev)):
                dst['salinity'][:,k,:,:] = nc_salt['salinity'][1:2,k,lat_idxs,lon_idxs]
                dst['temperature'][:,k,:,:] = nc_temp['temperature'][1:2,k,lat_idxs,lon_idxs]

        with Dataset(outdir / 'UV_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat 
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            #lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:2]
            #uvel
            dst.createVariable('water_u', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.0)
            dst['water_u'].long_name = "eastward_sea_water_velocity (m/s)"
            dst['water_u'].add_offset = 0.
            dst['water_u'].scale_factor = 0.001
            #vvel
            dst.createVariable('water_v', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.0)
            dst['water_v'].long_name = "northward_sea_water_velocity (m/s)"
            dst['water_v'].add_offset = 0.
            dst['water_v'].scale_factor = 0.001

            for k in np.arange(len(lev)):
                dst['water_u'][:,k,:,:] = nc_uvel['u'][1:2,k,lat_idxs,lon_idxs]
                dst['water_v'][:,k,:,:] = nc_vvel['v'][1:2,k,lat_idxs,lon_idxs]

       #symlink estaury.gr3 and *.in file 
       # os.symlink(estaury.gr3, './start/estuary.gr3')

class OpenBoundaryInventory():

    def __init__(self):

        pass

    def fetch_data(self, outdir: Union[str, os.PathLike], start_date, rnday, idx_min, idx_max, jdx_min, jdx_max):

        self.start_date = start_date
        self.rnday = rnday

        outdir = pathlib.Path(outdir)
        if outdir.name != '':
            outdir /= '3Dth'
        outdir.mkdir(exist_ok=True)

        #Go to yesterday's directory to get tomorrow's data
        dt = self.start_date - timedelta(days=1)
        nc_ssh = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_2ds_forecast_3hrly_diag')
        nc_salt = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_salt')
        nc_temp = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_temp')
        nc_uvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_uvel')
        nc_vvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_vvel')

        lon = nc_ssh['lon'][:]
        lat = nc_ssh['lat'][:]
        lev = nc_salt['lev'][:]

        xlon = lon[idx_min:idx_max+1]
        for ilon in range(xlon.size):
            if xlon[ilon] > 180.:
                xlon[ilon] = xlon[ilon] - 360.
        ylat = lat[jdx_min:jdx_max+1]

        with Dataset(outdir / 'SSH_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat 
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:rnday+2]
            #ssh
            dst.createVariable('surf_el', 'f4', ('time', 'ylat', 'xlon',), fill_value=-30000.0)
            dst['surf_el'].long_name = "sea_surface_elevation (m)"
            dst['surf_el'].add_offset = 0.
            dst['surf_el'].scale_factor = 0.001
            ssh = nc_ssh['ssh'][8:8*(rnday+1)+1:8,0,jdx_min:jdx_max+1,idx_min:idx_max+1]
            dst['surf_el'][:,:,:] = ssh

        with Dataset(outdir / 'TS_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat 
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            #lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:rnday+2]
            #salt 
            dst.createVariable('salinity', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.)
            dst['salinity'].long_name = "sea_water_salinity (psu)"
            dst['salinity'].add_offset = 20.
            dst['salinity'].scale_factor = 0.001
            sss = nc_salt['salinity'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]
            dst['salinity'][:,:,:,:] = sss
            #temp
            dst.createVariable('temperature', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.)
            dst['temperature'].long_name = "sea_water_potential_temperature (degc)"
            dst['temperature'].add_offset = 20.
            dst['temperature'].scale_factor = 0.001
            sst = nc_temp['temperature'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]
            dst['temperature'][:,:,:,:] = sst

        with Dataset(outdir / 'UV_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat 
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            #lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:rnday+2]
            #uvel
            dst.createVariable('water_u', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.)
            dst['water_u'].long_name = "eastward_sea_water_velocity (m/s)"
            dst['water_u'].add_offset = 0.
            dst['water_u'].scale_factor = 0.001
            uvel = nc_uvel['u'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]
            dst['water_u'][:,:,:,:] = uvel
            #vvel
            dst.createVariable('water_v', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.)
            dst['water_v'].long_name = "northward_sea_water_velocity (m/s)"
            dst['water_v'].add_offset = 0.
            dst['water_v'].scale_factor = 0.001
            vvel = nc_vvel['v'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]
            dst['water_v'][:,:,:,:] = vvel
