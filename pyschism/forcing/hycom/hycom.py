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
            dst.createDimension('lon', xlon.shape[0])
            dst.createDimension('lat', ylat.shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('lon', 'f4', ('lon',))
            dst['lon'].long_name = "Longitude"
            dst['lon'].standard_name = "longitude"
            dst['lon'].units = "degrees_east"
            dst['lon'][:] = xlon
            # lat 
            dst.createVariable('lat', 'f4', ('lat',))
            dst['lat'].long_name = "Latitude"
            dst['lat'].standard_name = "latitude"
            dst['lat'].units = "degrees_north"
            dst['lat'][:] = ylat
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].long_name = 'Time'
            dst['time'].standard_name = 'time'
            dst['time'][:] = nc_salt['time'][1:2]
            #ssh
            dst.createVariable('ssh', 'f4', ('time', 'lat', 'lon',), fill_value=True)
            dst['ssh'].long_name = "sea_surface_elevation (m)"
            dst['ssh'][:,:,:] = nc_ssh['ssh'][8:9,0,lat_idxs,lon_idxs]

        with Dataset(outdir / 'TS_1.nc', 'w', format='NETCDF3_CLASSIC') as dst: 
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('lon', xlon.shape[0])
            dst.createDimension('lat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('lon', 'f4', ('lon',))
            dst['lon'].long_name = "Longitude"
            dst['lon'].standard_name = "longitude"
            dst['lon'].units = "degrees_east"
            dst['lon'][:] = xlon
            # lat 
            dst.createVariable('lat', 'f4', ('lat',))
            dst['lat'].long_name = "Latitude"
            dst['lat'].standard_name = "latitude"
            dst['lat'].units = "degrees_north"
            dst['lat'][:] = ylat
            #lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].long_name = 'Time'
            dst['time'].standard_name = 'time'
            dst['time'][:] = nc_salt['time'][1:2]
            #salt 
            dst.createVariable('salinity', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
            dst['salinity'].long_name = "sea_water_salinity (psu)" 
            #temp
            dst.createVariable('temperature', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
            dst['temperature'].long_name = "sea_water_potential_temperature (degc)"
        
            for k in np.arange(len(lev)):
                dst['salinity'][:,k,:,:] = nc_salt['salinity'][1:2,k,lat_idxs,lon_idxs]
                dst['temperature'][:,k,:,:] = nc_temp['temperature'][1:2,k,lat_idxs,lon_idxs]

        with Dataset(outdir / 'UV_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('lon', xlon.shape[0])
            dst.createDimension('lat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('lon', 'f4', ('lon',))
            dst['lon'].long_name = "Longitude"
            dst['lon'].standard_name = "longitude"
            dst['lon'].units = "degrees_east"
            dst['lon'][:] = xlon
            # lat 
            dst.createVariable('lat', 'f4', ('lat',))
            dst['lat'].long_name = "Latitude"
            dst['lat'].standard_name = "latitude"
            dst['lat'].units = "degrees_north"
            dst['lat'][:] = ylat
            #lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].long_name = 'Time'
            dst['time'].standard_name = 'time'
            dst['time'][:] = nc_salt['time'][1:2]
            #uvel
            dst.createVariable('u', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
            dst['u'].long_name = "eastward_sea_water_velocity (m/s)"
            #vvel
            dst.createVariable('v', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
            dst['v'].long_name = "northward_sea_water_velocity (m/s)"

            for k in np.arange(len(lev)):
                dst['u'][:,k,:,:] = nc_uvel['u'][1:2,k,lat_idxs,lon_idxs]
                dst['v'][:,k,:,:] = nc_vvel['v'][1:2,k,lat_idxs,lon_idxs]

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
            dst.createDimension('lon', xlon.shape[0])
            dst.createDimension('lat', ylat.shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('lon', 'f4', ('lon',))
            dst['lon'].long_name = "Longitude"
            dst['lon'].standard_name = "longitude"
            dst['lon'].units = "degrees_east"
            dst['lon'][:] = xlon
            # lat 
            dst.createVariable('lat', 'f4', ('lat',))
            dst['lat'].long_name = "Latitude"
            dst['lat'].standard_name = "latitude"
            dst['lat'].units = "degrees_north"
            dst['lat'][:] = ylat
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].long_name = 'Time'
            dst['time'].standard_name = 'time'
            dst['time'][:] = nc_salt['time'][1:rnday+2]
            #ssh
            dst.createVariable('ssh', 'f4', ('time', 'lat', 'lon',), fill_value=True)
            dst['ssh'].long_name = "sea_surface_elevation (m)"
            dst['ssh'][:,:,:] = nc_ssh['ssh'][8:8*(rnday+1)+1:8,0,jdx_min:jdx_max+1,idx_min:idx_max+1]

        with Dataset(outdir / 'TS_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('lon', xlon.shape[0])
            dst.createDimension('lat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('lon', 'f4', ('lon',))
            dst['lon'].long_name = "Longitude"
            dst['lon'].standard_name = "longitude"
            dst['lon'].units = "degrees_east"
            dst['lon'][:] = xlon
            # lat 
            dst.createVariable('lat', 'f4', ('lat',))
            dst['lat'].long_name = "Latitude"
            dst['lat'].standard_name = "latitude"
            dst['lat'].units = "degrees_north"
            dst['lat'][:] = ylat
            #lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].long_name = 'Time'
            dst['time'].standard_name = 'time'
            dst['time'][:] = nc_salt['time'][1:rnday+2]
            #salt 
            dst.createVariable('salinity', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
            dst['salinity'].long_name = "sea_water_salinity (psu)"
            dst['salinity'][:,:,:,:] = nc_salt['salinity'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]
            #temp
            dst.createVariable('temperature', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
            dst['temperature'].long_name = "sea_water_potential_temperature (degc)"
            dst['temperature'][:,:,:,:] = nc_temp['temperature'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]

        with Dataset(outdir / 'UV_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('lon', xlon.shape[0])
            dst.createDimension('lat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('lon', 'f4', ('lon',))
            dst['lon'].long_name = "Longitude"
            dst['lon'].standard_name = "longitude"
            dst['lon'].units = "degrees_east"
            dst['lon'][:] = xlon
            # lat 
            dst.createVariable('lat', 'f4', ('lat',))
            dst['lat'].long_name = "Latitude"
            dst['lat'].standard_name = "latitude"
            dst['lat'].units = "degrees_north"
            dst['lat'][:] = ylat
            #lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].long_name = 'Time'
            dst['time'].standard_name = 'time'
            dst['time'][:] = nc_salt['time'][1:rnday+2]
            #uvel
            dst.createVariable('u', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
            dst['u'].long_name = "eastward_sea_water_velocity (m/s)"
            dst['u'][:,:,:,:] = nc_uvel['u'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]
            #vvel
            dst.createVariable('v', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
            dst['v'].long_name = "northward_sea_water_velocity (m/s)"
            dst['v'][:,:,:,:] = nc_vvel['v'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]
