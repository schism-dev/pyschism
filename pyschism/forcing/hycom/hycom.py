import numpy as np
from datetime import datetime
from netCDF4 import Dataset
from matplotlib.transforms import Bbox
import logging
import pathlib

from pyschism.mesh import Hgrid
from pyschism.dates import localize_datetime, nearest_cycle_date, pivot_time

logger = logging.getLogger(__name__)

class RTOFS(self, start_date, hgrid):

    logger.info('Fetching RTOFS data')
    self._hgrid = hgrid

    bbox = self._hgrid.get_bbox(crs='epsg:4326', output_type='bbox')
    #Bbox(x0=-98.0058874, y0=8.5344214, x1=-60.0400016, y1=45.83143077) 
    
    self.start_date = start_date
    nc_ssh = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
        + start_date.strftime(%Y%m%d)+f'/rtofs_glo_2ds_forecast_3hrly_diag')
    nc_salt = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
        + start_date.strftime(%Y%m%d)+f'/rtofs_glo_3dz_forecast_daily_salt')
    nc_temp = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
        + start_date.strftime(%Y%m%d)+f'/rtofs_glo_3dz_forecast_daily_temp')
    nc_uvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
        + start_date.strftime(%Y%m%d)+f'/rtofs_glo_3dz_forecast_daily_uvel')
    nc_vvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
        + start_date.strftime(%Y%m%d)+f'/rtofs_glo_3dz_forecast_daily_vvel')

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
        if xlon[ilon] > 180:

    xlon = lon[lon_idxs]
    for ilon in range(lon_idxs.size):
        if xlon[ilon] > 180.:
            xlon[ilon] = xlon[ilon] - 360.
    ylat = lat[lat_idxs]

    with Dataset('./SSH_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
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
        dst['time'][:] = nc_salt['time'][8:17:8]
        #ssh
        dst.createVariable('ssh', 'f4', ('time', 'lat', 'lon',), fill_value=True)
        dst['ssh'].long_name = "sea_surface_elevation (m)"
        dst['ssh'][:,:,:] = nc_ssh['ssh'][8:17:8,0,lat_idxs,lon_idxs]
    
    with Dataset('./TS_1.nc', 'w', format='NETCDF3_CLASSIC') as dst: 
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
        dst['time'][:] = nc_salt['time'][1:3]
        #salt 
        dst.createVariable('salinity', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
        dst['salinity'].long_name = "sea_water_salinity (psu)" 
        #temp
        dst.createVariable('temperature', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
        dst['temperature'].long_name = "sea_water_potential_temperature (degc)"
        
        for it in np.arange(2):
            for k in np.arange(len(lev)):
                dst['salinity'][it,k,:,:] = nc_salt['salinity'][it,k,lat_idxs,lon_idxs]
                dst['temperature'][it,k,:,:] = nc_temp['temperature'][it,k,lat_idxs,lon_idxs]

    with Dataset('./UV_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
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
        dst['time'][:] = nc_salt['time'][1:3]
        #uvel
        dst.createVariable('u', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
        dst['u'].long_name = "eastward_sea_water_velocity (m/s)"
        #vvel
        dst.createVariable('v', 'f4', ('time', 'lev', 'lat', 'lon',), fill_value=True)
        dst['v'].long_name = "northward_sea_water_velocity (m/s)"

        for it in np.arange(2):
            for k in np.arange(len(lev)):
                dst['u'][it,k,:,:] = nc_salt['u'][it,k,lat_idxs,lon_idxs]
                dst['v'][it,k,:,:] = nc_temp['v'][it,k,lat_idxs,lon_idxs]
