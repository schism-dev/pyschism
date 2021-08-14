from datetime import datetime, timedelta
import tempfile
import pathlib
from typing import Union
import logging
from multiprocessing import Pool

import numpy as np
import cdsapi
from netCDF4 import Dataset
import pandas as pd

from pyschism.forcing.atmosphere.nws.nws2.sflux import (
    SfluxDataset,
    AirComponent,
    PrcComponent,
    RadComponent,
)

logger = logging.getLogger(__name__)

class ERA5DataInventory:

    def __init__(self, start_date=None, rnday: Union[float, timedelta] = 4, bbox=None):

        self.start_date=start_date
        self.rnday=rnday
        self.enddate=self.start_date+timedelta(self.rnday)
        self._files = {_: None for _ in np.arange(
            self.start_date,
            self.enddate+timedelta(days=1),
            np.timedelta64(1, 'D'),
            dtype='datetime64')}
        self._bbox=bbox
        self.client=cdsapi.Client()

        for requested_time, _ in self._files.items():
            logger.info(f'Requesting ERA5 data for time {requested_time}')
            self._files[requested_time] = self.request_data(requested_time)
        #print(self._files)

    def request_data(self, requested_time):

        rt=pd.to_datetime(str(requested_time))
        r = self.client.retrieve(
            'reanalysis-era5-single-levels',
            {
            'variable':[
                '10m_u_component_of_wind','10m_v_component_of_wind','mean_sea_level_pressure',
                '2m_dewpoint_temperature','2m_temperature','mean_total_precipitation_rate',
                'mean_surface_downward_long_wave_radiation_flux','mean_surface_downward_short_wave_radiation_flux'
                ],
            'product_type':'reanalysis',
            'date':f"{rt.strftime('%Y-%m-%d')}/{rt.strftime('%Y-%m-%d')}",
            'time':[
                '00:00','01:00','02:00','03:00','04:00','05:00',
                '06:00','07:00','08:00','09:00','10:00','11:00',
                '12:00','13:00','14:00','15:00','16:00','17:00',
                '18:00','19:00','20:00','21:00','22:00','23:00'
                ],
            'area': [self._bbox.ymax, self._bbox.xmin, self._bbox.ymin, self._bbox.xmax], # North, West, South, East. Default: global
            'format':'netcdf'
            })
 
        filename = self.tmpdir / f"era5_{rt.strftime('%Y%m%d')}.nc"
        #r.download(self.tmpdir / f"era5_{requested_time.strftime('%Y%m%d')}.nc")
        r.download(filename)
        #self.nc=Dataset(self.tmpdir / f"era5_{start_date.strftime('%Y%m%d')}.nc")
        return filename 
        
    @property
    def tmpdir(self):
        if not hasattr(self, '_tmpdir'):
            self._tmpdir = tempfile.TemporaryDirectory()
        return pathlib.Path(self._tmpdir.name)

    @property
    def files(self):
        return sorted(list(self.tmpdir.glob('**/*.nc')))
        #return self._files

    @property
    def lon(self):
        if not hasattr(self, '_lon'):
            #nc = self._files[list(self._files.keys())[0]]
            self._lon = Dataset(self.files[0]).variables['longitude'][:]
            if not hasattr(self, '_lat'):
                self._lat = Dataset(self.files[0]).variables['latitude'][:]
        return self._lon

    @property
    def lat(self):
        if not hasattr(self, '_lat'):
            self._lat = Dataset(self.files[0]).variables['latitude'][:]
            if not hasattr(self, '_lon'):
                self._lon = Dataset(self.files[0]).variables['longitude'][:]
        return self._lat

    def _modified_bbox(self, bbox=None):
        if bbox is None:
            return Bbox.from_extents(0, -90, 360, 90)
        else:
            xmin = bbox.xmin + 360 if bbox.xmin < 0 else bbox.xmin
            xmax = bbox.xmax + 360 if bbox.xmax < 0 else bbox.xmax
            return Bbox.from_extents(xmin, bbox.ymin, xmax, bbox.ymax)

    def _modified_bbox_indexes(self):
        lat_idxs = np.where((self.lat >= self._bbox.ymin)
                            & (self.lat <= self._bbox.ymax))[0]
        lon_idxs = np.where((self.lon >= self._bbox.xmin)
                            & (self.lon <= self._bbox.xmax))[0]
        return lon_idxs, lat_idxs

    def xy_grid(self):
        lon_idxs, lat_idxs = self._modified_bbox_indexes()
        lon = []
        for x in self.lon[lon_idxs]:
            if x > 180:
                lon.append(x-360)
            else:
                lon.append(x)
        return np.meshgrid(np.array(lon), self.lat[lat_idxs])

def put_sflux_fields(iday, file, nx_grid, ny_grid, lon_idxs, lat_idxs, OUTDIR):
    #print(iday)
    #print(file)
    requested_date=file[0]
    filename=file[1]
    nc=Dataset(filename)
    rt=pd.to_datetime(str(requested_date))
    times=[i/24 for i in np.arange(24)]
    with Dataset(OUTDIR / "sflux_air_1.{:04}.nc".format(iday+1), 'w', format='NETCDF3_CLASSIC') as dst:
        dst.setncatts({"Conventions": "CF-1.0"})
        # dimensions
        dst.createDimension('nx_grid', nx_grid.shape[1])
        dst.createDimension('ny_grid', ny_grid.shape[0])
        dst.createDimension('time', None)
        # variables
        # lon
        dst.createVariable('lon', 'f4', ('ny_grid', 'nx_grid'))
        dst['lon'].long_name = "Longitude"
        dst['lon'].standard_name = "longitude"
        dst['lon'].units = "degrees_east"
        dst['lon'][:] = nx_grid
        # lat
        dst.createVariable('lat', 'f4', ('ny_grid', 'nx_grid'))
        dst['lat'].long_name = "Latitude"
        dst['lat'].standard_name = "latitude"
        dst['lat'].units = "degrees_north"
        dst['lat'][:] = ny_grid
        # time
        dst.createVariable('time', 'f4', ('time',))
        dst['time'].long_name = 'Time'
        dst['time'].standard_name = 'time'
        dst['time'].units = f'days since {rt.year}-{rt.month}'\
                            f'-{rt.day} 00:00 UTC'
        dst['time'].base_date = (rt.year, rt.month, rt.day, 0)
        dst['time'][:] = times

        # prmsl
        dst.createVariable('prmsl', 'f4', ('time', 'ny_grid', 'nx_grid'))
        dst['prmsl'].long_name = "Pressure reduced to MSL"
        dst['prmsl'].standard_name = "air_pressure_at_sea_level"
        dst['prmsl'].units = "Pa"
        dst['prmsl'][:,:,:]=nc['msl'][:,lat_idxs,lon_idxs]

        # spfh
        dst.createVariable('spfh', 'f4', ('time', 'ny_grid', 'nx_grid'))
        dst['spfh'].long_name = "Surface Specific Humidity "\
                                "(2m AGL)"
        dst['spfh'].standard_name = "specific_humidity"
        dst['spfh'].units = "1"
        dst['spfh'][:,:,:]=nc['d2m'][:,lat_idxs,lon_idxs]

        # stmp
        dst.createVariable('stmp', 'f4', ('time', 'ny_grid', 'nx_grid'))
        dst['stmp'].long_name = "Surface Air Temperature (2m AGL)"
        dst['stmp'].standard_name = "air_temperature"
        dst['stmp'].units = "K"
        dst['stmp'][:,:,:]=nc['t2m'][:,lat_idxs,lon_idxs]

        # uwind
        dst.createVariable('uwind', 'f4', ('time', 'ny_grid', 'nx_grid'))
        dst['uwind'].long_name = "Surface Eastward Air Velocity "\
            "(10m AGL)"
        dst['uwind'].standard_name = "eastward_wind"
        dst['uwind'].units = "m/s"
        dst['uwind'][:,:,:]=nc['u10'][:,lat_idxs,lon_idxs]

        # vwind
        dst.createVariable('vwind', 'f4', ('time', 'ny_grid', 'nx_grid'))
        dst['vwind'].long_name = "Surface Northward Air Velocity "\
            "(10m AGL)"
        dst['vwind'].standard_name = "northward_wind"
        dst['vwind'].units = "m/s"
        dst['vwind'][:,:,:]=nc['v10'][:,lat_idxs,lon_idxs]

    with Dataset(OUTDIR / "sflux_prc_1.{:04}.nc".format(iday+1), 'w', format='NETCDF3_CLASSIC') as dst:
        dst.setncatts({"Conventions": "CF-1.0"})
        # dimensions
        dst.createDimension('nx_grid', nx_grid.shape[1])
        dst.createDimension('ny_grid', ny_grid.shape[0])
        dst.createDimension('time', None)
        # variables
        # lon
        dst.createVariable('lon', 'f4', ('ny_grid', 'nx_grid'))
        dst['lon'].long_name = "Longitude"
        dst['lon'].standard_name = "longitude"
        dst['lon'].units = "degrees_east"
        dst['lon'][:] = nx_grid
        # lat
        dst.createVariable('lat', 'f4', ('ny_grid', 'nx_grid'))
        dst['lat'].long_name = "Latitude"
        dst['lat'].standard_name = "latitude"
        dst['lat'].units = "degrees_north"
        dst['lat'][:] = ny_grid
        # time
        dst.createVariable('time', 'f4', ('time',))
        dst['time'].long_name = 'Time'
        dst['time'].standard_name = 'time'
        dst['time'].units = f'days since {rt.year}-{rt.month}'\
                            f'-{rt.day} 00:00 UTC'
        dst['time'].base_date = (rt.year, rt.month, rt.day, 0)
        dst['time'][:] = times

        # prate
        dst.createVariable('prate', 'f4', ('time', 'ny_grid', 'nx_grid'))
        dst['prate'].long_name = "Surface Precipitation Rate"
        dst['prate'].standard_name = "air_pressure_at_sea_level"
        dst['prate'].units = "kg/m^2/s"
        dst['prate'][:,:,:]=nc['mtpr'][:,lat_idxs,lon_idxs]

    with Dataset(OUTDIR / "sflux_rad_1.{:04}.nc".format(iday+1), 'w', format='NETCDF3_CLASSIC') as dst:
        dst.setncatts({"Conventions": "CF-1.0"})
        # dimensions
        dst.createDimension('nx_grid', nx_grid.shape[1])
        dst.createDimension('ny_grid', ny_grid.shape[0])
        dst.createDimension('time', None)
        # variables
        # lon
        dst.createVariable('lon', 'f4', ('ny_grid', 'nx_grid'))
        dst['lon'].long_name = "Longitude"
        dst['lon'].standard_name = "longitude"
        dst['lon'].units = "degrees_east"
        dst['lon'][:] = nx_grid
        # lat
        dst.createVariable('lat', 'f4', ('ny_grid', 'nx_grid'))
        dst['lat'].long_name = "Latitude"
        dst['lat'].standard_name = "latitude"
        dst['lat'].units = "degrees_north"
        dst['lat'][:] = ny_grid
        # time
        dst.createVariable('time', 'f4', ('time',))
        dst['time'].long_name = 'Time'
        dst['time'].standard_name = 'time'
        dst['time'].units = f'days since {rt.year}-{rt.month}'\
                            f'-{rt.day} 00:00 UTC'
        dst['time'].base_date = (rt.year, rt.month, rt.day, 0)
        dst['time'][:] = times

        # dlwrf
        dst.createVariable('dlwrf', 'f4', ('time', 'ny_grid', 'nx_grid'))
        dst['dlwrf'].long_name = "Downward Long Wave Radiation Flux"
        dst['dlwrf'].standard_name = "surface_downwelling_longwave_flux_in_air"
        dst['dlwrf'].units = "W/m^2"
        dst['dlwrf'][:,:,:]=nc['msdwlwrf'][:,lat_idxs,lon_idxs]

        # dwrf
        dst.createVariable('dswrf', 'f4', ('time', 'ny_grid', 'nx_grid'))
        dst['dswrf'].long_name = "Downward Long Wave Radiation Flux"
        dst['dswrf'].standard_name = "surface_downwelling_shortwave_flux_in_air"
        dst['dswrf'].units = "W/m^2"
        dst['dswrf'][:,:,:]=nc['msdwswrf'][:,lat_idxs,lon_idxs]


class ERA5(SfluxDataset):

    def __init__(self, **kwargs):

        self.prmsl_name = 'prmslmsl'
        self.spfh_name = 'spfh2m'
        self.stmp_name = 'tmpsfc'
        self.uwind_name = 'ugrd10m'
        self.vwind_name = 'vgrd10m'
        self.prate_name = 'pratesfc'
        self.dlwrf_name = 'dlwrfsfc'
        self.dswrf_name = 'dswrfsfc'    
        self.air = None
        self.prc = None
        self.rad = None

    def gen_sflux(
        self,
        start_date: datetime = None,
        rnday: Union[float, timedelta] = 4,
        air: bool = True,
        prc: bool = True,
        rad: bool = True,
        bbox = None,
        outdir = None,
        nprocs=32,
    ):
        self.start_date=start_date
        self.rnday=rnday

        self.inventory = ERA5DataInventory(
            self.start_date,
            self.rnday,
            bbox,
        )
    
        nx_grid, ny_grid = self.inventory.xy_grid()
        lon_idxs, lat_idxs = self.inventory._modified_bbox_indexes()

        #outdir=pathlib.Path('./ERA5')
        #print(outdir)
        with Pool(processes=nprocs) as pool:
            pool.starmap(put_sflux_fields, [(iday, file, nx_grid, ny_grid, lon_idxs, lat_idxs, outdir)
                #for requested_date,file in self.inventory._files.items()])
                for iday,file in enumerate(self.inventory._files.items())])

