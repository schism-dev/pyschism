from datetime import datetime, timedelta
import tempfile
import pathlib
from typing import Union
import logging

import appdirs
import numpy as np
import cdsapi
import netCDF4 as nc4
from netCDF4 import Dataset
import pandas as pd
from matplotlib.transforms import Bbox

from pyschism.forcing.nws.nws2.sflux import (
    SfluxDataset,
    AirComponent,
    PrcComponent,
    RadComponent,
)

logger = logging.getLogger(__name__)

class ERA5DataInventory:

    def __init__(self, start_date=None, rnday: Union[float, timedelta] = 4, bbox=None):

        self.start_date = start_date
        self.rnday = rnday
        self.end_date = self.start_date+timedelta(self.rnday + 1)
        self.client=cdsapi.Client()
        self._bbox = bbox

        r = self.client.retrieve(
            'reanalysis-era5-single-levels',
            {
            'variable':[
                '10m_u_component_of_wind','10m_v_component_of_wind','mean_sea_level_pressure',
                '2m_dewpoint_temperature','2m_temperature','mean_total_precipitation_rate',
                'mean_surface_downward_long_wave_radiation_flux','mean_surface_downward_short_wave_radiation_flux'
                ],
            'product_type':'reanalysis',
            'date':f"{self.start_date.strftime('%Y-%m-%d')}/{self.end_date.strftime('%Y-%m-%d')}",
            'time':[
                '00:00','01:00','02:00','03:00','04:00','05:00',
                '06:00','07:00','08:00','09:00','10:00','11:00',
                '12:00','13:00','14:00','15:00','16:00','17:00',
                '18:00','19:00','20:00','21:00','22:00','23:00'
                ],
            'area': [self._bbox.ymax+0.5, self._bbox.xmin-0.5, self._bbox.ymin-0.5, self._bbox.xmax+0.5], # North, West, South, East. Default: global
            'format':'netcdf'
            })
 
        filename = self.tmpdir / f"era5_{self.start_date.strftime('%Y%m%d')}.nc"
        r.download(filename)
        
    @property
    def tmpdir(self):
        if not hasattr(self, '_tmpdir'):
            self._tmpdir = tempfile.TemporaryDirectory()
        return pathlib.Path(self._tmpdir.name)

    @property
    def files(self):
        return sorted(list(self.tmpdir.glob('**/era5_*.nc')))

    @property
    def lon(self):
        if not hasattr(self, '_lon'):
            self._lon = Dataset(self.files[0]).variables['longitude'][:]
            if not hasattr(self, '_lat'):
                self._lat = Dataset(self.files[0]).variables['latitude'][::-1]
        return self._lon

    @property
    def lat(self):
        if not hasattr(self, '_lat'):
            self._lat = Dataset(self.files[0]).variables['latitude'][::-1]
            if not hasattr(self, '_lon'):
                self._lon = Dataset(self.files[0]).variables['longitude'][:]
        return self._lat

    def xy_grid(self):
        lon = []
        for x in self.lon:
            if x > 180:
                lon.append(x-360)
            else:
                lon.append(x)
        return np.meshgrid(np.array(lon), self.lat)

def put_sflux_fields(iday, date, timevector, ds, nx_grid, ny_grid, air, rad, prc, output_interval, OUTDIR):
    rt=pd.to_datetime(str(date))
    idx=np.where(rt == timevector)[0].item()
    times=[i/24 for i in np.arange(0, 25, output_interval)]

    if air is True:
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
            dst['prmsl'][:,:,:]=ds['msl'][idx:idx+25:output_interval, ::-1, :]

            # spfh
            dst.createVariable('spfh', 'f4', ('time', 'ny_grid', 'nx_grid'))
            dst['spfh'].long_name = "Surface Specific Humidity "\
                                    "(2m AGL)"
            dst['spfh'].standard_name = "specific_humidity"
            dst['spfh'].units = "1"
            #convert dewpoint to specific humidity
            d2m = ds['d2m'][idx:idx+25:output_interval, ::-1, :]
            msl = ds['msl'][idx:idx+25:output_interval, ::-1, :]
            Td = d2m - 273.15
            e1 = 6.112*np.exp((17.67*Td)/(Td + 243.5))
            spfh = (0.622*e1)/(msl*0.01 - (0.378*e1))
            dst['spfh'][:,:,:]=spfh

            # stmp
            dst.createVariable('stmp', 'f4', ('time', 'ny_grid', 'nx_grid'))
            dst['stmp'].long_name = "Surface Air Temperature (2m AGL)"
            dst['stmp'].standard_name = "air_temperature"
            dst['stmp'].units = "K"
            dst['stmp'][:,:,:]=ds['t2m'][idx:idx+25:output_interval, ::-1, :]

            # uwind
            dst.createVariable('uwind', 'f4', ('time', 'ny_grid', 'nx_grid'))
            dst['uwind'].long_name = "Surface Eastward Air Velocity "\
                "(10m AGL)"
            dst['uwind'].standard_name = "eastward_wind"
            dst['uwind'].units = "m/s"
            dst['uwind'][:,:,:]=ds['u10'][idx:idx+25:output_interval, ::-1, :]

            # vwind
            dst.createVariable('vwind', 'f4', ('time', 'ny_grid', 'nx_grid'))
            dst['vwind'].long_name = "Surface Northward Air Velocity "\
                "(10m AGL)"
            dst['vwind'].standard_name = "northward_wind"
            dst['vwind'].units = "m/s"
            dst['vwind'][:,:,:]=ds['v10'][idx:idx+25:output_interval, ::-1, :]

    if prc is True:
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
            dst['prate'][:,:,:]=ds['mtpr'][idx:idx+25:output_interval, ::-1, :]

    if rad is True:
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
            dst['dlwrf'][:,:,:]=ds['msdwlwrf'][idx:idx+25:output_interval, ::-1, :]

            # dwrf
            dst.createVariable('dswrf', 'f4', ('time', 'ny_grid', 'nx_grid'))
            dst['dswrf'].long_name = "Downward Long Wave Radiation Flux"
            dst['dswrf'].standard_name = "surface_downwelling_shortwave_flux_in_air"
            dst['dswrf'].units = "W/m^2"
            dst['dswrf'][:,:,:]=ds['msdwswrf'][idx:idx+25:output_interval, ::-1, :]


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

    def write(
        self,
        outdir,
        start_date: datetime = None,
        rnday: Union[float, timedelta] = 4,
        air: bool = True,
        prc: bool = True,
        rad: bool = True,
        bbox: Bbox = None,
        overwrite: bool=False,
        output_interval: int = 1,
    ):
        self.start_date=start_date
        self.rnday=rnday
        self.end_date=self.start_date+timedelta(self.rnday + 1)

        dates = {_: None for _ in np.arange(
            self.start_date,
            self.end_date,
            np.timedelta64(1, 'D'),
            dtype='datetime64')}

        #write sflux_inputs.txt
        with open(f'{outdir}/sflux_inputs.txt', 'w') as f:
            f.write('&sflux_inputs\n/\n')

        logger.info('Start downloading ERA5')
        self.inventory = ERA5DataInventory(
            self.start_date,
            self.rnday,
            bbox,
        )

        logger.info('Finished downloading ERA5')
    
        nx_grid, ny_grid = self.inventory.xy_grid()

        ds=Dataset(self.inventory.files[0])
        time1=ds['time']
        times=nc4.num2date(time1,units=time1.units,only_use_cftime_datetimes=False)

        for iday, date in enumerate(dates):
            put_sflux_fields(iday, date, times, ds, nx_grid, ny_grid, air=air, rad=rad, prc=prc, output_interval=output_interval, OUTDIR=outdir)

