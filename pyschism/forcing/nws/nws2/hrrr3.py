from datetime import datetime, timedelta
import pathlib
import tempfile
import logging
from time import time
from typing import Union
import os
import glob
import multiprocessing as mp

#from appdirs import user_data_dir
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import numpy as np
import pandas as pd
import xarray as xr

from pyschism.dates import nearest_cycle

logger = logging.getLogger(__name__)

class AWSGrib2Inventory:

    def __init__(
            self,
            start_date: datetime = None,
            record = 2,
            pscr = None, #tmpdir to save grib files
            product='conus',
    ):
        self.start_date = nearest_cycle() if start_date is None else start_date
        self.pscr = pscr
        self.product = product

        self.forecast_cycle = self.start_date #nearest_cycle()

        paginator=self.s3.get_paginator('list_objects_v2')
        pages=paginator.paginate(Bucket=self.bucket, 
                Prefix=f'hrrr.{self.forecast_cycle.strftime("%Y%m%d")}'
                       f'/{self.product}/')
        data=[]
        for page in pages:
            for obj in page['Contents']:
                data.append(obj) 

        self.cycle=self.forecast_cycle.hour
        tz='t{:02d}z'.format(self.cycle)
        self.file_metadata = list(sorted([
            _['Key'] for _ in data if 'wrfsfcf' in _['Key'] and tz in _['Key'] and not 'idx' in _['Key']
        ]))

        for key in self.file_metadata[1:record*24+1]:
            filename = pathlib.Path(self.tmpdir) / key
            filename.parent.mkdir(parents=True, exist_ok=True)

            with open(filename, 'wb') as f:
                logger.info(f'Downloading file {key}, ')
                try:
                    self.s3.download_fileobj(self.bucket, key, f)
                except:
                    logger.info(f'file {key} is not available')
        #return filename

    @property
    def bucket(self):
        return 'noaa-hrrr-bdp-pds'

    @property
    def output_interval(self) -> timedelta:
        return {
            'conus': timedelta(hours=1)
        }[self.product]

    @property
    def s3(self):
        try:
            return self._s3
        except AttributeError:
            self._s3 = boto3.client(
                's3', config=Config(signature_version=UNSIGNED))
            return self._s3

    @property
    def tmpdir(self):
        if not hasattr(self, "_tmpdir"):
            self._tmpdir = tempfile.TemporaryDirectory(dir=self.pscr)
        return pathlib.Path(self._tmpdir.name)

    @property
    def files(self):
        grbfiles=glob.glob(f'{self.tmpdir}/hrrr.{self.forecast_cycle.strftime("%Y%m%d")}/conus/hrrr.t{self.cycle:02d}z.wrfsfcf*.grib2')
        grbfiles.sort()
        return grbfiles

class HRRR:

    def __init__(self, start_date=None, rnday=None, pscr=None, record=2, bbox=None):

        start_date = nearest_cycle() if start_date is None else start_date 
        self.bbox = bbox

        end_date = start_date + timedelta(days=rnday)
        logger.info(f'start_date is {start_date}, end_date is {end_date}')
        
        datevector = np.arange(start_date, end_date, np.timedelta64(1, 'D'),
                dtype='datetime64')
        datevector = pd.to_datetime(datevector)

        npool = len(datevector) if len(datevector) < mp.cpu_count() else mp.cpu_count()
        logger.info(f'npool is {npool}')
        pool = mp.Pool(npool)

        pool.starmap(self.gen_sflux, [(date, record, pscr) for date in datevector])

        pool.close()
        
    def gen_sflux(self, date, record, pscr):

        inventory = AWSGrib2Inventory(date, record, pscr)
        grbfiles = inventory.files
        cycle = date.hour

        path = pathlib.Path(date.strftime("%Y%m%d"))
        path.mkdir(parents=True, exist_ok=True)

        stmp = [] 
        spfh = []
        uwind = []
        vwind = []
        prmsl = []
        prate = [] 
        dlwrf = [] 
        dswrf = [] 

        Vars = {
            'group1': {'sh2': ['174096', spfh], 't2m': ['167', stmp],'u10': ['165', uwind], 'v10': ['166', vwind]},
            'group2': {'mslma': ['meanSea', prmsl]},
            'group3': {'prate': ['surface', prate], 'dlwrf': ['surface', dlwrf], 'dswrf': ['surface', dswrf]}
        }

        #Get lon/lat
        lon, lat, idx_ymin, idx_ymax, idx_xmin, idx_xmax = self.modified_latlon(grbfiles[0])

        for ifile, file in enumerate(grbfiles):
            logger.info(f'file {ifile} is {file}')

            for key, value in Vars.items():
                if key == 'group1':
                    for key2, value2 in value.items():
                        ds=xr.open_dataset(file,
                            engine='cfgrib',
                            backend_kwargs=dict(filter_by_keys={'paramId':int(value2[0])}))
                        value2[1].append(ds[key2][idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32'))
                        ds.close()

                elif key == 'group2':
                    ds=xr.open_dataset(file,
                        engine='cfgrib',
                        backend_kwargs=dict(filter_by_keys={'typeOfLevel':'meanSea'}))
                    for key2, value2 in value.items():
                        value2[1].append(ds[key2][idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32'))
                    ds.close()

                else:
                    ds=xr.open_dataset(file,
                        engine='cfgrib',
                        backend_kwargs=dict(filter_by_keys={'stepType': 'instant','typeOfLevel': 'surface'}))
                    for key2, value2 in value.items():
                        value2[1].append(ds[key2][idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32'))
                    ds.close()    

        #write netcdf
        fout = xr.Dataset({
            'stmp': (['time', 'ny_grid', 'nx_grid'], np.array(stmp)),
            'spfh': (['time', 'ny_grid', 'nx_grid'], np.array(spfh)),
            'uwind': (['time', 'ny_grid', 'nx_grid'], np.array(uwind)),
            'vwind': (['time', 'ny_grid', 'nx_grid'], np.array(vwind)),
            'prmsl': (['time', 'ny_grid', 'nx_grid'], np.array(prmsl)),
            'prate': (['time', 'ny_grid', 'nx_grid'], np.array(prate)),
            'dlwrf': (['time', 'ny_grid', 'nx_grid'], np.array(dlwrf)),
            'dswrf': (['time', 'ny_grid', 'nx_grid'], np.array(dswrf)),
            },
            coords={
                'time': np.round(np.arange(1, len(grbfiles)+1)/24, 4).astype('float32'),
                'lon': (['ny_grid', 'nx_grid'], lon),
                'lat': (['ny_grid', 'nx_grid'], lat)})

        bdate = date.strftime('%Y %m %d %H').split(' ')
        bdate = [int(q) for q in bdate[:4]] + [0]

        fout.time.attrs = {
            'long_name': 'Time',
            'standard_name': 'time',
            'base_date': bdate,
            'units': f"days since {date.year}-{date.month}-{date.day} {cycle:02d}:00 UTC"
        }

        fout.lat.attrs = {
            'units': 'degrees_north',
            'long_name': 'Latitude',
            'standard_name': 'latitude',
        }

        fout.lon.attrs = {
            'units': 'degrees_east',
            'long_name': 'Longitude',
            'standard_name': 'longitude',
        }

        fout.uwind.attrs={
            'units': 'm/s',
            'long_name': '10m_above_ground/UGRD',
            'standard_name':'eastward_wind'
        }
      
        fout.vwind.attrs={
            'units': 'm/s',
            'long_name': '10m_above_ground/UGRD',
            'standard_name':'northward_wind'
        }

        fout.spfh.attrs={
            'units': 'kg kg-1',
            'long_name': '2m_above_ground/Specific Humidity',
            'standard_name':'specific_humidity'
        }

        fout.prmsl.attrs = {
            'units': 'Pa',
            'long_name': 'Pressure reduced to MSL',
            'standard_name': 'air_pressure_at_sea_level'
        }

        fout.stmp.attrs={
            'units': 'K',
            'long_name': '2m_above_ground/Temperature',
        }

        fout.prate.attrs={
            'units': 'kg m-2 s-1',
            'long_name': 'Precipitation rate'
        }
        fout.dlwrf.attrs = {
            'units': 'W m-2',
            'long_name': 'Downward short-wave radiation flux'
        }

        fout.dswrf.attrs = {
            'units': 'W m-2',
            'long_name': 'Downward long-wave radiation flux'
        }
                         
        fout.to_netcdf(path / f'hrrr_{date.strftime("%Y%m%d")}{cycle:02d}.nc','w', 'NETCDF3_CLASSIC', unlimited_dims='time')

    def modified_latlon(self, grbfile):

        xmin, xmax, ymin, ymax = self.bbox.xmin, self.bbox.xmax, self.bbox.ymin, self.bbox.ymax
        xmin = xmin + 360 if xmin < 0 else xmin
        xmax = xmax + 360 if xmax < 0 else xmax
       
        ds=xr.open_dataset(grbfile,
            engine='cfgrib',
            backend_kwargs=dict(filter_by_keys={'stepType': 'instant','typeOfLevel': 'surface'}))
        lon=ds.longitude.astype('float32')
        lon_idxs=(lon.values >= xmin) & (lon.values <= xmax)
        lat=ds.latitude.astype('float32')
        lat_idxs=(lat.values >= ymin) & (lat.values <= ymax + 2.0)
        idxs = lon_idxs & lat_idxs
        idxs = np.argwhere(idxs)
        idx_ymin = np.min(idxs[:,0])
        idx_ymax = np.max(idxs[:,0])
        idx_xmin = np.min(idxs[:,1])
        idx_xmax = np.max(idxs[:,1])
     
        lon2 = lon.values[idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1]
        lat2 = lat.values[idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1]
        idxs = np.where(lon2 > 180)
        lon2[idxs] -= 360
        logger.info(f'idx_ymin is {idx_ymin}, idx_ymax is {idx_ymax}, idx_xmin is {idx_xmin}, idx_xmax is {idx_xmax}')
        ds.close()

        return lon2, lat2, idx_ymin, idx_ymax, idx_xmin, idx_xmax
