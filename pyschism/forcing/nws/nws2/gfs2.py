import os
from datetime import datetime, timedelta
import logging
import pathlib
import tempfile
from time import time
import glob
import multiprocessing as mp

import boto3
from botocore import UNSIGNED
from botocore.config import Config
import numpy as np
import xarray as xr
import pandas as pd

from pyschism.dates import nearest_cycle

logger = logging.getLogger(__name__)

class AWSGrib2Inventory:

    def __init__(
            self,
            start_date: datetime = None,
            record = 1,
            pscr = None,
            product='atmos',
    ):
        """
        Download GFS data from AWS.
        This dataset GFS V16.3.0 starts on Feb 26, 2021.
        """
        self.start_date = nearest_cycle() if start_date is None else start_date
        self.pscr = pscr
        self.product = product

        self.forecast_cycle = self.start_date

        timevector = np.arange(
            self.start_date,
            self.start_date + timedelta(hours=record*24+1),
            np.timedelta64(1, 'h')
        ).astype(datetime)

        file_metadata = self.get_file_namelist(timevector)

        for dt in timevector:
            
            outfile_name = f"gfs.{self.start_date.strftime('%Y%m%d')}/gfs.pgrb2.0p25.{dt.strftime('%Y%m%d%H')}.grib2"
            filename = pathlib.Path(self.tmpdir) / outfile_name
            filename.parent.mkdir(parents=True, exist_ok=True)

            with open(filename, 'wb') as f:
                while(file_metadata[dt]):
                    try:
                        key = file_metadata[dt].pop(0)
                        logger.info(f"Downloading file {key} for {dt}") 
                        self.s3.download_fileobj(self.bucket, key, f)
                        logger.info("Success!")
                        break
                    except:
                        if not file_metadata[dt]:
                            logger.info(f'No file for {dt}')
                            if os.path.exists(filename):
                                os.remove(filename)
                        else:
                            logger.info(f'file {key} is not available, try next file')
                            continue

    def get_file_namelist(self, requested_dates):

        file_metadata = {}
        hours = (datetime.utcnow() - self.start_date).days * 24 + (datetime.utcnow() - self.start_date).seconds // 3600
        n_cycles = hours // 6 if hours < 25 else 4
        cycle_index = (int(self.start_date.hour) - 1) // 6
        dt = requested_dates[0]
        for it, dt in enumerate(requested_dates[:n_cycles*6+1]):
            levels = 3
            i = 0
            fhour = int(dt.hour)
            cycle_index = (fhour - 1) // 6
            date2 = (dt - timedelta(days=1)).strftime('%Y%m%d') if dt.hour == 0 else dt.strftime('%Y%m%d')

            while (levels):

                cycle = self.fcst_cycles[cycle_index - i]
                fhour2 = fhour + i * 6 if cycle_index == 0 else fhour - cycle_index * 6 + i * 6
           
                file_metadata.setdefault(dt, []).append(f"gfs.{date2}/{cycle}/{self.product}/gfs.t{cycle}z.pgrb2.0p25.f{fhour2:03d}")
                levels -= 1
                i += 1

        if it < 25:
            date2 = (dt - timedelta(days=1)).strftime('%Y%m%d') if dt.hour == 0 else dt.strftime('%Y%m%d')
            for it, dt in enumerate(requested_dates[n_cycles*6+1:]):
                levels = 3
                i = 0
                hours = (dt - self.forecast_cycle).days * 24 + (dt - self.forecast_cycle).seconds // 3600
                
                while (levels):
                    #starting from the last cycle
                    cycle = self.fcst_cycles[cycle_index - i]
                    fhour2 = (hours - (n_cycles - 1) * 6) + i * 6
                    file_metadata.setdefault(dt, []).append(f"gfs.{date2}/{cycle}/{self.product}/gfs.t{cycle}z.pgrb2.0p25.f{fhour2:03d}")
                    levels -= 1
                    i += 1

        return file_metadata

    @property
    def bucket(self):
        return 'noaa-gfs-bdp-pds'

    @property
    def fcst_cycles(self):
        return ['00', '06', '12', '18']

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
        grbfiles=glob.glob(f'{self.tmpdir}/gfs.{self.forecast_cycle.strftime("%Y%m%d")}/gfs.pgrb2.0p25.*.grib2')
        grbfiles.sort()
        return grbfiles

class GFS:
    def __init__(self, level=1, bbox=None, pscr=None, outdir=None): 
        self.level = level
        self.bbox = bbox
        self.pscr = pscr
        self.outdir = outdir
        self.record = 1

    def write(self, start_date, rnday, air: bool=True, prc: bool=True, rad: bool=True):

        start_date = nearest_cycle() if start_date is None else start_date 

        if (start_date + timedelta(days=rnday)) > datetime.utcnow():
            logger.info(f'End date is beyond the current time, set rnday to 1 day and record to 5 days')
            rnday = 1
            self.record = 5 #days

        end_date = start_date + timedelta(hours=rnday * self.record + 1)
        logger.info(f'start time is {start_date}, end time is {end_date}')

        if self.outdir is None:
            #self.outdir = pathlib.Path(start_date.strftime("%Y%m%d"))
            self.outdir = pathlib.Path("sflux")
            self.outdir.mkdir(parents=True, exist_ok=True)
        
        datevector = np.arange(
            start_date, 
            start_date + timedelta(days=rnday), 
            np.timedelta64(1, 'D'),
            dtype='datetime64',
         )
        datevector = pd.to_datetime(datevector)

        npool = len(datevector) if len(datevector) < mp.cpu_count()/2 else mp.cpu_count()/2
        logger.info(f'npool is {npool}')
        pool = mp.Pool(int(npool))

        pool.starmap(self.gen_sflux, [(istack+1, date, air, prc, rad) for istack, date in enumerate(datevector)])
        pool.close()
        #self.gen_sflux(1, datevector[0], air, prc, rad)

    def gen_sflux(
            self, 
            istack,
            date,
            air: bool = True,
            prc: bool = True,
            rad: bool = True,
        ):

        inventory = AWSGrib2Inventory(date, self.record, self.pscr)
        grbfiles = inventory.files
        #cycle = date.hour
        
        prate = []
        dlwrf = []
        dswrf = []
        stmp = []
        spfh = []
        uwind = []
        vwind = []
        prmsl = []

        Vars = {
            'group1': {'sh2': ['174096', spfh], 't2m': ['167', stmp],
            'u10': ['165', uwind], 'v10': ['166', vwind]},
            'group2': {'prmsl': ['meanSea', prmsl]},
            'group3': {'prate': ['surface', prate]},
            'group4': {'dlwrf': ['surface', dlwrf], 'dswrf': ['surface', dswrf]}
        }

        #Get lon/lat
        lon, lat, idx_ymin, idx_ymax, idx_xmin, idx_xmax = self.modified_latlon(grbfiles[0])

        times = []
        for ifile, file in enumerate(grbfiles):
            logger.info(f'file {ifile} is {file}')

            for key, value in Vars.items():
                if key == 'group1':
                    for key2, value2 in value.items():
                        ds = xr.open_dataset(file,
                                engine='cfgrib',
                                backend_kwargs=dict(filter_by_keys={'paramId':int(value2[0])}))
                        tmp = ds[key2][idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32')
                        value2[1].append(tmp[::-1,:])
                        ds.close()

                elif key == 'group2':
                    ds = xr.open_dataset(file,
                            engine='cfgrib',
                            backend_kwargs=dict(filter_by_keys={'typeOfLevel':'meanSea'}))
                    for key2, value2 in value.items():
                        tmp = ds[key2][idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32')
                        value2[1].append(tmp[::-1,:])
                    times.append(ds.valid_time.values)
                    ds.close()

                elif key == 'group3':
                    ds = xr.open_dataset(file,
                            engine='cfgrib',
                            backend_kwargs=dict(filter_by_keys={'stepType': 'instant','typeOfLevel': 'surface'}))
                    for key2, value2 in value.items():
                        tmp = ds[key2][idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32')
                        value2[1].append(tmp[::-1, :])
                    ds.close()

                else:
                    ds = xr.open_dataset(file,
                            engine='cfgrib',
                            backend_kwargs=dict(filter_by_keys={'stepType': 'avg','typeOfLevel': 'surface'}))
                    for key2, value2 in value.items():
                        tmp = ds[key2][idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32')
                        value2[1].append(tmp[::-1, :])
                    ds.close()

        #write to netcdf
        bdate = date.strftime('%Y %m %d %H').split(' ')
        bdate = [int(q) for q in bdate[:4]] + [0]

        if air:
            ds = xr.Dataset({
                'stmp': (['time', 'ny_grid', 'nx_grid'], np.array(stmp)),
                'spfh': (['time', 'ny_grid', 'nx_grid'], np.array(spfh)),
                'uwind': (['time', 'ny_grid', 'nx_grid'], np.array(uwind)),
                'vwind': (['time', 'ny_grid', 'nx_grid'], np.array(vwind)),
                'prmsl': (['time', 'ny_grid', 'nx_grid'], np.array(prmsl)),
            },
                coords={
                    'time': np.round((times - date.to_datetime64()) / np.timedelta64(1, 'D'), 5).astype('float32'),
                    'lon': (['ny_grid', 'nx_grid'], lon),
                    'lat': (['ny_grid', 'nx_grid'], lat)
            })

            ds.time.attrs = { 
                'long_name': 'Time',
                'standard_name': 'time',
                'base_date': bdate,
                'units': f"days since {date.strftime('%Y-%m-%d %H:00')} UTC",
            }

            ds.lat.attrs = {
                'units': 'degrees_north',
                'long_name': 'Latitude',
                'standard_name': 'latitude',
            }

            ds.lon.attrs = {
                'units': 'degrees_east',
                'long_name': 'Longitude',
                'standard_name': 'longitude',
            }

            ds.uwind.attrs={
                'units': 'm/s',
                'long_name': '10m_above_ground/UGRD',
                'standard_name':'eastward_wind'
            }

            ds.vwind.attrs={
                'units': 'm/s',
                'long_name': '10m_above_ground/VGRD',
                'standard_name':'northward_wind'
            }
            
            ds.spfh.attrs={
                'units': 'kg kg-1',
                'long_name': '2m_above_ground/Specific Humidity',
                'standard_name':'specific_humidity'
            }

            ds.prmsl.attrs = {
                'units': 'Pa',
                'long_name': 'Pressure reduced to MSL',
                'standard_name': 'air_pressure_at_sea_level'
            }

            ds.stmp.attrs={
                'units': 'K',
                'long_name': '2m_above_ground/Temperature',
            }

            ds.to_netcdf(f'{self.outdir}/sflux_air_{self.level}.{istack:04d}.nc','w', 'NETCDF3_CLASSIC', unlimited_dims='time')
            ds.close()

        if prc:
            ds = xr.Dataset({
                'prate': (['time', 'ny_grid', 'nx_grid'], np.array(prate)),
            },
                coords={
                    'time': np.round((times - date.to_datetime64()) / np.timedelta64(1, 'D'), 5).astype('float32'),
                    'lon': (['ny_grid', 'nx_grid'], lon),
                    'lat': (['ny_grid', 'nx_grid'], lat)
            })

            ds.time.attrs = { 
                'long_name': 'Time',
                'standard_name': 'time',
                'base_date': bdate,
                'units': f"days since {date.strftime('%Y-%m-%d %H:00')} UTC"
            }

            ds.lat.attrs = {
                'units': 'degrees_north',
                'long_name': 'Latitude',
                'standard_name': 'latitude',
            }

            ds.lon.attrs = {
                'units': 'degrees_east',
                'long_name': 'Longitude',
                'standard_name': 'longitude',
            }

            ds.prate.attrs={
                'units': 'kg m-2 s-1',
                'long_name': 'Precipitation rate'
            }

            ds.to_netcdf(f'{self.outdir}/sflux_prc_{self.level}.{istack:04d}.nc','w', 'NETCDF3_CLASSIC', unlimited_dims='time')
            ds.close()

        if rad:
            ds = xr.Dataset({
                'dlwrf': (['time', 'ny_grid', 'nx_grid'], np.array(dlwrf)),
                'dswrf': (['time', 'ny_grid', 'nx_grid'], np.array(dswrf)),
            },
                coords={
                    'time': np.round((times - date.to_datetime64()) / np.timedelta64(1, 'D'), 5).astype('float32'),
                    'lon': (['ny_grid', 'nx_grid'], lon),
                    'lat': (['ny_grid', 'nx_grid'], lat)
            })

            ds.time.attrs = { 
                'long_name': 'Time',
                'standard_name': 'time',
                'base_date': bdate,
                'units': f"days since {date.strftime('%Y-%m-%d %H:00')} UTC"
            }

            ds.lat.attrs = {
                'units': 'degrees_north',
                'long_name': 'Latitude',
                'standard_name': 'latitude',
            }

            ds.lon.attrs = {
                'units': 'degrees_east',
                'long_name': 'Longitude',
                'standard_name': 'longitude',
            }

            ds.dlwrf.attrs = {
                'units': 'W m-2',
                'long_name': 'Downward short-wave radiation flux'
            }

            ds.dswrf.attrs = {
                'units': 'W m-2',
                'long_name': 'Downward long-wave radiation flux'
            }
                             
            ds.to_netcdf(f'{self.outdir}/sflux_rad_{self.level}.{istack:04d}.nc','w', 'NETCDF3_CLASSIC', unlimited_dims='time')
            ds.close()

    def modified_latlon(self, grbfile):
        xmin, xmax, ymin, ymax = self.bbox.xmin, self.bbox.xmax, self.bbox.ymin, self.bbox.ymax
        xmin = xmin + 360 if xmin < 0 else xmin
        xmax = xmax + 360 if xmax < 0 else xmax

        ds=xr.open_dataset(grbfile,
                engine='cfgrib',
                backend_kwargs=dict(filter_by_keys={'stepType': 'instant','typeOfLevel': 'surface'}))
        lon=ds.longitude.astype('float32')
        lat=ds.latitude.astype('float32')
        lon_idxs=np.where((lon.values >= xmin-1.0) & (lon.values <= xmax+1.0))[0]
        lat_idxs=np.where((lat.values >= ymin-1.0) & (lat.values <= ymax+1.0))[0]
        idx_ymin = lat_idxs[0]
        idx_ymax = lat_idxs[-1]
        idx_xmin = lon_idxs[0]
        idx_xmax = lon_idxs[-1]
        lon2 = lon[lon_idxs]
        lat2 = lat[lat_idxs]
        idxs = np.where(lon2 > 180)
        lon2[idxs] -= 360
        logger.info(f'idx_ymin is {idx_ymin}, idx_ymax is {idx_ymax}, idx_xmin is {idx_xmin}, idx_xmax is {idx_xmax}')
        #make sure lat is in ascending order
        nx_grid, ny_grid=np.meshgrid(lon2, lat2[::-1])

        ds.close()

        return nx_grid, ny_grid, idx_ymin, idx_ymax, idx_xmin, idx_xmax
