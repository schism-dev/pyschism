from datetime import datetime, timedelta
import pathlib
import os
from typing import Union
import logging

import numpy as np
import pandas as pd
import xarray as xr
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import fsspec
import zarr
from matplotlib.transforms import Bbox

from pyschism.dates import nearest_cycle

logger = logging.getLogger(__name__)

class AWSZarrInventory:
    '''HRRR AWS Zarr bucket
    Create sflux from HRRr Zarr files.

    Notes
    -------
    DSWRF has missing data:
    https://mesowest.utah.edu/html/hrrr/zarr_documentation/html/zarr_variables.html
    If only need air/prc, this method is faster than generating sflux from grib2 files (hrrr3).
    '''

    def __init__(self, bbox: Union[list, Bbox]=None):
        '''
        Parameters
        -------
        bbox: matplotlib.transforms.Bbox object or a list [xmin, ymin, xmax, ymax]
        '''

        if bbox is None:
            raise ValueError('Bbox is needed!')

        if not isinstance(bbox, Bbox) and not isinstance(bbox, list):
            raise TypeError('bbox is a list or matplotlib.transforms.Bbox object!')

        if isinstance(bbox, list):
            xmin, ymin, xmax, ymax = [v for v in bbox]
            self.bbox =  Bbox.from_extents(xmin, ymin, xmax, ymax)
        else:
            self.bbox = bbox

    def gen_sflux(self, outdir: Union[str, os.PathLike], start_date=None, stack=None, air=True, prc=True, rad=True): 
        self.start_date = nearest_cycle() if start_date is None else start_date
        airVars = {'stmp': ['TMP', '2m_above_ground'],
                   'spfh': ['SPFH', '2m_above_ground'],
                   'uwind': ['UGRD', '10m_above_ground'],
                   'vwind': ['VGRD', '10m_above_ground'],
                   'prmsl': ['MSLMA', 'mean_sea_level']}
        prcVars = {'prate': ['PRATE', 'surface']}
        radVars = {'dlwrf': ['DLWRF', 'surface'],
                   'dswrf': ['DSWRF', 'surface']}

        lon, lat, idx_ymin, idx_ymax, idx_xmin, idx_xmax = self.modified_latlon() 
        path = pathlib.Path(outdir)
        path.mkdir(parents=True, exist_ok=True)

        self.cycle = start_date.hour
        base_url = f's3://hrrrzarr/sfc/{self.start_date.strftime("%Y%m%d")}/' \
            + f'{self.start_date.strftime("%Y%m%d")}_{self.cycle:02d}z_fcst.zarr'

        if air:
            xdat = []
            xdat.append(self.lon_DataArray(lon1D=lon))
            xdat.append(self.lat_DataArray(lat1D=lat))
            for key, value in airVars.items():
                url = f'{base_url}/{value[1]}/{value[0]}/{value[1]}/'
                logger.info(url)
                ds = xr.open_zarr(fsspec.get_mapper(url, anon=True), consolidated=False)
                ds = ds.rename({'projection_y_coordinate': 'ny_grid'})
                ds = ds.rename({'projection_x_coordinate': 'nx_grid'})
                val = ds[value[0]][:, idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32')
                xdat.append(xr.DataArray(val, dims=("time", "ny_grid", "nx_grid",), name=key))

                #add variable time
                times = ds['time'].values.astype('float32')
               
                #zarr data already discard t00z, which has all junk values
                times = (times+1)/24
                xdat.append(self.time_DataArray(time1D=times))

            fout = xr.merge(xdat)
            fout = fout.rename_vars({'time2': 'time'})
            #fout.to_netcdf(path / f'hrrr_{start_date.strftime("%Y%m%d")}{cycle:02d}.nc', 'w', unlimited_dims='time')
            fout.to_netcdf(path / f'sflux_air_2.{stack:04d}.nc', 'w', unlimited_dims='time')
            fout.close()

        if prc:
            xdat = []
            xdat.append(self.lon_DataArray(lon1D=lon))
            xdat.append(self.lat_DataArray(lat1D=lat))
            for key, value in prcVars.items():
                url = f'{base_url}/{value[1]}/{value[0]}/{value[1]}/'
                logger.info(url)
                ds = xr.open_zarr(fsspec.get_mapper(url, anon=True), consolidated=False)
                ds = ds.rename({'projection_y_coordinate': 'ny_grid'})
                ds = ds.rename({'projection_x_coordinate': 'nx_grid'})
                val = ds[value[0]][:, idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32')
                xdat.append(xr.DataArray(val, dims=("time", "ny_grid", "nx_grid",), name=key))

                #add variable time
                times = ds['time'].values.astype('float32')
                times = (times+1)/24
                #bdate = pd.to_datetime(start_date + timedelta(hours=cycle)).strftime('%Y %m %d %H').split(' ')
                #bdate = [int(q) for q in bdate[:4]] + [0]
                xdat.append(self.time_DataArray(time1D=times))
            fout = xr.merge(xdat)
            fout = fout.rename_vars({'time2': 'time'})
            #fout.to_netcdf(path / f'hrrr_{start_date.strftime("%Y%m%d")}{cycle:02d}.nc', 'w', unlimited_dims='time')
            fout.to_netcdf(path / f'sflux_prc_2.{stack:04d}.nc', 'w', unlimited_dims='time')
            fout.close()

        if rad:
            xdat = []
            xdat.append(self.lon_DataArray(lon1D=lon))
            xdat.append(self.lat_DataArray(lat1D=lat))
            for key, value in radVars.items():
                url = f'{base_url}/{value[1]}/{value[0]}/{value[1]}/'
                logger.info(url)
                ds = xr.open_zarr(fsspec.get_mapper(url, anon=True), consolidated=False)
                ds = ds.rename({'projection_y_coordinate': 'ny_grid'})
                ds = ds.rename({'projection_x_coordinate': 'nx_grid'})
                val = ds[value[0]][:, idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32')
                #DSWRF are all junk values in this dataset
                #if np.sum(np.isnan(val)) > 0:
                #    val[:] = 0.0
                xdat.append(xr.DataArray(val, dims=("time", "ny_grid", "nx_grid",), name=key))

                #add variable time
                times = ds['time'].values.astype('float32')
                times = (times+1)/24
                #bdate = pd.to_datetime(start_date + timedelta(hours=cycle)).strftime('%Y %m %d %H').split(' ')
                #bdate = [int(q) for q in bdate[:4]] + [0]
                xdat.append(self.time_DataArray(time1D=times))
            fout = xr.merge(xdat)
            fout = fout.rename_vars({'time2': 'time'})
            #fout.to_netcdf(path / f'hrrr_{start_date.strftime("%Y%m%d")}{cycle:02d}.nc', 'w', unlimited_dims='time')
            fout.to_netcdf(path / f'sflux_rad_2.{stack:04d}.nc', 'w', unlimited_dims='time')
            fout.close()
        
        #fout = xr.merge(xdat)
        #fout=fout.rename_vars({'time2': 'time'})
        #fout.to_netcdf(path / f'hrrr_{start_date.strftime("%Y%m%d")}{cycle:02d}.nc', 'w', unlimited_dims='time')

    def write(
        self,
        outdir,
        start_date: datetime = None,
        rnday: Union[float, timedelta] = 2,
        air: bool = True,
        prc: bool = False,
        rad: bool = False,
    ):

        start_date = nearest_cycle() if start_date is None else start_date

        if not isinstance(rnday, timedelta):
            try:
                rnday = timedelta(days=rnday)
            except:
                raise TypeError("rnday must be either float or datetime.timedelta!")

        timevector = np.arange(
            start_date,
            rnday,
            np.timedelta64(1, "D"),
            dtype="datetime64",
        ).astype(datetime)

        for i, date in enumerate(timevector):
            self.gen_sflux(outdir, date, i+1, air=air, prc=prc, rad=rad)

    def modified_latlon(self):
        ds = xr.open_dataset(self.grid_file)
        lon = ds['longitude'].astype("float32")
        lat = ds['latitude'].astype("float32")
        xmin = self.bbox.xmin
        xmax = self.bbox.xmax
        ymin = self.bbox.ymin
        ymax = self.bbox.ymax
        lon_idxs=(lon.values >= xmin) & (lon.values <= xmax)
        lat_idxs=(lat.values >= ymin) & (lat.values <= ymax+2.0)
        idxs = lon_idxs & lat_idxs
        idxs = np.argwhere(idxs)
        idx_ymin = np.min(idxs[:,0])
        idx_ymax = np.max(idxs[:,0])
        idx_xmin = np.min(idxs[:,1])
        idx_xmax = np.max(idxs[:,1])
        #print(f'idx_ymin is {idx_ymin}, idx_ymax is {idx_ymax}, idx_xmin is {idx_xmin}, idx_xmax is {idx_xmax}')
        return lon[idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1], lat[idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1], idx_ymin, idx_ymax, idx_xmin, idx_xmax


    def lon_DataArray(self, lon1D):
        return xr.DataArray(
            data=lon1D, 
            dims=("ny_grid", "nx_grid",), 
            name="lon",
            attrs=dict(
                long_name = "Longitude",
                standard_name = "longitude",
                units = "degrees_east"
            )
        )

    def lat_DataArray(self, lat1D):
        return xr.DataArray(
            data=lat1D, 
            dims=("ny_grid", "nx_grid",), 
            name="lat",
            attrs=dict(
                long_name = "Latitude",
                standard_name = "latitude",
                units = "degrees_north"
            )
        )

    def time_DataArray(self, time1D):
        #bdate = pd.to_datetime(self.start_date + timedelta(hours=self.cycle)).strftime('%Y %m %d %H').split(' ')
        bdate = pd.to_datetime(self.start_date).strftime('%Y %m %d %H').split(' ')
        bdate = [int(q) for q in bdate[:4]] + [0]
        
        return xr.DataArray(
            time1D, 
            dims=('time'), 
            name='time2',
            attrs=dict(
                long_name="Time",
                standard_name = "time",
                base_date = bdate,
                units = f"days since {self.start_date.year}-{self.start_date.month}-{self.start_date.day} {self.cycle:02d}:00 UTC"
            )
        )

    @property
    def s3(self):
        try:
            return self._s3
        except AttributeError:
            self._s3 = boto3.client(
                "s3", config = Config(signature_version=UNSIGNED))
            return self._s3

    @property
    def bucket(self):
        return "hrrrzarr"

    @property
    def grid_file(self):
        #if not hasattr(self, "_grid_file"):
        self._grid_file = "HRRR_latlon.h5"
        with open(self._grid_file, "wb") as f:
            self.s3.download_fileobj(f"{self.bucket}", "grid/HRRR_latlon.h5", f)
        return self._grid_file
