from datetime import datetime, timedelta
import pathlib
import os
from typing import Union

import numpy as np
import pandas as pd
import xarray as xr
import fsspec
import zarr

from pyschism.dates import nearest_cycle

class AWSZarrInventory:
    '''
    12/21/2021, L. Cui
    There is an issue with DSWRF in this dataset. After communicating with 
    MesoWest group, it is confirmed that they are unable to successfully 
    generate zarr-formatted output for this vaiable for the zarr-formatted
    forecast files. At this point, I'll keep this script in the repo. 
    '''

    def __init__(self, bbox=None):
        self.bbox = bbox

    def gen_sflux(self, outdir: Union[str, os.PathLike], start_date=None, air=True, prc=True, rad=True):

        start_date = nearest_cycle() if start_date is None else start_date
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

        #for cycle in range(0, 24, 6):

        cycle = start_date.hour
        base_url = f's3://hrrrzarr/sfc/{start_date.strftime("%Y%m%d")}/' \
            + f'{start_date.strftime("%Y%m%d")}_{cycle:02d}z_fcst.zarr'

        xdat = []
        lon2 = xr.DataArray(
            data=lon, 
            dims=("ny_grid", "nx_grid",), 
            name="lon",
            attrs=dict(
                long_name = "Longitude",
                standard_name = "longitude",
                units = "degrees_east"
            )
        )
        xdat.append(lon2)
        lat2 = xr.DataArray(
            data=lat, 
            dims=("ny_grid", "nx_grid",), 
            name="lat",
            attrs=dict(
                long_name = "Latitude",
                standard_name = "latitude",
                units = "degrees_north"
            )
        )
        xdat.append(lat2)

        if air:
            for key, value in airVars.items():
                url = f'{base_url}/{value[1]}/{value[0]}/{value[1]}/'
                print(url)
                ds = xr.open_zarr(fsspec.get_mapper(url, anon=True))
                ds = ds.rename({'projection_y_coordinate': 'ny_grid'})
                ds = ds.rename({'projection_x_coordinate': 'nx_grid'})
                val = ds[value[0]][:, idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32')
                xdat.append(xr.DataArray(val, dims=("time", "ny_grid", "nx_grid",), name=key))

        if prc:
            for key, value in prcVars.items():
                url = f'{base_url}/{value[1]}/{value[0]}/{value[1]}/'
                print(url)
                ds = xr.open_zarr(fsspec.get_mapper(url, anon=True))
                ds = ds.rename({'projection_y_coordinate': 'ny_grid'})
                ds = ds.rename({'projection_x_coordinate': 'nx_grid'})
                val = ds[value[0]][:, idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32')
                xdat.append(xr.DataArray(val, dims=("time", "ny_grid", "nx_grid",), name=key))

        if rad:
            for key, value in radVars.items():
                url = f'{base_url}/{value[1]}/{value[0]}/{value[1]}/'
                print(url)
                ds = xr.open_zarr(fsspec.get_mapper(url, anon=True))
                ds = ds.rename({'projection_y_coordinate': 'ny_grid'})
                ds = ds.rename({'projection_x_coordinate': 'nx_grid'})
                val = ds[value[0]][:, idx_ymin:idx_ymax+1, idx_xmin:idx_xmax+1].astype('float32')
                #DSWRF are all junk values in this dataset
                #if np.sum(np.isnan(val)) > 0:
                #    val[:] = 0.0
                xdat.append(xr.DataArray(val, dims=("time", "ny_grid", "nx_grid",), name=key))
        
        time = ds['time'].values.astype('float32')
        time = (time+1)/24
        bdate = pd.to_datetime(start_date + timedelta(hours=cycle)).strftime('%Y %m %d %H').split(' ')
        bdate = [int(q) for q in bdate[:4]] + [0]
        xdat.append(
            xr.DataArray(
                time, dims=('time'), 
                name='time2',
                attrs=dict(
                    long_name="Time",
                    standard_name = "time",
                    base_date = bdate,
                    units = f"days since {start_date.year}-{start_date.month}-{start_date.day} {cycle:02d}:00 UTC"
                )))
        fout = xr.merge(xdat)
        fout=fout.rename_vars({'time2': 'time'})
        fout.to_netcdf(path / f'hrrr_{start_date.strftime("%Y%m%d")}{cycle:02d}.nc', 'w', unlimited_dims='time')

    def modified_latlon(self):
        latlon_h5_file = "HRRR_latlon.h5"
        ds = xr.open_dataset(latlon_h5_file)
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
