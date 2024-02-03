import os
import sys
from datetime import datetime,timedelta
import logging
import pathlib
import tempfile
import subprocess
import shutil
from typing import Union
from time import time

import numpy as np
import scipy as sp
#from numba import jit, prange
import netCDF4 as nc
from netCDF4 import Dataset
from matplotlib.transforms import Bbox
import seawater as sw
import xarray as xr

from pyschism.mesh.base import Nodes, Elements
from pyschism.mesh.vgrid import Vgrid

logger = logging.getLogger(__name__)

def convert_longitude(ds, bbox):
#https://stackoverflow.com/questions/53345442/about-changing-longitude-array-from-0-360-to-180-to-180-with-python-xarray
#Light_B's solution didn't generate the correct result
#Michael's solution works, but it takes significantly longer to write nc file (~30 mins compared with 5 mins)
#TODO: figure out why it takes much longer with the second method
    #lon_attr = ds.coords['lon'].attrs
    if bbox.xmin < 0:
        logger.info(f'Convert HYCOM longitude from [0, 360) to [-180, 180):')
        #ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
        ds['_lon_adjusted'] = xr.where(ds['lon'] > 180, ds['lon'] - 360, ds['lon'])
    elif bbox.xmin > 0:
        logger.info(f'Convert HYCOM longitude from [-180, 180) to [0, 360): ')
        #ds.coords['lon'] = (ds.coords['lon'] + 360) % 360 - 180
        ds['_lon_adjusted'] = xr.where(ds['lon'] < 0, ds['lon'] + 360, ds['lon'])

    t0 = time()
    ds = (
        ds.swap_dims({'lon': '_lon_adjusted'})
        .sel(**{'_lon_adjusted': sorted(ds._lon_adjusted)})
        .drop('lon')
    )
    ds = ds.rename({'_lon_adjusted': 'lon'})
    #ds = ds.sortby(ds.lon)
    #ds.coords['lon'].attrs = lon_attr
    logger.info(f'swap dims took {time()-t0} seconds!')

    #make sure it is clipped to the bbox
    ds = ds.sel(lon=slice(bbox.xmin - 0.5, bbox.xmax + 0.5))

    return ds

def get_database(date, Bbox=None):
    if date >= datetime(2018, 12, 4):
        database = f'GLBy0.08/expt_93.0'
    elif date >= datetime(2018, 1, 1) and date < datetime(2018, 12, 4):
        database = f'GLBv0.08/expt_93.0'
    elif date >= datetime(2017, 10, 1) and date < datetime(2018, 1, 1):
        database = f'GLBv0.08/expt_92.9'
    elif date >= datetime(2017, 6, 1) and date < datetime(2017, 10, 1):
        database = f'GLBv0.08/expt_57.7'
    elif date >= datetime(2017, 2, 1) and date < datetime(2017, 6, 1):
        database = f'GLBv0.08/expt_92.8'
    elif date >= datetime(2016, 5, 1) and date < datetime(2017, 2, 1):
        database = f'GLBv0.08/expt_57.2'
    elif date >= datetime(2016, 1, 1) and date < datetime(2016, 5, 1):
        database = f'GLBv0.08/expt_56.3'
    elif date >= datetime(1994, 1, 1) and date < datetime(2016, 1, 1):
        database = f'GLBv0.08/expt_53.X/data/{date.year}'
    else:
        raise ValueError(f'No data fro {date}!') 
    return database

def get_idxs(date, database, bbox, lonc=None, latc=None):

    if date >= datetime.utcnow():
        date2 = datetime.utcnow() - timedelta(days=1)
        baseurl = f'https://tds.hycom.org/thredds/dodsC/{database}/FMRC/runs/GLBy0.08_930_FMRC_RUN_{date2.strftime("%Y-%m-%dT12:00:00Z")}?depth[0:1:-1],lat[0:1:-1],lon[0:1:-1],time[0:1:-1]'
    else:
        baseurl=f'https://tds.hycom.org/thredds/dodsC/{database}?lat[0:1:-1],lon[0:1:-1],time[0:1:-1],depth[0:1:-1]'

    ds=Dataset(baseurl)
    time1=ds['time']
    times=nc.num2date(time1,units=time1.units,only_use_cftime_datetimes=False)
    
    lon=ds['lon'][:]
    lat=ds['lat'][:]
    dep=ds['depth'][:]
    
    #check if hycom's lon is the same range as schism's
    same = True
    if not (bbox.xmin >= lon.min() and bbox.xmax <= lon.max()):
        same = False
        if lon.min() >= 0:
            logger.info(f'Convert HYCOM longitude from [0, 360) to [-180, 180):')
            idxs = lon>=180
            lon[idxs] = lon[idxs]-360
        elif lon.min() <= 0:
            logger.info(f'Convert HYCOM longitude from [-180, 180) to [0, 360):')
            idxs = lon<=0
            lon[idxs] = lon[idxs]+360

    lat_idxs=np.where((lat>=bbox.ymin-0.5)&(lat<=bbox.ymax+0.5))[0]
    lon_idxs=np.where((lon>=bbox.xmin-0.5) & (lon<=bbox.xmax+0.5))[0]
    lon=lon[lon_idxs]
    lat=lat[lat_idxs]
    #logger.info(lon_idxs)
    #logger.info(lat_idxs)
    lon_idx1=lon_idxs[0].item()
    lon_idx2=lon_idxs[-1].item()
    #logger.info(f'lon_idx1 is {lon_idx1}, lon_idx2 is {lon_idx2}')
    lat_idx1=lat_idxs[0].item()
    lat_idx2=lat_idxs[-1].item()
    #logger.info(f'lat_idx1 is {lat_idx1}, lat_idx2 is {lat_idx2}')
    
    if lonc is None:
        lonc = lon.mean()
    #logger.info(f'lonc is {lonc}')
    if latc is None:
        latc = lat.mean()
    #logger.info(f'latc is {latc}')
    x2, y2=transform_ll_to_cpp(lon, lat, lonc, latc)

    idxs=np.where( date == times)[0]
    #check if time_idx is empty
    if len(idxs) == 0:
        #If there is missing data, use the data from the next days, the maximum searching days is 3. Otherwise, stop.
        for i in np.arange(0,3):
            date_before=(date + timedelta(days=int(i)+1)) #.astype(datetime)
            logger.info(f'Try replacing the missing data from {date_before}')
            idxs=np.where(date_before == times)[0]
            if len(idxs) == 0:
                continue
            else:
                break
    if len(idxs) ==0:
        logger.info(f'No date for date {date}')
        sys.exit()
    time_idx=idxs.item()  

    ds.close()

    return time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2, same

def transform_ll_to_cpp(lon, lat, lonc=-77.07, latc=24.0):
    longitude=lon/180*np.pi
    latitude=lat/180*np.pi
    radius=6378206.4
    loncc=lonc/180*np.pi
    latcc=latc/180*np.pi
    lon_new=[radius*(longitude[i]-loncc)*np.cos(latcc) for i in np.arange(len(longitude))]
    lat_new=[radius*latitude[i] for i in np.arange(len(latitude))]

    return np.array(lon_new), np.array(lat_new)

def interp_to_points_3d(dep, y2, x2, bxyz, val):
    idxs = np.where(abs(val) > 10000)
    val[idxs] = float('nan')

    if not np.all(x2[:-1] <= x2[1:]):
        logger.info('x2 is not in stricitly ascending order! Sorting x2 and val')
        idxs = np.argsort(x2)
        x2 = x2[idxs]
        val = val[:, :, idxs]

    val_fd = sp.interpolate.RegularGridInterpolator((dep,y2,x2),np.squeeze(val),'linear', bounds_error=False, fill_value = float('nan'))
    val_int = val_fd(bxyz)
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        val_int[idxs] = sp.interpolate.griddata(bxyz[~idxs,:], val_int[~idxs], bxyz[idxs,:],'nearest')
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        logger.info(f'There is still missing value for {val}')
        sys.exit()
    return val_int

def interp_to_points_2d(y2, x2, bxy, val):
    idxs = np.where(abs(val) > 10000)
    val[idxs] = float('nan')

    if not np.all(x2[:-1] <= x2[1:]):
        logger.info('x2 is not in stricitly ascending order! Sorting x2 and val')
        idxs = np.argsort(x2)
        x2 = x2[idxs]
        val = val[:, idxs]

    val_fd = sp.interpolate.RegularGridInterpolator((y2,x2),np.squeeze(val),'linear', bounds_error=False, fill_value = float('nan'))
    val_int = val_fd(bxy)
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        val_int[idxs] = sp.interpolate.griddata(bxy[~idxs,:], val_int[~idxs], bxy[idxs,:],'nearest')
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        logger.info(f'There is still missing value for {val}')
        sys.exit()
    return val_int

def ConvertTemp(salt, temp, dep):
    nz = temp.shape[0]
    ny = temp.shape[1]
    nx = temp.shape[2]
    pr = np.ones(temp.shape)
    pre = pr*dep[:,None, None]
    Pr = np.zeros(temp.shape)
    ptemp = sw.ptmp(salt, temp, pre, Pr)*1.00024
    return ptemp

class OpenBoundaryInventory:

    def __init__(self, hgrid, vgrid=None):
        self.hgrid = hgrid
        self.vgrid = Vgrid.default() if vgrid is None else vgrid

    def fetch_data(self, outdir: Union[str, os.PathLike], start_date, rnday, ocean_bnd_ids = [0], elev2D=True, TS=True, UV=True, restart=False, adjust2D=False, lats=None, msl_shifts=None): 
        outdir = pathlib.Path(outdir)

        self.start_date = start_date
        self.rnday=rnday
        self.timevector=np.arange(
            self.start_date,
            self.start_date + timedelta(days=self.rnday+1),
            timedelta(days=1)).astype(datetime)

        #Get open boundary 
        gdf=self.hgrid.boundaries.open.copy()
        opbd=[]
        #for boundary in gdf.itertuples():
        #    opbd.extend(list(boundary.indexes))
        for ibnd in ocean_bnd_ids:
            opbd.extend(list(gdf.iloc[ibnd].indexes))
        blon = self.hgrid.coords[opbd,0]
        blat = self.hgrid.coords[opbd,1]
        #logger.info(f'blon min {np.min(blon)}, max {np.max(blon)}')
        NOP = len(blon)

        #calculate zcor for 3D
        if TS or UV:
            vd=Vgrid.open(self.vgrid)
            sigma=vd.sigma

            #get bathymetry
            depth = self.hgrid.values

            #compute zcor
            zcor = depth[:,None]*sigma
            nvrt=zcor.shape[1]

            #zcor2=zcor[opbd,:]
            #idxs=np.where(zcor2 > 5000)
            #zcor2[idxs]=5000.0-1.0e-6

            #construct schism grid
            #x2i=np.tile(xi,[nvrt,1]).T
            #y2i=np.tile(yi,[nvrt,1]).T
            #bxyz=np.c_[zcor2.reshape(np.size(zcor2)),y2i.reshape(np.size(y2i)),x2i.reshape(np.size(x2i))]
            #logger.info('Computing SCHISM zcor is done!')

        #create netcdf
        ntimes=self.rnday+1
        nComp1=1
        nComp2=2
        one=1
        #ndt=np.zeros([ntimes])

        if elev2D and restart == False:
            #timeseries_el=np.zeros([ntimes,NOP,nComp1])
            #create netcdf 
            dst_elev = Dataset(outdir / 'elev2D.th.nc', 'w', format='NETCDF4')
            #dimensions
            dst_elev.createDimension('nOpenBndNodes', NOP)
            dst_elev.createDimension('one', one)
            dst_elev.createDimension('time', None)
            dst_elev.createDimension('nLevels', one)
            dst_elev.createDimension('nComponents', nComp1)

            #variables
            dst_elev.createVariable('time_step', 'f', ('one',))
            dst_elev['time_step'][:] = 86400

            dst_elev.createVariable('time', 'f', ('time',))
            #dst_elev['time'][:] = ndt

            dst_elev.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
            #dst_elev['time_series'][:,:,:,:] = timeseries_el
        elif elev2D and restart:
            dst_elev = Dataset(outdir / 'elev2D.th.nc', 'a', format='NETCDF4')
            time_idx_restart = dst_elev['time'][:].shape[0]

        if TS and restart == False:
            #timeseries_s=np.zeros([ntimes,NOP,nvrt,nComp1])
            dst_salt = Dataset(outdir / 'SAL_3D.th.nc', 'w', format='NETCDF4')
            #dimensions
            dst_salt.createDimension('nOpenBndNodes', NOP)
            dst_salt.createDimension('one', one)
            dst_salt.createDimension('time', None)
            dst_salt.createDimension('nLevels', nvrt)
            dst_salt.createDimension('nComponents', nComp1)
            #variables
            dst_salt.createVariable('time_step', 'f', ('one',))
            dst_salt['time_step'][:] = 86400

            dst_salt.createVariable('time', 'f', ('time',))
            #dst_salt['time'][:] = ndt

            dst_salt.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))

            #temp
            #timeseries_t=np.zeros([ntimes,NOP,nvrt,nComp1])

            dst_temp = Dataset(outdir / 'TEM_3D.th.nc', 'w', format='NETCDF4')
            #dimensions
            dst_temp.createDimension('nOpenBndNodes', NOP)
            dst_temp.createDimension('one', one)
            dst_temp.createDimension('time', None)
            dst_temp.createDimension('nLevels', nvrt)
            dst_temp.createDimension('nComponents', nComp1)
            #variables
            dst_temp.createVariable('time_step', 'f', ('one',))
            dst_temp['time_step'][:] = 86400

            dst_temp.createVariable('time', 'f', ('time',))
            #dst_temp['time'][:] = ndt

            dst_temp.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
            #dst_temp['time_series'][:,:,:,:] = timeseries_t
        elif TS and restart:
            dst_salt = Dataset(outdir / 'SAL_3D.th.nc', 'a', format='NETCDF4')
            dst_temp = Dataset(outdir / 'TEM_3D.th.nc', 'a', format='NETCDF4')
            time_idx_restart = dst_salt['time'][:].shape[0]

        if UV and restart == False:
            #timeseries_uv=np.zeros([ntimes,NOP,nvrt,nComp2])
            dst_uv = Dataset(outdir / 'uv3D.th.nc', 'w', format='NETCDF4')
            #dimensions
            dst_uv.createDimension('nOpenBndNodes', NOP)
            dst_uv.createDimension('one', one)
            dst_uv.createDimension('time', None)
            dst_uv.createDimension('nLevels', nvrt)
            dst_uv.createDimension('nComponents', nComp2)
            #variables
            dst_uv.createVariable('time_step', 'f', ('one',))
            dst_uv['time_step'][:] = 86400

            dst_uv.createVariable('time', 'f', ('time',))
            #dst_uv['time'][:] = ndt
            
            dst_uv.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
            #dst_uv['time_series'][:,:,:,:] = timeseries_uv

        elif UV and restart:
            dst_uv = Dataset(outdir / 'uv3D.th.nc', 'a', format='NETCDF4')
            time_idx_restart = dst_uv['time'][:].shape[0]

        logger.info('**** Accessing GOFS data*****')
        t0=time()

        if restart == False: 
            timevector = self.timevector
            it0 = 0
        elif restart:
            #restart from one day earlier
            timevector = self.timevector[time_idx_restart-1:]
            it0 = time_idx_restart-1

        for it1, date in enumerate(timevector):
            
            it = it0 + it1
        
            database=get_database(date)
            logger.info(f'Fetching data for {date} from database {database}')

            #loop over each open boundary
            ind1 = 0
            ind2 = 0
            #for boundary in gdf.itertuples():
            for ibnd in ocean_bnd_ids:
                #opbd = list(boundary.indexes)
                opbd = list(gdf.iloc[ibnd].indexes)
                ind1 = ind2
                ind2 = ind1 + len(opbd)
                #logger.info(f'ind1 = {ind1}, ind2 = {ind2}')
                blon = self.hgrid.coords[opbd,0]
                blat = self.hgrid.coords[opbd,1]
                blonc = blon.mean()
                blatc = blat.mean()
                #logger.info(f'blonc = {blon.mean()}, blatc = {blat.mean()}')
                xi,yi = transform_ll_to_cpp(blon, blat, blonc, blatc)
                bxy = np.c_[yi, xi]

                if TS or UV:
                    zcor2=zcor[opbd,:]
                    idxs=np.where(zcor2 > 5000)
                    zcor2[idxs]=5000.0-1.0e-6

                    #construct schism grid
                    x2i=np.tile(xi,[nvrt,1]).T
                    y2i=np.tile(yi,[nvrt,1]).T
                    bxyz=np.c_[zcor2.reshape(np.size(zcor2)),y2i.reshape(np.size(y2i)),x2i.reshape(np.size(x2i))]

                xmin, xmax = np.min(blon), np.max(blon)
                ymin, ymax = np.min(blat), np.max(blat)
                bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)

                time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2, _ = get_idxs(date, database, bbox, lonc=blonc, latc=blatc)

                if date >= datetime.utcnow():
                    date2 = datetime.utcnow() - timedelta(days=1)
                    url = f'https://tds.hycom.org/thredds/dodsC/{database}/FMRC/runs/GLBy0.08_930_FMRC_RUN_' + \
                        f'{date2.strftime("%Y-%m-%dT12:00:00Z")}?depth[0:1:-1],lat[{lat_idx1}:1:{lat_idx2}],' + \
                        f'lon[{lon_idx1}:1:{lon_idx2}],time[{time_idx}],' + \
                        f'surf_el[{time_idx}][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                        f'water_temp[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                        f'salinity[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                        f'water_u[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                        f'water_v[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}]'
                        
                else:
                    url=f'https://tds.hycom.org/thredds/dodsC/{database}?lat[{lat_idx1}:1:{lat_idx2}],' + \
                        f'lon[{lon_idx1}:1:{lon_idx2}],depth[0:1:-1],time[{time_idx}],' + \
                        f'surf_el[{time_idx}][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                        f'water_temp[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                        f'salinity[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                        f'water_u[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                        f'water_v[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}]'
                #logger.info(url)
                 
                ds=Dataset(url)
                dep=ds['depth'][:]

                logger.info(f'****Interpolation starts for boundary {ibnd}****')

                #ndt[it]=it*24*3600.

                if elev2D:
                    #ssh
                    ssh=np.squeeze(ds['surf_el'][:,:])

                    ssh_int = interp_to_points_2d(y2, x2, bxy, ssh)
                    dst_elev['time'][it] = it*24*3600.
                    if adjust2D:
                        elev_adjust = np.interp(blat, lats, msl_shifts)
                        dst_elev['time_series'][it,ind1:ind2,0,0] = ssh_int + elev_adjust
                    else:
                        dst_elev['time_series'][it,ind1:ind2,0,0] = ssh_int 

                if TS:
                    #salt
                    salt = np.squeeze(ds['salinity'][:,:,:])

                    salt_int = interp_to_points_3d(dep, y2, x2, bxyz, salt)
                    salt_int = salt_int.reshape(zcor2.shape)
                    #timeseries_s[it,:,:,0]=salt_int
                    dst_salt['time'][it] = it*24*3600.
                    dst_salt['time_series'][it,ind1:ind2,:,0] = salt_int

                    #temp
                    temp = np.squeeze(ds['water_temp'][:,:,:])

                    #Convert temp to potential temp
                    ptemp = ConvertTemp(salt, temp, dep)

                    temp_int = interp_to_points_3d(dep, y2, x2, bxyz, ptemp)
                    temp_int = temp_int.reshape(zcor2.shape)
                    #timeseries_t[it,:,:,0]=temp_int
                    dst_temp['time'][it] = it*24*3600.
                    dst_temp['time_series'][it,ind1:ind2,:,0] = temp_int

                if UV:
                    uvel=np.squeeze(ds['water_u'][:,:,:])
                    vvel=np.squeeze(ds['water_v'][:,:,:])

                    dst_uv['time'][it] = it*24*3600.
                    #uvel
                    uvel_int = interp_to_points_3d(dep, y2, x2, bxyz, uvel)
                    uvel_int = uvel_int.reshape(zcor2.shape)
                    dst_uv['time_series'][it,ind1:ind2,:,0] = uvel_int

                    #vvel
                    vvel_int = interp_to_points_3d(dep, y2, x2, bxyz, vvel)
                    vvel_int = vvel_int.reshape(zcor2.shape)
                    dst_uv['time_series'][it,ind1:ind2,:,1] = vvel_int
                    #timeseries_uv[it,:,:,1]=vvel_int

                ds.close()
        logger.info(f'Writing *th.nc takes {time()-t0} seconds')

class Nudge:

    def __init__(self, hgrid=None, ocean_bnd_ids=None):

        if hgrid is None:
            raise ValueError('No hgrid information!') 
        else:
            self.hgrid = hgrid

        
        if ocean_bnd_ids is None:
            raise ValueError('Please specify indexes for ocean boundaries!') 
        else:
            self.ocean_bnd_ids = ocean_bnd_ids


    def gen_nudge(self, outdir: Union[str, os.PathLike], rlmax = 1.5, rnu_day=0.25):
        """
        set up nudge zone within rlmax distance from the ocean boundary;
        modify the nudging zone width rlmax.
        rlmax can be a uniform value, e.g., rlmax = 1.5 (degree if hgrid is lon/lat)
        """
          
        outdir = pathlib.Path(outdir)

        rnu_max = 1.0 / rnu_day / 86400.0

        #get nudge zone
        lon = self.hgrid.coords[:,0]
        lat = self.hgrid.coords[:,1]
        gdf = self.hgrid.boundaries.open.copy()
        elnode = self.hgrid.elements.array
        NE, NP = elnode.shape[0],len(lon)
        nudge_coeff = np.zeros(NP, dtype=float)

        global_idxs = {}

        t0 = time()
        nudge_coeff = np.zeros(NP, dtype=float)
        for i in self.ocean_bnd_ids:
            print(f'boundary {i}')
            bnd_idxs = gdf.iloc[i].indexes

            dis = abs((lon + 1j*lat)[:, None] - (lon[bnd_idxs] + 1j*lat[bnd_idxs])[None, :]).min(axis=1)
            out = (1-dis/rlmax)*rnu_max
            out[out<0] = 0
            out[out>rnu_max] = rnu_max
            fp = out>0
            nudge_coeff[fp] = np.maximum(out[fp], nudge_coeff[fp])

            idxs_nudge=np.zeros(NP, dtype=int)
            idxs=np.where(out > 0)[0]
            idxs_nudge[idxs]=1
  
            #expand nudging marker to neighbor nodes
            i34 = self.hgrid.elements.i34
            fp = i34==3
            idxs=np.where(np.max(out[elnode[fp, 0:3]], axis=1) > 0)[0]
            idxs_nudge[elnode[fp,0:3][idxs,:]]=1
            idxs=np.where(np.max(out[elnode[~fp, :]], axis=1) > 0)[0]
            idxs_nudge[elnode[~fp,:][idxs,:]]=1

            idxs=np.where(idxs_nudge == 1)[0]
            global_idxs[i] = idxs


        #logger.info(f'len of nudge idxs is {len(idxs)}')
        logger.info(f'It took {time() -t0} sencods to calcuate nudge coefficient')

        nudge = [f"rlmax={rlmax}, rnu_day={rnu_day}"]
        nudge.extend("\n")
        nudge.append(f"{NE} {NP}")
        nudge.extend("\n")
        hgrid = self.hgrid.to_dict()
        nodes = hgrid['nodes']
        elements = hgrid['elements']
        for idn, (coords, values) in nodes.items():
            line = [f"{idn}"]
            line.extend([f"{x:<.7e}" for x in coords])
            line.extend([f"{nudge_coeff[int(idn)-1]:<.7e}"])
            line.extend("\n")
            nudge.append(" ".join(line))

        for id, element in elements.items():
            line = [f"{id}"]
            line.append(f"{len(element)}")
            line.extend([f"{e}" for e in element])
            line.extend("\n")
            nudge.append(" ".join(line))

        with open(outdir / 'TEM_nudge.gr3','w+') as fid:
            fid.writelines(nudge)

        shutil.copy2(outdir / 'TEM_nudge.gr3', outdir / 'SAL_nudge.gr3')

        return global_idxs

    def fetch_data(self, outdir: Union[str, os.PathLike], vgrid, start_date, rnday, restart = False, rlmax = None, rnu_day=None):
        """
        fetch data from the database and generate nudge file
        see gen_nudge for the meaning of rlmax, rnu_day
        """

        outdir = pathlib.Path(outdir)

        self.start_date = start_date
        self.rnday=rnday
        self.timevector=np.arange(
            self.start_date,
            self.start_date + timedelta(days=self.rnday+1),
            timedelta(days=1)).astype(datetime)

        vd=Vgrid.open(vgrid)
        sigma=vd.sigma

        #define nudge zone and strength 
        rlmax = 1.5 if rlmax is None else rlmax
        rnu_day = 0.25 if rnu_day is None else rnu_day
        logger.info(f'Max relax distance is {rlmax} degree, max relax strengh is {rnu_day} days.')
        #Get the index for nudge
        global_idxs = self.gen_nudge(outdir, rlmax = rlmax, rnu_day=rnu_day)

        #get bathymetry
        depth = self.hgrid.values

        #compute zcor
        zcor = depth[:,None]*sigma
        nvrt=zcor.shape[1]

        #allocate output variables
        include = global_idxs[0]
        nbnd = len(self.ocean_bnd_ids)
        if nbnd > 1:
            include = np.concatenate((global_idxs[0], global_idxs[1],))
            if nbnd > 2:
                for i in self.ocean_bnd_ids[2:]:
                    include = np.concatenate((include, global_idxs[i]))
        else:
            include = global_idxs[0]

        nNode = include.shape[0]
        one = 1
        ntimes = self.rnday+1

        #timeseries_s=np.zeros([ntimes,nNode,nvrt,one])
        #timeseries_t=np.zeros([ntimes,nNode,nvrt,one])
        #ndt=np.zeros([ntimes])
        if restart:
            dst_temp = Dataset(outdir / 'TEM_nu.nc', 'a', format='NETCDF4')
            dst_salt = Dataset(outdir / 'SAL_nu.nc', 'a', format='NETCDF4')
            time_idx_restart = dst_temp['time'][:].shape[0]
        else:
            dst_temp = Dataset(outdir / 'TEM_nu.nc', 'w', format='NETCDF4')
            #dimensions
            dst_temp.createDimension('node', nNode)
            dst_temp.createDimension('nLevels', nvrt)
            dst_temp.createDimension('one', one)
            dst_temp.createDimension('time', None)
            #variables
            dst_temp.createVariable('time', 'f', ('time',))
            #dst_temp['time'][:] = ndt

            dst_temp.createVariable('map_to_global_node', 'i4', ('node',))
            dst_temp['map_to_global_node'][:] = include+1

            dst_temp.createVariable('tracer_concentration', 'f', ('time', 'node', 'nLevels', 'one'))
            #dst_temp['tracer_concentration'][:,:,:,:] = timeseries_t

            #salinity
            dst_salt = Dataset(outdir / 'SAL_nu.nc', 'w', format='NETCDF4')
            #dimensions
            dst_salt.createDimension('node', nNode)
            dst_salt.createDimension('nLevels', nvrt)
            dst_salt.createDimension('one', one)
            dst_salt.createDimension('time', None)
            #variables
            dst_salt.createVariable('time', 'f', ('time',))
            #dst_salt['time'][:] = ndt

            dst_salt.createVariable('map_to_global_node', 'i4', ('node',))
            dst_salt['map_to_global_node'][:] = include+1

            dst_salt.createVariable('tracer_concentration', 'f', ('time', 'node', 'nLevels', 'one'))
            #dst_salt['tracer_concentration'][:,:,:,:] = timeseries_s

        logger.info('**** Accessing GOFS data*****')
        if restart:
            #restart from one day earlier to make sure all files consistant
            timevector = self.timevector[time_idx_restart-1:]
            it0 = time_idx_restart-1
        else:
            timevector = self.timevector
            it0 = 0

        t0=time()
        for it1, date in enumerate(timevector):

            it = it0 + it1

            database=get_database(date)
            logger.info(f'Fetching data for {date} from database {database}')

            ind1 = 0
            ind2 = 0
            for ibnd in self.ocean_bnd_ids:
                include = global_idxs[ibnd]

                ind1 = ind2
                ind2 = ind1 + include.shape[0]
                #Get open nudge array 
                nlon = self.hgrid.coords[include, 0]
                nlat = self.hgrid.coords[include, 1]
                nlonc = nlon.mean()
                nlatc = nlat.mean()
                xi,yi = transform_ll_to_cpp(nlon, nlat, nlonc, nlatc)
                bxy = np.c_[yi, xi]

                zcor2=zcor[include,:]
                idxs=np.where(zcor2 > 5000)
                zcor2[idxs]=5000.0-1.0e-6

                #construct schism grid
                x2i=np.tile(xi,[nvrt,1]).T
                y2i=np.tile(yi,[nvrt,1]).T
                bxyz=np.c_[zcor2.reshape(np.size(zcor2)),y2i.reshape(np.size(y2i)),x2i.reshape(np.size(x2i))]
                logger.info('Computing SCHISM zcor is done!')

                xmin, xmax = np.min(nlon), np.max(nlon)
                ymin, ymax = np.min(nlat), np.max(nlat)
                bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)
                logger.info(f'bbox for nudge is {bbox}')

                time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2, _ = get_idxs(date, database, bbox, lonc=nlonc, latc=nlatc)

                if date >= datetime.utcnow():
                    date2 = datetime.utcnow() - timedelta(days=1)
                    url = f'https://tds.hycom.org/thredds/dodsC/{database}/FMRC/runs/GLBy0.08_930_FMRC_RUN_' + \
                        f'{date2.strftime("%Y-%m-%dT12:00:00Z")}?depth[0:1:-1],lat[{lat_idx1}:1:{lat_idx2}],' + \
                        f'lon[{lon_idx1}:1:{lon_idx2}],time[{time_idx}],' + \
                        f'water_temp[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                        f'salinity[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}]'

                else:
                    url=f'https://tds.hycom.org/thredds/dodsC/{database}?lat[{lat_idx1}:1:{lat_idx2}],' + \
                        f'lon[{lon_idx1}:1:{lon_idx2}],depth[0:1:-1],time[{time_idx}],' + \
                        f'water_temp[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                        f'salinity[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}]'

                ds=Dataset(url)
                salt=np.squeeze(ds['salinity'][:,:,:])
                temp=np.squeeze(ds['water_temp'][:,:,:])

                #Convert temp to potential temp
                dep=ds['depth'][:]
                ptemp = ConvertTemp(salt, temp, dep)

                logger.info('****Interpolation starts****')

                #ndt[it]=it
                #salt
                dst_salt['time'][it] = it
                salt_int = interp_to_points_3d(dep, y2, x2, bxyz, salt)
                salt_int = salt_int.reshape(zcor2.shape)
                #timeseries_s[it,:,:,0]=salt_int
                dst_salt['tracer_concentration'][it,ind1:ind2,:,0] = salt_int

                #temp
                dst_temp['time'][it] = it
                temp_int = interp_to_points_3d(dep, y2, x2, bxyz, ptemp)
                temp_int = temp_int.reshape(zcor2.shape)
                #timeseries_t[it,:,:,0]=temp_int
                dst_temp['tracer_concentration'][it,ind1:ind2,:,0] = temp_int

                ds.close()
 
        #dst_temp.close()
        #dst_salt.close()
        
        logger.info(f'Writing *_nu.nc takes {time()-t0} seconds')

class DownloadHycom:

    def __init__(self, hgrid=None, bbox=None):

        if hgrid is None and bbox is None:
            raise ValueError('Either hgrid or bbox must be provided!') 

        if hgrid is not None:        
            xmin, xmax = hgrid.coords[:, 0].min(), hgrid.coords[:, 0].max()
            ymin, ymax = hgrid.coords[:, 1].min(), hgrid.coords[:, 1].max()
            self.bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)
        elif bbox is not None:
            self.bbox = bbox

    def fetch_data(self, start_date, rnday=1, fmt='schism', bnd=False, nudge=False, sub_sample=1, outdir=None):
        '''
        start_date: datetime.datetime
        rnday: integer
        fmt: 'schism' - for Fortran code; 'hycom' - raw netCDF from HYCOM
        bnd: file names are SSH_*.nc, TS_*.nc, UV_*.nc used in gen_hot_3Dth_from_hycom.f90
        nudge: file name is TS_*.nc used in gen_nudge_from_hycom.f90
        outdir: directory for output files
        '''
        if rnday == 1:
            timevector = [start_date]
        else:
            timevector = np.arange(
                start_date, start_date + timedelta(days=rnday+1), timedelta(days=1)
            ).astype(datetime)

        for i, date in enumerate(timevector):
            database=get_database(date)
            logger.info(f'Fetching data for {date} from database {database}')

            time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2, isLonSame = get_idxs(date, database, self.bbox)

            url_ssh = f'https://tds.hycom.org/thredds/dodsC/{database}?lat[{lat_idx1}:{sub_sample}:{lat_idx2}],' + \
                f'lon[{lon_idx1}:{sub_sample}:{lon_idx2}],depth[0:1:-1],time[{time_idx}],' + \
                f'surf_el[{time_idx}][{lat_idx1}:{sub_sample}:{lat_idx2}][{lon_idx1}:{sub_sample}:{lon_idx2}],' + \
                f'water_u[{time_idx}][0:1:39][{lat_idx1}:{sub_sample}:{lat_idx2}][{lon_idx1}:{sub_sample}:{lon_idx2}],' + \
                f'water_v[{time_idx}][0:1:39][{lat_idx1}:{sub_sample}:{lat_idx2}][{lon_idx1}:{sub_sample}:{lon_idx2}],' + \
                f'water_temp[{time_idx}][0:1:39][{lat_idx1}:{sub_sample}:{lat_idx2}][{lon_idx1}:{sub_sample}:{lon_idx2}],' + \
                f'salinity[{time_idx}][0:1:39][{lat_idx1}:{sub_sample}:{lat_idx2}][{lon_idx1}:{sub_sample}:{lon_idx2}]'

            foutname = f'hycom_{date.strftime("%Y%m%d")}.nc'
            if fmt == 'schism':
                #foutname = f'TS_{i+1}.nc'
                logger.info(f'filename is {foutname}')
                ds = xr.open_dataset(url_ssh)

                #convert in-situ temperature to potential temperature
                temp = ds.water_temp.values
                salt = ds.salinity.values
                dep = ds.depth.values

                ptemp = ConvertTemp(salt, temp, dep)
                #drop water_temp variable and add new temperature variable
                ds = ds.drop('water_temp')
                ds['temperature']=(['time','depth','lat','lon'], ptemp)
                ds.temperature.attrs = {
                    'long_name': 'Sea water potential temperature',
                    'standard_name': 'sea_water_potential_temperature',
                    'units': 'degC'
                }

                if not isLonSame:
                    logger.info('Lon is not the same!')
                    ds = convert_longitude(ds, self.bbox)

                ds = ds.rename_dims({'lon':'xlon'})
                ds = ds.rename_dims({'lat':'ylat'})
                ds = ds.rename_vars({'lat':'ylat'})
                ds = ds.rename_vars({'lon':'xlon'})

                t0 =  time()
                logger.info(f'Start writing nc file!')
                ds.to_netcdf(foutname, 'w', unlimited_dims='time', encoding={'temperature':{'dtype': 'h', '_FillValue': -30000.,'scale_factor': 0.001, 'add_offset': 20., 'missing_value': -30000.}})
                ds.close()
                logger.info(f'It took {time()-t0} seconds to write nc file!')

                #link output file to Fortran required files
                dir_path = os.path.abspath(outdir)
                logger.info(f'current dir is {dir_path}')
                src = f'{dir_path}/{foutname}'
                if bnd:
                    names = ['SSH', 'TS', 'UV']
                    for name in names:
                        dst = f'{dir_path}/{name}_{i+1}.nc'
                        os.symlink(src, dst)
                elif nudge:
                        dst = f'{dir_path}/TS_{i+1}.nc'
                        os.symlink(src, dst)

            elif fmt == 'hycom':
                #url=f'https://tds.hycom.org/thredds/dodsC/{database}?lat[{lat_idx1}:1:{lat_idx2}],' + \
                #    f'lon[{lon_idx1}:1:{lon_idx2}],depth[0:1:-1],time[{time_idx}],' + \
                #    f'surf_el[{time_idx}][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                #    f'water_temp[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                #    f'salinity[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                #    f'water_u[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
                #    f'water_v[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}]'

                #foutname = f'hycom_{date.strftime("%Y%m%d")}.nc' 
                ds = xr.open_dataset(url)
                ds.to_netcdf(foutname, 'w')

                ds.close()
