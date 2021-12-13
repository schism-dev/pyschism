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
from numba import jit, prange
import netCDF4 as nc
from netCDF4 import Dataset
from matplotlib.transforms import Bbox
import seawater as sw

from pyschism.mesh.base import Nodes, Elements
from pyschism.mesh.vgrid import Vgrid
from pyschism.forcing.hycom.hycom2schism import Nudge

logger = logging.getLogger(__name__)

def get_idxs(date, ds, bbox):

    time1=ds['time']
    times=nc.num2date(time1,units=time1.units,only_use_cftime_datetimes=False)
    
    lon=ds['lon'][:]
    lat=ds['lat'][:]
    lat_idxs=np.where((lat>=bbox.ymin-2.0)&(lat<=bbox.ymax+2.0))[0]
    lon_idxs=np.where((lon>=bbox.xmin-2.0) & (lon<=bbox.xmax+2.0))[0]
    lon=lon[lon_idxs]
    lat=lat[lat_idxs]
    #print(lon_idxs)
    #print(lat_idxs)
    lon_idx1=lon_idxs[0].item()
    lon_idx2=lon_idxs[-1].item()
    #print(f'lon_idx1 is {lon_idx1}, lon_idx2 is {lon_idx2}')
    lat_idx1=lat_idxs[0].item()
    lat_idx2=lat_idxs[-1].item()
    #print(f'lat_idx1 is {lat_idx1}, lat_idx2 is {lat_idx2}')
    
    for ilon in np.arange(len(lon)):
        if lon[ilon] > 180:
            lon[ilon] = lon[ilon]-360.
    #lonc=(np.max(lon)+np.min(lon))/2.0
    #print(f'lonc is {lonc}')
    #latc=(np.max(lat)+np.min(lat))/2.0
    #print(f'latc is {latc}')
    x2, y2=transform_ll_to_cpp(lon, lat)

    idxs=np.where( date == times)[0]
    #check if time_idx is empty
    if len(idxs) == 0:
        #If there is missing data, use the data from the next days, the maximum searching days is 3. Otherwise, stop.
        for i in np.arange(0,3):
            date_before=(date + timedelta(days=int(i)+1)) #.astype(datetime)
            print(f'Try replacing the missing data from {date_before}')
            idxs=np.where(date_before == times)[0]
            if len(idxs) == 0:
                continue
            else:
                break
    if len(idxs) ==0:
        print(f'No date for date {date}')
        sys.exit()
    time_idx=idxs.item()  
    print(f'time_idx is {time_idx}')

    return time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2

def transform_ll_to_cpp(lon, lat, lonc=-77.07, latc=24.0):
    #lonc=(np.max(lon)+np.min(lon))/2.0
    #print(f'lonc is {lonc}')
    #latc=(np.max(lat)+np.min(lat))/2.0
    #print(f'latc is {latc}')
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

    val_fd = sp.interpolate.RegularGridInterpolator((dep,y2,x2),np.squeeze(val),'linear', bounds_error=False, fill_value = float('nan'))
    val_int = val_fd(bxyz)
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        val_int[idxs] = sp.interpolate.griddata(bxyz[~idxs,:], val_int[~idxs], bxyz[idxs,:],'nearest')
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        print(f'There is still missing value for {val}')
        sys.exit()
    return val_int

def interp_to_points_2d(y2, x2, bxy, val):
    idxs = np.where(abs(val) > 10000)
    val[idxs] = float('nan')

    val_fd = sp.interpolate.RegularGridInterpolator((y2,x2),np.squeeze(val),'linear', bounds_error=False, fill_value = float('nan'))
    val_int = val_fd(bxy)
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        val_int[idxs] = sp.interpolate.griddata(bxy[~idxs,:], val_int[~idxs], bxy[idxs,:],'nearest')
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        print(f'There is still missing value for {val}')
        sys.exit()
    return val_int

class OpenBoundaryInventory:

    def __init__(self, hgrid, vgrid=None):
        self.hgrid = hgrid
        self.vgrid = Vgrid.default() if vgrid is None else vgrid

    def fetch_data(self, outdir: Union[str, os.PathLike], start_date, rnday, elev2D=True, TS=True, UV=True, adjust2D=False, lats=None, msl_shifts=None): 
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
        for boundary in gdf.itertuples():
            opbd.extend(list(boundary.indexes))
        blon = self.hgrid.coords[opbd,0]
        blat = self.hgrid.coords[opbd,1]
        #print(f'blon min {np.min(blon)}, max {np.max(blon)}')
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
            #print('Computing SCHISM zcor is done!')

        #create netcdf
        ntimes=self.rnday+1
        nComp1=1
        nComp2=2
        one=1
        ndt=np.zeros([ntimes])

        if elev2D:
            timeseries_el=np.zeros([ntimes,NOP,nComp1])
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
            dst_elev['time'][:] = ndt

            dst_elev.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
            dst_elev['time_series'][:,:,:,:] = timeseries_el

        if TS:
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
            dst_salt['time'][:] = ndt

            dst_salt.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))

            #temp
            timeseries_t=np.zeros([ntimes,NOP,nvrt,nComp1])

            dst_temp =  Dataset(outdir / 'TEM_3D.th.nc', 'w', format='NETCDF4')
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
            dst_temp['time'][:] = ndt

            dst_temp.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
            dst_temp['time_series'][:,:,:,:] = timeseries_t

        if UV:
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

        print('**** Accessing RTOFS data*****')
        baseurl = f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
        for it, date in enumerate(self.timevector):
            t0=time()
            print(f'Fetching data for {date}')

            #loop over each open boundary
            ind1 = 0
            ind2 = 0
            for boundary in gdf.itertuples():

                opbd = list(boundary.indexes)
                ind1 = ind2
                ind2 = ind1 + len(opbd)
                #print(f'ind1 = {ind1}, ind2 = {ind2}')
                blon = self.hgrid.coords[opbd,0]
                blat = self.hgrid.coords[opbd,1]
                xi,yi = transform_ll_to_cpp(blon, blat)
                bxy = np.c_[yi, xi]

                if TS or UV:
                    zcor2=zcor[opbd,:]
                    idxs=np.where(zcor2 > 5500)
                    zcor2[idxs]=5500.0-1.0e-6
                    print(f'zcor2[200,:] is {zcor2[200,:]}')

                    #construct schism grid
                    x2i=np.tile(xi,[nvrt,1]).T
                    y2i=np.tile(yi,[nvrt,1]).T
                    bxyz=np.c_[zcor2.reshape(np.size(zcor2)),y2i.reshape(np.size(y2i)),x2i.reshape(np.size(x2i))]

                xmin, xmax = np.min(blon), np.max(blon)
                ymin, ymax = np.min(blat), np.max(blat)

                # convert hgrid lon [180, 180) to [0, 360)
                xmin = xmin + 360. if xmin < 0 else xmin
                xmax = xmax + 360. if xmax < 0 else xmax
                bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)
                print(f'xmin is {xmin}, xmax is {xmax}')

                ssh_url = f'{baseurl}{date.strftime("%Y%m%d")}/rtofs_glo_2ds_forecast_3hrly_diag'
                salt_url = f'{baseurl}{date.strftime("%Y%m%d")}/rtofs_glo_3dz_forecast_daily_salt'
                temp_url = f'{baseurl}{date.strftime("%Y%m%d")}/rtofs_glo_3dz_forecast_daily_temp'
                uvel_url = f'{baseurl}{date.strftime("%Y%m%d")}/rtofs_glo_3dz_forecast_daily_uvel'
                vvel_url = f'{baseurl}{date.strftime("%Y%m%d")}/rtofs_glo_3dz_forecast_daily_vvel'

                print('****Interpolation starts****')

                #ndt[it]=it*24*3600.

                if elev2D:
                    #ssh
                    ds=Dataset(ssh_url)
                    print(f'ssh_url is {ssh_url}')
                    time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2 = get_idxs(date, ds, bbox)
                    ssh=np.squeeze(ds['ssh'][time_idx,0,lat_idx1:lat_idx2+1,lon_idx1:lon_idx2+1])

                    ssh_int = interp_to_points_2d(y2, x2, bxy, ssh)
                    dst_elev['time'][:] = it*24*3600.
                    if adjust2D:
                        elev_adjust = np.interp(blat, lats, msl_shifts)
                        dst_elev['time_series'][it,ind1:ind2,0,0] = ssh_int + elev_adjust
                    else:
                        dst_elev['time_series'][it,ind1:ind2,0,0] = ssh_int 
                    ds.close()

                if TS:
                    #salt
                    ds = Dataset(salt_url)
                    print(f'salt_url is {salt_url}')
                    time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2 = get_idxs(date, ds, bbox)
                    dep = ds['lev'][:]
                    salt = np.squeeze(ds['salinity'][time_idx,:,lat_idx1:lat_idx2+1,lon_idx1:lon_idx2+1])

                    salt_int = interp_to_points_3d(dep, y2, x2, bxyz, salt)
                    salt_int = salt_int.reshape(zcor2.shape)
                    #timeseries_s[it,:,:,0]=salt_int
                    dst_salt['time'][:] = it*24*3600.
                    dst_salt['time_series'][it,ind1:ind2,:,0] = salt_int

                    ds.close()

                    #temp
                    ds = Dataset(temp_url)
                    print(f'temp_url is {temp_url}')
                    dep = ds['lev'][:]
                    temp = np.squeeze(ds['temperature'][time_idx,:,lat_idx1:lat_idx2+1,lon_idx1:lon_idx2+1])

                    temp_int = interp_to_points_3d(dep, y2, x2, bxyz, temp)
                    temp_int = temp_int.reshape(zcor2.shape)
                    #timeseries_t[it,:,:,0]=temp_int
                    dst_temp['time'][:] = it*24*3600.
                    dst_temp['time_series'][it,ind1:ind2,:,0] = temp_int

                    ds.close()
                if UV:
                    ds = Dataset(uvel_url)
                    print(f'uvel_url is {uvel_url}')
                    time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2 = get_idxs(date, ds, bbox)
                    dep = ds['lev'][:]
                    uvel = np.squeeze(ds['u'][time_idx,:,lat_idx1:lat_idx2+1,lon_idx1:lon_idx2+1])

                    dst_uv['time'][:] = it*24*3600.
                    #uvel
                    uvel_int = interp_to_points_3d(dep, y2, x2, bxyz, uvel)
                    uvel_int = uvel_int.reshape(zcor2.shape)
                    dst_uv['time_series'][it,ind1:ind2,:,0] = uvel_int
                    ds.close()

                    #vvel
                    ds = Dataset(vvel_url)
                    print(f'vvel_url is {vvel_url}')
                    dep = ds['lev'][:]
                    vvel = np.squeeze(ds['v'][time_idx,:,lat_idx1:lat_idx2+1,lon_idx1:lon_idx2+1])

                    vvel_int = interp_to_points_3d(dep, y2, x2, bxyz, vvel)
                    vvel_int = vvel_int.reshape(zcor2.shape)
                    dst_uv['time_series'][it,ind1:ind2,:,1] = vvel_int
                    #timeseries_uv[it,:,:,1]=vvel_int

        print(f'Writing *th.nc takes {time()-t0} seconds')

class NudgeTS:

    def __init__(self):
        pass

    def fetch_data(self, outdir: Union[str, os.PathLike], hgrid, vgrid, start_date, rnday, include):

        outdir = pathlib.Path(outdir)

        timevector = np.arange(start_date, start_date + timedelta(days=rnday+1), \
            timedelta(days=1)).astype(datetime)

        vd = Vgrid.open(vgrid)
        sigma = vd.sigma

        #get coords of SCHISM
        loni=hgrid.nodes.coords[:,0]
        lati=hgrid.nodes.coords[:,1]

        #get bathymetry
        depth = hgrid.values

        #compute zcor
        zcor = depth[:,None]*sigma
        nvrt=zcor.shape[1]
        #print(f'zcor at node 1098677 is {zcor[1098676,:]}')

        #Get open nudge array 
        print(len(include))
        nlon = hgrid.coords[include, 0]
        nlat = hgrid.coords[include, 1]
        xi,yi = transform_ll_to_cpp(nlon, nlat)
        bxy = np.c_[yi, xi]

        zcor2=zcor[include,:]
        idxs=np.where(zcor2 > 5500)
        zcor2[idxs]=5500.0-1.0e-6
        #print(f'zcor2 at node 200 is {zcor2[199,:]}')

        #construct schism grid
        x2i=np.tile(xi,[nvrt,1]).T
        y2i=np.tile(yi,[nvrt,1]).T
        bxyz=np.c_[zcor2.reshape(np.size(zcor2)),y2i.reshape(np.size(y2i)),x2i.reshape(np.size(x2i))]
        print('Computing SCHISM zcor is done!')

        #allocate output variables
        nNode=len(include)
        one=1
        ntimes=rnday+1

        timeseries_s=np.zeros([ntimes,nNode,nvrt,one])
        timeseries_t=np.zeros([ntimes,nNode,nvrt,one])
        ndt=np.zeros([ntimes])

        print('**** Accessing RTOFS data*****')

        xmin, xmax = np.min(nlon), np.max(nlon)
        ymin, ymax = np.min(nlat), np.max(nlat)

        # convert hgrid lon [180, 180) to [0, 360)
        xmin = xmin + 360. if xmin < 0 else xmin
        xmax = xmax + 360. if xmax < 0 else xmax
        bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)
        print(f'xmin is {xmin}, xmax is {xmax}')

        baseurl = f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
        for it, date in enumerate(timevector):
            t0=time()

            ndt[it]=it
             
            ssh_url = f'{baseurl}{date.strftime("%Y%m%d")}/rtofs_glo_2ds_forecast_3hrly_diag'
            salt_url = f'{baseurl}{date.strftime("%Y%m%d")}/rtofs_glo_3dz_forecast_daily_salt'
            temp_url = f'{baseurl}{date.strftime("%Y%m%d")}/rtofs_glo_3dz_forecast_daily_temp'
            uvel_url = f'{baseurl}{date.strftime("%Y%m%d")}/rtofs_glo_3dz_forecast_daily_uvel'
            vvel_url = f'{baseurl}{date.strftime("%Y%m%d")}/rtofs_glo_3dz_forecast_daily_vvel'

            #salt
            ds = Dataset(salt_url)
            time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2 = get_idxs(date, ds, bbox)
            dep = ds['lev'][:]
            print(f'max depth is {np.max(dep)}')
            salt = np.squeeze(ds['salinity'][time_idx,:,lat_idx1:lat_idx2+1,lon_idx1:lon_idx2+1])

            salt_int = interp_to_points_3d(dep, y2, x2, bxyz, salt)
            salt_int = salt_int.reshape(zcor2.shape)
            #timeseries_s[it,:,:,0]=salt_int
            timeseries_s[it,:,:,0] = salt_int

            ds.close()

            #temp
            ds = Dataset(temp_url)
            dep = ds['lev'][:]
            temp = np.squeeze(ds['temperature'][time_idx,:,lat_idx1:lat_idx2+1,lon_idx1:lon_idx2+1])

            temp_int = interp_to_points_3d(dep, y2, x2, bxyz, temp)
            temp_int = temp_int.reshape(zcor2.shape)
            #timeseries_t[it,:,:,0]=temp_int
            timeseries_s[it,:,:,0] = temp_int
            ds.close()

        with Dataset(outdir / 'SAL_nu.nc', 'w', format='NETCDF4') as dst:
        #dimensions
            dst.createDimension('node', nNode)
            dst.createDimension('nLevels', nvrt)
            dst.createDimension('one', one)
            dst.createDimension('time', None)
        #variables
            dst.createVariable('time', 'f', ('time',))
            dst['time'][:] = ndt

            dst.createVariable('map_to_global_node', 'i4', ('node',))
            dst['map_to_global_node'][:] = include+1

            dst.createVariable('tracer_concentration', 'f', ('time', 'node', 'nLevels', 'one'))
            dst['tracer_concentration'][:,:,:,:] = timeseries_s

        with Dataset(outdir / 'TEM_nu.nc', 'w', format='NETCDF4') as dst:
        #dimensions
            dst.createDimension('node', nNode)
            dst.createDimension('nLevels', nvrt)
            dst.createDimension('one', one)
            dst.createDimension('time', None)
        #variables
            dst.createVariable('time', 'f', ('time',))
            dst['time'][:] = ndt

            dst.createVariable('map_to_global_node', 'i4', ('node',))
            dst['map_to_global_node'][:] = include+1

            dst.createVariable('tracer_concentration', 'f', ('time', 'node', 'nLevels', 'one'))
            dst['tracer_concentration'][:,:,:,:] = timeseries_t

        print(f'Writing *_nu.nc takes {time()-t0} seconds')
