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
from pyschism.mesh.gridgr3 import Gr3Field


logger = logging.getLogger(__name__)

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
        print('No data for {date}')
    return database

def get_idxs(date, database, bbox):

    if date.strftime("%Y-%m-%d") >= datetime.now().strftime("%Y-%m-%d"):
        date2 = datetime.now() - timedelta(days=1)
        baseurl = f'https://tds.hycom.org/thredds/dodsC/{database}/FMRC/runs/GLBy0.08_930_FMRC_RUN_{date2.strftime("%Y-%m-%dT12:00:00Z")}?depth[0:1:-1],lat[0:1:-1],lon[0:1:-1],time[0:1:-1]'
    else:
        baseurl=f'https://tds.hycom.org/thredds/dodsC/{database}?lat[0:1:-1],lon[0:1:-1],time[0:1:-1],depth[0:1:-1]'

    ds=Dataset(baseurl)
    time1=ds['time']
    times=nc.num2date(time1,units=time1.units,only_use_cftime_datetimes=False)
    
    lon=ds['lon'][:]
    lat=ds['lat'][:]
    dep=ds['depth'][:]
    lat_idxs=np.where((lat>=bbox.ymin)&(lat<=bbox.ymax))[0]
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
    x2, y2=transform_ll_to_cpp(lon, lat)

    idxs=np.where( date == times)[0]
    #check if time_idx is empty
    if len(idxs) == 0:
        #If there is missing data, use the data from the previous days, the maximum searching days is 3. Otherwise, stop.
        for i in np.arange(0,3):
            date_before=(date - timedelta(days=int(i)+1)) #.astype(datetime)
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

class OpenBoundaryInventory():

    def __init__(self):

        pass

    def fetch_data(self, outdir: Union[str, os.PathLike], hgrid, vgrid, start_date, rnday):

        outdir = pathlib.Path(outdir)

        self.start_date = start_date
        self.rnday=rnday
        self.timevector=np.arange(
            self.start_date,
            self.start_date + timedelta(days=self.rnday+1),
            timedelta(days=1)).astype(datetime)

        vd=Vgrid.open(vgrid)
        sigma=vd.sigma

        #get bathymetry
        depth = hgrid.values

        #compute zcor
        zcor = depth[:,None]*sigma
        nvrt=zcor.shape[1]
        #print(f'zcor at node 1098677 is {zcor[1098676,:]}')

        #Get open boundary 
        gdf=hgrid.boundaries.open.copy()
        opbd=[]
        for boundary in gdf.itertuples():
            opbd.extend(list(boundary.indexes))
        blon = hgrid.coords[opbd,0]
        blat = hgrid.coords[opbd,1]
        
        NOP = len(blon)
        xi,yi = transform_ll_to_cpp(blon, blat)
        bxy = np.c_[yi, xi]

        zcor2=zcor[opbd,:]
        idxs=np.where(zcor2 > 5000)
        #print(idxs)
        zcor2[idxs]=5000.0-1.0e-6
        print(f'zcor2 at node 200 is {zcor2[199,:]}')

        xmin, xmax = np.min(blon), np.max(blon)
        ymin, ymax = np.min(blat), np.max(blat)

        #construct schism grid
        x2i=np.tile(xi,[nvrt,1]).T
        y2i=np.tile(yi,[nvrt,1]).T
        bxyz=np.c_[zcor2.reshape(np.size(zcor2)),y2i.reshape(np.size(y2i)),x2i.reshape(np.size(x2i))]
        print('Computing SCHISM zcor is done!')

        ntimes=self.rnday+1
        nComp1=1
        nComp2=2
        one=1
        timeseries_s=np.zeros([ntimes,NOP,nvrt,nComp1])
        timeseries_t=np.zeros([ntimes,NOP,nvrt,nComp1])
        timeseries_el=np.zeros([ntimes,NOP,nComp1])
        timeseries_uv=np.zeros([ntimes,NOP,nvrt,nComp2])
        ndt=np.zeros([ntimes])

        print('**** Accessing GOFS data*****')
        for it, date in enumerate(self.timevector):
            t0=time()

            database=get_database(date)
            print(f'Fetching data for {date} from database {database}')

            if date.strftime("%Y-%m-%d") >= datetime(2017, 10, 1).strftime("%Y-%m-%d"):
                xmin = xmin + 360. if xmin < 0 else xmin
                xmax = xmax + 360. if xmax < 0 else xmax
                bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)
            else:
                bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)

            time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2 = get_idxs(date, database, bbox)

            if date.strftime("%Y-%m-%d") >= datetime.now().strftime("%Y-%m-%d"):
                date2 = datetime.now() - timedelta(days=1)
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
            print(url)
             
            ds=Dataset(url)
            salt=np.squeeze(ds['salinity'][:,:,:])
            temp=np.squeeze(ds['water_temp'][:,:,:])
            uvel=np.squeeze(ds['water_u'][:,:,:])
            vvel=np.squeeze(ds['water_v'][:,:,:])
            ssh=np.squeeze(ds['surf_el'][:,:])
            #print(f'The shape of temp is {temp.shape}')

            #Convert temp to potential temp
            nz=temp.shape[0]
            ny=temp.shape[1]
            nx=temp.shape[2]
            dep=ds['depth'][:]
            pr=np.ones(temp.shape)
            pre=pr*dep[:,None, None]
            Pr=np.zeros(temp.shape)
            ptemp=sw.ptmp(salt, temp, pre, Pr)*1.00024

            print('****Interpolation starts****')

            ndt[it]=it*24*3600.
            #salt
            salt_int = interp_to_points_3d(dep, y2, x2, bxyz, salt)
            salt_int = salt_int.reshape(zcor2.shape)
            timeseries_s[it,:,:,0]=salt_int

            #temp
            temp_int = interp_to_points_3d(dep, y2, x2, bxyz, ptemp)
            temp_int = temp_int.reshape(zcor2.shape)
            timeseries_t[it,:,:,0]=temp_int

            #uvel
            uvel_int = interp_to_points_3d(dep, y2, x2, bxyz, uvel)
            uvel_int = uvel_int.reshape(zcor2.shape)
            timeseries_uv[it,:,:,0]=uvel_int

            #vvel
            vvel_int = interp_to_points_3d(dep, y2, x2, bxyz, vvel)
            vvel_int = vvel_int.reshape(zcor2.shape)
            timeseries_uv[it,:,:,1]=vvel_int

            #ssh
            ssh_int = interp_to_points_2d(y2, x2, bxy, ssh)
            timeseries_el[it,:,0]=ssh_int

            print(f'Interpolation takes {time()-t0} seconds')

        #Write to files
        t0=time()

        with Dataset(outdir / 'SAL_3D.th.nc', 'w', format='NETCDF4') as dst:
        #dimensions
            dst.createDimension('nOpenBndNodes', NOP)
            dst.createDimension('one', one)
            dst.createDimension('time', None)
            dst.createDimension('nLevels', nvrt)
            dst.createDimension('nComponents', nComp1)
        #variables
            dst.createVariable('time_step', 'f', ('one',))
            dst['time_step'][:] = 86400

            dst.createVariable('time', 'f', ('time',))
            dst['time'][:] = ndt

            dst.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
            dst['time_series'][:,:,:,:] = timeseries_s

        with Dataset(outdir / 'TEM_3D.th.nc', 'w', format='NETCDF4') as dst:
            #dimensions
            dst.createDimension('nOpenBndNodes', NOP)
            dst.createDimension('one', one)
            dst.createDimension('time', None)
            dst.createDimension('nLevels', nvrt)
            dst.createDimension('nComponents', nComp1)
            #variables
            dst.createVariable('time_step', 'f', ('one',))
            dst['time_step'][:] = 86400

            dst.createVariable('time', 'f', ('time',))
            dst['time'][:] = ndt

            dst.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
            dst['time_series'][:,:,:,:] = timeseries_t

        with Dataset(outdir / 'uv3D.th.nc', 'w', format='NETCDF4') as dst:
            #dimensions
            dst.createDimension('nOpenBndNodes', NOP)
            dst.createDimension('one', one)
            dst.createDimension('time', None)
            dst.createDimension('nLevels', nvrt)
            dst.createDimension('nComponents', nComp2)
            #variables
            dst.createVariable('time_step', 'f', ('one',))
            dst['time_step'][:] = 86400

            dst.createVariable('time', 'f', ('time',))
            dst['time'][:] = ndt
            
            dst.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
            dst['time_series'][:,:,:,:] = timeseries_uv

        with Dataset(outdir / 'elev2D.th.nc', 'w', format='NETCDF4') as dst:
            #dimensions
            dst.createDimension('nOpenBndNodes', NOP)
            dst.createDimension('one', one)
            dst.createDimension('time', None)
            dst.createDimension('nLevels', one)
            dst.createDimension('nComponents', nComp1)
            #variables
            dst.createVariable('time_step', 'f', ('one',))
            dst['time_step'][:] = 86400

            dst.createVariable('time', 'f', ('time',))
            dst['time'][:] = ndt

            dst.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
            dst['time_series'][:,:,:,:] = timeseries_el

        print(f'Writing *th.nc takes {time()-t0} seconds')

class Nudge:

    def __init__(self):

        self.include = None


    def gen_nudge(self, outdir: Union[str, os.PathLike], hgrid, rlmax = 1.5, rnu_day=0.25):

        @jit(nopython=True, parallel=True)
        def compute_nudge(lon, lat, nnode, opbd, out):
          
            rnu_max = 1.0 / rnu_day / 86400.0
            rnu = 0

            for idn in prange(nnode):
                if idn in opbd:
                    rnu = rnu_max
                    distmin = 0.
                else:
                    distmin = np.finfo(np.float64).max
                    for j in opbd:
                        tmp = np.square(lon[idn]-lon[j]) + np.square(lat[idn]-lat[j])
                        rl2 = np.sqrt(tmp)
                        if rl2 < distmin:
                            distmin=rl2
                rnu = 0.
                if distmin <= rlmax:
                    rnu = (1-distmin/rlmax)*rnu_max
                    #idxs_nudge[idn]=1 #idn
                out[idn] = rnu

        outdir = pathlib.Path(outdir)

        #get nudge zone
        lon=hgrid.coords[:,0]
        lat=hgrid.coords[:,1]

        #Get open boundary 
        gdf=hgrid.boundaries.open.copy()
        opbd=[]
        for boundary in gdf.itertuples():
            opbd.extend(list(boundary.indexes))
        opbd = np.array(opbd)

        elnode=hgrid.elements.array
        NE, NP = len(elnode), len(lon)

        out = np.zeros([NP])
        idxs_nudge=np.zeros(NP, dtype=int)
        t0 = time()
        #compute_nudge(lon, lat, NP, opbd2, idxs_nudge, out)
        compute_nudge(lon, lat, NP, opbd, out)

        idxs=np.where(out > 0)[0]
        idxs_nudge[idxs]=1
        #expand nudging marker to neighbor nodes
        idxs=np.where(np.max(out[elnode], axis=1) > 0)[0]
        fp=elnode[idxs,-1] < 0
        idxs_nudge[elnode[idxs[fp],:3]]=1
        idxs_nudge[elnode[idxs[~fp],:]]=1

        #idxs_nudge=np.delete(idxs_nudge, np.where(idxs_nudge == -99))
        idxs=np.where(idxs_nudge == 1)[0]
        self.include=idxs
        print(f'len of nudge idxs is {len(idxs)}')
        print(f'It took {time() -t0} sencods to calcuate nudge coefficient')

        nudge = [f"{rlmax}, {rnu_day}"]
        nudge.extend("\n")
        nudge.append(f"{NE} {NP}")
        nudge.extend("\n")
        hgrid = hgrid.to_dict()
        nodes = hgrid['nodes']
        elements = hgrid['elements']
        for idn, (coords, values) in nodes.items():
            line = [f"{idn}"]
            line.extend([f"{x:<.7e}" for x in coords])
            line.extend([f"{out[int(idn)-1]:<.7e}"])
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

        return self.include

    def fetch_data(self, outdir: Union[str, os.PathLike], hgrid, vgrid, start_date, rnday):

        outdir = pathlib.Path(outdir)

        self.start_date = start_date
        self.rnday=rnday
        self.timevector=np.arange(
            self.start_date,
            self.start_date + timedelta(days=self.rnday+1),
            timedelta(days=1)).astype(datetime)

        vd=Vgrid.open(vgrid)
        sigma=vd.sigma

        #Get the index for nudge
        include = self.gen_nudge(outdir,hgrid)

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
        nlon = hgrid.coords[include, 0]
        nlat = hgrid.coords[include, 1]
        xi,yi = transform_ll_to_cpp(nlon, nlat)
        bxy = np.c_[yi, xi]

        zcor2=zcor[include,:]
        idxs=np.where(zcor2 > 5000)
        #print(idxs)
        zcor2[idxs]=5000.0-1.0e-6
        print(f'zcor2 at node 200 is {zcor2[199,:]}')

        xmin, xmax = np.min(nlon), np.max(nlon)
        ymin, ymax = np.min(nlat), np.max(nlat)

        #construct schism grid
        x2i=np.tile(xi,[nvrt,1]).T
        y2i=np.tile(yi,[nvrt,1]).T
        bxyz=np.c_[zcor2.reshape(np.size(zcor2)),y2i.reshape(np.size(y2i)),x2i.reshape(np.size(x2i))]
        print('Computing SCHISM zcor is done!')

        #allocate output variables
        nNode=len(include)
        one=1
        ntimes=self.rnday+1

        timeseries_s=np.zeros([ntimes,nNode,nvrt,one])
        timeseries_t=np.zeros([ntimes,nNode,nvrt,one])
        ndt=np.zeros([ntimes])

        print('**** Accessing GOFS data*****')
        for it, date in enumerate(self.timevector):
            t0=time()

            database=get_database(date)
            print(database)

            if date.strftime("%Y-%m-%d") >= datetime(2017, 10, 1).strftime("%Y-%m-%d"):
                xmin = xmin + 360. if xmin < 0 else xmin
                xmax = xmax + 360. if xmax < 0 else xmax
                bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)
            else:
                bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)

            time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2 = get_idxs(date, database, bbox)

            if date.strftime("%Y-%m-%d") >= datetime.now().strftime("%Y-%m-%d"):
                date2 = datetime.now() - timedelta(days=1)
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
            print(url)

            ds=Dataset(url)
            salt=np.squeeze(ds['salinity'][:,:,:])
            temp=np.squeeze(ds['water_temp'][:,:,:])
            print(f'The shape of temp is {temp.shape}')

            #Convert temp to potential temp
            nz=temp.shape[0]
            ny=temp.shape[1]
            nx=temp.shape[2]
            dep=ds['depth'][:]
            pr=np.ones(temp.shape)
            pre=pr*dep[:,None, None]
            Pr=np.zeros(temp.shape)
            ptemp=sw.ptmp(salt, temp, pre, Pr)*1.00024

            print('****Interpolation starts****')

            ndt[it]=it
            #salt
            salt_int = interp_to_points_3d(dep, y2, x2, bxyz, salt)
            salt_int = salt_int.reshape(zcor2.shape)
            timeseries_s[it,:,:,0]=salt_int

            #temp
            temp_int = interp_to_points_3d(dep, y2, x2, bxyz, ptemp)
            temp_int = salt_int.reshape(zcor2.shape)
            timeseries_t[it,:,:,0]=temp_int
 
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
