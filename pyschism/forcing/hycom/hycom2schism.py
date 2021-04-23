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
from netCDF4 import Dataset
from matplotlib.transforms import Bbox

from pyschism.mesh.base import Nodes, Elements
from pyschism.mesh.vgrid import Vgrid

from pyschism.dates import localize_datetime, nearest_cycle_date, pivot_time

logger = logging.getLogger(__name__)

class Nudge:

    def __init__(self):
  
        pass

    def gen_nudge(self, outdir: Union[str, os.PathLike], hgrid):

        @jit(nopython=True, parallel=True)
        def compute_nudge(lon, lat, nnode, opbd2, out):

            for idn in prange(nnode):
                if idn in opbd2:
                    rnu = rnu_max
                    distmin = 0.
                else:
                    distmin = np.finfo(np.float64).max
                    for j in opbd2:
                        tmp = np.square(lon[idn]-lon[j-1]) + np.square(lat[idn]-lat[j-1])
                        rl2 = np.sqrt(tmp)
                        if rl2 < distmin:
                            distmin=rl2
                rnu = 0.
                if distmin <= rlmax:
                    rnu = (1-distmin/rlmax)*rnu_max
                out[idn] = rnu
        
        #self.hgrid = hgrid
        outdir = pathlib.Path(outdir)

        hgrid = hgrid.to_dict()
        nodes = hgrid['nodes']
        elements = hgrid['elements']
        NE, NP = len(elements), len(nodes)
        lon = []
        lat = []
        for id, (coords, values) in nodes.items():
            lon.append(coords[0])
            lat.append(coords[1])

        bnd = hgrid['boundaries']
        opbd = bnd[None][0]['indexes']
        opbd2 = []
        for idn in opbd:
            opbd2.append(int(idn))

        #Max relax distance in degr
        rlmax = 1.5
        #Max relax strength in days
        rnu_day = 0.25
        rnu = 0
        rnu_max = 1./rnu_day/86400.

        out = np.zeros([NP])
        t0 = time()
        compute_nudge(lon, lat, NP, opbd2, out)
        print(f'It took {time() -t0}')

        nudge = [f"{rlmax}, {rnu_day}"]
        nudge.extend("\n")
        nudge.append(f"{NE} {NP}")
        nudge.extend("\n")
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

        with open(outdir / 'nudge_pyschism.gr3','w+') as fid:
            fid.writelines(nudge)

class InitialTS():

    def __init__(self):

        pass

    logger.info('Fetching RTOFS data')
    '''
    The tmpdir is not enough for saving executed files, which will
    cause segmentation fault.
    '''

    def fetch_data(self, outdir: Union[str, os.PathLike], hgrid, vgrid, start_date):

        outdir = pathlib.Path(outdir)

        self.start_date = start_date
        date = self.start_date - timedelta(days=1)
    
        vd=Vgrid()
        sigma=vd.read_vgrid(vgrid)
        
        self.hgrid = hgrid
        bbox = self.hgrid.get_bbox(crs='epsg:4326', output_type='bbox')
        xmin = bbox.x0 + 360. if bbox.x0 < 0 else bbox.x0
        xmax = bbox.x1 + 360. if bbox.x1 < 0 else bbox.x1
        bbox = Bbox.from_extents(xmin, bbox.y0, xmax, bbox.y1)

        hgrid = self.hgrid.to_dict()
        nodes = Nodes(hgrid['nodes'])

        #get coords of SCHISM
        loni = nodes.coords[:,0]
        lati = nodes.coords[:,1]
        bxy = np.c_[lati, loni]

        #get bathymetry
        depth = nodes.values
        idxs = np.where(depth < 0.11)
        depth[idxs] = 0.11
        #compute zcor
        zcor = depth[:,None]*sigma
        nvrt=zcor.shape[1]

        #Get properties of elements (tris, quads, side)
        eles = hgrid['elements']
        elements = Elements(nodes, eles)

        tris = elements.triangles
        tidxs = elements.tri_idxs
        #print(tris[0,:])
        #tnd = tris.shape[1]

        quads = elements.quads
        qidxs = elements.qua_idxs
        #qnd = quads.shape[1]

        side = elements.side
        nside = side.shape[0]

        #construct schism grid
        lon2i=np.tile(loni,[nvrt,1]).T
        lat2i=np.tile(lati,[nvrt,1]).T
        bxyz=np.c_[zcor.reshape(np.size(zcor)),lat2i.reshape(np.size(lat2i)),lon2i.reshape(np.size(lon2i))]
        print('Computing SCHISM zcor is done!')

        #Get hycom data
        planA=True
        planB=False
        planC=False

        #Plan A: ROTFS 
        if planA:
            try:
                print('**** Accessing RTOFS data*****')
                base_url='http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
                nc_salt=Dataset(f'{base_url}'
                    f'{date.strftime("%Y%m%d")}/rtofs_glo_3dz_nowcast_daily_salt')
                nc_temp=Dataset(f'{base_url}'
                    f'{date.strftime("%Y%m%d")}/rtofs_glo_3dz_nowcast_daily_temp')
                lon=nc_salt['lon'][:]
                lat=nc_salt['lat'][:]
                dep=nc_salt['lev'][:]

                lat_idxs=np.where((lat>=bbox.ymin) & (lat<=bbox.ymax))[0]
                lon_idxs=np.where((lon>=bbox.xmin-1.0) & (lon<=bbox.xmax+1.0))[0]

                salt=np.squeeze(nc_salt['salinity'][1,:,lat_idxs,lon_idxs])
                temp=np.squeeze(nc_temp['temperature'][1,:,lat_idxs,lon_idxs])

            except Exception as e:
                print(e)
                planB = True

        if planB:
            try:
                '''
                Here temperature is in-situ temperature.
                TODO: convert in-situ temperature into potential temperature.
                '''
                print('**** Accessing GOFS data*****') 
                base_url='https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/FMRC/runs/GLBy0.08_930_FMRC_RUN_'
                nc=Dataset(f'{base_url}'
                    f'{date.strftime("%Y-%m-%dT12:00:00Z")}')
                lon=nc['lon'][:]
                lat=nc['lat'][:]
                dep=nc['depth'][:]

                lat_idxs=np.where((lat>=bbox.ymin)&(lat<=bbox.ymax))[0]
                lon_idxs=np.where((lon>=bbox.xmin-1.0) & (lon<=bbox.xmax+1.0))[0]

                salt=np.squeeze(nc['salinity'][4,:,lat_idxs,lon_idxs])
                temp=np.squeeze(nc['water_temp'][4,:,lat_idxs,lon_idxs])
            except Exception as e:
                print(e)
                planC=True

        if planC:
            print('Use pre-downloaded data!')

        print('****Interpolation starts****')

        #Interpolate HYCOM to SCHISM
        lon=lon[lon_idxs]
        for ilon in range(lon_idxs.size):
            if lon[ilon] > 180:
                lon[ilon] = lon[ilon]-360.
        lat=lat[lat_idxs]

        #change missing value to nan
        idxs = np.where(salt > 30000)
        salt[idxs]=float('nan')
        idxs = np.where(temp > 30000)
        temp[idxs]=float('nan')

        t0 = time()
        salt_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),salt,'nearest', bounds_error=False, fill_value = float('nan'))
        salt_int = salt_fd(bxyz)
        idxs = np.isnan(salt_int)
        if np.sum(idxs)!=0:
            salt_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], salt_int[~idxs], bxyz[idxs,:],'nearest')
        idxs = np.isnan(salt_int)
        if np.sum(idxs)!=0:
            print(f'There is still missing value for salinity!')
            sys.exit()
        salt_int = salt_int.reshape(zcor.shape)
        #print(salt_int.shape)

        #temperature
        temp_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),temp,'nearest', bounds_error=False, fill_value=float('nan'))
        temp_int = temp_fd(bxyz)
        idxs = np.isnan(temp_int)
        if np.sum(idxs)!=0:
            temp_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], temp_int[~idxs], bxyz[idxs,:],'nearest')
        idxs = np.isnan(temp_int)
        if np.sum(idxs)!=0:
            print(f'There is still missing value for temperature!')
            sys.exit()
        temp_int = temp_int.reshape(zcor.shape)
        #print(temp_int.shape)

        print(f'Interpolation takes {time()-t0} seconds')

        #Create hotstart.nc
        t0=time()
        NP = len(loni)
        NE = len(eles)
        ntracers = 2
        one = 1

        #Compute tr_nd and tr_el
        tr_nd = np.zeros([NP, nvrt, 2])
        tr_nd[:,:,0] = salt_int
        tr_nd[:,:,1] = temp_int
        tr_el = np.zeros([NE, nvrt, 2])

        for k in np.arange(nvrt):
            #triangles
            tr_el[tidxs,k,0] = np.mean(salt_int[tris[:,:],k], axis=1)
            tr_el[tidxs,k,1] = np.mean(temp_int[tris[:,:],k], axis=1)
            #quads 
            tr_el[qidxs,k,0] = np.mean(salt_int[quads[:,:],k], axis=1)
            tr_el[qidxs,k,1] = np.mean(temp_int[quads[:,:],k], axis=1)

        with Dataset(outdir / 'hotstart.nc', 'w', format='NETCDF4') as dst:
            #dimensions
            dst.createDimension('node', NP)
            dst.createDimension('elem', NE)
            dst.createDimension('side', nside)
            dst.createDimension('nVert', nvrt)
            dst.createDimension('ntracers', ntracers)
            dst.createDimension('one', one)

            #variables
            dst.createVariable('time', 'd', ('one',))
            dst['time'][:] = 0.

            dst.createVariable('iths', 'i4', ('one',))
            dst['iths'][:] = 0

            dst.createVariable('ifile', 'i4', ('one',))
            dst['ifile'][:] = 1

            dst.createVariable('idry_e', 'i4', ('elem',))
            dst['idry_e'][:] = np.zeros(NE).astype('int32')

            dst.createVariable('idry_s', 'i4', ('side',))
            dst['idry_s'][:] = np.zeros(nside).astype('int32')

            dst.createVariable('idry', 'i4', ('node',))
            dst['idry'][:] = np.zeros(NP).astype('int32')

            dst.createVariable('eta2', 'd', ('node',))
            dst['eta2'][:] = np.zeros(NP)

            dst.createVariable('we', 'd', ('elem', 'nVert'))
            dst['we'][:,:] = np.zeros([NE,nvrt])

            dst.createVariable('tr_el', 'd', ('elem', 'nVert', 'ntracers'))
            dst['tr_el'][:,:,:] = tr_el

            dst.createVariable('tr_nd', 'd', ('node', 'nVert', 'ntracers'))
            dst['tr_nd'][:,:,:] = tr_nd

            dst.createVariable('tr_nd0', 'd', ('node', 'nVert', 'ntracers'))
            dst['tr_nd0'][:,:,:] = tr_nd

            dst.createVariable('su2', 'd', ('side', 'nVert'))
            dst['su2'][:,:] = np.zeros([nside,nvrt])

            dst.createVariable('sv2', 'd', ('side', 'nVert'))
            dst['sv2'][:,:] = np.zeros([nside,nvrt])
 
            dst.createVariable('q2', 'd', ('node', 'nVert'))
            dst['q2'][:,:] = np.zeros([NP,nvrt])

            dst.createVariable('xl', 'd', ('node', 'nVert'))
            dst['xl'][:,:] = np.zeros([NP,nvrt])

            dst.createVariable('dfv', 'd', ('node', 'nVert'))
            dst['dfv'][:,:] = np.zeros([NP,nvrt])

            dst.createVariable('dfh', 'd', ('node', 'nVert'))
            dst['dfh'][:,:] = np.zeros([NP,nvrt])

            dst.createVariable('dfq1', 'd', ('node', 'nVert'))
            dst['dfq1'][:,:] = np.zeros([NP,nvrt])

            dst.createVariable('dfq2', 'd', ('node', 'nVert'))
            dst['dfq2'][:,:] = np.zeros([NP,nvrt])

        print(f'Writing hotstart.nc takes {time()-t0} seconds')

class OpenBoundaryInventory():

    def __init__(self):

        pass

    def fetch_data(self, outdir: Union[str, os.PathLike], hgrid, vgrid, start_date, rnday):

        outdir = pathlib.Path(outdir)

        self.start_date = start_date
        self.rnday=rnday
        date = self.start_date - timedelta(days=1)

        vd=Vgrid()
        sigma=vd.read_vgrid(vgrid)

        hgrid = hgrid.to_dict()
        nodes = Nodes(hgrid['nodes']) 
        #get coords of SCHISM
        loni = nodes.coords[:,0]
        lati = nodes.coords[:,1]

        #get bathymetry
        depth = nodes.values
        idxs = np.where(depth < 0.11)
        depth[idxs] = 0.11
        #compute zcor
        zcor = depth[:,None]*sigma
        nvrt=zcor.shape[1]

        #Get open boundary
        opbd = hgrid['boundaries'][None][0]['indexes']
        blon = []
        blat = []
        idxs = []
        for idn in opbd:
            blon.append(loni[int(idn)-1])
            blat.append(lati[int(idn)-1])
            idxs.append(int(idn)-1)
        NOP = len(blon)
        bxy = np.c_[blat, blon]

        zcor2=zcor[idxs,:]

        xmin, xmax = np.min(blon), np.max(blon)
        ymin, ymax = np.min(blat), np.max(blat)

        xmin = xmin + 360. if xmin < 0 else xmin
        xmax = xmax + 360. if xmax < 0 else xmax
        bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)

        #construct schism grid
        lon2i=np.tile(blon,[nvrt,1]).T
        lat2i=np.tile(blat,[nvrt,1]).T
        bxyz=np.c_[zcor2.reshape(np.size(zcor2)),lat2i.reshape(np.size(lat2i)),lon2i.reshape(np.size(lon2i))]
        print('Computing SCHISM zcor is done!')

        #Get hycom data
        baseurl_rtofs='http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
        baseurl_gofs='https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/FMRC/runs/GLBy0.08_930_FMRC_RUN_'

        planA=True
        planB=False
        planC=False

        #Plan A: ROTFS 
        if planA:
            try:
                print('**** Accessing RTOFS data*****')
                nc_salt=Dataset(f'{baseurl_rtofs}'
                    f'{date.strftime("%Y%m%d")}/rtofs_glo_3dz_nowcast_daily_salt')
                nc_temp=Dataset(f'{baseurl_rtofs}'
                    f'{date.strftime("%Y%m%d")}/rtofs_glo_3dz_nowcast_daily_temp')
                nc_uvel=Dataset(f'{baseurl_rtofs}'
                    f'{date.strftime("%Y%m%d")}/rtofs_glo_3dz_nowcast_daily_uvel')
                nc_vvel=Dataset(f'{baseurl_rtofs}'
                    f'{date.strftime("%Y%m%d")}/rtofs_glo_3dz_nowcast_daily_vvel')
                nc_ssh=Dataset(f'{baseurl_rtofs}'
                    f'{date.strftime("%Y%m%d")}/rtofs_glo_2ds_forecast_3hrly_diag')
                lon=nc_salt['lon'][:]
                lat=nc_salt['lat'][:]
                dep=nc_salt['lev'][:]

                lat_idxs=np.where((lat>=bbox.ymin) & (lat<=bbox.ymax))[0]
                lon_idxs=np.where((lon>=bbox.xmin-2.0) & (lon<=bbox.xmax+2.0))[0]
                print(f'The lengh is lon {len(lon_idxs)}, lat {len(lat_idxs)}')

                salt=nc_salt['salinity'][1:,:,lat_idxs,lon_idxs]
                temp=nc_temp['temperature'][1:,:,lat_idxs,lon_idxs]
                uvel=nc_uvel['u'][1:,:,lat_idxs,lon_idxs]
                vvel=nc_vvel['v'][1:,:,lat_idxs,lon_idxs]
                ssh=nc_ssh['ssh'][8::8,lat_idxs,lon_idxs]
                print(f'Salinity shape is {salt.shape}')

            except Exception as e:
                print(e)
                planB = True

        if planB:
            try:
                print('**** Accessing GOFS data*****')
                nc=Dataset(f'{baseurl_gofs}'
                    f'{date.strftime("%Y-%m-%dT12:00:00Z")}')
                lon=nc['lon'][:]
                lat=nc['lat'][:]
                dep=nc['depth'][:]

                lat_idxs=np.where((lat>=bbox.ymin)&(lat<=bbox.ymax))[0]
                lon_idxs=np.where((lon>=bbox.xmin-1.0) & (lon<=bbox.xmax+1.0))[0]

                salt=nc['salinity'][4::8,:,lat_idxs,lon_idxs]
                temp=nc['water_temp'][4::8,:,lat_idxs,lon_idxs]
                uvel=nc['water_u'][4::8,:,lat_idxs,lon_idxs]
                vvel=nc['water_v'][4::8,:,lat_idxs,lon_idxs]
                ssh=nc['surf_el'][4::8,lat_idxs,lon_idxs]

            except Exception as e:
                print(e)
                planC=True

        if planC:
            try:
                print('**** Accessing GOFS best time series dataset*****')
                nc=Dataset('https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/FMRC/GLBy0.08_930_FMRC_best.ncd')
                lon=nc['lon'][:]
                lat=nc['lat'][:]
                dep=nc['depth'][:]

                lat_idxs=np.where((lat>=bbox.ymin)&(lat<=bbox.ymax))[0]
                lon_idxs=np.where((lon>=bbox.xmin-1.0) & (lon<=bbox.xmax+1.0))[0]

                salt=nc['salinity'][4::8,:,lat_idxs,lon_idxs]
                temp=nc['water_temp'][4::8,:,lat_idxs,lon_idxs]
                uvel=nc['water_u'][4::8,:,lat_idxs,lon_idxs]
                vvel=nc['water_v'][4::8,:,lat_idxs,lon_idxs]
                ssh=nc['surf_el'][4::8,lat_idxs,lon_idxs]

            except Exception as e:
                print(e)

        print('****Interpolation starts****')
        #Interpolate HYCOM to SCHISM
        lon=lon[lon_idxs]
        for ilon in range(lon_idxs.size):
            if lon[ilon] > 180:
                lon[ilon] = lon[ilon]-360.
        lat=lat[lat_idxs]

        #change missing value to nan
        #vars={'salt', 'temp', 'uvel', 'vvel', 'ssh'}
        idxs = np.where(salt > 30000)
        salt[idxs]=float('nan')
        idxs = np.where(temp > 30000)
        temp[idxs]=float('nan')
        idxs = np.where(uvel > 30000)
        uvel[idxs]=float('nan')
        idxs = np.where(vvel > 30000)
        vvel[idxs]=float('nan')
        idxs = np.where(ssh > 30000)
        ssh[idxs]=float('nan')
        #print(f'The shape of salt is {salt.shape}')

        ntimes=self.rnday+1
        nComp1=1
        nComp2=2
        one=1
        timeseries_s=np.zeros([ntimes,NOP,nvrt,nComp1])
        timeseries_t=np.zeros([ntimes,NOP,nvrt,nComp1])
        timeseries_el=np.zeros([ntimes,NOP,nComp1])
        timeseries_uv=np.zeros([ntimes,NOP,nvrt,nComp2])
        ndt=np.zeros([ntimes])
        t0 = time()

        for it in np.arange(ntimes):
            ndt[it]=it*24*3600.
            #salt
            salt_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(salt[it,:,:,:]),'nearest', bounds_error=False, fill_value = float('nan'))
            salt_int = salt_fd(bxyz)
            idxs = np.isnan(salt_int)
            if np.sum(idxs)!=0:
                salt_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], salt_int[~idxs], bxyz[idxs,:],'nearest')
            idxs = np.isnan(salt_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for salinity!')
                sys.exit()
            salt_int = salt_int.reshape(zcor2.shape)
            timeseries_s[it,:,:,0]=salt_int

            #temp
            temp_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(temp[it,:,:,:]),'nearest', bounds_error=False, fill_value = float('nan'))
            temp_int = temp_fd(bxyz)
            idxs = np.isnan(temp_int)
            if np.sum(idxs)!=0:
                temp_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], temp_int[~idxs], bxyz[idxs,:],'nearest')
            idxs = np.isnan(temp_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for temperature!')
                sys.exit()
            temp_int = temp_int.reshape(zcor2.shape)
            timeseries_t[it,:,:,0]=temp_int

            #uvel
            uvel_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(uvel[it,:,:,:]),'nearest', bounds_error=False, fill_value = float('nan'))
            uvel_int = uvel_fd(bxyz)
            idxs = np.isnan(uvel_int)
            if np.sum(idxs)!=0:
                uvel_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], uvel_int[~idxs], bxyz[idxs,:],'nearest')
            idxs = np.isnan(uvel_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for u-velocity!')
                sys.exit()
            uvel_int = uvel_int.reshape(zcor2.shape)
            timeseries_uv[it,:,:,0]=uvel_int

            #vvel
            vvel_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(vvel[it,:,:,:]),'nearest', bounds_error=False, fill_value = float('nan'))
            vvel_int = vvel_fd(bxyz)
            idxs = np.isnan(vvel_int)
            if np.sum(idxs)!=0:
                vvel_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], vvel_int[~idxs], bxyz[idxs,:],'nearest')
            idxs = np.isnan(vvel_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for v-velocity!')
                sys.exit()
            vvel_int = vvel_int.reshape(zcor2.shape)
            timeseries_uv[it,:,:,1]=vvel_int

            #ssh
            ssh_fd=sp.interpolate.RegularGridInterpolator((lat,lon),np.squeeze(ssh[it,:,:]),'nearest', bounds_error=False, fill_value = float('nan'))
            ssh_int = ssh_fd(bxy)
            idxs = np.isnan(ssh_int)
            if np.sum(idxs)!=0:
                ssh_int[idxs]=sp.interpolate.griddata(bxy[~idxs,:], ssh_int[~idxs], bxy[idxs,:],'nearest')
            idxs = np.isnan(ssh_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for ssh!')
                sys.exit()
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
                dst.createDimension('nComponents', nComp1)
                #variables
                dst.createVariable('time_step', 'f', ('one',))
                dst['time_step'][:] = 86400

                dst.createVariable('time', 'f', ('time',))
                dst['time'][:] = ndt

                dst.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nComponents'))
                dst['time_series'][:,:,:] = timeseries_el

            print(f'Writing hotstart.nc takes {time()-t0} seconds')