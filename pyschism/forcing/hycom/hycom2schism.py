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

class HotStartInventory():

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
                    f'{date.strftime("%Y-%m-%dT12:00:00Z")}/')
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

    def fetch_data(self, outdir: Union[str, os.PathLike], hgrid, start_date, rnday, idx_min, idx_max, jdx_min, jdx_max):

        pass
