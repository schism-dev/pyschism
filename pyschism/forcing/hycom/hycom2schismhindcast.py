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
#from metpy.units import units
#from metpy.calc import height_to_pressure_std
import seawater as sw

from pyschism.mesh.base import Nodes, Elements
from pyschism.mesh.vgrid import Vgrid


logger = logging.getLogger(__name__)

def get_database(date, Bbox=None):
    if date >= datetime(2018, 12, 4):
        database = f'GLBy0.08/expt_93.0'
    elif date >= datetime(2018, 1, 1) and date < datetime(2020, 2, 18):
        database = f'GLBv0.08/expt_93.0'
    elif date >= datetime(2017, 10, 1) and date < datetime(2017, 12, 31):
        database = f'GLBv0.08/expt_92.9'
    elif date >= datetime(2017, 6, 1) and date < datetime(2017, 9, 30):
        database = f'GLBv0.08/expt_57.7'
    elif date >= datetime(2017, 2, 1) and date < datetime(2017, 5, 31):
        database = f'GLBv0.08/expt_92.8'
    elif date >= datetime(2016, 5, 1) and date < datetime(2017, 1, 31):
        database = f'GLBv0.08/expt_57.2'
    elif date >= datetime(2016, 1, 1) and date < datetime(2016, 4, 30):
        database = f'GLBv0.08/expt_56.3'
    elif date >= datetime(1994, 1, 1) and date < datetime(2015, 12, 31):
        database = f'GLBv0.08/expt_53.X/data/{self.start_date.year}'
    else:
        print('No data for {self.start_date}')
    return database

def transform_ll_to_cpp(lon, lat, lonc=-77.07, latc=24.0):
    #lonc=(np.max(lon)+np.min(lon))/2.0
    print(f'lonc is {lonc}')
    #latc=(np.max(lat)+np.min(lat))/2.0
    print(f'latc is {latc}')
    longitude=lon/180*np.pi
    latitude=lat/180*np.pi
    radius=6378206.4
    loncc=lonc/180*np.pi
    latcc=latc/180*np.pi
    lon_new=[radius*(longitude[i]-loncc)*np.cos(latcc) for i in np.arange(len(longitude))]
    lat_new=[radius*latitude[i] for i in np.arange(len(latitude))]

    return np.array(lon_new), np.array(lat_new)

class Nudge:

    def __init__(self):

        self.include = None
  

    def gen_nudge(self, outdir: Union[str, os.PathLike], hgrid):

        @jit(nopython=True, parallel=True)
        def compute_nudge(lon, lat, nnode, opbd2, out):

            for idn in prange(nnode):
                if idn+1 in opbd2:
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
                    #idxs_nudge[idn]=1 #idn
                out[idn] = rnu
        
        #self.hgrid = hgrid
        outdir = pathlib.Path(outdir)

        elnode=hgrid.elements.array
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
        idxs_nudge=np.zeros(NP, dtype=int)
        t0 = time()
        #compute_nudge(lon, lat, NP, opbd2, idxs_nudge, out)
        compute_nudge(lon, lat, NP, opbd2, out)
        
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
        print(f'It took {time() -t0} sencods to calcuate nudge coefficient')

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

        with open(outdir / 'TEM_nudge.gr3','w+') as fid:
            fid.writelines(nudge)

        shutil.copy2(outdir / 'TEM_nudge.gr3', outdir / 'SAL_nudge.gr3')
 
        return self.include

    def fetch_data(self, outdir: Union[str, os.PathLike], hgrid, vgrid, start_date, rnday):

        outdir = pathlib.Path(outdir)

        self.start_date = start_date
        self.rnday=rnday
        #date = self.start_date - timedelta(days=1)
        self.timevector=np.arange(
            self.start_date,
            self.start_date + timedelta(days=self.rnday+1),
            timedelta(days=1))

        vd=Vgrid()
        sigma=vd.read_vgrid(vgrid)

        #Get the index for nudge
        include = self.gen_nudge(outdir,hgrid)
        #print(include[:10])

        hgrid=hgrid.to_dict()
        nodes=Nodes(hgrid['nodes'])
        #get coords of SCHISM
        loni=nodes.coords[:,0]
        lati=nodes.coords[:,1]

        #get bathymetry
        depth=nodes.values
        #idxs=np.where(depth > 0.11)
        #depth[idxs]=0.11
        #compute zcor
        zcor=-depth[:,None]*sigma
        zcor2=zcor[include,::-1]
        idxs=np.where(zcor2 > 5000)
        zcor2[idxs]=5000.0-1.0e-6
        print(f'zcor2 at node 200 is {zcor2[199,:]}')

        nvrt=zcor.shape[1]

        loni=np.array(loni)
        lati=np.array(lati)
        nlon=loni[include] 
        nlat=lati[include]
        xi, yi = transform_ll_to_cpp(nlon, nlat)

        xmin, xmax = np.min(nlon), np.max(nlon)
        ymin, ymax = np.min(nlat), np.max(nlat)

        #xmin = xmin + 360. if xmin < 0 else xmin
        #xmax = xmax + 360. if xmax < 0 else xmax
        bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)

        #construct schism grid
        #lon2i=np.tile(nlon,[nvrt,1]).T
        #lat2i=np.tile(nlat,[nvrt,1]).T
        #bxyz=np.c_[zcor2.reshape(np.size(zcor2)),lat2i.reshape(np.size(lat2i)),lon2i.reshape(np.size(lon2i))]
        x2i=np.tile(xi,[nvrt,1]).T
        y2i=np.tile(yi,[nvrt,1]).T
        bxyz=np.c_[zcor2.reshape(np.size(zcor2)),y2i.reshape(np.size(y2i)),x2i.reshape(np.size(x2i))]
        print('Computing SCHISM zcor is done!')

        #Get hycom data
        if self.start_date >= datetime(2018, 12, 4):
            database = f'GLBy0.08/expt_93.0'
            xmin = xmin + 360. if xmin < 0 else xmin
            xmax = xmax + 360. if xmax < 0 else xmax
            bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)
        elif self.start_date >= datetime(2018, 1, 1) and self.start_date < datetime(2020, 2, 18):
            database = f'GLBv0.08/expt_93.0'
            xmin = xmin + 360. if xmin < 0 else xmin
            xmax = xmax + 360. if xmax < 0 else xmax
            bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)
        elif self.start_date >= datetime(2017, 10, 1) and self.start_date < datetime(2017, 12, 31):
            database = f'GLBv0.08/expt_92.9'
        elif self.start_date >= datetime(2017, 6, 1) and self.start_date < datetime(2017, 9, 30):
            database = f'GLBv0.08/expt_57.7'
        elif self.start_date >= datetime(2017, 2, 1) and self.start_date < datetime(2017, 5, 31):
            database = f'GLBv0.08/expt_92.8'
        elif self.start_date >= datetime(2016, 5, 1) and self.start_date < datetime(2017, 1, 31):
            database = f'GLBv0.08/expt_57.2'
        elif self.start_date >= datetime(2016, 1, 1) and self.start_date < datetime(2016, 4, 30):
            database = f'GLBv0.08/expt_56.3'
        elif self.start_date >= datetime(1994, 1, 1) and self.start_date < datetime(2015, 12, 31):
            database = f'GLBv0.08/expt_53.X/data/{self.start_date.year}'
        else:
            print('No data for {self.start_date}')

        print(database)
        baseurl=f'https://tds.hycom.org/thredds/dodsC/{database}?lat[0:1:-1],lon[0:1:-1],time[0:1:-1],depth[0:1:39]'
        ds=Dataset(baseurl)
        time1=ds['time']
        times=nc.num2date(time1,units=time1.units,only_use_cftime_datetimes=False)

        lon=ds['lon'][:]
        lat=ds['lat'][:]
        dep=ds['depth'][:]
        lat_idxs=np.where((lat>=bbox.ymin-2.0)&(lat<=bbox.ymax+2.0))[0]
        lon_idxs=np.where((lon>=bbox.xmin-2.0) & (lon<=bbox.xmax+2.0))[0]
        lon=lon[lon_idxs]
        lat=lat[lat_idxs]
        #print(lon_idxs)
        #print(lat_idxs)
        lon_idx1=lon_idxs[0].item()
        lon_idx2=lon_idxs[-1].item()
        print(f'lon_idx1 is {lon_idx1}, lon_idx2 is {lon_idx2}')
        lat_idx1=lat_idxs[0].item()
        lat_idx2=lat_idxs[-1].item()
        print(f'lat_idx1 is {lat_idx1}, lat_idx2 is {lat_idx2}')
    
        for ilon in np.arange(len(lon)):
            if lon[ilon] > 180:
                lon[ilon] = lon[ilon]-360.
        x2, y2=transform_ll_to_cpp(lon, lat)

        #allocate output variables
        nNode=len(include)
        one=1
        ntimes=self.rnday+1

        timeseries_s=np.zeros([ntimes,nNode,nvrt,one])
        timeseries_t=np.zeros([ntimes,nNode,nvrt,one])
        ndt=np.zeros([ntimes])

        print('**** Accessing GOFS data*****')
        t0 = time()
        for it, date in enumerate(self.timevector):

            time_idx=np.where( date == times)[0].item()
            print(time_idx)
            url=f'https://tds.hycom.org/thredds/dodsC/{database}?lat[{lat_idx1}:1:{lat_idx2}],' + \
            f'lon[{lon_idx1}:1:{lon_idx2}],depth[0:1:39],time[{time_idx}],' + \
            f'water_temp[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
            f'salinity[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}]'
            #f'water_u[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}],' + \
            #f'water_v[{time_idx}][0:1:39][{lat_idx1}:1:{lat_idx2}][{lon_idx1}:1:{lon_idx2}]'
            print(url)

            ds=Dataset(url)
            salt=np.squeeze(ds['salinity'][:,:,:])
            temp=np.squeeze(ds['water_temp'][:,:,:])

            print(f'The shape of temp is {temp.shape}')
            nz=temp.shape[0]
            ny=temp.shape[1]
            nx=temp.shape[2]

            #Convert temp to potential temp
            dep=ds['depth'][:]
            #In the ocean at 1m depth is approximately 1 dbar
            pre=np.tile(dep, ny*nx).reshape(nz, ny, nx)
            Pr=np.zeros(temp.shape)
            ptemp=sw.ptmp(salt, temp, pre, Pr)*1.00024

            print('****Interpolation starts****')

            #change missing value to nan
            #vars={'salt', 'temp'}
            idxs = np.where(abs(salt) > 10000)
            salt[idxs]=float('nan')
            idxs = np.where(abs(ptemp) > 10000)
            ptemp[idxs]=float('nan')
            #print(f'The shape of salt is {salt.shape}')

            ndt[it]=it #*24*3600.
            #salt
            #salt_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(salt[it,:,:,:]),'linear', bounds_error=False, fill_value = float('nan'))
            salt_fd=sp.interpolate.RegularGridInterpolator((dep,y2,x2),np.squeeze(salt),'nearest', bounds_error=False, fill_value = float('nan'))
            salt_int = salt_fd(bxyz)
            idxs = np.isnan(salt_int)
            if np.sum(idxs)!=0:
                salt_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], salt_int[~idxs], bxyz[idxs,:],'nearest')
            idxs = np.isnan(salt_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for salinity!')
                sys.exit()
            salt_int = salt_int.reshape(zcor2.shape)
            timeseries_s[it,:,:,0]=salt_int[:,::-1]

            #temp
            #temp_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(potentialT[it,:,:,:]),'linear', bounds_error=False, fill_value = float('nan'))
            temp_fd=sp.interpolate.RegularGridInterpolator((dep,y2,x2),np.squeeze(ptemp),'nearest', bounds_error=False, fill_value = float('nan'))
            temp_int = temp_fd(bxyz)
            idxs = np.isnan(temp_int)
            if np.sum(idxs)!=0:
                temp_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], temp_int[~idxs], bxyz[idxs,:],'nearest')
            idxs = np.isnan(temp_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for temperature!')
                sys.exit()
            temp_int = temp_int.reshape(zcor2.shape)
            timeseries_t[it,:,:,0]=temp_int[:,::-1]

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

class InitialTS():

    def __init__(self):

        pass

    def fetch_data(
            self, 
            outdir: Union[str, os.PathLike], 
            hgrid, 
            vgrid, 
            start_date, 
            include_eluv=True):

        outdir = pathlib.Path(outdir)

        self.start_date = start_date
        date = self.start_date - timedelta(days=1)
    
        vd=Vgrid()
        sigma=vd.read_vgrid(vgrid)
        
        self.hgrid = hgrid
        #bbox = self.hgrid.get_bbox(crs='epsg:4326', output_type='bbox')
        #xmin = bbox.x0 + 360. if bbox.x0 < 0 else bbox.x0
        #xmax = bbox.x1 + 360. if bbox.x1 < 0 else bbox.x1
        #bbox = Bbox.from_extents(xmin, bbox.y0, xmax, bbox.y1)

        hgrid = self.hgrid.to_dict()
        nodes = Nodes(hgrid['nodes'])

        #get coords of SCHISM
        loni = nodes.coords[:,0]
        lati = nodes.coords[:,1]
        xi,yi = transform_ll_to_cpp(np.array(loni), np.array(lati))
        bxy = np.c_[lati, loni]

        #get bathymetry
        depth = nodes.values
        print(f'maximum depth is {np.max(depth)}')
        idxs = np.where(depth < 0.11)
        depth[idxs] = 0.11
        print(f'maximum depth is {np.max(depth)}')
        #compute zcor
        zcor = -depth[:,None]*sigma
        zcor=zcor[:,::-1]
        idxs=np.where(zcor > 5000)
        zcor[idxs]=5000.0-1.0e-6
        #print(f'zcor at 647885 is {zcor[647884,:]}')
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

        side = elements.sides
        nside = side.shape[0]
        print(f'# of side is {nside}')

        #construct schism grid
        #lon2i=np.tile(loni,[nvrt,1]).T
        #lat2i=np.tile(lati,[nvrt,1]).T
        #bxyz=np.c_[zcor.reshape(np.size(zcor)),lat2i.reshape(np.size(lat2i)),lon2i.reshape(np.size(lon2i))]
        x2i=np.tile(xi,[nvrt,1]).T
        y2i=np.tile(yi,[nvrt,1]).T
        bxyz=np.c_[zcor.reshape(np.size(zcor)),x2i.reshape(np.size(x2i)),y2i.reshape(np.size(y2i))]
        print('Computing SCHISM zcor is done!')

        #Get hycom data
        '''
        Here temperature is in-situ temperature.
        TODO: convert in-situ temperature into potential temperature.
        '''
        print('**** Accessing GOFS data*****') 
        baseurl_gofs=f'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0?lat[2200:1:3162],' + \
            f'lon[3267:1:3758],depth[0:1:39],time[6671],surf_el[6671][2200:1:3162][3267:1:3758],' + \
            f'depth[0:1:39],water_temp[6671][0:1:39][2200:1:3162][3267:1:3758],' + \
            f'salinity[6671][0:1:39][2200:1:3162][3267:1:3758],' + \
            f'water_u[6671][0:1:39][2200:1:3162][3267:1:3758],' + \
            f'water_v[6671][0:1:39][2200:1:3162][3267:1:3758]'
        ds=Dataset(f'{baseurl_gofs}')
        lon=ds['lon'][:]
        lat=ds['lat'][:]
        dep=ds['depth'][:]

        #convert lon [0, 360) to [-180, 180)
        for ilon in range(len(lon)):
            if lon[ilon] > 180:
                lon[ilon] = lon[ilon]-360.
        x2, y2=transform_ll_to_cpp(lon, lat)

        #lat_idxs=np.where((lat>=bbox.ymin)&(lat<=bbox.ymax))[0]
        #lon_idxs=np.where((lon>=bbox.xmin-2.0) & (lon<=bbox.xmax+2.0))[0]

        salt=np.squeeze(ds['salinity'][:,:,:])
        temp=np.squeeze(ds['water_temp'][:,:,:])

        #Convert temp to potential temp
        pre=np.tile(dep, ny*nx).reshape(nz, ny, nx)
        Pr=np.zeros(temp.shape)
        ptemp=sw.ptmp(salt, temp, pre, Pr)*1.00024

        if include_eluv:
            uvel=nc['water_u'][:,:,:]
            vvel=nc['water_v'][:,:,:]
            ssh=nc['surf_el'][:,:]

        #change missing value to nan
        idxs = np.where(abs(salt) > 10000)
        salt[idxs]=float('nan')

        idxs = np.where(abs(ptemp) > 10000)
        ptemp[idxs]=float('nan')

        if include_eluv:
            idxs = np.where(abs(uvel) > 10000)
            uvel[idxs]=float('nan')
            idxs = np.where(abs(vvel) > 10000)
            vvel[idxs]=float('nan')
            idxs = np.where(abs(ssh) > 10000)
            ssh[idxs]=float('nan')

        t0 = time()
        print('Interpolate salinity')
        #salt_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),salt,'linear', bounds_error=False, fill_value = float('nan'))
        salt_fd=sp.interpolate.RegularGridInterpolator((dep,y2,x2),salt,'linear', bounds_error=False, fill_value = float('nan'))
        print('Done interpolator')
        salt_int = salt_fd(bxyz)
        print('Done interpolation on schism grid')
        idxs = np.isnan(salt_int)
        print(np.sum(idxs))
        if np.sum(idxs)!=0:
            #idxs2=np.where(~np.isnan(salt_int))
            salt_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], salt_int[~idxs], bxyz[idxs,:],'linear')
            #salt_fd=sp.interpolate.NearestNDInterpolator(bxyz[idxs2], salt_int[idxs2])
            #salt_int[idxs]=salt_fd(bxyz[idxs])
        print('Done interpolation on missing values')
        idxs = np.isnan(salt_int)
        if np.sum(idxs)!=0:
            print(f'There is still missing value for salinity!')
            sys.exit()
        salt_int = salt_int.reshape(zcor.shape)
        #print(salt_int.shape)

        #temperature
        print('Interpolate temperature')
        #temp_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon), temp, 'linear', bounds_error=False, fill_value=float('nan'))
        temp_fd=sp.interpolate.RegularGridInterpolator((dep,y2,x2), ptemp, 'linear', bounds_error=False, fill_value=float('nan'))
        temp_int = temp_fd(bxyz)
        idxs = np.isnan(temp_int)
        if np.sum(idxs)!=0:
            temp_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], temp_int[~idxs], bxyz[idxs,:],'linear')
        idxs = np.isnan(temp_int)
        if np.sum(idxs)!=0:
            print(f'There is still missing value for temperature!')
            sys.exit()
        temp_int = temp_int.reshape(zcor.shape)
        #print(temp_int.shape)

        if include_eluv:
            #uvel
            print('Interpolate uvel')
            #uvel_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),uvel,'linear', bounds_error=False, fill_value = float('nan'))
            uvel_fd=sp.interpolate.RegularGridInterpolator((dep,y2,x2),uvel,'linear', bounds_error=False, fill_value = float('nan'))
            uvel_int = uvel_fd(bxyz)
            idxs = np.isnan(uvel_int)
            if np.sum(idxs)!=0:
                uvel_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], uvel_int[~idxs], bxyz[idxs,:],'linear')
            idxs = np.isnan(uvel_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for u-velocity!')
                sys.exit()
            uvel_int = uvel_int.reshape(zcor.shape)
            print(f'uvel_int size is {uvel_int.shape}')

            #vvel
            print('Interpolate vvel')
            #vvel_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),vvel,'linear', bounds_error=False, fill_value = float('nan'))
            vvel_fd=sp.interpolate.RegularGridInterpolator((dep,y2,x2),vvel,'linear', bounds_error=False, fill_value = float('nan'))
            vvel_int = vvel_fd(bxyz)
            idxs = np.isnan(vvel_int)
            if np.sum(idxs)!=0:
                vvel_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], vvel_int[~idxs], bxyz[idxs,:],'linear')
            idxs = np.isnan(vvel_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for v-velocity!')
                sys.exit()
            vvel_int = vvel_int.reshape(zcor.shape)

            #ssh
            print('Interpolate ssh')
            ssh_fd=sp.interpolate.RegularGridInterpolator((lat,lon),ssh,'linear', bounds_error=False, fill_value = float('nan'))
            ssh_int = ssh_fd(bxy)
            idxs = np.isnan(ssh_int)
            if np.sum(idxs)!=0:
                ssh_int[idxs]=sp.interpolate.griddata(bxy[~idxs,:], ssh_int[~idxs], bxy[idxs,:],'linear')

            idxs = np.isnan(ssh_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for ssh!')
                sys.exit()
            print(f'ssh_int size is {ssh_int.shape}')

        print(f'Interpolation takes {time()-t0} seconds')

        #Create hotstart.nc
        t0=time()
        NP = len(loni)
        NE = len(eles)
        ntracers = 2
        one = 1

        #Compute tr_nd and tr_el
        tr_nd = np.zeros([NP, nvrt, 2])
        tr_nd[:,:,0] = temp_int[:,::-1]
        tr_nd[:,:,1] = salt_int[:,::-1]
        tr_el = np.zeros([NE, nvrt, 2])

        for k in np.arange(nvrt):
            #triangles
            tr_el[tidxs,k,0] = np.mean(temp_int[tris[:,:],k], axis=1)
            tr_el[tidxs,k,1] = np.mean(salt_int[tris[:,:],k], axis=1)
            #quads 
            tr_el[qidxs,k,0] = np.mean(temp_int[quads[:,:],k], axis=1)
            tr_el[qidxs,k,1] = np.mean(salt_int[quads[:,:],k], axis=1)

        #Compute su2 and sv2
        su2 = np.zeros([nside, nvrt])
        sv2 = np.zeros([nside, nvrt])
        if include_eluv:
            for iside in np.arange(nside):
                id1=side[iside,0]
                id2=side[iside,1]
                su2[iside,:]=(uvel_int[id1,::-1]+uvel_int[id2,::-1])/2.0
                sv2[iside,:]=(vvel_int[id1,::-1]+vvel_int[id2,::-1])/2.0
        print(f'su2 size is {su2.shape}')

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

            dst.createVariable('nsteps_from_cold', 'i4', ('one',))
            dst['nsteps_from_cold'][:] = 0

            dst.createVariable('idry_e', 'i4', ('elem',))
            dst['idry_e'][:] = np.zeros(NE).astype('int32')

            dst.createVariable('idry_s', 'i4', ('side',))
            dst['idry_s'][:] = np.zeros(nside).astype('int32')

            dst.createVariable('idry', 'i4', ('node',))
            dst['idry'][:] = np.zeros(NP).astype('int32')

            dst.createVariable('eta2', 'd', ('node',))
            if include_eluv:
                dst['eta2'][:] = ssh_int
            else:
                dst['eta2'][:] = np.zeros(NP)

            dst.createVariable('cumsum_eta', 'd', ('node',))
            if include_eluv:
                dst['cumsum_eta'][:] = ssh_int
            else:
                dst['cumsum_eta'][:] = np.zeros(NP)
    
            dst.createVariable('we', 'd', ('elem', 'nVert'))
            dst['we'][:,:] = np.zeros([NE,nvrt])

            dst.createVariable('tr_el', 'd', ('elem', 'nVert', 'ntracers'))
            dst['tr_el'][:,:,:] = tr_el

            dst.createVariable('tr_nd', 'd', ('node', 'nVert', 'ntracers'))
            dst['tr_nd'][:,:,:] = tr_nd

            dst.createVariable('tr_nd0', 'd', ('node', 'nVert', 'ntracers'))
            dst['tr_nd0'][:,:,:] = tr_nd

            dst.createVariable('su2', 'd', ('side', 'nVert'))
            dst['su2'][:,:] = su2 #np.zeros([nside,nvrt])

            dst.createVariable('sv2', 'd', ('side', 'nVert'))
            dst['sv2'][:,:] = sv2 #np.zeros([nside,nvrt])
 
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
        self.timevector=np.arange(
            self.start_date,
            self.start_date + timedelta(days=self.rnday+1),
            timedelta(days=1))

        vd=Vgrid()
        sigma=vd.read_vgrid(vgrid)

        hgrid = hgrid.to_dict()
        nodes = Nodes(hgrid['nodes']) 
        #get coords of SCHISM
        loni = nodes.coords[:,0]
        lati = nodes.coords[:,1]

        #get bathymetry
        depth = nodes.values
        #idxs = np.where(depth > 0.11)
        #depth[idxs] = 0.11
        #compute zcor
        zcor = -depth[:,None]*sigma
        nvrt=zcor.shape[1]
        #print(f'zcor at node 1098677 is {zcor[1098676,:]}')

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
        xi,yi = transform_ll_to_cpp(np.array(blon), np.array(blat))
        #bxy = np.c_[xi, yi]
        #blon=np.array(blon)*1.e5
        #blat=np.array(blat)*1.e5
        bxy = np.c_[blat, blon]

        zcor2=zcor[idxs,::-1]
        idxs=np.where(zcor2 > 5000)
        #print(idxs)
        zcor2[idxs]=5000.0-1.0e-6
        #print(f'zcor2 at node 200 is {zcor2[199,:]}')

        xmin, xmax = np.min(blon), np.max(blon)
        ymin, ymax = np.min(blat), np.max(blat)

        #xmin = xmin + 360. if xmin < 0 else xmin
        #xmax = xmax + 360. if xmax < 0 else xmax
        bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)

        #construct schism grid
        #lon2i=np.tile(blon,[nvrt,1]).T
        #lat2i=np.tile(blat,[nvrt,1]).T
        #bxyz=np.c_[zcor2.reshape(np.size(zcor2)),lat2i.reshape(np.size(lat2i)),lon2i.reshape(np.size(lon2i))]
        x2i=np.tile(xi,[nvrt,1]).T
        y2i=np.tile(yi,[nvrt,1]).T
        bxyz=np.c_[zcor2.reshape(np.size(zcor2)),y2i.reshape(np.size(y2i)),x2i.reshape(np.size(x2i))]
        print('Computing SCHISM zcor is done!')

        #Get hycom data
        if self.start_date >= datetime(2018, 12, 4):
            database = f'GLBy0.08/expt_93.0'
            xmin = xmin + 360. if xmin < 0 else xmin
            xmax = xmax + 360. if xmax < 0 else xmax
            bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)
        elif self.start_date >= datetime(2018, 1, 1) and self.start_date < datetime(2020, 2, 18):
            database = f'GLBv0.08/expt_93.0'
            xmin = xmin + 360. if xmin < 0 else xmin
            xmax = xmax + 360. if xmax < 0 else xmax
            bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)
        elif self.start_date >= datetime(2017, 10, 1) and self.start_date < datetime(2017, 12, 31):
            database = f'GLBv0.08/expt_92.9'
        elif self.start_date >= datetime(2017, 6, 1) and self.start_date < datetime(2017, 9, 30):
            database = f'GLBv0.08/expt_57.7'
        elif self.start_date >= datetime(2017, 2, 1) and self.start_date < datetime(2017, 5, 31):
            database = f'GLBv0.08/expt_92.8'
        elif self.start_date >= datetime(2016, 5, 1) and self.start_date < datetime(2017, 1, 31):
            database = f'GLBv0.08/expt_57.2'
        elif self.start_date >= datetime(2016, 1, 1) and self.start_date < datetime(2016, 4, 30):
            database = f'GLBv0.08/expt_56.3'
        elif self.start_date >= datetime(1994, 1, 1) and self.start_date < datetime(2015, 12, 31):
            database = f'GLBv0.08/expt_53.X/data/{self.start_date.year}'
        else:
            print('No data for {self.start_date}')

        #baseurl=f'https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0?lat[0:1:3250],lon[0:1:4499],time[0:1:6127],depth[0:1:39]'
        print(database)
        baseurl=f'https://tds.hycom.org/thredds/dodsC/{database}?lat[0:1:-1],lon[0:1:-1],time[0:1:-1],depth[0:1:39]'
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
        print(lat_idxs)
        lon_idx1=lon_idxs[0].item()
        lon_idx2=lon_idxs[-1].item()
        print(f'lon_idx1 is {lon_idx1}, lon_idx2 is {lon_idx2}')
        lat_idx1=lat_idxs[0].item()
        lat_idx2=lat_idxs[-1].item()
        print(f'lat_idx1 is {lat_idx1}, lat_idx2 is {lat_idx2}')

        for ilon in np.arange(len(lon)):
            if lon[ilon] > 180:
                lon[ilon] = lon[ilon]-360.
        x2, y2=transform_ll_to_cpp(lon, lat)

        #baseurl_gofs=f'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0?lat[2209:1:3156],' + \
        #    f'lon[3270:1:3758],time[6671:8:7786],surf_el[6671:8:7786][2209:1:3156][3270:1:3758],' + \
        #    f'depth[0:1:39],water_temp[6671:8:7786][0:1:39][2209:1:3156][3270:1:3758],' + \
        #    f'salinity[6671:8:7786][0:1:39][2209:1:3156][3270:1:3758],' + \
        #    f'water_u[6671:8:7786][0:1:39][2209:1:3156][3270:1:3758],' + \
        #    f'water_v[6671:8:7786][0:1:39][2209:1:3156][3270:1:3758]'

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

            time_idx=np.where( date == times)[0].item()
            print(time_idx)
            url=f'https://tds.hycom.org/thredds/dodsC/{database}?lat[{lat_idx1}:1:{lat_idx2}],' + \
            f'lon[{lon_idx1}:1:{lon_idx2}],depth[0:1:39],time[{time_idx}],' + \
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
            print(f'The shape of temp is {temp.shape}')

            #Convert temp to potential temp
            nz=temp.shape[0]
            ny=temp.shape[1]
            nx=temp.shape[2]
            dep=ds['depth'][:]
            pre=np.tile(dep, ny*nx).reshape(nz, ny, nx)
            Pr=np.zeros(temp.shape)
            ptemp=sw.ptmp(salt, temp, pre, Pr)*1.00024

            print('****Interpolation starts****')

            #change missing value to nan
            #vars={'salt', 'temp', 'uvel', 'vvel', 'ssh'}
            idxs = np.where(abs(salt) > 10000)
            salt[idxs]=float('nan')
            #idxs = np.where(abs(temp) > 10000)
            #temp[idxs]=float('nan')
            idxs = np.where(abs(ptemp) > 10000)
            ptemp[idxs]=float('nan')
            idxs = np.where(abs(uvel) > 10000)
            uvel[idxs]=float('nan')
            idxs = np.where(abs(vvel) > 10000)
            vvel[idxs]=float('nan')
            idxs = np.where(abs(ssh) > 10000)
            ssh[idxs]=float('nan')
            #print(f'The shape of salt is {salt.shape}')

            ndt[it]=it*24*3600.
            #salt
            #salt_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(salt),'nearest', bounds_error=False, fill_value = float('nan'))
            #salt_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(salt),'linear', bounds_error=False)
            salt_fd=sp.interpolate.RegularGridInterpolator((dep,y2,x2),np.squeeze(salt),'linear', bounds_error=False)
            salt_int = salt_fd(bxyz)
            idxs = np.isnan(salt_int)
            if np.sum(idxs)!=0:
                salt_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], salt_int[~idxs], bxyz[idxs,:],'nearest')
            idxs = np.isnan(salt_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for salinity!')
                sys.exit()
            salt_int = salt_int.reshape(zcor2.shape)
            timeseries_s[it,:,:,0]=salt_int[:,::-1]

            #temp
            #temp_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(potentialT),'linear', bounds_error=False, fill_value = float('nan'))
            #temp_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(temp),'linear', bounds_error=False)
            temp_fd=sp.interpolate.RegularGridInterpolator((dep,y2,x2),np.squeeze(ptemp),'linear', bounds_error=False)
            temp_int = temp_fd(bxyz)
            idxs = np.isnan(temp_int)
            if np.sum(idxs)!=0:
                temp_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], temp_int[~idxs], bxyz[idxs,:],'nearest')
            idxs = np.isnan(temp_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for temperature!')
                sys.exit()
            temp_int = temp_int.reshape(zcor2.shape)
            timeseries_t[it,:,:,0]=temp_int[:,::-1]

            #uvel
            #uvel_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(uvel),'linear', bounds_error=False, fill_value = float('nan'))
            uvel_fd=sp.interpolate.RegularGridInterpolator((dep,y2,x2),np.squeeze(uvel),'linear', bounds_error=False, fill_value = float('nan'))
            #uvel_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(uvel),'linear', bounds_error=False)
            uvel_int = uvel_fd(bxyz)
           
            idxs = np.isnan(uvel_int)
            if np.sum(idxs)!=0:
                uvel_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], uvel_int[~idxs], bxyz[idxs,:],'nearest')
            idxs = np.isnan(uvel_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for u-velocity!')
                sys.exit()
            uvel_int = uvel_int.reshape(zcor2.shape)
            timeseries_uv[it,:,:,0]=uvel_int[:,::-1]

            #vvel
            #vvel_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(vvel),'linear', bounds_error=False, fill_value = float('nan'))
            vvel_fd=sp.interpolate.RegularGridInterpolator((dep,y2,x2),np.squeeze(vvel),'linear', bounds_error=False, fill_value = float('nan'))
            #vvel_fd=sp.interpolate.RegularGridInterpolator((dep,lat,lon),np.squeeze(vvel),'linear', bounds_error=False)
            vvel_int = vvel_fd(bxyz)
            idxs = np.isnan(vvel_int)
            if np.sum(idxs)!=0:
                vvel_int[idxs]=sp.interpolate.griddata(bxyz[~idxs,:], vvel_int[~idxs], bxyz[idxs,:],'nearest')
            idxs = np.isnan(vvel_int)
            if np.sum(idxs)!=0:
                print(f'There is still missing value for v-velocity!')
                sys.exit()
            vvel_int = vvel_int.reshape(zcor2.shape)
            timeseries_uv[it,:,:,1]=vvel_int[:,::-1]

            #ssh
            ssh_fd=sp.interpolate.RegularGridInterpolator((lat,lon),np.squeeze(ssh),'linear', bounds_error=False, fill_value = float('nan'))
            #ssh_fd=sp.interpolate.RegularGridInterpolator((y2,x2),np.squeeze(ssh),'linear', bounds_error=False, fill_value = float('nan'))
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
