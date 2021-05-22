import os
import numpy as np
from datetime import datetime,timedelta
import logging
import pathlib
import tempfile
import subprocess
from typing import Union
from numba import jit, prange


from netCDF4 import Dataset
from matplotlib.transforms import Bbox


logger = logging.getLogger(__name__)


class Nudge:

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
                        rl2 = np.sqrt(
                            np.square(lon[idn] - lon[j-1]) +
                            np.square(lat[idn] - lat[j-1]))
                        if rl2 < distmin:
                            distmin = rl2
                rnu = 0.
                if distmin <= rlmax:
                    rnu = (1-distmin/rlmax)*rnu_max
                out[idn] = rnu

        # self.hgrid = hgrid
        outdir = pathlib.Path(outdir)

        hgrid = hgrid.to_dict()
        nodes = hgrid['nodes']
        elements = hgrid['elements']
        _, NP = len(elements), len(nodes)
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

        # Max relax distance in degr
        rlmax = 1.5

        # Max relax strength in days
        rnu_day = 0.25

        # rnu = 0
        rnu_max = 1./rnu_day/86400.

        out = np.zeros([NP])
        # t0 = time()
        compute_nudge(lon, lat, NP, opbd2, out)
        # print(f'It took {time() -t0}')

        # nudge = [f"{rlmax}, {rnu_day}"]
        # nudge.extend("\n")
        # nudge.append(f"{NE} {NP}")
        # nudge.extend("\n")
        # for idn, (coords, values) in nodes.items():
        #     line = [f"{idn}"]
        #     line.extend([f"{x:<.7e}" for x in coords])
        #     line.extend([f"{out[int(idn)-1]:<.7e}"])
        #     line.extend("\n")
        #     nudge.append(" ".join(line))

        # for id, element in elements.items():
        #     line = [f"{id}"]
        #     line.append(f"{len(element)}")
        #     line.extend([f"{e}" for e in element])
        #     line.extend("\n")
        #     nudge.append(" ".join(line))

        # with open(outdir / 'nudge_pyschism.gr3','w+') as fid:
        #     fid.writelines(nudge)


class HotStartInventory():

    def fetch_data(self, outdir: Union[str, os.PathLike], hgrid, start_date):

        logger.info('Fetching RTOFS data')

        self.hgrid = hgrid
        self.start_date = start_date

        outdir = pathlib.Path(outdir)
        if outdir.name != '':
            outdir /= 'hotstart'
        outdir.mkdir(exist_ok=True)

        # _tmpdir = tempfile.TemporaryDirectory()
        # tmpdir = pathlib.Path(_tmpdir.name)

        bbox = self.hgrid.get_bbox(crs='epsg:4326', output_type='bbox')

        # Go to yesterday's directory to get tomorrow's data
        dt = self.start_date - timedelta(days=1)
        nc_ssh = Dataset(
            'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
            + dt.strftime('%Y%m%d')+'/rtofs_glo_2ds_forecast_3hrly_diag')
        nc_salt = Dataset(
            'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
            + dt.strftime('%Y%m%d')+'/rtofs_glo_3dz_forecast_daily_salt')
        nc_temp = Dataset(
            'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
            + dt.strftime('%Y%m%d')+'/rtofs_glo_3dz_forecast_daily_temp')
        nc_uvel = Dataset(
            'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
            + dt.strftime('%Y%m%d')+'/rtofs_glo_3dz_forecast_daily_uvel')
        nc_vvel = Dataset(
            'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
            + dt.strftime('%Y%m%d')+'/rtofs_glo_3dz_forecast_daily_vvel')

        lon = nc_ssh['lon'][:]
        lat = nc_ssh['lat'][:]
        lev = nc_salt['lev'][:]

        xmin = bbox.x0 + 360. if bbox.x0 < 0 else bbox.x0
        xmax = bbox.x1 + 360. if bbox.x1 < 0 else bbox.x1
        bbox = Bbox.from_extents(xmin, bbox.y0, xmax, bbox.y1)

        lat_idxs = np.where((lat >= bbox.ymin) & (lat <= bbox.ymax))[0]
        lon_idxs = np.where((lon >= bbox.xmin) & (lon <= bbox.xmax))[0]

        xlon = lon[lon_idxs]
        for ilon in range(lon_idxs.size):
            if xlon[ilon] > 180.:
                xlon[ilon] = xlon[ilon] - 360.
        ylat = lat[lat_idxs]

        with Dataset(outdir / 'SSH_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            # time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:2]
            # ssh
            dst.createVariable('surf_el', 'f4', ('time', 'ylat', 'xlon',),
                               fill_value=-30000.0)
            dst['surf_el'].long_name = "sea_surface_elevation (m)"
            dst['surf_el'].add_offset = 0.
            dst['surf_el'].scale_factor = 0.001
            dst['surf_el'][:, :, :] = nc_ssh['ssh'][8:9, 0, lat_idxs, lon_idxs]

        with Dataset(outdir / 'TS_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            # dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            # lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            # time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:2]
            # salt
            dst.createVariable(
                'salinity', 'f4', ('time', 'lev', 'ylat', 'xlon',),
                fill_value=-30000.0)
            dst['salinity'].long_name = "sea_water_salinity (psu)"
            dst['salinity'].add_offset = 20.
            dst['salinity'].scale_factor = 0.001
            # temp
            dst.createVariable(
                'temperature', 'f4', ('time', 'lev', 'ylat', 'xlon',),
                fill_value=-30000.0)
            dst['temperature'].long_name = "sea_water_potential_temperature (degc)"
            dst['temperature'].add_offset = 20.
            dst['temperature'].scale_factor = 0.001

            for k in np.arange(len(lev)):
                dst['salinity'][:, k, :, :] = nc_salt[
                    'salinity'][1:2, k, lat_idxs, lon_idxs]
                dst['temperature'][:, k, :, :] = nc_temp[
                    'temperature'][1:2, k, lat_idxs, lon_idxs]

        with Dataset(outdir / 'UV_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            # dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            # lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            # time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:2]
            # uvel
            dst.createVariable(
                'water_u', 'f4', ('time', 'lev', 'ylat', 'xlon',),
                fill_value=-30000.0)
            dst['water_u'].long_name = "eastward_sea_water_velocity (m/s)"
            dst['water_u'].add_offset = 0.
            dst['water_u'].scale_factor = 0.001
            # vvel
            dst.createVariable(
                'water_v', 'f4', ('time', 'lev', 'ylat', 'xlon',),
                fill_value=-30000.0)
            dst['water_v'].long_name = "northward_sea_water_velocity (m/s)"
            dst['water_v'].add_offset = 0.
            dst['water_v'].scale_factor = 0.001

            for k in np.arange(len(lev)):
                dst['water_u'][:, k, :, :] = nc_uvel['u'][
                    1:2, k, lat_idxs, lon_idxs]
                dst['water_v'][:, k, :, :] = nc_vvel['v'][
                    1:2, k, lat_idxs, lon_idxs]

        # symlink estuary.gr3 and *.in file
        hgrid.write(outdir / 'hgrid.gr3')
        hgrid.write(outdir / 'hgrid.ll')
        src = '/sciclone/data10/lcui01/ICOGS3D_dep/estuary.gr3'
        dst = outdir / 'estuary.gr3'
        os.symlink(src, dst)
        src2 = '/sciclone/data10/lcui01/ICOGS3D_dep/gen_hot_3Dth_from_nc.in'
        dst2 = outdir / 'gen_hot_3Dth_from_nc.in'
        os.symlink(src2, dst2)
        src3 = '/sciclone/data10/lcui01/ICOGS3D_dep/vgrid.in'
        dst3 = outdir / 'vgrid.in'
        os.symlink(src3, dst3)
        cmd = 'gen_hot_3Dth_from_hycom.exe < gen_hot_3Dth_from_nc.in'
        subprocess.check_call(cmd, shell=True, cwd=outdir)
        # shutil.move(outdir / 'hotstart.nc')


class OpenBoundaryInventory():

    def fetch_data(self, outdir: Union[str, os.PathLike], hgrid, start_date, rnday, idx_min, idx_max, jdx_min, jdx_max):

        self.start_date = start_date
        self.rnday = rnday

        outdir = pathlib.Path(outdir)
        if outdir.name != '':
            outdir /= '3Dth'
        outdir.mkdir(exist_ok=True)

        _tmpdir = tempfile.TemporaryDirectory()
        tmpdir = pathlib.Path(_tmpdir.name)

        #Go to yesterday's directory to get tomorrow's data
        dt = self.start_date - timedelta(days=1)
        nc_ssh = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_2ds_forecast_3hrly_diag')
        nc_salt = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_salt')
        nc_temp = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_temp')
        nc_uvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_uvel')
        nc_vvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global' \
            + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_vvel')

        lon = nc_ssh['lon'][:]
        lat = nc_ssh['lat'][:]
        lev = nc_salt['lev'][:]

        xlon = lon[idx_min:idx_max+1]
        for ilon in range(xlon.size):
            if xlon[ilon] > 180.:
                xlon[ilon] = xlon[ilon] - 360.
        ylat = lat[jdx_min:jdx_max+1]

        #xlon = lon[idx_min:idx_max+1]
        #for ilon in range(xlon.size):
        #    if xlon[ilon] > 180.:
        #        xlon[ilon] = xlon[ilon] - 360.
        #ylat = lat[jdx_min:jdx_max+1]

        with Dataset(outdir / 'SSH_1.nc', 'w', format='NETCDF3_64BIT_OFFSET') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat 
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:rnday+2]
            #ssh
            dst.createVariable('surf_el', 'f4', ('time', 'ylat', 'xlon',), fill_value=-30000.0)
            dst['surf_el'].long_name = "sea_surface_elevation (m)"
            dst['surf_el'].add_offset = 0.
            dst['surf_el'].scale_factor = 0.001
            ssh = nc_ssh['ssh'][8:8*(rnday+1)+1:8,0,jdx_min:jdx_max+1,idx_min:idx_max+1]
            dst['surf_el'][:,:,:] = ssh

        with Dataset(outdir / 'TS_1.nc', 'w', format='NETCDF3_64BIT_OFFSET') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat 
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            #lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:rnday+2]
            #salt 
            dst.createVariable('salinity', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.)
            dst['salinity'].long_name = "sea_water_salinity (psu)"
            dst['salinity'].add_offset = 20.
            dst['salinity'].scale_factor = 0.001
            sss = nc_salt['salinity'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]
            dst['salinity'][:,:,:,:] = sss
            #temp
            dst.createVariable('temperature', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.)
            dst['temperature'].long_name = "sea_water_potential_temperature (degc)"
            dst['temperature'].add_offset = 20.
            dst['temperature'].scale_factor = 0.001
            sst = nc_temp['temperature'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]
            dst['temperature'][:,:,:,:] = sst

        with Dataset(outdir / 'UV_1.nc', 'w', format='NETCDF3_64BIT_OFFSET') as dst:
            dst.setncatts({"Conventions": "cf-1.0"})
            #dimensions
            dst.createDimension('xlon', xlon.shape[0])
            dst.createDimension('ylat', ylat.shape[0])
            dst.createDimension('lev', nc_salt['lev'].shape[0])
            dst.createDimension('time', None)

            # variables
            # lon
            dst.createVariable('xlon', 'f4', ('xlon',))
            dst['xlon'].long_name = "Longitude"
            dst['xlon'].standard_name = "longitude"
            dst['xlon'].units = "degrees_east"
            dst['xlon'][:] = xlon
            # lat 
            dst.createVariable('ylat', 'f4', ('ylat',))
            dst['ylat'].long_name = "Latitude"
            dst['ylat'].standard_name = "latitude"
            dst['ylat'].units = "degrees_north"
            dst['ylat'][:] = ylat
            #lev
            dst.createVariable('lev', 'f4', ('lev',))
            dst['lev'].long_name = "altitude"
            dst['lev'].units = "millibar"
            dst['lev'][:] = nc_salt['lev'][:]
            #time
            dst.createVariable('time', 'f4', ('time',))
            dst['time'].standard_name = "time"
            dst['time'].units = "days since 1-1-1 00:00:0.0"
            dst['time'][:] = nc_salt['time'][1:rnday+2]
            #uvel
            dst.createVariable('water_u', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.)
            dst['water_u'].long_name = "eastward_sea_water_velocity (m/s)"
            dst['water_u'].add_offset = 0.
            dst['water_u'].scale_factor = 0.001
            uvel = nc_uvel['u'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]
            dst['water_u'][:,:,:,:] = uvel
            #vvel
            dst.createVariable('water_v', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.)
            dst['water_v'].long_name = "northward_sea_water_velocity (m/s)"
            dst['water_v'].add_offset = 0.
            dst['water_v'].scale_factor = 0.001
            vvel = nc_vvel['v'][1:rnday+2,:,jdx_min:jdx_max+1,idx_min:idx_max+1]
            dst['water_v'][:,:,:,:] = vvel

        #symlink estaury.gr3 and *.in file 
        hgrid.write(outdir / 'hgrid.gr3')
        hgrid.write(outdir / 'hgrid.ll')
        src = '/sciclone/data10/lcui01/ICOGS3D_dep/estuary.gr3'
        dst = outdir / 'estuary.gr3'
        os.symlink(src, dst)
        src2 = '/sciclone/data10/lcui01/ICOGS3D_dep/gen_hot_3Dth_from_nc.in'
        dst2 = outdir / 'gen_hot_3Dth_from_nc.in'
        os.symlink(src2, dst2)
        src3 = '/sciclone/data10/lcui01/ICOGS3D_dep/vgrid.in'
        dst3 = outdir / 'vgrid.in'
        os.symlink(src3, dst3)
        cmd = 'gen_hot_3Dth_from_hycom.exe < gen_hot_3Dth_from_nc.in'
        subprocess.check_call(cmd, shell=True, cwd=outdir)
        #shutil.move(tmpdir /'elev2D.th.nc', outdir / 'elev2D.th.nc')
        #shutil.move(tmpdir /'SAL_3D.th.nc', outdir / 'SAL_3D.th.nc')
        #shutil.move(tmpdir /'TEM_3D.th.nc', outdir / 'TEM_3D.th.nc')
        #shutil.move(tmpdir /'uv3D.th.nc', outdir / 'uv3D.th.nc')

        #Generate nudging input
        subprocess.check_call(['gen_nudge2'], cwd=outdir)
        subprocess.check_call(['gen_nudge_from_hycom'], cwd=tmpdir)
        #shutil.move(tmpdir /'SAL_nu.nc', outdir)
        #shutil.move(tmpdir /'TEM_nu.nc', outdir)
