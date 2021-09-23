# from datetime import datetime, timedelta
# import logging
# from typing import Dict
# import warnings

# from matplotlib.transforms import Bbox
# import numpy as np
# from scipy.interpolate import griddata
# import xarray


# import tqdm
# import tqdm_logging_wrapper

# from pyschism import dates
# # from pyschism.forcing.baroclinic.base import BaroclinicForcing, BaroclinicComponent
# from pyschism.mesh.base import Gr3

# logger = logging.getLogger(__name__)


# class RTOFSBaroclinicComponent:

#     def interpolate(self, gr3: Gr3, time: datetime, lev_index: int = 0):
#         bbox = self._modified_bbox(
#             gr3.get_bbox(crs='epsg:4326', output_type='bbox'))

#         zi = None
#         # for pivot_date, dataset in self.datasets.items():
#         #     check_date = np.datetime64(dates.nearest_zulu(time).replace(
#         #                 tzinfo=None))
#         #     idx = np.where(
#         #         ~np.isnan(dataset.time.where(dataset.time == check_date)))[0]
#         #     if len(idx) > 0:
#         #         if idx[0] == 0:
#         #             continue
#         #         else:
#         #             logging.info(
#         #                 f'Downloading RTOFS data for variable {self.ncvar} '
#         #                 f'for requested date {time}, rounded as {check_date}.')
#         #             lon_idxs, lat_idxs = self._modified_bbox_indexes(bbox)
#         #             zi = dataset[self.ncvar][idx[0], lev_index,
#         #                                      lat_idxs, lon_idxs].values
#         # if zi is None:
#         for pivot_date, dataset in self.datasets.items():
#             check_date = np.datetime64(dates.nearest_zulu(time).replace(
#                     tzinfo=None))
#             idx = np.where(
#                 ~np.isnan(dataset.time.where(dataset.time == check_date)))[0]
#             if len(idx) > 0:
#                 if idx[0] == 0:
#                     continue
#                 else:
#                     logging.info(
#                         f'Downloading RTOFS data for variable {self.ncvar} '
#                         f'for requested date {time}, rounded as {check_date}.')

#                     # items = [1, 2, 3]
#                     # items_iter = tqdm.tqdm(items)
#                     # logger.info(f"Items: {items}")
#                     # with tqdm_logging_wrapper.wrap_logging_for_tqdm(items_iter), items_iter:
#                     #     for item in items_iter:
#                     #         logger.info(f"Item: {item}")

#                     lon_idxs, lat_idxs = self._modified_bbox_indexes(
#                         bbox,
#                         pixel_buffer=100
#                     )
#                     items_iter = tqdm.tqdm(lat_idxs)
#                     zi0 = np.full((len(lat_idxs), len(lon_idxs)), np.nan)
#                     with tqdm_logging_wrapper.wrap_logging_for_tqdm(
#                             items_iter), items_iter:
#                         for i, row in enumerate(items_iter):
#                             zi0[i, :] = dataset[self.ncvar][
#                                 idx[0], lev_index, row, lon_idxs].values

#                     lon_idxs, lat_idxs = self._modified_bbox_indexes(
#                         bbox,
#                         # pixel_buffer=100
#                     )
#                     items_iter = tqdm.tqdm(lat_idxs)
#                     zi1 = np.full((len(lat_idxs), len(lon_idxs)), np.nan)
#                     with tqdm_logging_wrapper.wrap_logging_for_tqdm(
#                             items_iter), items_iter:
#                         for i, row in enumerate(items_iter):
#                             zi1[i, :] = dataset[self.ncvar][
#                                 idx[0], lev_index, row, lon_idxs].values

#                     breakpoint()

#         if zi is None:
#             raise ValueError(f'No RTOFS data for requested date {time}.')

#         xi, yi = self._xy_grid(lon_idxs, lat_idxs)
#         xi, yi, zi = xi.flatten(), yi.flatten(), zi.flatten()
#         rtofs_idxs = np.where(~np.isnan(zi))
#         values = griddata(
#             (xi[rtofs_idxs], yi[rtofs_idxs]),
#             zi[rtofs_idxs],
#             (gr3.x, gr3.y),
#             method='linear',
#             fill_value=self.fill_value
#         )

#         nan_idxs = np.where(np.isnan(values))

#         values[nan_idxs] = griddata(
#             (xi[rtofs_idxs], yi[rtofs_idxs]),
#             zi[rtofs_idxs],
#             (gr3.x[nan_idxs], gr3.y[nan_idxs]),
#             method='nearest',
#         )

#         return values

#     @property
#     def base_url(self) -> str:
#         return 'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'

#     @property
#     def output_interval(self):
#         return timedelta(days=1)

#     @property
#     def datasets(self) -> Dict[datetime, xarray.Dataset]:
#         if not hasattr(self, '_datasets'):
#             self._datasets = {}
#             for pivot_date in [dates.nearest_zulu() - timedelta(days=i)
#                                for i in range(3)]:
#                 logger.info(
#                     'Checking for available RTOFS data for pivot date: '
#                     f'{pivot_date}.')
#                 try:
#                     with warnings.catch_warnings():
#                         warnings.filterwarnings(
#                             "ignore", category=xarray.SerializationWarning)
#                         self._datasets.setdefault(
#                             pivot_date,
#                             xarray.open_dataset(
#                                 f'{self.base_url}'
#                                 f'{pivot_date.strftime("%Y%m%d")}'
#                                 f'/{self.product}')
#                             )
#                     logger.info('Success!')
#                 except OSError as e:
#                     if e.errno == -70:
#                         print()
#                         continue
#                 if len(self._datasets) == 2:
#                     break
#         return self._datasets

#     @property
#     def base_date(self):
#         return self._base_date

#     def _xy_grid(self, lon_idxs, lat_idxs):
#         dataset = list(self.datasets.values())[-1]
#         lon = []
#         for x in dataset.lon[lon_idxs]:
#             if not (x >= dataset.lon.min()
#                     and x < 180.):
#                 lon.append(x-360.)
#             else:
#                 lon.append(x)
#         return np.meshgrid(np.array(lon), dataset.lat[lat_idxs].values)

#     def _modified_bbox(self, bbox=None):
#         dataset = list(self.datasets.values())[-1]
#         if bbox is None:
#             return Bbox.from_extents(
#                 dataset.lon.min(),
#                 dataset.lat.min(),
#                 dataset.lon.max(),
#                 dataset.lat.max()
#             )
#         else:
#             xmin = bbox.xmin + 360. if not (
#                 bbox.xmin >= dataset.lon.min()
#                 and bbox.xmin < 180.) else bbox.xmin

#             xmax = bbox.xmax + 360. if not (
#                 bbox.xmax >= dataset.lon.min()
#                 and bbox.xmax < 180.) else bbox.xmax

#             return Bbox.from_extents(
#                 np.min([xmin, xmax]),
#                 bbox.ymin,
#                 np.max([xmin, xmax]),
#                 bbox.ymax
#             )

#     def _modified_bbox_indexes(
#             self,
#             bbox,
#             pixel_buffer=0
#     ):
#         dataset = list(self.datasets.values())[-1]
#         lat_idxs = np.where((dataset.lat >= bbox.ymin)
#                             & (dataset.lat <= bbox.ymax))[0]
#         lon_idxs = np.where((dataset.lon >= bbox.xmin)
#                             & (dataset.lon <= bbox.xmax))[0]
#         lon_idxs = lon_idxs.tolist()
#         lat_idxs = lat_idxs.tolist()
#         for i in range(pixel_buffer):
#             lon_idxs.insert(lon_idxs[0] - 1, 0)
#             lon_idxs.append(lon_idxs[-1] + 1)
#             lat_idxs.insert(lat_idxs[0] - 1, 0)
#             lat_idxs.append(lat_idxs[-1] + 1)
#         return lon_idxs, lat_idxs


# class RTOFSTemperature(RTOFSBaroclinicComponent):

#     @property
#     def product(self) -> str:
#         return 'rtofs_glo_3dz_forecast_daily_temp'

#     @property
#     def ncvar(self):
#         return 'temperature'

#     @property
#     def fill_value(self):
#         return np.nan

#     @property
#     def nowcast_varname(self):
#         return 'temp'


# class RTOFSSalinity(RTOFSBaroclinicComponent):

#     @property
#     def product(self) -> str:
#         return 'rtofs_glo_3dz_forecast_daily_salt'

#     @property
#     def ncvar(self):
#         return 'salinity'

#     @property
#     def fill_value(self):
#         return 0.

#     @property
#     def nowcast_varname(self):
#         return 'salt'

from pyschism.forcing.hycom import Hycom


class RTOFS(Hycom):
    '''Public API for assigning ROTFS forcing to model.'''

    # def __init__(self):
    #     self._temperature = RTOFSTemperature()
    #     self._salinity = RTOFSSalinity()











# # import os
# # import numpy as np
# # from datetime import datetime, timedelta
# # import logging
# # import pathlib
# # from typing import Union

# # from netCDF4 import Dataset
# # from matplotlib.transforms import Bbox

# # from pyschism.mesh import Hgrid
# # from pyschism.dates import localize_datetime, nearest_cycle_date, pivot_time

# # logger = logging.getLogger(__name__)


# # class HotStartInventory():

# #     def __init__(self):

# #         pass

# #     logger.info('Fetching RTOFS data')

# #     def fetch_data(self, outdir: Union[str, os.PathLike], hgrid, start_date):

# #         self.hgrid = hgrid
# #         self.start_date = start_date

# #         outdir = pathlib.Path(outdir)
# #         if outdir.name != '':
# #             outdir /= 'hotstart'
# #         outdir.mkdir(exist_ok=True)

# #         bbox = self.hgrid.get_bbox(crs='epsg:4326', output_type='bbox')

# #         # Go to yesterday's directory to get tomorrow's data
# #         dt = self.start_date - timedelta(days=1)
# #         nc_ssh = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
# #                          + dt.strftime('%Y%m%d')+f'/rtofs_glo_2ds_forecast_3hrly_diag')
# #         nc_salt = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
# #                           + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_salt')
# #         nc_temp = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
# #                           + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_temp')
# #         nc_uvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
# #                           + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_uvel')
# #         nc_vvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
# #                           + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_vvel')

# #         lon = nc_ssh['lon'][:]
# #         lat = nc_ssh['lat'][:]
# #         lev = nc_salt['lev'][:]

# #         xmin = bbox.x0 + 360. if bbox.x0 < 0 else bbox.x0
# #         xmax = bbox.x1 + 360. if bbox.x1 < 0 else bbox.x1
# #         bbox = Bbox.from_extents(xmin, bbox.y0, xmax, bbox.y1)

# #         lat_idxs = np.where((lat >= bbox.ymin) & (lat <= bbox.ymax))[0]
# #         lon_idxs = np.where((lon >= bbox.xmin) & (lon <= bbox.xmax))[0]

# #         xlon = lon[lon_idxs]
# #         for ilon in range(lon_idxs.size):
# #             if xlon[ilon] > 180.:
# #                 xlon[ilon] = xlon[ilon] - 360.
# #         ylat = lat[lat_idxs]

# #         with Dataset(outdir / 'SSH_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
# #             dst.setncatts({"Conventions": "cf-1.0"})
# #             # dimensions
# #             dst.createDimension('xlon', xlon.shape[0])
# #             dst.createDimension('ylat', ylat.shape[0])
# #             dst.createDimension('time', None)

# #             # variables
# #             # lon
# #             dst.createVariable('xlon', 'f4', ('xlon',))
# #             dst['xlon'].long_name = "Longitude"
# #             dst['xlon'].standard_name = "longitude"
# #             dst['xlon'].units = "degrees_east"
# #             dst['xlon'][:] = xlon
# #             # lat
# #             dst.createVariable('ylat', 'f4', ('ylat',))
# #             dst['ylat'].long_name = "Latitude"
# #             dst['ylat'].standard_name = "latitude"
# #             dst['ylat'].units = "degrees_north"
# #             dst['ylat'][:] = ylat
# #             # time
# #             dst.createVariable('time', 'f4', ('time',))
# #             dst['time'].standard_name = "time"
# #             dst['time'].units = "days since 1-1-1 00:00:0.0"
# #             dst['time'][:] = nc_salt['time'][1:2]
# #             # ssh
# #             dst.createVariable(
# #                 'surf_el', 'f4', ('time', 'ylat', 'xlon',), fill_value=-30000.0)
# #             dst['surf_el'].long_name = "sea_surface_elevation (m)"
# #             dst['surf_el'].add_offset = 0.
# #             dst['surf_el'].scale_factor = 0.001
# #             dst['surf_el'][:, :, :] = nc_ssh['ssh'][8:9, 0, lat_idxs, lon_idxs]

# #         with Dataset(outdir / 'TS_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
# #             dst.setncatts({"Conventions": "cf-1.0"})
# #             # dimensions
# #             dst.createDimension('xlon', xlon.shape[0])
# #             dst.createDimension('ylat', ylat.shape[0])
# #             dst.createDimension('lev', nc_salt['lev'].shape[0])
# #             dst.createDimension('time', None)

# #             # variables
# #             # lon
# #             dst.createVariable('xlon', 'f4', ('xlon',))
# #             dst['xlon'].long_name = "Longitude"
# #             dst['xlon'].standard_name = "longitude"
# #             dst['xlon'].units = "degrees_east"
# #             dst['xlon'][:] = xlon
# #             # lat
# #             dst.createVariable('ylat', 'f4', ('ylat',))
# #             dst['ylat'].long_name = "Latitude"
# #             dst['ylat'].standard_name = "latitude"
# #             dst['ylat'].units = "degrees_north"
# #             dst['ylat'][:] = ylat
# #             # lev
# #             dst.createVariable('lev', 'f4', ('lev',))
# #             dst['lev'].long_name = "altitude"
# #             dst['lev'].units = "millibar"
# #             dst['lev'][:] = nc_salt['lev'][:]
# #             # time
# #             dst.createVariable('time', 'f4', ('time',))
# #             dst['time'].standard_name = "time"
# #             dst['time'].units = "days since 1-1-1 00:00:0.0"
# #             dst['time'][:] = nc_salt['time'][1:2]
# #             # salt
# #             dst.createVariable(
# #                 'salinity', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.0)
# #             dst['salinity'].long_name = "sea_water_salinity (psu)"
# #             dst['salinity'].add_offset = 20.
# #             dst['salinity'].scale_factor = 0.001
# #             # temp
# #             dst.createVariable('temperature', 'f4', ('time',
# #                                                      'lev', 'ylat', 'xlon',), fill_value=-30000.0)
# #             dst['temperature'].long_name = "sea_water_potential_temperature (degc)"
# #             dst['temperature'].add_offset = 20.
# #             dst['temperature'].scale_factor = 0.001

# #             for k in np.arange(len(lev)):
# #                 dst['salinity'][:, k, :, :] = nc_salt['salinity'][1:2,
# #                                                                   k, lat_idxs, lon_idxs]
# #                 dst['temperature'][:, k, :,
# #                                    :] = nc_temp['temperature'][1:2, k, lat_idxs, lon_idxs]

# #         with Dataset(outdir / 'UV_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
# #             dst.setncatts({"Conventions": "cf-1.0"})
# #             # dimensions
# #             dst.createDimension('xlon', xlon.shape[0])
# #             dst.createDimension('ylat', ylat.shape[0])
# #             dst.createDimension('lev', nc_salt['lev'].shape[0])
# #             dst.createDimension('time', None)

# #             # variables
# #             # lon
# #             dst.createVariable('xlon', 'f4', ('xlon',))
# #             dst['xlon'].long_name = "Longitude"
# #             dst['xlon'].standard_name = "longitude"
# #             dst['xlon'].units = "degrees_east"
# #             dst['xlon'][:] = xlon
# #             # lat
# #             dst.createVariable('ylat', 'f4', ('ylat',))
# #             dst['ylat'].long_name = "Latitude"
# #             dst['ylat'].standard_name = "latitude"
# #             dst['ylat'].units = "degrees_north"
# #             dst['ylat'][:] = ylat
# #             # lev
# #             dst.createVariable('lev', 'f4', ('lev',))
# #             dst['lev'].long_name = "altitude"
# #             dst['lev'].units = "millibar"
# #             dst['lev'][:] = nc_salt['lev'][:]
# #             # time
# #             dst.createVariable('time', 'f4', ('time',))
# #             dst['time'].standard_name = "time"
# #             dst['time'].units = "days since 1-1-1 00:00:0.0"
# #             dst['time'][:] = nc_salt['time'][1:2]
# #             # uvel
# #             dst.createVariable(
# #                 'water_u', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.0)
# #             dst['water_u'].long_name = "eastward_sea_water_velocity (m/s)"
# #             dst['water_u'].add_offset = 0.
# #             dst['water_u'].scale_factor = 0.001
# #             # vvel
# #             dst.createVariable(
# #                 'water_v', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.0)
# #             dst['water_v'].long_name = "northward_sea_water_velocity (m/s)"
# #             dst['water_v'].add_offset = 0.
# #             dst['water_v'].scale_factor = 0.001

# #             for k in np.arange(len(lev)):
# #                 dst['water_u'][:, k, :, :] = nc_uvel['u'][1:2, k, lat_idxs, lon_idxs]
# #                 dst['water_v'][:, k, :, :] = nc_vvel['v'][1:2, k, lat_idxs, lon_idxs]

# #        # symlink estaury.gr3 and *.in file
# #        # os.symlink(estaury.gr3, './start/estuary.gr3')


# # class OpenBoundaryInventory():

# #     def __init__(self):

# #         pass

# #     def fetch_data(self, outdir: Union[str, os.PathLike], start_date, rnday, idx_min, idx_max, jdx_min, jdx_max):

# #         self.start_date = start_date
# #         self.rnday = rnday

# #         outdir = pathlib.Path(outdir)
# #         if outdir.name != '':
# #             outdir /= '3Dth'
# #         outdir.mkdir(exist_ok=True)

# #         # Go to yesterday's directory to get tomorrow's data
# #         dt = self.start_date - timedelta(days=1)
# #         nc_ssh = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
# #                          + dt.strftime('%Y%m%d')+f'/rtofs_glo_2ds_forecast_3hrly_diag')
# #         nc_salt = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
# #                           + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_salt')
# #         nc_temp = Dataset(f'http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
# #                           + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_temp')
# #         nc_uvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
# #                           + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_uvel')
# #         nc_vvel = Dataset('http://nomads.ncep.noaa.gov:80/dods/rtofs/rtofs_global'
# #                           + dt.strftime('%Y%m%d')+f'/rtofs_glo_3dz_forecast_daily_vvel')

# #         lon = nc_ssh['lon'][:]
# #         lat = nc_ssh['lat'][:]
# #         lev = nc_salt['lev'][:]

# #         xlon = lon[idx_min:idx_max+1]
# #         for ilon in range(xlon.size):
# #             if xlon[ilon] > 180.:
# #                 xlon[ilon] = xlon[ilon] - 360.
# #         ylat = lat[jdx_min:jdx_max+1]

# #         with Dataset(outdir / 'SSH_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
# #             dst.setncatts({"Conventions": "cf-1.0"})
# #             # dimensions
# #             dst.createDimension('xlon', xlon.shape[0])
# #             dst.createDimension('ylat', ylat.shape[0])
# #             dst.createDimension('time', None)

# #             # variables
# #             # lon
# #             dst.createVariable('xlon', 'f4', ('xlon',))
# #             dst['xlon'].long_name = "Longitude"
# #             dst['xlon'].standard_name = "longitude"
# #             dst['xlon'].units = "degrees_east"
# #             dst['xlon'][:] = xlon
# #             # lat
# #             dst.createVariable('ylat', 'f4', ('ylat',))
# #             dst['ylat'].long_name = "Latitude"
# #             dst['ylat'].standard_name = "latitude"
# #             dst['ylat'].units = "degrees_north"
# #             dst['ylat'][:] = ylat
# #             # time
# #             dst.createVariable('time', 'f4', ('time',))
# #             dst['time'].standard_name = "time"
# #             dst['time'].units = "days since 1-1-1 00:00:0.0"
# #             dst['time'][:] = nc_salt['time'][1:rnday+2]
# #             # ssh
# #             dst.createVariable(
# #                 'surf_el', 'f4', ('time', 'ylat', 'xlon',), fill_value=-30000.0)
# #             dst['surf_el'].long_name = "sea_surface_elevation (m)"
# #             dst['surf_el'].add_offset = 0.
# #             dst['surf_el'].scale_factor = 0.001
# #             ssh = nc_ssh['ssh'][8:8*(rnday+1)+1:8, 0,
# #                                 jdx_min:jdx_max+1, idx_min:idx_max+1]
# #             dst['surf_el'][:, :, :] = ssh

# #         with Dataset(outdir / 'TS_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
# #             dst.setncatts({"Conventions": "cf-1.0"})
# #             # dimensions
# #             dst.createDimension('xlon', xlon.shape[0])
# #             dst.createDimension('ylat', ylat.shape[0])
# #             dst.createDimension('lev', nc_salt['lev'].shape[0])
# #             dst.createDimension('time', None)

# #             # variables
# #             # lon
# #             dst.createVariable('xlon', 'f4', ('xlon',))
# #             dst['xlon'].long_name = "Longitude"
# #             dst['xlon'].standard_name = "longitude"
# #             dst['xlon'].units = "degrees_east"
# #             dst['xlon'][:] = xlon
# #             # lat
# #             dst.createVariable('ylat', 'f4', ('ylat',))
# #             dst['ylat'].long_name = "Latitude"
# #             dst['ylat'].standard_name = "latitude"
# #             dst['ylat'].units = "degrees_north"
# #             dst['ylat'][:] = ylat
# #             # lev
# #             dst.createVariable('lev', 'f4', ('lev',))
# #             dst['lev'].long_name = "altitude"
# #             dst['lev'].units = "millibar"
# #             dst['lev'][:] = nc_salt['lev'][:]
# #             # time
# #             dst.createVariable('time', 'f4', ('time',))
# #             dst['time'].standard_name = "time"
# #             dst['time'].units = "days since 1-1-1 00:00:0.0"
# #             dst['time'][:] = nc_salt['time'][1:rnday+2]
# #             # salt
# #             dst.createVariable(
# #                 'salinity', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.)
# #             dst['salinity'].long_name = "sea_water_salinity (psu)"
# #             dst['salinity'].add_offset = 20.
# #             dst['salinity'].scale_factor = 0.001
# #             sss = nc_salt['salinity'][1:rnday+2, :,
# #                                       jdx_min:jdx_max+1, idx_min:idx_max+1]
# #             dst['salinity'][:, :, :, :] = sss
# #             # temp
# #             dst.createVariable('temperature', 'f4', ('time',
# #                                                      'lev', 'ylat', 'xlon',), fill_value=-30000.)
# #             dst['temperature'].long_name = "sea_water_potential_temperature (degc)"
# #             dst['temperature'].add_offset = 20.
# #             dst['temperature'].scale_factor = 0.001
# #             sst = nc_temp['temperature'][1:rnday+2, :,
# #                                          jdx_min:jdx_max+1, idx_min:idx_max+1]
# #             dst['temperature'][:, :, :, :] = sst

# #         with Dataset(outdir / 'UV_1.nc', 'w', format='NETCDF3_CLASSIC') as dst:
# #             dst.setncatts({"Conventions": "cf-1.0"})
# #             # dimensions
# #             dst.createDimension('xlon', xlon.shape[0])
# #             dst.createDimension('ylat', ylat.shape[0])
# #             dst.createDimension('lev', nc_salt['lev'].shape[0])
# #             dst.createDimension('time', None)

# #             # variables
# #             # lon
# #             dst.createVariable('xlon', 'f4', ('xlon',))
# #             dst['xlon'].long_name = "Longitude"
# #             dst['xlon'].standard_name = "longitude"
# #             dst['xlon'].units = "degrees_east"
# #             dst['xlon'][:] = xlon
# #             # lat
# #             dst.createVariable('ylat', 'f4', ('ylat',))
# #             dst['ylat'].long_name = "Latitude"
# #             dst['ylat'].standard_name = "latitude"
# #             dst['ylat'].units = "degrees_north"
# #             dst['ylat'][:] = ylat
# #             # lev
# #             dst.createVariable('lev', 'f4', ('lev',))
# #             dst['lev'].long_name = "altitude"
# #             dst['lev'].units = "millibar"
# #             dst['lev'][:] = nc_salt['lev'][:]
# #             # time
# #             dst.createVariable('time', 'f4', ('time',))
# #             dst['time'].standard_name = "time"
# #             dst['time'].units = "days since 1-1-1 00:00:0.0"
# #             dst['time'][:] = nc_salt['time'][1:rnday+2]
# #             # uvel
# #             dst.createVariable(
# #                 'water_u', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.)
# #             dst['water_u'].long_name = "eastward_sea_water_velocity (m/s)"
# #             dst['water_u'].add_offset = 0.
# #             dst['water_u'].scale_factor = 0.001
# #             uvel = nc_uvel['u'][1:rnday+2, :,
# #                                 jdx_min:jdx_max+1, idx_min:idx_max+1]
# #             dst['water_u'][:, :, :, :] = uvel
# #             # vvel
# #             dst.createVariable(
# #                 'water_v', 'f4', ('time', 'lev', 'ylat', 'xlon',), fill_value=-30000.)
# #             dst['water_v'].long_name = "northward_sea_water_velocity (m/s)"
# #             dst['water_v'].add_offset = 0.
# #             dst['water_v'].scale_factor = 0.001
# #             vvel = nc_vvel['v'][1:rnday+2, :,
# #                                 jdx_min:jdx_max+1, idx_min:idx_max+1]
# #             dst['water_v'][:, :, :, :] = vvel
