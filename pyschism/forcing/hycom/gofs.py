from abc import ABC, abstractmethod
from datetime import datetime, timedelta
from functools import lru_cache
import logging
from typing import Dict, Union

from matplotlib.transforms import Bbox
import seawater as sw
from netCDF4 import Dataset
import numpy as np
import requests
from scipy.interpolate import griddata, RegularGridInterpolator, NearestNDInterpolator
import tqdm
import tqdm_logging_wrapper
import xmltodict

from pyschism import dates
from pyschism.forcing.hycom.base import Hycom, HycomComponent

logger = logging.getLogger(__name__)


class GofsDatasetCollection(ABC):

    start_date = dates.StartDate()
    end_date = dates.EndDate()

    def __init__(
            self,
            start_date: datetime = None,
            end_date: Union[datetime, timedelta] = None,
            output_interval: timedelta = None
    ):
        self.start_date = start_date
        self.end_date = self.maximum_end_date - self.sampling_interval if \
            end_date is None else end_date 
        self.output_interval = self.sampling_interval if output_interval is None \
            else output_interval
        if (self.output_interval % self.sampling_interval).total_seconds() != 0:
            raise ValueError(
                'Argument output_interval must be modulus of '
                f'{self.sampling_interval} hours.')

    @property
    @abstractmethod
    def datasets(self):
        '''Get the datasets based on the set start_date and end_date dates.'''

    @property
    def sampling_frequency(self):
        return 1./self.sampling_interval.total_seconds()

    @property
    def sampling_interval(self):
        return timedelta(hours=3)

    @property
    def required_datevector(self):
        required_start_date = dates.nearest_cycle(self.start_date)
        required_end_date = dates.nearest_cycle(
            self.end_date + self.output_interval,
            period=int(self.output_interval/timedelta(hours=1))
        )
        return np.arange(
            required_start_date,
            required_end_date + self.output_interval,
            self.output_interval
        ).astype(datetime)


class GofsForecastDatasets(GofsDatasetCollection):

    catalog_url = 'https://tds.hycom.org/thredds/catalog/GLBy0.08/expt_93.0/FMRC/runs/catalog.xml'
    base_url = 'https://tds.hycom.org/thredds/dodsC/'

    @property
    def datasets(self):
        datasets = {}
        for required_date in self.required_datevector:
            for i, datevector in enumerate(self.datevectors):
                if required_date in datevector:
                    opendap_url = self.base_url + self.xmlcatalog[
                        'catalog']['dataset']['dataset'][i]['@urlPath']
                    datasets.setdefault(
                        required_date,
                        Dataset(opendap_url)
                    )
                    break
        return datasets
        
    @property
    def datevectors(self):
        for file in self.xmlcatalog['catalog']['dataset']['dataset']:
            start = datetime.strptime(
                file['timeCoverage']['start'], '%Y-%m-%dT%H:%M:%SZ')
            end = datetime.strptime(
                file['timeCoverage']['end'], '%Y-%m-%dT%H:%M:%SZ')
            steps = int(
                (end - start).total_seconds() * self.sampling_frequency) + 1
            datevec = [start + x*timedelta(
                seconds=1/self.sampling_frequency) for x in range(steps)]
            yield datevec

    @property
    def xmlcatalog(self):
        return xmltodict.parse(requests.get(self.catalog_url).content)

    @property
    def minimum_datetime(self):
        start_dates = []
        for rec in self.xmlcatalog['catalog']['dataset']['dataset']:
            start_dates.append(dates.localize_datetime(datetime.strptime(
                rec['timeCoverage']['start'], '%Y-%m-%dT%H:%M:%SZ')))
        return np.min(start_dates)

    @property
    def maximum_datetime(self):
        end_dates = []
        for rec in self.xmlcatalog['catalog']['dataset']['dataset']:
            end_dates.append(dates.localize_datetime(datetime.strptime(
                rec['timeCoverage']['end'], '%Y-%m-%dT%H:%M:%SZ')))
        return np.max(end_dates)

    @property
    def maximum_time_range(self):
        return self.maximum_datetime - self.minimum_datetime

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
    elif date >= datetime(1994, 1, 1) and date < datetime(2016, 1, 3):
        database = f'GLBv0.08/expt_53.X/data/{date.year}'
    else:
        print('No data for {date}')
    return database

class GofsHindcastDatasets(GofsDatasetCollection):

    @property
    def datasets(self):
        datasets = {}
        for required_date in self.required_datevector:
            database = get_database(required_date)
            print(f'Database for {required_date} is {database}')
            #baseurl = f'https://tds.hycom.org/thredds/dodsC/{database}?lat[0:1:3250],lon[0:1:4499],' + \
            #    f'time[0:1:6127],depth[0:1:39]' 
            #ds=Dataset(baseurl)
            #time1=ds['time']
            #times=nc4.num2date(time1,units=time1.units,only_use_cftime_datetimes=False)
            #time_idx=np.where( required_date == times)[0].item()
             
            opendap_url = f'https://tds.hycom.org/thredds/dodsC/{database}?lat[0:-1],lon[0:-1],' + \
                f'time[0:-1],surf_el[0:-1][0:-1][0:-1],depth[0:-1],' + \
                f'water_temp[0:-1][0:-1][0:-1][0:1:-1],' + \
                f'salinity[0:-1][0:-1][0:-1][0:-1],' + \
                f'water_u[0:-1][0:-1][0:-1][0:-1],' + \
                f'water_v[0:-1][0:-1][0:-1][0:-1]'
            datasets.setdefault(
                required_date,
                Dataset(opendap_url)
                    )
        return datasets
        #raise NotImplementedError('Need to return the datasets.')

        

    def pad_datasets(self, datasets):
        for dataset in datasets.values():
            if dataset is None:
                raise NotImplementedError(
                    'Hindcast data fetching is unavailable.')


class GofsDatasets:

    def __init__(self, start_date, run_days, output_interval):
        print(f'start_date is {start_date}')
        print(f'today is {datetime.now().strftime("%Y-%m-%d")}')
        self.start_date = start_date
        if start_date.strftime("%Y-%m-%d") < datetime.now().strftime("%Y-%m-%d"):
            self.hindcast = GofsHindcastDatasets(start_date, run_days, output_interval)
        else:
            self.forecast = GofsForecastDatasets(start_date, run_days, output_interval)

    @property
    def datasets(self):
        if self.start_date.strftime("%Y-%m-%d") < datetime.now().strftime("%Y-%m-%d"):
            datasets = self.hindcast.datasets
        else:
            datasets = self.forecast.datasets
        #self.hindcast.pad_datasets(datasets)
        return datasets


class GOFSComponent(HycomComponent):

    @lru_cache(maxsize=None)
    def get_datasets(
            self,
            start_date: datetime,
            run_days: Union[float, timedelta],
            output_interval=timedelta(days=1)
    ) -> Dict[datetime, Dataset]:
        return GofsDatasets(start_date, run_days, output_interval).datasets

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

class GOFSElevation(GOFSComponent):

    @property
    def ncvar(self):
        return 'surf_el'

    def put_boundary_ncdata(
            self,
            boundary,
            dst,
            start_date,
            run_days,
            overwrite=False,
            offset=0,
            output_interval=timedelta(hours=24),
            pixel_buffer=10,
            progress_bar=True
    ):
        for i, (time, dataset) in enumerate(
            self.get_datasets(
                start_date,
                run_days,
                output_interval
            ).items()
        ):
            if start_date.strftime("%Y-%m-%d") < datetime.now().strftime("%Y-%m-%d"):
                ds_base_date = datetime.strptime(
                    ''.join(dataset['time'].units.split()[2:]), 
                    '%Y-%m-%d%H:%M:%S')
                #time1 = dataset['time']
                #ds_timevector = nc4.num2date(
                #    time1,
                #    units=time1.units,
                #    only_use_cftime_datetimes=False) 
            else:
                ds_base_date = datetime.strptime(
                    ''.join(dataset['time'].units.split()[2:-1]),
                    '%Y-%m-%d%H:%M:%S.%f')
            ds_timevector = [ds_base_date + timedelta(hours=x)
                         for x in dataset['time'][:]]
            requested_date = dates.nearest_cycle(
                start_date + i*output_interval,
                period=3).replace(tzinfo=None)
            time_idx = ds_timevector.index(requested_date)
            logger.info(
                'Saving GOFS elev data for date: '
                f'{start_date+i*output_interval} '
                f'approximated as {ds_timevector[time_idx]} for '
                f'boundary id={boundary.id}'
                )
            bounds = boundary.geometry.bounds
            dx = (dataset['lon'][-1] - dataset['lon'][0]) / len(dataset['lon'])
            dy = (dataset['lat'][-1] - dataset['lat'][0]) / len(dataset['lat'])
            bounds = (
                bounds[0] - 2*dx,
                bounds[1] - 2*dy,
                bounds[2] + 2*dx,
                bounds[3] + 2*dy,
                )
            bbox = self._modified_bbox(
                dataset, Bbox.from_extents(*bounds))
            lon_idxs, lat_idxs = self._modified_bbox_indexes(
                    bbox,
                    dataset,
                    pixel_buffer
                )
            #zi = np.full((len(lat_idxs), len(lon_idxs)), np.nan)
            if progress_bar is True:
                items_iter = tqdm.tqdm(lat_idxs)
                with tqdm_logging_wrapper.wrap_logging_for_tqdm(
                        items_iter), items_iter:
                    #for k, lat_idx in enumerate(items_iter):
                    #    zi[k, :] = dataset[self.ncvar][
                    #                        time_idx, lat_idx, lon_idxs]
                    zi = dataset[self.ncvar][time_idx, lat_idxs, lon_idxs]
            else:
                #for k, lat_idx in enumerate(lat_idxs):
                #    zi[k, :] = dataset[self.ncvar][
                #                        time_idx, lat_idx, lon_idxs]
                zi = dataset[self.ncvar][time_idx, lat_idxs, lon_idxs]
            idxs = np.where(abs(zi) > 10000)
            zi[idxs] = float('nan')

            xi = dataset['lon'][lon_idxs]
            for idx in range(len(xi)):
                if xi[idx] > 180:
                    xi[idx] = xi[idx]-360.
            yi = dataset['lat'][lat_idxs]
            xi, yi = np.meshgrid(xi, yi)
            xi = xi.flatten()
            yi = yi.flatten()
            zi = zi.flatten()
            xyq = np.array(boundary.geometry.coords)
            zq = griddata(
                (xi, yi),
                zi,
                (xyq[:, 0], xyq[:, 1]),
                method='linear',
                fill_value=np.nan,
            )
            nan_idxs = np.where(np.isnan(zq))
            q_non_nan = np.where(~np.isnan(zi))
            zq[nan_idxs] = griddata(
                (xi[q_non_nan], yi[q_non_nan]),
                zi[q_non_nan],
                (xyq[nan_idxs, 0], xyq[nan_idxs, 1]),
                method='nearest',
            )
            if np.any(np.isnan(zq)):
                raise ValueError('Boundary contains NaNs.')
            print(f'the shape of zq is {len(zq)}, max zq is {np.max(zq)}, min zq is {np.min(zq)}')
            dst['time_series'][i, offset:offset+len(zq)] = zq


class GOFSVelocity(GOFSComponent):

    @property
    def ncvar(self):
        return 'water_u', 'water_v'

    def put_boundary_ncdata(
            self,
            hgrid,
            vgrid,
            boundary,
            dst,
            start_date,
            run_days,
            overwrite=False,
            offset=0,
            output_interval=timedelta(hours=24),
            pixel_buffer=10,
            progress_bar=True
    ):

        for i, (time, dataset) in enumerate(
            self.get_datasets(
                start_date,
                run_days,
                output_interval
            ).items()
        ):
            if start_date.strftime("%Y-%m-%d") < datetime.now().strftime("%Y-%m-%d"):
                ds_base_date = datetime.strptime(
                    ''.join(dataset['time'].units.split()[2:]),
                    '%Y-%m-%d%H:%M:%S')
                #time1 = dataset['time']
                #ds_timevector = nc4.num2date(
                #    time1,
                #    units=time1.units,
                #    only_use_cftime_datetimes=False) 
            else:
                ds_base_date = datetime.strptime(
                    ''.join(dataset['time'].units.split()[2:-1]),
                    '%Y-%m-%d%H:%M:%S.%f')
            ds_timevector = [ds_base_date + timedelta(hours=x)
                             for x in dataset['time'][:]]
            requested_date = dates.nearest_cycle(
                start_date + i*output_interval,
                period=3).replace(tzinfo=None)
            time_idx = ds_timevector.index(requested_date)
            logger.info(
                'Saving GOFS water_uv data for date: '
                f'{start_date+i*output_interval} '
                f'approximated as {ds_timevector[time_idx]} for '
                f'boundary id={boundary.id}'
                )
            bounds = boundary.geometry.bounds
            dx = (dataset['lon'][-1] - dataset['lon'][0]) / len(dataset['lon'])
            dy = (dataset['lat'][-1] - dataset['lat'][0]) / len(dataset['lat'])
            bounds = (
                bounds[0] - 2*dx,
                bounds[1] - 2*dy,
                bounds[2] + 2*dx,
                bounds[3] + 2*dy,
                )
            bbox = self._modified_bbox(
                dataset, Bbox.from_extents(*bounds))
            lon_idxs, lat_idxs = self._modified_bbox_indexes(
                    bbox,
                    dataset,
                    pixel_buffer
                )
            uvar, vvar = self.ncvar
            # z_ui_idxs = list(range(dataset['depth'].shape[0]))
            z_idxs = list(range(dataset['depth'].shape[0]))  # TODO: subset?
            #ui = np.full((len(z_idxs), len(lat_idxs), len(lon_idxs)), np.nan)
            #vi = np.full((len(z_idxs), len(lat_idxs), len(lon_idxs)), np.nan)
            if progress_bar is True:
                items_iter = tqdm.tqdm(lat_idxs)
                with tqdm_logging_wrapper.wrap_logging_for_tqdm(
                        items_iter), items_iter:
                    #for k, lat_idx in enumerate(items_iter):
                    #    ui[:, k, :] = dataset[uvar][time_idx, z_idxs, lat_idx, lon_idxs]
                    #    vi[:, k, :] = dataset[vvar][time_idx, z_idxs, lat_idx, lon_idxs]
                    ui = dataset[uvar][time_idx, :, lat_idxs, lon_idxs]
                    vi = dataset[vvar][time_idx, :, lat_idxs, lon_idxs]
            else:
                #for k, lat_idx in enumerate(lat_idxs):
                #    ui[:, k, :] = dataset[uvar][time_idx, z_idxs, lat_idx, lon_idxs]
                #    vi[:, k, :] = dataset[vvar][time_idx, z_idxs, lat_idx, lon_idxs]
                ui = dataset[uvar][time_idx, :, lat_idxs, lon_idxs]
                vi = dataset[vvar][time_idx, :, lat_idxs, lon_idxs]
            #change missing value to nan
            idxs = np.where(abs(ui) > 10000)
            ui[idxs] = float('nan')
            idxs = np.where(abs(vi) > 10000)
            vi[idxs] = float('nan')

            loni = dataset['lon'][lon_idxs]
            for idx in range(len(loni)):
                if loni[idx] > 180:
                    loni[idx] = loni[idx]-360.
            lati = dataset['lat'][lat_idxs]
            xi, yi = transform_ll_to_cpp(loni, lati)

            if vgrid.ivcor == 1:
                bz = (hgrid.values[:, None]*vgrid.sigma)[boundary.indexes, :]
                idxs = np.where(bz > 5000.0)
                bz[idxs] = 5000.0 - 1.0e-6 
            else:
                raise NotImplementedError('vgrid.ivcor!=1')

            xy = hgrid.get_xy(crs='epsg:4326')
            lonb = xy[boundary.indexes, 0]
            latb = xy[boundary.indexes, 1]
            xb, yb = transform_ll_to_cpp(lonb, latb)
            bx = np.tile(xb, [bz.shape[1],1]).T
            by = np.tile(yb, [bz.shape[1],1]).T
            bzyx = np.c_[bz.reshape(np.size(bz)), by.reshape(np.size(by)), bx.reshape(np.size(bx))]
            zi = dataset['depth'][z_idxs]

            # First try with RegularGridInterpolator
            ui_fd = RegularGridInterpolator(
                (zi, yi, xi),
                ui,
                method='linear',
                bounds_error=False,
                fill_value=np.nan
            )
            u_interp = ui_fd(bzyx)

            # the boundary and the data don't intersect
            if np.all(np.isnan(u_interp)):
                ui_idxs = np.where(~np.isnan(ui))
                xyzi = np.vstack([
                    np.tile(xi, ui.shape)[ui_idxs].flatten(),
                    np.tile(yi, ui.shape)[ui_idxs].flatten(),
                    np.tile(zi, ui.shape)[ui_idxs].flatten()
                ]).T
                ui_fd = NearestNDInterpolator(xyzi, ui[ui_idxs].flatten())
                u_interp = ui_fd(np.vstack([bx, by, -bz.flatten()]).T)
            # the boundary and the data partially intersect
            elif np.any(np.isnan(u_interp)):
                ui_idxs = np.where(~np.isnan(u_interp))
                ui_fd = NearestNDInterpolator(bzyx[ui_idxs], u_interp[ui_idxs])
                ui_idxs = np.where(np.isnan(u_interp))
                u_interp[ui_idxs] = ui_fd(bzyx[ui_idxs])

            if np.any(np.isnan(u_interp)):
                raise ValueError('No boundary  u velocity data for GOFS. '
                                 'Try increasing pixel_buffer argument.')

            vi_fd = RegularGridInterpolator(
                (zi, yi, xi),
                vi,
                method='linear',
                bounds_error=False,
                fill_value=np.nan
            )

            v_interp = vi_fd(bzyx)

            if np.all(np.isnan(v_interp)):
                vi_idxs = np.where(~np.isnan(vi))
                xi = np.tile(xi, vi.shape)[vi_idxs].flatten()
                yi = np.tile(yi, vi.shape)[vi_idxs].flatten()
                zi = np.tile(zi, vi.shape)[vi_idxs].flatten()
                xyzi = np.vstack([xi, yi, zi]).T
                vi_fd = NearestNDInterpolator(xyzi, vi[ui_idxs].flatten())
                v_interp = ui_fd(np.vstack([bx, by, -bz.flatten()]).T)

            elif np.any(np.isnan(v_interp)):
                vi_idxs = np.where(~np.isnan(v_interp))
                vi_fd = NearestNDInterpolator(bzyx[vi_idxs], v_interp[vi_idxs])
                vi_idxs = np.where(np.isnan(v_interp))
                v_interp[vi_idxs] = vi_fd(bzyx[vi_idxs])

            if np.any(np.isnan(v_interp)):
                raise ValueError('No boundary data for GOFS. Try increasing pixel_buffer argument.')

            dst['time_series'][i, offset:offset+bz.shape[0], :, 0] = u_interp.reshape(bz.shape)
            dst['time_series'][i, offset:offset+bz.shape[0], :, 1] = v_interp.reshape(bz.shape)


class GOFSTemperature(GOFSComponent):

    @property
    def ncvar(self):
        return 'water_temp'

    def put_boundary_ncdata(
            self,
            hgrid,
            vgrid,
            boundary,
            dst,
            start_date,
            run_days,
            overwrite=False,
            offset=0,
            output_interval=timedelta(hours=24),
            pixel_buffer=10,
            progress_bar=True
    ):

        for i, (time, dataset) in enumerate(
            self.get_datasets(
                start_date,
                run_days,
                output_interval
            ).items()
        ):
            if start_date.strftime("%Y-%m-%d") < datetime.now().strftime("%Y-%m-%d"):
                ds_base_date = datetime.strptime(
                    ''.join(dataset['time'].units.split()[2:]),
                    '%Y-%m-%d%H:%M:%S')
                #time1 = dataset['time']
                #ds_timevector = nc4.num2date(
                #    time1,
                #    units=time1.units,
                #    only_use_cftime_datetimes=False) 
            else:
                ds_base_date = datetime.strptime(
                    ''.join(dataset['time'].units.split()[2:-1]),
                    '%Y-%m-%d%H:%M:%S.%f')
            ds_timevector = [ds_base_date + timedelta(hours=x)
                             for x in dataset['time'][:]]
            requested_date = dates.nearest_cycle(
                start_date + i*output_interval,
                period=3).replace(tzinfo=None)
            time_idx = ds_timevector.index(requested_date)
            logger.info(
                'Saving GOFS temp data for date: '
                f'{start_date+i*output_interval} '
                f'approximated as {ds_timevector[time_idx]} for '
                f'boundary id={boundary.id}'
                )
            bounds = boundary.geometry.bounds
            dx = (dataset['lon'][-1] - dataset['lon'][0]) / len(dataset['lon'])
            dy = (dataset['lat'][-1] - dataset['lat'][0]) / len(dataset['lat'])
            bounds = (
                bounds[0] - 2*dx,
                bounds[1] - 2*dy,
                bounds[2] + 2*dx,
                bounds[3] + 2*dy,
                )
            bbox = self._modified_bbox(
                dataset, Bbox.from_extents(*bounds))
            lon_idxs, lat_idxs = self._modified_bbox_indexes(
                    bbox,
                    dataset,
                    pixel_buffer
                )
            # z_ui_idxs = list(range(dataset['depth'].shape[0]))
            z_idxs = list(range(dataset['depth'].shape[0]))  # TODO: subset?
            #temp = np.full((len(z_idxs), len(lat_idxs), len(lon_idxs)), np.nan)
            if progress_bar is True:
                items_iter = tqdm.tqdm(lat_idxs)
                with tqdm_logging_wrapper.wrap_logging_for_tqdm(
                        items_iter), items_iter:
                    #for k, lat_idx in enumerate(items_iter):
                    #    temp[:, k, :] = dataset[self.ncvar][time_idx, z_idxs, lat_idx, lon_idxs]
                    temp = dataset[self.ncvar][time_idx, :, lat_idxs, lon_idxs]
                    salt = dataset['salinity'][time_idx, :, lat_idxs, lon_idxs]
            else:
                #for k, lat_idx in enumerate(lat_idxs):
                #    temp[:, k, :] = dataset[self.ncvar][time_idx, z_idxs, lat_idx, lon_idxs]
                temp = np.squeeze(dataset[self.ncvar][time_idx, :, lat_idxs, lon_idxs])
                salt = np.squeeze(dataset['salinity'][time_idx, :, lat_idxs, lon_idxs])

            #convert in-situ temperature to potential temperature
            print(f'The shape of temp is {temp.shape}')
            nz = temp.shape[0]
            ny = temp.shape[1]
            nx = temp.shape[2]
            dep = dataset['depth'][:]
            pr=np.ones(temp.shape)
            pre = pr*dep[:,None, None]
            Pr = np.zeros(temp.shape)
            ptemp = sw.ptmp(salt, temp, pre, Pr)*1.00024

            #change missing value to nan
            idxs = np.where(abs(ptemp) > 10000)
            ptemp[idxs] = float('nan')

            loni = dataset['lon'][lon_idxs]
            for idx in range(len(loni)):
                if loni[idx] > 180:
                    loni[idx] = loni[idx]-360.
            lati = dataset['lat'][lat_idxs]
            xi, yi = transform_ll_to_cpp(loni, lati)

            if vgrid.ivcor == 1:
                bz = (hgrid.values[:, None]*vgrid.sigma)[boundary.indexes, :]
                idxs = np.where(bz > 5000.0)
                bz[idxs] = 5000.0 - 1.0e-6
                print(f'zcor at 200 is {bz[200,:]}')
            else:
                raise NotImplementedError('vgrid.ivcor!=1')

            xy = hgrid.get_xy(crs='epsg:4326')
            lonb = xy[boundary.indexes, 0]
            latb = xy[boundary.indexes, 1]
            xb, yb = transform_ll_to_cpp(lonb, latb)
            #bx = np.tile(xy[boundary.indexes, 0], [bz.shape[1],1])
            #by = np.tile(xy[boundary.indexes, 1], [bz.shape[1],1])
            #bzyx = np.vstack([bz.flatten(), by, bx]).T
            bx = np.tile(xb, [bz.shape[1],1]).T
            by = np.tile(yb, [bz.shape[1],1]).T
            bzyx = np.c_[bz.reshape(np.size(bz)), by.reshape(np.size(by)), bx.reshape(np.size(bx))]
            zi = dataset['depth'][z_idxs]

            # First try with RegularGridInterpolator
            ptemp_fd = RegularGridInterpolator(
                (zi, yi, xi),
                ptemp,
                method='linear',
                bounds_error=False,
                fill_value=np.nan
            )
            ptemp_interp = ptemp_fd(bzyx)

            # the boundary and the data don't intersect
            if np.all(np.isnan(ptemp_interp)):
                ptemp_idxs = np.where(~np.isnan(ptemp))
                xyzi = np.vstack([
                    np.tile(xi, ptemp.shape)[ptemp_idxs].flatten(),
                    np.tile(yi, ptemp.shape)[ptemp_idxs].flatten(),
                    np.tile(zi, ptemp.shape)[ptemp_idxs].flatten()
                ]).T
                ptemp_fd = NearestNDInterpolator(xyzi, ptemp[ptemp_idxs].flatten())
                ptemp_interp = ptemp_fd(np.vstack([bx, by, -bz.flatten()]).T)
            # the boundary and the data partially intersect
            elif np.any(np.isnan(ptemp_interp)):
                ptemp_idxs = np.where(~np.isnan(ptemp_interp))
                ptemp_fd = NearestNDInterpolator(bzyx[ptemp_idxs], ptemp_interp[ptemp_idxs])
                ptemp_idxs = np.where(np.isnan(ptemp_interp))
                ptemp_interp[ptemp_idxs] = ptemp_fd(bzyx[ptemp_idxs])

            if np.any(np.isnan(ptemp_interp)):
                raise ValueError('No boundary  temperature data for GOFS. '
                                 'Try increasing pixel_buffer argument.')
            dst['time_series'][i, offset:offset+bz.shape[0], :, :] = ptemp_interp.reshape(bz.shape)


class GOFSSalinity(GOFSComponent):

    @property
    def ncvar(self):
        return 'salinity'

    def put_boundary_ncdata(
            self,
            hgrid,
            vgrid,
            boundary,
            dst,
            start_date,
            run_days,
            overwrite=False,
            offset=0,
            output_interval=timedelta(hours=24),
            pixel_buffer=10,
            progress_bar=True
    ):

        for i, (time, dataset) in enumerate(
            self.get_datasets(
                start_date,
                run_days,
                output_interval
            ).items()
        ):
            if start_date.strftime("%Y-%m-%d") < datetime.now().strftime("%Y-%m-%d"):
                ds_base_date = datetime.strptime(
                    ''.join(dataset['time'].units.split()[2:]),
                    '%Y-%m-%d%H:%M:%S')
                #time1 = dataset['time']
                #ds_timevector = nc4.num2date(
                #    time1,
                #    units=time1.units,
                #    only_use_cftime_datetimes=False) 
            else:
                ds_base_date = datetime.strptime(
                    ''.join(dataset['time'].units.split()[2:-1]),
                    '%Y-%m-%d%H:%M:%S.%f')
            ds_timevector = [ds_base_date + timedelta(hours=x)
                             for x in dataset['time'][:]]
            requested_date = dates.nearest_cycle(
                start_date + i*output_interval,
                period=3).replace(tzinfo=None)
            time_idx = ds_timevector.index(requested_date)
            print(f'time_idx for {requested_date} is {time_idx}')
            logger.info(
                'Saving GOFS salinity data for date: '
                f'{start_date+i*output_interval} '
                f'approximated as {ds_timevector[time_idx]} for '
                f'boundary id={boundary.id}'
                )
            bounds = boundary.geometry.bounds
            dx = (dataset['lon'][-1] - dataset['lon'][0]) / len(dataset['lon'])
            dy = (dataset['lat'][-1] - dataset['lat'][0]) / len(dataset['lat'])
            bounds = (
                bounds[0] - 2*dx,
                bounds[1] - 2*dy,
                bounds[2] + 2*dx,
                bounds[3] + 2*dy,
                )
            bbox = self._modified_bbox(
                dataset, Bbox.from_extents(*bounds))
            lon_idxs, lat_idxs = self._modified_bbox_indexes(
                    bbox,
                    dataset,
                    pixel_buffer
                )
            print(f'lon_idxs is {lon_idxs}')
            print(f'lat_idxs is {lat_idxs}')

            # z_ui_idxs = list(range(dataset['depth'].shape[0]))
            z_idxs = list(range(dataset['depth'].shape[0]))  # TODO: subset?
            #salt = np.full((len(z_idxs), len(lat_idxs), len(lon_idxs)), np.nan)
            if progress_bar is True:
                items_iter = tqdm.tqdm(lat_idxs)
                with tqdm_logging_wrapper.wrap_logging_for_tqdm(
                        items_iter), items_iter:
                    #for k, lat_idx in enumerate(items_iter):
                    #    salt[:, k, :] = dataset[self.ncvar][time_idx, z_idxs, lat_idx, lon_idxs]
                    salt = dataset[self.ncvar][time_idx, :, lat_idxs, lon_idxs]
            else:
                #for k, lat_idx in enumerate(lat_idxs):
                #    salt[:, k, :] = dataset[self.ncvar][time_idx, z_idxs, lat_idx, lon_idxs]
                salt = dataset[self.ncvar][time_idx, :, lat_idxs, lon_idxs]
            #change missing value to nan
            idxs = np.where(abs(salt) > 10000)
            salt[idxs] = float('nan')

            loni = dataset['lon'][lon_idxs]
            for idx in range(len(loni)):
                if loni[idx] > 180:
                    loni[idx] = loni[idx]-360.
            lati = dataset['lat'][lat_idxs]
            xi, yi = transform_ll_to_cpp(loni, lati)

            if vgrid.ivcor == 1:
                bz = (hgrid.values[:, None]*vgrid.sigma)[boundary.indexes, :]
                print(f'zcor is {bz[200,:]}')
                idxs = np.where(bz > 5000.0)
                bz[idxs] = 5000.0 - 1.0e-6
            else:
                raise NotImplementedError('vgrid.ivcor!=1')

            xy = hgrid.get_xy(crs='epsg:4326')
            lonb = xy[boundary.indexes, 0]
            latb = xy[boundary.indexes, 1]
            xb, yb = transform_ll_to_cpp(lonb, latb)
            bx = np.tile(xb, [bz.shape[1],1]).T
            by = np.tile(yb, [bz.shape[1],1]).T
            bzyx = np.c_[bz.reshape(np.size(bz)), by.reshape(np.size(by)), bx.reshape(np.size(bx))]
            zi = dataset['depth'][z_idxs]

            # First try with RegularGridInterpolator
            salt_fd = RegularGridInterpolator(
                (zi, yi, xi),
                salt,
                method='linear',
                bounds_error=False,
                fill_value=np.nan
            )
            salt_interp = salt_fd(bzyx)

            # the boundary and the data don't intersect
            if np.all(np.isnan(salt_interp)):
                salt_idxs = np.where(~np.isnan(salt))
                xyzi = np.vstack([
                    np.tile(xi, salt.shape)[salt_idxs].flatten(),
                    np.tile(yi, salt.shape)[salt_idxs].flatten(),
                    np.tile(zi, salt.shape)[salt_idxs].flatten()
                ]).T
                salt_fd = NearestNDInterpolator(xyzi, salt[salt_idxs].flatten())
                salt_interp = salt_fd(np.vstack([bx, by, -bz.flatten()]).T)
            # the boundary and the data partially intersect
            elif np.any(np.isnan(salt_interp)):
                salt_idxs = np.where(~np.isnan(salt_interp))
                salt_fd = NearestNDInterpolator(bzyx[salt_idxs], salt_interp[salt_idxs])
                salt_idxs = np.where(np.isnan(salt_interp))
                salt_interp[salt_idxs] = salt_fd(bzyx[salt_idxs])

            if np.any(np.isnan(salt_interp)):
                raise ValueError('No boundary  salt data for GOFS. '
                                 'Try increasing pixel_buffer argument.')
            dst['time_series'][i, offset:offset + bz.shape[0], :, :] = salt_interp.reshape(bz.shape)


class GOFS(Hycom):
    '''Public interface for GOFS model forcings.'''

    def __init__(self):
        self._elevation = GOFSElevation()
        self._velocity = GOFSVelocity()
        self._temperature = GOFSTemperature()
        self._salinity = GOFSSalinity()

    @property
    def elevation(self) -> GOFSElevation:
        return self._elevation

    @property
    def velocity(self) -> GOFSVelocity:
        return self._velocity

    @property
    def temperature(self) -> GOFSTemperature:
        return self._temperature

    @property
    def salinity(self) -> GOFSSalinity:
        return self._salinity
