from datetime import datetime, timedelta
from functools import lru_cache
import logging
from typing import Dict, Union


from matplotlib.transforms import Bbox
from metpy.units import units
from metpy.calc import potential_temperature, height_to_pressure_std
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata, RegularGridInterpolator, NearestNDInterpolator
import tqdm
import tqdm_logging_wrapper


from pyschism import dates
from pyschism.forcing.baroclinic.base import (
    BaroclinicForcing, BaroclinicComponent)

logger = logging.getLogger(__name__)


class GofsForecastDatasets:

    def __init__(self):
        base_url = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/' +\
                   'expt_93.0/FMRC/runs/GLBy0.08_930_FMRC_RUN_'
        datasets = {}
        for pivot_date in [dates.nearest_zulu() - timedelta(days=i)
                           for i in range(8)]:
            pivot_date += timedelta(hours=12)
            logger.info(
                'Checking for available GOFS data for pivot date: '
                f'{pivot_date}.')
            try:
                datasets.setdefault(
                    pivot_date,
                    Dataset(
                        f'{base_url}' +
                        f'{pivot_date.strftime("%Y-%m-%dT%H:%M:%SZ")}'
                        )
                    )
                logger.info('Success!')
                # break
            except OSError as e:
                logger.info('Not available.')
                if e.errno == -70:
                    print()
                    continue
        self.datasets = datasets

    def datevectors(self):
        for ds in self.values:
            base_date = datetime.strptime(
                ''.join(ds['time'].units.split()[2:4]), '%Y-%m-%d%H:%M:%S.%f')
            yield [base_date+timedelta(hours=x) for x in ds['time'][:]]

    @property
    def values(self):
        return list(self.datasets.values())

    @property
    def start_dates(self):
        return list(self.datasets.keys())


class GOFSBaroclinicComponent(BaroclinicComponent):

    forecast_datasets = GofsForecastDatasets()
    # hindcast_datasets = GofsHindcastDatasets()

    @lru_cache(maxsize=None)
    def get_datasets(
            self,
            start_date: datetime,
            run_days: Union[float, timedelta],
            output_interval=None
    ) -> Dict[datetime, Dataset]:

        if not isinstance(start_date, datetime):
            raise TypeError(
                f'Argument start_date must be of type {datetime}, '
                f'not type {type(start_date)}.')

        if not isinstance(run_days, timedelta):
            run_days = timedelta(days=float(run_days))

        sampling_interval = output_interval if output_interval is not None else self.sampling_interval

        required_start_date = dates.nearest_zulu(start_date) + timedelta(hours=12.)
        required_end_date = dates.nearest_cycle(
            start_date + run_days + sampling_interval,
            period=int(sampling_interval/timedelta(hours=1))
        )
        required_datevector = np.arange(
            required_start_date,
            required_end_date + sampling_interval,
            sampling_interval
        ).astype(datetime)
        datasets = {required_date: None for required_date in required_datevector}
        for required_date in required_datevector:
            for i, datevector in enumerate(self.forecast_datasets.datevectors()):
                if required_date in datevector:
                    datasets[required_date] = self.forecast_datasets.values[i]
                    break

        missing_dates = []
        for date, ds in datasets.items():
            if ds is None:
                missing_dates.append(date)
        if len(missing_dates) > 0:
            raise ValueError(f'No data for dates {missing_dates}\n, got {datasets.keys()}.')

        return datasets

    @property
    def output_interval(self):
        return timedelta(hours=3)

    def _modified_bbox(self, dataset, bbox=None):
        # dataset = list(self.datasets.values())[-1]
        if bbox is None:
            return Bbox.from_extents(
                dataset['lon'][:].min(),
                dataset['lat'][:].min(),
                dataset['lon'][:].max(),
                dataset['lat'][:].max()
            )
        else:
            xmin = bbox.xmin + 360. if not (
                bbox.xmin >= dataset['lon'][:].min()
                and bbox.xmin < 180.) else bbox.xmin

            xmax = bbox.xmax + 360. if not (
                bbox.xmax >= dataset['lon'][:].min()
                and bbox.xmax < 180.) else bbox.xmax

            return Bbox.from_extents(
                np.min([xmin, xmax]),
                bbox.ymin,
                np.max([xmin, xmax]),
                bbox.ymax
            )

    def _modified_bbox_indexes(
            self,
            bbox,
            dataset,
            pixel_buffer=0
    ):
        # dataset = list(self.datasets.values())[-1]
        lat_idxs = np.where((dataset['lat'][:] >= bbox.ymin)
                            & (dataset['lat'][:] <= bbox.ymax))[0]
        lon_idxs = np.where((dataset['lon'][:] >= bbox.xmin)
                            & (dataset['lon'][:] <= bbox.xmax))[0]
        lon_idxs = lon_idxs.tolist()
        lat_idxs = lat_idxs.tolist()
        for i in range(pixel_buffer):
            lon_idxs.insert(0, lon_idxs[0] - 1)
            lon_idxs.append(lon_idxs[-1] + 1)
            lat_idxs.insert(0, lat_idxs[0] - 1)
            lat_idxs.append(lat_idxs[-1] + 1)
        return lon_idxs, lat_idxs


class GOFSElevation(GOFSBaroclinicComponent):

    @property
    def product(self) -> str:
        return 'rtofs_glo_2ds_forecast_3hrly_diag'

    @property
    def ncvar(self):
        return 'surf_el'

    @property
    def sampling_interval(self):
        return timedelta(hours=3.)

    def put_ncdata(self, boundary, dst, start_date, run_days, overwrite=False,
                   offset=0, output_interval=timedelta(hours=24),
                   pixel_buffer=10):
        for i, (time, dataset) in enumerate(
            self.get_datasets(
                start_date,
                run_days,
                output_interval
            ).items()
        ):
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
            zi = np.full((len(lat_idxs), len(lon_idxs)), np.nan)
            items_iter = tqdm.tqdm(lat_idxs)
            with tqdm_logging_wrapper.wrap_logging_for_tqdm(
                    items_iter), items_iter:
                for k, lat_idx in enumerate(items_iter):
                    zi[k, :] = dataset[self.ncvar][
                                        time_idx, lat_idx, lon_idxs]
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
            dst['time_series'][time_idx, offset:len(zq)] = zq


class GOFSVelocity(GOFSBaroclinicComponent):

    @property
    def product(self) -> str:
        return 'rtofs_glo_2ds_forecast_3hrly_diag'

    @property
    def ncvar(self):
        return 'water_u', 'water_v'

    @property
    def sampling_interval(self):
        return timedelta(hours=3.)

    def put_ncdata(self, hgrid, vgrid, boundary, dst, start_date, run_days,
                   overwrite=False, offset=0,
                   output_interval=timedelta(hours=24), pixel_buffer=10):

        for i, (time, dataset) in enumerate(
            self.get_datasets(
                start_date,
                run_days,
                output_interval
            ).items()
        ):
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
            ui = np.full((len(z_idxs), len(lat_idxs), len(lon_idxs)), np.nan)
            vi = np.full((len(z_idxs), len(lat_idxs), len(lon_idxs)), np.nan)
            items_iter = tqdm.tqdm(lat_idxs)
            with tqdm_logging_wrapper.wrap_logging_for_tqdm(
                    items_iter), items_iter:
                for k, lat_idx in enumerate(items_iter):
                    ui[:, k, :] = dataset[uvar][time_idx, z_idxs, lat_idx, lon_idxs]
                    vi[:, k, :] = dataset[vvar][time_idx, z_idxs, lat_idx, lon_idxs]

            xi = dataset['lon'][lon_idxs]
            for idx in range(len(xi)):
                if xi[idx] > 180:
                    xi[idx] = xi[idx]-360.
            yi = dataset['lat'][lat_idxs]

            if vgrid.ivcor == 1:
                bz = (hgrid.values[:, None]*vgrid.sigma)[boundary.indexes, :]
            else:
                raise NotImplementedError('vgrid.ivcor!=1')

            xy = hgrid.get_xy(crs='epsg:4326')
            bx = np.tile(xy[boundary.indexes, 0], (bz.shape[1],))
            by = np.tile(xy[boundary.indexes, 1], (bz.shape[1],))
            bzyx = np.vstack([-bz.flatten(), by, bx]).T
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


class GOFSTemperature(GOFSBaroclinicComponent):

    def put_ncdata(self, hgrid, vgrid, boundary, dst, start_date, run_days,
                   overwrite=False, offset=0,
                   output_interval=timedelta(hours=24), pixel_buffer=10):

        for i, (time, dataset) in enumerate(
            self.get_datasets(
                start_date,
                run_days,
                output_interval
            ).items()
        ):
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
            temp = np.full((len(z_idxs), len(lat_idxs), len(lon_idxs)), np.nan)
            items_iter = tqdm.tqdm(lat_idxs)
            with tqdm_logging_wrapper.wrap_logging_for_tqdm(
                    items_iter), items_iter:
                for k, lat_idx in enumerate(items_iter):
                    temp[:, k, :] = dataset[self.ncvar][time_idx, z_idxs, lat_idx, lon_idxs]

            xi = dataset['lon'][lon_idxs]
            for idx in range(len(xi)):
                if xi[idx] > 180:
                    xi[idx] = xi[idx]-360.
            yi = dataset['lat'][lat_idxs]

            if vgrid.ivcor == 1:
                bz = (hgrid.values[:, None]*vgrid.sigma)[boundary.indexes, :]
            else:
                raise NotImplementedError('vgrid.ivcor!=1')

            xy = hgrid.get_xy(crs='epsg:4326')
            bx = np.tile(xy[boundary.indexes, 0], (bz.shape[1],))
            by = np.tile(xy[boundary.indexes, 1], (bz.shape[1],))
            bzyx = np.vstack([-bz.flatten(), by, bx]).T
            zi = dataset['depth'][z_idxs]

            # First try with RegularGridInterpolator
            temp_fd = RegularGridInterpolator(
                (zi, yi, xi),
                temp,
                method='linear',
                bounds_error=False,
                fill_value=np.nan
            )
            temp_interp = temp_fd(bzyx)

            # the boundary and the data don't intersect
            if np.all(np.isnan(temp_interp)):
                temp_idxs = np.where(~np.isnan(temp))
                xyzi = np.vstack([
                    np.tile(xi, temp.shape)[temp_idxs].flatten(),
                    np.tile(yi, temp.shape)[temp_idxs].flatten(),
                    np.tile(zi, temp.shape)[temp_idxs].flatten()
                ]).T
                temp_fd = NearestNDInterpolator(xyzi, temp[temp_idxs].flatten())
                temp_interp = temp_fd(np.vstack([bx, by, -bz.flatten()]).T)
            # the boundary and the data partially intersect
            elif np.any(np.isnan(temp_interp)):
                temp_idxs = np.where(~np.isnan(temp_interp))
                temp_fd = NearestNDInterpolator(bzyx[temp_idxs], temp_interp[temp_idxs])
                temp_idxs = np.where(np.isnan(temp_interp))
                temp_interp[temp_idxs] = temp_fd(bzyx[temp_idxs])

            if np.any(np.isnan(temp_interp)):
                raise ValueError('No boundary  temperature data for GOFS. '
                                 'Try increasing pixel_buffer argument.')
            pressure = height_to_pressure_std(units('meter')*bz.flatten())
            temp_interp = units('degC')*temp_interp
            temp_interp = potential_temperature(pressure, temp_interp).to('degC')
            dst['time_series'][i, offset:offset+bz.shape[0], :, :] = temp_interp.reshape(bz.shape)

    @property
    def product(self) -> str:
        return 'rtofs_glo_3dz_forecast_daily_temp'

    @property
    def ncvar(self):
        return 'water_temp'

    @property
    def fill_value(self):
        return np.nan

    @property
    def nowcast_varname(self):
        return 'temp'

    @property
    def name(self):
        return 'temperature'


class GOFSSalinity(GOFSBaroclinicComponent):

    def put_ncdata(self, hgrid, vgrid, boundary, dst, start_date, run_days,
                   overwrite=False, offset=0,
                   output_interval=timedelta(hours=24), pixel_buffer=10):

        for i, (time, dataset) in enumerate(
            self.get_datasets(
                start_date,
                run_days,
                output_interval
            ).items()
        ):
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
            # z_ui_idxs = list(range(dataset['depth'].shape[0]))
            z_idxs = list(range(dataset['depth'].shape[0]))  # TODO: subset?
            salt = np.full((len(z_idxs), len(lat_idxs), len(lon_idxs)), np.nan)
            items_iter = tqdm.tqdm(lat_idxs)
            with tqdm_logging_wrapper.wrap_logging_for_tqdm(
                    items_iter), items_iter:
                for k, lat_idx in enumerate(items_iter):
                    salt[:, k, :] = dataset[self.ncvar][time_idx, z_idxs, lat_idx, lon_idxs]

            xi = dataset['lon'][lon_idxs]
            for idx in range(len(xi)):
                if xi[idx] > 180:
                    xi[idx] = xi[idx]-360.
            yi = dataset['lat'][lat_idxs]

            if vgrid.ivcor == 1:
                bz = (hgrid.values[:, None]*vgrid.sigma)[boundary.indexes, :]
            else:
                raise NotImplementedError('vgrid.ivcor!=1')

            xy = hgrid.get_xy(crs='epsg:4326')
            bx = np.tile(xy[boundary.indexes, 0], (bz.shape[1],))
            by = np.tile(xy[boundary.indexes, 1], (bz.shape[1],))
            bzyx = np.vstack([-bz.flatten(), by, bx]).T
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
            dst['time_series'][i, offset:offset+bz.shape[0], :, :] = salt_interp.reshape(bz.shape)

    @property
    def product(self) -> str:
        return 'rtofs_glo_3dz_forecast_daily_salt'

    @property
    def ncvar(self):
        return 'salinity'

    @property
    def fill_value(self):
        return 0.

    @property
    def nowcast_varname(self):
        return 'salt'


class GOFS(BaroclinicForcing):

    def __init__(self):
        self.elevation = GOFSElevation()
        self.velocity = GOFSVelocity()
        self.temperature = GOFSTemperature()
        self.salinity = GOFSSalinity()
