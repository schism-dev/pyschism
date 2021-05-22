from datetime import datetime, timedelta
from functools import lru_cache
from typing import Dict, Union

import logging
import pathlib


from matplotlib.transforms import Bbox
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata


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

        for date, ds in datasets.items():
            if ds is None:
                raise ValueError(f'No data for date {date}.')

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

    def write(self, path, hgrid, start_date, run_days, overwrite=False):

        path = pathlib.Path(path)
        if path.exists() and overwrite is not True:
            raise IOError('File exists and overwrite is not True.')

        nOpenBndNodes = 0
        for boundary in hgrid.boundaries.ocean().itertuples():
            nOpenBndNodes += len(boundary.indexes)
        output_interval = timedelta(days=1)

        with Dataset(path, 'w', format='NETCDF4') as dst:

            # dimensions
            dst.createDimension('nOpenBndNodes', nOpenBndNodes)
            dst.createDimension('one', 1)
            dst.createDimension('time', None)
            dst.createDimension('nComponents', 1)

            # variables
            dst.createVariable('time', 'f', ('time',))
            dst.createVariable('time_series', 'f',
                               ('time', 'nOpenBndNodes', 'nComponents'))
            dst.createVariable('time_step', 'f', ('one',))
            dst['time_step'][:] = int(output_interval.total_seconds())
            for i, (time, dataset) in enumerate(
                self.get_datasets(
                    start_date, run_days, output_interval).items()
            ):
                logger.info(
                    'Saving GOFS surf_el data for date: '
                    f'{start_date+i*output_interval} '
                    f'approximated as {time}.'
                    )
                ds_base_date = datetime.strptime(
                    ''.join(dataset['time'].units.split()[2:-1]),
                    '%Y-%m-%d%H:%M:%S.%f')
                ds_timevector = [ds_base_date + timedelta(hours=x)
                                 for x in dataset['time'][:]]
                time_idx = ds_timevector.index(time)
                offset = 0
                for j, boundary in enumerate(hgrid.boundaries.ocean().to_crs(
                        'epsg:4326').itertuples()):
                    bbox = self._modified_bbox(
                        dataset, Bbox.from_extents(*boundary.geometry.bounds))
                    lon_idxs, lat_idxs = self._modified_bbox_indexes(
                            bbox,
                            dataset,
                            pixel_buffer=2
                        )
                    zi = np.full((len(lat_idxs), len(lon_idxs)), np.nan)
                    for k, lat_idx in enumerate(lat_idxs):
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
                    if np.any(np.isnan(zq)):
                        raise ValueError('Boundary contains NaNs.')
                    dst['time_series'][time_idx, offset+j:len(zq)] = zq
                    offset += len(zq)


def read_vgrid(fname):
    fid = open(fname, 'r')
    lines = fid.readlines()
    fid.close()

    ivcor = int(lines[0].strip().split()[0])
    nvrt = int(lines[1].strip().split()[0])

    if ivcor == 1:
        lines = lines[2:]
        kbp = np.array([int(i.split()[1])-1 for i in lines])
        NP = len(kbp)
        print(NP)
        sigma = -np.ones([NP, nvrt])
        for i, line in enumerate(lines):
            sigma[i, kbp[i]:] = np.array(line.strip().split()[2:]).astype('float')
    elif ivcor == 2:
        kz, h_s = lines[1].strip().split()[1:3]
        kz = int(kz)
        h_s = float(kz)

        # read z grid
        ztot = []
        irec = 2
        for i in np.arange(kz):
            irec = irec+1
            ztot.append(lines[irec].strip().split()[1])
        ztot = np.array(ztot).astype('float')

        # read s grid
        sigma = []
        irec = irec+2
        nsigma = nvrt-kz+1
        h_c, theta_b, theta_f = np.array(lines[irec].strip().split()[:3]).astype('float')
        for i in np.arange(nsigma):
            irec = irec+1
            sigma.append(lines[irec].strip().split()[1])
        sigma = np.array(sigma).astype('float')

    return sigma


class GOFSVelocity(GOFSBaroclinicComponent):

    def write(self, path, hgrid, vgrid, start_date, run_days, overwrite=False):

        path = pathlib.Path(path)
        if path.exists() and overwrite is not True:
            raise IOError('File exists and overwrite is not True.')

        nOpenBndNodes = 0
        for boundary in hgrid.boundaries.ocean().itertuples():
            nOpenBndNodes += len(boundary.indexes)
        output_interval = timedelta(days=1)
        import tempfile
        tmpdir = tempfile.TemporaryDirectory()
        vgrid.write(tmpdir.name+'/vgrid.in')
        sigma = read_vgrid(tmpdir.name+'/vgrid.in')
        depth = hgrid.values
        idxs = np.where(depth < 0.11)
        depth[idxs] = 0.11
        breakpoint()
        zcor = depth[:, None] * vgrid.sigma
        nvrt = zcor.shape[1]
        print(zcor)
        print(nvrt)
        exit()
        with Dataset(path, 'w', format='NETCDF4') as dst:

            # dimensions
            dst.createDimension('nOpenBndNodes', nOpenBndNodes)
            dst.createDimension('one', 1)
            dst.createDimension('time', None)
            dst.createDimension('nComponents', 1)
            dst.createDimension('nLevels', nvrt)

            # variables
            dst.createVariable('time', 'f', ('time',))
            dst.createVariable('time_series', 'f',
                               ('time', 'nOpenBndNodes', 'nComponents'))
            dst.createVariable('time_step', 'f', ('one',))
            dst['time_step'][:] = int(output_interval.total_seconds())
            for i, (time, dataset) in enumerate(
                self.get_datasets(
                    start_date, run_days, output_interval).items()
            ):
                logger.info(
                    'Saving GOFS surf_el data for date: '
                    f'{start_date+i*output_interval} '
                    f'approximated as {time}.'
                    )
                ds_base_date = datetime.strptime(
                    ''.join(dataset['time'].units.split()[2:-1]),
                    '%Y-%m-%d%H:%M:%S.%f')
                ds_timevector = [ds_base_date + timedelta(hours=x)
                                 for x in dataset['time'][:]]
                time_idx = ds_timevector.index(time)
                offset = 0
                for j, boundary in enumerate(hgrid.boundaries.ocean().to_crs(
                        'epsg:4326').itertuples()):
                    bbox = self._modified_bbox(
                        dataset, Bbox.from_extents(*boundary.geometry.bounds))
                    lon_idxs, lat_idxs = self._modified_bbox_indexes(
                            bbox,
                            dataset,
                            pixel_buffer=2
                        )
                    zi = np.full((len(lat_idxs), len(lon_idxs)), np.nan)
                    for k, lat_idx in enumerate(lat_idxs):
                        zi[k, :] = dataset['surf_el'][
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
                    if np.any(np.isnan(zq)):
                        raise ValueError('Boundary contains NaNs.')
                    dst['time_series'][time_idx, offset+j:len(zq)] = zq
                    offset += len(zq)


class GOFSTemperature(GOFSBaroclinicComponent):

    @property
    def product(self) -> str:
        return 'rtofs_glo_3dz_forecast_daily_temp'

    @property
    def ncvar(self):
        return 'temperature'

    @property
    def fill_value(self):
        return np.nan

    @property
    def nowcast_varname(self):
        return 'temp'


class GOFSSalinity(GOFSBaroclinicComponent):

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

    # def write(
    #         self,
    #         path,
    #         outdir,
    #         start_date,
    #         run_days,

    # ):

