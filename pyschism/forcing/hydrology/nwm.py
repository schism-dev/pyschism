from collections import defaultdict
from datetime import datetime, timedelta
import logging
import multiprocessing
import pathlib
import tempfile
import warnings
from typing import Union

from appdirs import user_data_dir
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import geopandas as gpd
from netCDF4 import Dataset
import numpy as np
from psutil import cpu_count
import pytz
from scipy.spatial import cKDTree
from shapely import ops
from shapely.geometry import (
    LinearRing, Point, MultiPoint, LineString)
# from tqdm import tqdm
import wget

from pyschism.mesh import Hgrid, Gr3
from pyschism.forcing.hydrology.base import Hydrology, Sources, Sinks


DATADIR = pathlib.Path(user_data_dir()) / 'nwm'
DATADIR.mkdir(exist_ok=True, parents=True)
NWM_FILE = DATADIR / 'nwm_v12.gdb.zip'

_logger = logging.getLogger(__name__)


class NWMGeoDataFrame:

    def __get__(self, obj, val):
        gdf = obj.__dict__.get('gdf')
        if gdf is None:
            bbox = obj._hgrid.bbox
            gdf = gpd.read_file(
                NWM_FILE,
                bbox=(bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymax),
                # bbox=(-75.889435, 38.895308, -74.604034, 39.477546),  # Delaware Bay, debug  # noqa: E501
                )
            obj.__dict__['gdf'] = gdf
        return gdf

    def __delete__(self, obj):
        del obj.__dict__['gdf']


class NWMElementPairings:

    _gdf = NWMGeoDataFrame()

    def __init__(self, hgrid):
        _logger.info('Initiliaze NWMElementPairings')
        _logger.debug('This debug message should also appear.')
        self._hgrid = hgrid

        # An STR-Index returns the reaches that are near the boundaries of the
        # mesh. This subsamples the NWM network, but also is not the exact
        # result. This is used to speed-up computations by filtering the input
        # data.
        nwm_r_index = self._gdf.sindex

        # The r-index is used to find intersections between mesh boundary edges
        # and NWM reaches (approximate results)
        possible_indexes = set()
        for edge in hgrid.hull.edges().itertuples():
            for index in list(nwm_r_index.intersection(edge.geometry.bounds)):
                possible_indexes.add(index)
        possible_matches = self._gdf.iloc[list(possible_indexes)]
        del possible_indexes
        del nwm_r_index

        # The hull rings itersections is used to find the exact NWM reaches
        # that intersect the mesh's hull.
        exact_indexes = set()
        for pm in possible_matches.itertuples():
            if hgrid.hull.rings().geometry.intersects(pm.geometry).any():
                exact_indexes.add(pm.Index)
        reaches = self._gdf.iloc[list(exact_indexes)]

        # release some memory
        del possible_matches
        del exact_indexes
        del self._gdf

        # Pair each reach with corresponding element.
        # 1) Find reach-hull intersections.
        data = []
        intersections: gpd.GeoDataFrame
        for i, reach in enumerate(reaches.itertuples()):
            for ring in hgrid.hull.rings().itertuples():
                if ring.geometry.intersects(reach.geometry):
                    data.append({
                        "geometry": ring.geometry.intersection(reach.geometry),
                        "reachIndex": i})
                    break

        intersections = gpd.GeoDataFrame(data, crs=hgrid.crs)
        del data

        # 2) Generate element centroid KDTree
        centroids = []
        for element in hgrid.elements().values():
            cent = LinearRing(hgrid.nodes.coord()[list(
                map(hgrid.nodes.get_index_by_id, element))]
            ).centroid
            centroids.append((cent.x, cent.y))
        tree = cKDTree(centroids)
        del centroids

        # 3) Match reach/boundary intersection to nearest element centroid
        element_index = {}
        for inters in intersections.itertuples():
            geom = inters.geometry
            if isinstance(geom, Point):
                geom = MultiPoint(geom.coords)
            element_index[inters.reachIndex] = []
            for point in geom:
                _, idx = tree.query(point)
                element_index[inters.reachIndex].append(idx)
        del tree

        elements = hgrid.elements.geodataframe()
        del self._hgrid  # release

        sources = defaultdict(set)
        sinks = defaultdict(set)
        for reach_index, paired_elements_idxs in element_index.items():
            reach = reaches.iloc[reach_index]
            point_of_intersection = intersections.loc[
                intersections["reachIndex"] == reach_index]
            for element_idx in paired_elements_idxs:
                element = elements.iloc[element_idx]
                if not isinstance(reach.geometry, LineString):
                    geom = ops.linemerge(reach.geometry)
                else:
                    geom = reach.geometry
                for segment in map(LineString, zip(geom.coords[:-1],
                                                   geom.coords[1:])):
                    if segment.intersects(
                            point_of_intersection.iloc[0].geometry.buffer(
                                np.finfo(np.float32).eps)):
                        downstream = segment.coords[-1]
                        if element.geometry.contains(Point(downstream)):
                            sources[element.id].add(reach.featureID)
                        else:
                            sinks[element.id].add(reach.featureID)
                        break
        self._sources = sources
        self._sinks = sinks

    @property
    def sources(self):
        return self._sources

    @property
    def sinks(self):
        return self._sinks

    @property
    def _hgrid(self):
        return self.__hgrid

    @_hgrid.setter
    def _hgrid(self, hgrid: Gr3):
        hgrid = Hgrid(**hgrid.to_dict())
        hgrid.transform_to(gpd.read_file(NWM_FILE, rows=1).crs)
        self.__hgrid = hgrid

    @_hgrid.deleter
    def _hgrid(self):
        del self.__hgrid


class NWMDataGetter:

    def __init__(
            self,
            pairings: NWMElementPairings,
            start_date: datetime,
            rnday: timedelta
    ):
        self._pairings = pairings
        # naïve-datetime
        if start_date.tzinfo is None or \
                start_date.tzinfo.utcoffset(start_date) is None:
            start_date = pytz.timezone('utm').localize(start_date)
        self._start_date = start_date
        self._rnday = rnday


def streamflow_lookup(file, pairings):
    nc = Dataset(file)
    sources = []
    # TODO: read scaling factor directly from netcdf file?
    for element_id, features in pairings.sources.items():
        sources.append(0.01*np.sum(nc['streamflow'][
            np.where(np.isin(nc['feature_id'], list(features)))]))
    sinks = []
    for element_id, features in pairings.sinks.items():
        sinks.append(-0.01*np.sum(nc['streamflow'][
            np.where(np.isin(nc['feature_id'], list(features)))]))
    return (sources, sinks)


def pivot_time(input_datetime=None, period=6):
    if input_datetime is None:
        input_datetime = pytz.timezone('UTC').localize(datetime.utcnow())
    current_cycle = int(period * np.floor(input_datetime.hour / period))
    return pytz.timezone('UTC').localize(
        datetime(input_datetime.year, input_datetime.month,
                 input_datetime.day, current_cycle))


class AWSDataInventory:

    def __init__(
            self,
            start_date: datetime,
            rnday: Union[int, float, timedelta] = None,
            product='medium_range_mem1',
            verbose=False,
            fallback=True,
    ):
        """This will download the latest National Water Model data.

        NetCDF files are saved to the system's temporary directory.
        The AWS data goes back 30 days. For requesting hindcast data from
        before we need a different data source
        """
        _logger.info('Initialize AWSDataInventory')
        self.start_date = start_date
        self.rnday = rnday

        requested_time = pivot_time(start_date) - timedelta(days=1)

        while requested_time <= self.pivot_time:
            self._get_data_from_bucket(requested_time, product, fallback)
            requested_time += timedelta(days=1)

    @property
    def Bucket(self):
        return 'noaa-nwm-pds'

    @property
    def pivot_time(self):
        return pivot_time()

    @property
    def s3(self):
        try:
            return self._s3
        except AttributeError:
            self._s3 = boto3.client(
                's3', config=Config(signature_version=UNSIGNED))
            return self._s3

    @property
    def tmpdir(self):
        try:
            return self._tmpdir
        except AttributeError:
            self._tmpdir = tempfile.TemporaryDirectory()
            return self._tmpdir

    def _get_data_from_bucket(self, requested_time, product, fallback):
        end_date = self.start_date + self.rnday
        nwm_time = requested_time.strftime('%Y%m%d')
        res = self.s3.list_objects_v2(Bucket=self.Bucket, Delimiter='/',
                                      Prefix=f'nwm.{nwm_time}/{product}/')

        if requested_time == pivot_time() and res['KeyCount'] == 0:
            if fallback is True:
                warnings.warn('NWM data is late, defaulting to previous.')
                requested_time -= timedelta(days=1)
                nwm_time = requested_time.strftime('%Y%m%d')
                res = self.s3.list_objects_v2(
                    Bucket=self.Bucket,
                    Delimiter='/',
                    Prefix=f'nwm.{nwm_time}/{product}/'
                )
            else:
                raise IOError(
                    'NWM data is "late", no NWM-data for current cycle.')
        for data in res['Contents']:
            if 'channel' in data['Key']:
                filename = pathlib.Path(
                    self.tmpdir.name) / data['Key'].split('/')[-1]
                with open(filename, 'wb') as f:
                    self.s3.download_fileobj(self.Bucket, data['Key'], f)
                nc = Dataset(filename)
                # TODO: Will crash if the user localized with timedelta object.
                if end_date.tzinfo.localize(
                    datetime.strptime(nc.model_output_valid_time,
                                      '%Y-%m-%d_%H:%M:%S')) > end_date:
                    break

    def __call__(self, pairings: NWMElementPairings, h0=1e-1, nprocs=-1):
        _logger.info('Will pair NWM data to elements...')
        from time import time as _time
        start = _time()
        files = sorted(list(pathlib.Path(self.tmpdir.name).glob('*')))
        with multiprocessing.Pool(
                processes=cpu_count() if nprocs == -1 else nprocs
        ) as pool:
            res = pool.starmap(
                streamflow_lookup, [(file, pairings) for file in files])
        pool.join()
        sources = Sources()
        sinks = Sinks()
        for i, file in enumerate(files):
            nc = Dataset(file)
            time = pytz.timezone('UTC').localize(
                datetime.strptime(
                    nc.model_output_valid_time,
                    "%Y-%m-%d_%H:%M:%S"))
            # TODO: This is slow, it might change if add_data is vectorized.
            for j, element_id in enumerate(pairings.sources.keys()):
                sources.add_data(time, element_id, res[i][0][j], -9999, 0.)
            for k, element_id in enumerate(pairings.sinks.keys()):
                sinks.add_data(time, element_id, res[i][1][k])
        _logger.info(f'Done pairing, took {_time()-start} seconds...')
        return sources, sinks


class AWSDataGetter(NWMDataGetter):

    def __call__(self, product: str = 'medium_range_mem1',
                 verbose: bool = False, h0=1e-1, nprocs=-1):
        """
        We just picked up a server_config. This part must be slurm-aware and
        potentially
        """
        return AWSDataInventory(
            start_date=self._start_date,
            rnday=self._rnday,
            product=product,
            verbose=False
        )(self._pairings, h0, nprocs)


class FTPDataGetter(NWMDataGetter):

    def __call__(self):
        raise NotImplementedError('FTPDataGetter')


class NationalWaterModel(Hydrology):

    def __init__(self):
        self._nwm_file = NWM_FILE
        if not self._nwm_file.exists():
            _logger.warning(
                "Downloading National Water Model stream network file to "
                "the pyschism cache...")
            wget.download(
                "https://www.dropbox.com/s/3w8i46uumbcs49v/nwm_v12.gdb.zip"
                "?dl=1", out=str(self._nwm_file))
        super().__init__()

    def __call__(self, model_driver, nramp_ss: bool = False, dramp_ss=None,
                 h0=1e-1, nprocs=-1):
        """Initializes the NWM data.
        Used by :class:`pyschism.driver.ModelDriver`

        Will pick the "best" data source based on start date and rnday.
        The are fringe cases not yet covered, for example when the data spans
        more than 1 data source.
        """
        super().__call__(model_driver)
        _logger.info('NationalWaterModel.__call__')
        pairings = NWMElementPairings(model_driver.model_domain.hgrid)
        start_date = model_driver.param.opt.start_date
        rnday = model_driver.param.core.rnday
        if start_date >= pivot_time() - timedelta(days=30):
            _logger.info('Fetching NWM data.')
            AWSData = AWSDataGetter(pairings, start_date, rnday)
            sources, sinks = AWSData(h0=h0, nprocs=nprocs)

        else:
            raise NotImplementedError(
                'start_date is less than pivot_time - 30 days ')

        # turn 'on' the source/sink system in SCHISM.
        model_driver.param.opt.if_source = 1

        # set the ramps if applicable, no ramp by default.
        if int(nramp_ss) != 0:
            # nramp_ss = 1 # needed if if_source=1; ramp-up flag for
            # source/sinks
            model_driver.param.opt.nramp_ss = nramp_ss
            # dramp_ss = 2 # needed if if_source=1; ramp-up period in days
        if int(nramp_ss) != 0 and dramp_ss is not None:
            model_driver.param.opt.dramp_ss = dramp_ss

    @staticmethod
    def from_files(msource, vsource, vsink):
        raise NotImplementedError
        # TODO: File parsers.
        nwm = NationalWaterModel()
        nwm._msource = None
        nwm._vsource = None
        nwm._vsink = None
        return nwm
