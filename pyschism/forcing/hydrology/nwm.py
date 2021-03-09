from collections import defaultdict
from datetime import datetime, timedelta
import logging
import multiprocessing
import pathlib
import tempfile
import warnings
import tarfile
import tempfile
from time import time
from typing import Union
import urllib

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
from shapely.geometry import LinearRing, Point, MultiPoint, LineString
import wget

from pyschism.dates import localize_datetime, nearest_cycle_date, pivot_time
from pyschism.mesh import Hgrid, Gr3
from pyschism.forcing.hydrology.base import Hydrology, Sources, Sinks


DATADIR = pathlib.Path(user_data_dir('nwm'))
DATADIR.mkdir(exist_ok=True, parents=True)
NWM_FILE = DATADIR / 'NWM_channel_hydrofabric.tar.gz'

logger = logging.getLogger(__name__)


class NWMElementPairings:

    def __init__(self, hgrid):

        logger.info('Computing NWMElementPairings...')
        self._hgrid = hgrid

        # An STR-Index returns the reaches that are near the boundaries of the
        # mesh. This subsamples the NWM network, but also is not the exact
        # result. This is used to speed-up computations by filtering the input
        # data.
        logger.info('Computing r_index.')
        start = time()
        nwm_r_index = self.gdf.sindex
        logger.info(f'Computing r_index took {time() - start}.')

        # The r-index is used to find intersections between mesh boundary edges
        # and NWM reaches (approximate results)
        logger.info('Use r_index to filter features.')
        start = time()
        possible_indexes = set()
        for edge in hgrid.hull.edges().itertuples():
            for index in list(nwm_r_index.intersection(edge.geometry.bounds)):
                possible_indexes.add(index)
        possible_matches = self.gdf.iloc[list(possible_indexes)]
        logger.info(f'Filtering features took {time()-start}.')
        del possible_indexes
        del nwm_r_index

        # The hull rings itersections is used to find the exact NWM reaches
        # that intersect the mesh's hull.
        logger.info('Finding exact features intersections.')
        start = time()
        exact_indexes = set()
        for pm in possible_matches.itertuples():
            if hgrid.hull.rings().geometry.intersects(pm.geometry).any():
                exact_indexes.add(pm.Index)
        reaches = self.gdf.iloc[list(exact_indexes)]
        logger.info(f'Finding exact features took {time()-start}.')

        # release some memory
        del possible_matches
        del exact_indexes
        del self._gdf

        logger.info('Pairing features to corresponding element.')

        # Pair each reach with corresponding element.
        # 1) Find reach-hull intersections.
        start = time()
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
                _, idx = tree.query(point, workers=-1)
                element_index[inters.reachIndex].append(idx)
        del tree
        logger.info(
            'Pairing features to corresponding element took '
            f'{time()-start}.')

        start = time()
        elements = hgrid.elements.geodataframe()
        logger.info(
            'Generating mesh-element geodataframe took: '
            f'{time()-start}.')
        del self._hgrid  # release

        start = time()
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
                            sources[element.id].add(reach.feature_id)
                        else:
                            sinks[element.id].add(reach.feature_id)
                        break
        logger.info(
            'Sorting features into sources and sinks took: '
            f'{time()-start}.')
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
        hgrid.transform_to(gpd.read_file(self.nwm_file, rows=1, layer=0).crs)
        self.__hgrid = hgrid

    @_hgrid.deleter
    def _hgrid(self):
        del self.__hgrid

    @property
    def nwm_file(self):
        if not hasattr(self, '_nwm_file'):
            self._tmpdir = tempfile.TemporaryDirectory()
            with tarfile.open(NWM_FILE, "r:gz") as src:
                src.extractall(self._tmpdir.name)
            self._nwm_file = (
                pathlib.Path(self._tmpdir.name) /
                'NWM_v2.0_channel_hydrofabric/nwm_v2_0_hydrofabric.gdb'
            )
        return self._nwm_file

    @property
    def gdf(self):
        if not hasattr(self, '_gdf'):
            bbox = self._hgrid.bbox
            self._gdf = gpd.read_file(
                self.nwm_file,
                bbox=(bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymax),
                # bbox=(-75.889435, 38.895308, -74.604034, 39.477546),  # Delaware Bay, debug  # noqa: E501
                layer=0,
                )
        return self._gdf


class NWMDataGetter:

    def __init__(
            self,
            pairings: NWMElementPairings,
            start_date: datetime,
            rnday: timedelta
    ):
        self._pairings = pairings
        # naïve-datetime
        self._start_date = localize_datetime(start_date).astimezone(pytz.utc)
        self._rnday = rnday


def streamflow_lookup(file, pairings):
    nc = Dataset(file)
    streamflow = nc['streamflow'][:]
    feature_id = nc['feature_id'][:]
    sources = []
    # TODO: read scaling factor directly from netcdf file?
    for features in pairings.sources.values():
        in_file = np.in1d(feature_id, list(features), assume_unique=True)
        sources.append(0.01*np.sum(streamflow[np.where(in_file)]))
    sinks = []
    for features in pairings.sinks.values():
        in_file = np.in1d(feature_id, list(features), assume_unique=True)
        sinks.append(-0.01*np.sum(streamflow[np.where(in_file)]))
    return sources, sinks


class AWSDataInventory:

    def __init__(
            self,
            start_date: datetime = None,
            rnday: Union[int, float, timedelta] = timedelta(days=5.),
            product='medium_range_mem1',
            verbose=False,
            fallback=True,
    ):
        """This will download the latest National Water Model data.

        NetCDF files are saved to the system's temporary directory.
        The AWS data goes back 30 days. For requesting hindcast data from
        before we need a different data source
        """
        self.start_date = nearest_cycle_date() if start_date is None \
            else localize_datetime(start_date).astimezone(pytz.utc)
        self.rnday = rnday if isinstance(rnday, timedelta) \
            else timedelta(days=rnday)
        self.product = product
        self.fallback = fallback

        # if the model start_date aligns to a  "zero" with the NWM data, then
        # fetching the data is trivial
        self._files = {dt: None for dt in self.timevector}
        if self.start_date == nearest_cycle_date(self.start_date):
            self._fetch_data()

        # if they don't align then we need to inject a "zero" entry at the
        # beginning. My suggestion is to repeat the first value. Program will
        # raise if that is the case, so we can address this special case later.
        # TODO: Put a "zero" entry if start_date and NWM dates do not align.
        else:
            raise NotImplementedError(
                f'Model start_date={str(self.start_date)} is not a "pivot" '
                'time.')

    @property
    def bucket(self):
        return 'noaa-nwm-pds'

    @property
    def nearest_cycle_date(self) -> datetime:
        return nearest_cycle_date(self.start_date)

    @property
    def output_interval(self) -> timedelta:
        return {
            'medium_range_mem1': timedelta(hours=3)
        }[self.product]

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

    @property
    def timevector(self):
        return np.arange(
                    self.start_date,
                    self.start_date + self.rnday + self.output_interval,
                    self.output_interval
                ).astype(datetime)

    def _fetch_data(self):

        # this needs to be checked if there is no "zero" alignment
        requested_time = pivot_time(self.start_date)

        # This download here could be more robust. Right now it tries the
        # "nowcast" first and since it is expected to fail (NWM is 6 hours
        # offset from nowcast) then immediately downloads the previous nowcast.
        # This has been kept like this as a reminder that the "nowcast" never
        # exists for NWM, but this is not true for other models with shorter
        # lags (for example GFS).
        nwm_time = requested_time.strftime('%Y%m%d')
        res = self.s3.list_objects_v2(
            Bucket=self.bucket,
            Delimiter='/',
            Prefix=f'nwm.{nwm_time}/{self.product}/'
                   )

        if 'Contents' in res:
            # contents will be empty when t00z is called.
            data = list(reversed(sorted([
                data['Key'] for data in res['Contents'] if 'channel' in data['Key']
                ])))
        else:
            data = []

        # In reality the previous will always fail, so we need to get the pivot
        # time
        nwm_time = (requested_time - timedelta(days=1)).strftime('%Y%m%d')
        res = self.s3.list_objects_v2(
            Bucket=self.bucket,
            Delimiter='/',
            Prefix=f'nwm.{nwm_time}/{self.product}/'
                   )

        data.extend(list(reversed(sorted([
            data['Key'] for data in res['Contents'] if 'channel' in data['Key']
            ]))))

        nearest_cycle = int(6 * np.floor(requested_time.hour/6))
        previous_cycle = (nearest_cycle - 6) % 24

        if f't{nearest_cycle:02d}z' not in data[0] \
                and f't{previous_cycle:02d}z' in data[0]:
            if self.fallback is True:
                warnings.warn(
                        f'NWM data for cycle t{nearest_cycle:02d}z is not yet '
                        'on the server, defaulting to previous cycle.')
            else:
                raise IOError(
                    'Unknown error while fetching NWM data.')

        for d in data:

            base_date_str = f'{d.split("/")[0].split(".")[-1]}'
            timedelta_str = d.split(
                'channel_rt_1.')[-1].split('.')[0].strip('f')
            file_datetime = datetime.strptime(base_date_str, '%Y%m%d') \
                + timedelta(hours=int(d.split('.')[2].strip('tz'))) \
                + timedelta(hours=float(timedelta_str))

            if file_datetime in self._files:
                file = self._files[file_datetime]
                if file is None:
                    filename = pathlib.Path(self.tmpdir.name) / d
                    filename.parent.mkdir(parents=True, exist_ok=True)
                    with open(filename, 'wb') as f:
                        logger.info(f'Downloading file {d}.')
                        self.s3.download_fileobj(self.bucket, d, f)
                    self._files[file_datetime] = filename

        for dt, data in self._files.items():
            if data is None:
                raise IOError(f'No NWM data for time {str(dt)}.')

    def __call__(self, pairings: NWMElementPairings, h0=1e-1, nprocs=-1):
        logger.info('Launching streamflow lookup...')
        start = time()
        with multiprocessing.Pool(
                processes=cpu_count() if nprocs == -1 else nprocs
        ) as pool:
            res = pool.starmap(
                streamflow_lookup, [(file, pairings) for file in self.files])
        pool.join()
        logger.info(f'streamflow lookup took {time()-start}...')

        logger.info('Adding streamflow data as sources and sinks...')
        start = time()
        sources = Sources()
        sinks = Sinks()
        for i, file in enumerate(self.files):
            nc = Dataset(file)
            _time = pytz.timezone('UTC').localize(
                datetime.strptime(
                    nc.model_output_valid_time,
                    "%Y-%m-%d_%H:%M:%S"))
            # TODO: This is slow, it might change if add_data is vectorized.
            for j, element_id in enumerate(pairings.sources.keys()):
                sources.add_data(_time, element_id, res[i][0][j], -9999, 0.)
            for k, element_id in enumerate(pairings.sinks.keys()):
                sinks.add_data(_time, element_id, res[i][1][k])
        logger.info(f'Adding source/sink data took {time()-start}.')
        return sources, sinks

    @property
    def files(self):
        return sorted(list(pathlib.Path(self.tmpdir.name).glob('**/*.nc')))

    @property
    def requested_product(self):
        return {
            'medium_range_mem1': 'medium_range.channel_rt_1'
        }[self.product]


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
            logger.warning(
                "Downloading National Water Model stream network file to "
                "the pyschism cache...")
            try:
                wget.download(
                    'http://www.nohrsc.noaa.gov/pub/staff/keicher/NWM_live/'
                    'web/data_tools/NWM_channel_hydrofabric.tar.gz',
                    out=str(self._nwm_file))
            except urllib.error.HTTPError as e:
                logger.fatal(
                    'Could not download NWM_channel_hydrofabric.tar.gz')
                raise e

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
        logger.info('NationalWaterModel.__call__')
        pairings = NWMElementPairings(model_driver.model_domain.hgrid)
        start_date = model_driver.param.opt.start_date
        rnday = model_driver.param.core.rnday

        # forecast
        if start_date >= pivot_time() - timedelta(days=30): 
            logger.info('Fetching NWM data.')
            AWSData = AWSDataGetter(pairings, start_date, rnday)
            sources, sinks = AWSData(h0=h0, nprocs=nprocs)

        # hindcast
        else:
            raise NotImplementedError('Hindcast is not implemented 30 days.')

        # turn 'on' the source/sink system in SCHISM.
        model_driver.param.opt.if_source = 1

        # set the ramps if applicable, no ramp by default.
        if int(nramp_ss) != 0:
            # nramp_ss = 1 # needed if if_source=1; ramp-up flag for
            # source/sinks
            model_driver.param.opt.nramp_ss = nramp_ss
            # dramp_ss = 2 # needed if if_source=1; ramp-up period in days
            if dramp_ss is not None:
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
