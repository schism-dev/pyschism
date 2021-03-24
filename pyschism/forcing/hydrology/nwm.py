from collections import defaultdict
from datetime import datetime, timedelta
import logging
from multiprocessing import Pool
import pathlib
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
from shapely.geometry import LinearRing, Point, MultiPoint, LineString, box
import wget

from pyschism.dates import localize_datetime, nearest_cycle_date, pivot_time
from pyschism.mesh import Hgrid, Gr3
from pyschism.forcing.hydrology.base import Hydrology


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
        intersection: gpd.GeoDataFrame
        for i, reach in enumerate(reaches.itertuples()):
            for ring in hgrid.hull.rings().itertuples():
                if ring.geometry.intersects(reach.geometry):
                    _intersections = ring.geometry.intersection(reach.geometry)
                    if isinstance(_intersections, MultiPoint):
                        # features with even number of intersections are
                        # discarded.
                        if len(_intersections.geoms) % 2 == 0:
                            continue
                        # if a feature has an odd number of intersections > 1,
                        # we only take the last point.
                        _intersections = _intersections.geoms[-1]
                    data.append({
                        "geometry": _intersections,
                        "reachIndex": i})
                    break

        intersection = gpd.GeoDataFrame(data, crs=hgrid.crs)
        del data

        # 2) Generate element centroid KDTree
        centroids = []
        for element in hgrid.elements.elements.values():
            cent = LinearRing(hgrid.nodes.coord[list(
                map(hgrid.nodes.get_index_by_id, element))]
            ).centroid
            centroids.append((cent.x, cent.y))
        tree = cKDTree(centroids)
        del centroids

        # 3) Match reach/boundary intersection to nearest element centroid

        coords = [np.array(inters.geometry.coords)
                  for inters in intersection.itertuples()]
        _, idxs = tree.query(np.vstack(coords), workers=-1)

        element_index = defaultdict(list)
        for i, idx in enumerate(idxs):
            element_index[intersection.iloc[i].reachIndex].append(idx)
        del tree
        logger.info(
            'Pairing features to corresponding element took '
            f'{time()-start}.')

        hull = hgrid.hull.multipolygon()
        # del self._hgrid  # release

        start = time()
        sources = defaultdict(set)
        sinks = defaultdict(set)
        for reach_index, paired_elements_idxs in element_index.items():
            reach = reaches.iloc[reach_index]
            point_of_intersection = intersection.loc[
                intersection["reachIndex"] == reach_index]
            for element_idx in paired_elements_idxs:
                element = hgrid.elements.gdf.iloc[element_idx]
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
                        if box(*segment.bounds).intersection(hull).intersects(
                                Point(downstream)):
                            sources[element.id].add(reach.feature_id)
                        else:
                            sinks[element.id].add(reach.feature_id)

        logger.info(
            'Sorting features into sources and sinks took: '
            f'{time()-start}.')

        # verification plot
        # data = []
        # for eid, features in sources.items():
        #     # for feat in features:
        #     eidx = self._hgrid.elements.get_index_by_id(eid)
        #     data.append({
        #             'geometry': elements.iloc[eidx].geometry
        #         })
        # src_gdf = gpd.GeoDataFrame(data)
        # data = []
        # for eid, features in sinks.items():
        #     # for feat in features:
        #     eidx = self._hgrid.elements.get_index_by_id(eid)
        #     data.append({
        #             'geometry': elements.iloc[eidx].geometry
        #         })
        # snk_gdf = gpd.GeoDataFrame(data)
        # import matplotlib.pyplot as plt
        # ax = elements.plot(facecolor="none",  edgecolor='black', lw=0.7)
        # src_gdf.plot(color='red', ax=ax, alpha=0.5)
        # snk_gdf.plot(color='blue', ax=ax, alpha=0.5)
        # plt.show()

        self.sources = sources
        self.sinks = sinks

    @property
    def sources_gdf(self):
        if not hasattr(self, '_sources_gdf'):
            data = []
            for eid, features in self.sources.items():
                eidx = self.hgrid.elements.get_index_by_id(eid)
                data.append(
                    {
                        'element_id': eid,
                        'geometry':
                            self.hgrid.elements.gdf.loc[eidx].geometry,
                        'features': list(features),
                    })
            self._sources_gdf = gpd.GeoDataFrame(data)
        return self._sources_gdf

    @property
    def sinks_gdf(self):
        if not hasattr(self, '_sinks_gdf'):
            data = []
            for eid, features in self.sinks.items():
                eidx = self.hgrid.elements.get_index_by_id(eid)
                data.append(
                    {
                        'element_id': eid,
                        'geometry':
                            self.hgrid.elements.gdf.loc[eidx].geometry,
                        'features': list(features),

                    })
            self._sinks_gdf = gpd.GeoDataFrame(data)
        return self._sinks_gdf

    @property
    def hgrid(self):
        return self._hgrid

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
            bbox = self._hgrid.get_bbox(crs='EPSG:4326', output_type='bbox')
            self._gdf = gpd.read_file(
                self.nwm_file,
                bbox=(bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymax),
                # bbox=(-75.889435, 38.895308, -74.604034, 39.477546),  # Delaware Bay, debug  # noqa: E501
                layer=0,
                )
        return self._gdf


def streamflow_lookup(file, indexes):
    nc = Dataset(file)
    streamflow = nc['streamflow'][:]
    data = []
    # TODO: read scaling factor directly from netcdf file?
    for indxs in indexes:
        data.append(0.01*np.sum(streamflow[indxs]))
    return data


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

    def get_nc_pairing_indexes(self, pairings: NWMElementPairings):
        nc_feature_id = Dataset(self.files[0])['feature_id'][:]

        def get_aggregated_features(features):
            aggregated_features = []
            for source_feats in features:
                aggregated_features.extend(list(source_feats))
            in_file = np.where(
                np.in1d(nc_feature_id, aggregated_features,
                        assume_unique=True))[0]
            in_file_2 = []
            sidx = 0
            for source_feats in features:
                eidx = sidx + len(source_feats)
                in_file_2.append(in_file[sidx:eidx].tolist())
                sidx = eidx
            return in_file_2

        sources = get_aggregated_features(pairings.sources.values())
        sinks = get_aggregated_features(pairings.sinks.values())
        return sources, sinks

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
                data['Key'] for data in res['Contents'] if 'channel'
                in data['Key']])))
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

    @property
    def files(self):
        return sorted(list(pathlib.Path(self.tmpdir.name).glob('**/*.nc')))

    @property
    def requested_product(self):
        return {
            'medium_range_mem1': 'medium_range.channel_rt_1'
        }[self.product]


class NationalWaterModel(Hydrology):

    def __init__(self, aggregation_radius=None):
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
        self.aggregation_radius = aggregation_radius

    def __call__(self, model_driver, nramp_ss: bool = False, dramp_ss=None,
                 nprocs=-1):
        """Initializes the NWM data.
        Used by :class:`pyschism.driver.ModelDriver`

        Will pick the "best" data source based on start date and rnday.
        The are fringe cases not yet covered, for example when the data spans
        more than 1 data source.
        """
        super().__init__(
            model_driver.param.opt.start_date,
            model_driver.param.core.rnday
        )
        logger.info('Initializing NationalWaterModel data.')
        pairings = NWMElementPairings(model_driver.model_domain.hgrid)

        # start_date = model_driver.param.opt.start_date
        # rnday = model_driver.param.core.rnday

        # forecast
        if self.start_date >= pivot_time() - timedelta(days=30):
            logger.info('Fetching NWM data.')
            # (self._pairings, h0, nprocs)
            inventory = AWSDataInventory(
                start_date=self.start_date,
                rnday=self.rnday,
                product='medium_range_mem1',
                verbose=False
            )

            logger.info('Launching streamflow lookup...')
            source_indexes, sinks_indexes = inventory.get_nc_pairing_indexes(
                pairings)
            start = time()
            with Pool(processes=cpu_count()) as pool:
                sources = pool.starmap(
                    streamflow_lookup,
                    [(file, source_indexes) for file in inventory.files])
                sinks = pool.starmap(
                    streamflow_lookup,
                    [(file, sinks_indexes) for file in inventory.files])
            pool.join()
            logger.info(f'streamflow lookup took {time()-start}...')

            # Pass NWM data to Hydrology class.
            logging.info('Generating per-element hydrologic timeseries...')
            start = time()
            hydro = Hydrology(
                self.start_date,
                self.rnday
            )
            for i, file in enumerate(inventory.files):
                nc = Dataset(file)
                _time = localize_datetime(datetime.strptime(
                    nc.model_output_valid_time,
                    "%Y-%m-%d_%H:%M:%S"))
                for j, element_id in enumerate(pairings.sources.keys()):
                    hydro.add_data(_time, element_id, sources[i][j], -9999, 0.)
                for k, element_id in enumerate(pairings.sinks.keys()):
                    hydro.add_data(_time, element_id, -sinks[i][k])
            logging.info(
                'Generating per-element hydrologic timeseries took '
                f'{time() - start}.')

        # hindcast
        else:
            raise NotImplementedError('Hindcast is not implemented 30 days.')

        # aggregate timeseries
        if self.aggregation_radius is not None:
            aggregation_radius = float(self.aggregation_radius)
            logging.info(
                'Aggregating hydrology timeseries/elements using a '
                f'radius of {aggregation_radius} meters.')
            start = time()
            hydro.aggregate_by_radius(model_driver.model_domain.hgrid,
                                      aggregation_radius)
            logging.info(
                f'Aggregating NWM elements took {time() - start}.')

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
