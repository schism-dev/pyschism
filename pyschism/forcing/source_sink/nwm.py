from collections import defaultdict
from datetime import datetime, timedelta
import logging
from multiprocessing import Pool, cpu_count
import pathlib

import tarfile
import tempfile
from time import time
from typing import Union
import urllib
import appdirs
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import fiona
import geopandas as gpd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from shapely import ops
from shapely.geometry import LinearRing, Point, MultiPoint, LineString, box
import wget

from pyschism import dates
from pyschism.mesh.hgrid import Hgrid
from pyschism.mesh.base import Gr3

from pyschism.forcing.source_sink.base import SourceSink

DATADIR = pathlib.Path(appdirs.user_data_dir("pyschism/nwm"))
DATADIR.mkdir(exist_ok=True, parents=True)

logger = logging.getLogger(__name__)


class NWMElementPairings:
    def __init__(self, hgrid, nwm_file=None):

        self._nwm_file = nwm_file

        logger.info("Computing NWMElementPairings...")
        self._hgrid = hgrid

        # An STR-Index returns the reaches that are near the boundaries of the
        # mesh. This subsamples the NWM network, but also is not the exact
        # result. This is used to speed-up computations by filtering the input
        # data.
        logger.info("Computing r_index.")
        start = time()
        nwm_r_index = self.gdf.sindex
        logger.info(f"Computing r_index took {time() - start}.")

        # The r-index is used to find intersections between mesh boundary edges
        # and NWM reaches (approximate results)
        logger.info("Use r_index to filter features.")
        start = time()
        possible_indexes = set()
        for edge in hgrid.hull.edges().itertuples():
            for index in list(nwm_r_index.intersection(edge.geometry.bounds)):
                possible_indexes.add(index)
        possible_matches = self.gdf.iloc[list(possible_indexes)]
        logger.info(f"Filtering features took {time()-start}.")
        del possible_indexes
        del nwm_r_index

        # The hull rings itersections is used to find the exact NWM reaches
        # that intersect the mesh's hull.
        logger.info("Finding exact features intersections.")
        start = time()
        exact_indexes = set()
        for pm in possible_matches.itertuples():
            if hgrid.hull.rings().geometry.intersects(pm.geometry).any():
                exact_indexes.add(pm.Index)
        reaches = self.gdf.iloc[list(exact_indexes)]
        logger.info(f"Finding exact features took {time()-start}.")

        # release some memory
        del possible_matches
        del exact_indexes
        del self._gdf

        logger.info("Pairing features to corresponding element.")

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
                    data.append({"geometry": _intersections, "reachIndex": i})
                    break

        if len(data) == 0:
            # TODO: change for warning in future.
            raise IOError("No National Water model intersections found on the mesh.")
        intersection = gpd.GeoDataFrame(data, crs=hgrid.crs)
        del data

        # 2) Generate element centroid KDTree
        centroids = []
        for element in hgrid.elements.elements.values():
            cent = LinearRing(
                hgrid.nodes.coord[list(map(hgrid.nodes.get_index_by_id, element))]
            ).centroid
            centroids.append((cent.x, cent.y))
        tree = cKDTree(centroids)
        del centroids

        # 3) Match reach/boundary intersection to nearest element centroid

        coords = [
            np.array(inters.geometry.coords) for inters in intersection.itertuples()
        ]
        _, idxs = tree.query(np.vstack(coords), workers=-1)

        element_index = defaultdict(list)
        for i, idx in enumerate(idxs):
            element_index[intersection.iloc[i].reachIndex].append(idx)
        del tree
        logger.info(
            "Pairing features to corresponding element took " f"{time()-start}."
        )

        hull = hgrid.hull.multipolygon()
        # del self._hgrid  # release

        start = time()
        sources = defaultdict(set)
        sinks = defaultdict(set)
        for reach_index, paired_elements_idxs in element_index.items():
            reach = reaches.iloc[reach_index]
            point_of_intersection = intersection.loc[
                intersection["reachIndex"] == reach_index
            ]
            for element_idx in paired_elements_idxs:
                element = hgrid.elements.gdf.iloc[element_idx]
                if not isinstance(reach.geometry, LineString):
                    geom = ops.linemerge(reach.geometry)
                else:
                    geom = reach.geometry
                for segment in map(LineString, zip(geom.coords[:-1], geom.coords[1:])):
                    if segment.intersects(
                        point_of_intersection.iloc[0].geometry.buffer(
                            np.finfo(np.float32).eps
                        )
                    ):
                        downstream = segment.coords[-1]
                        if (
                            box(*segment.bounds)
                            .intersection(hull)
                            .intersects(Point(downstream))
                        ):
                            sources[element.id].add(reach.feature_id)
                        else:
                            sinks[element.id].add(reach.feature_id)

        logger.info("Sorting features into sources and sinks took: " f"{time()-start}.")

        self.sources = sources
        self.sinks = sinks

    def make_plot(self):
        # verification plot
        data = []
        egdf = self.hgrid.elements.geodataframe()
        for eid, features in self.sources.items():
            # for feat in features:
            eidx = self.hgrid.elements.get_index_by_id(eid)
            data.append({
                    'geometry': egdf.iloc[eidx].geometry
                })
        src_gdf = gpd.GeoDataFrame(data)
        data = []
        for eid, features in self.sinks.items():
            # for feat in features:
            eidx = self._hgrid.elements.get_index_by_id(eid)
            data.append({
                    'geometry': egdf.iloc[eidx].geometry
                })
        snk_gdf = gpd.GeoDataFrame(data)

        ax = egdf.plot(facecolor="none",  edgecolor='black', lw=0.7)
        src_gdf.plot(color='red', ax=ax, alpha=0.5)
        snk_gdf.plot(color='blue', ax=ax, alpha=0.5)
        plt.show()

    @property
    def sources_gdf(self):
        if not hasattr(self, "_sources_gdf"):
            data = []
            for eid, features in self.sources.items():
                eidx = self.hgrid.elements.get_index_by_id(eid)
                data.append(
                    {
                        "element_id": eid,
                        "geometry": self.hgrid.elements.gdf.loc[eidx].geometry,
                        "features": list(features),
                    }
                )
            self._sources_gdf = gpd.GeoDataFrame(data)
        return self._sources_gdf

    @property
    def sinks_gdf(self):
        if not hasattr(self, "_sinks_gdf"):
            data = []
            for eid, features in self.sinks.items():
                eidx = self.hgrid.elements.get_index_by_id(eid)
                data.append(
                    {
                        "element_id": eid,
                        "geometry": self.hgrid.elements.gdf.loc[eidx].geometry,
                        "features": list(features),
                    }
                )
            self._sinks_gdf = gpd.GeoDataFrame(data)
        return self._sinks_gdf

    @property
    def hgrid(self):
        return self._hgrid

    # @property
    # def _hgrid(self):
    #     return self.__hgrid

    # @_hgrid.setter
    # def _hgrid(self, hgrid: Gr3):
    #     hgrid = Hgrid(**hgrid.to_dict())
    #     hgrid.transform_to(
    #         gpd.read_file(self.nwm_file, rows=1, layer=0).crs
    #         # 'epsg:4269',
    #         )
    #     self.__hgrid = hgrid

    # @_hgrid.deleter
    # def _hgrid(self):
    #     del self.__hgrid

    # @property
    # def nwm_file(self):

    #     # NOTE: There's a tradeoff between time spent decompressing the file vs. the memory it takes to store it.

    #     if not hasattr(self, "_nwm_file"):
    #         self._tmpdir = tempfile.TemporaryDirectory()
    #         with tarfile.open(NWM_FILE, "r:gz") as src:
    #             src.extractall(self._tmpdir.name)
    #         self._nwm_file = (
    #             list(pathlib.Path(self._tmpdir.name).glob('**/*.gdb'))[0]
    #         )
    #     return self._nwm_file

    @property
    def gdf(self):
        if not hasattr(self, "_gdf"):
            gdf_coll = []
            for reach_layer in [reach_layer for reach_layer in fiona.listlayers(self.nwm_file) if'reaches' in reach_layer]:
                layer_crs = gpd.read_file(self.nwm_file, rows=1, layer=reach_layer).crs
                bbox = self.hgrid.get_bbox(crs=layer_crs)
                gdf_coll.append(
                    gpd.read_file(
                        self.nwm_file,
                        bbox=(bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymax),
                        layer=reach_layer
                        )
                )
            self._gdf = pd.concat(gdf_coll)
        return self._gdf

    @property
    def nwm_file(self):
        return self._nwm_file

    @property
    def _nwm_file(self):
        return self.__nwm_file

    @_nwm_file.setter
    def _nwm_file(self, nwm_file):
        nwm_file = list(DATADIR.glob('**/*hydrofabric*.gdb')) if nwm_file is None else nwm_file
        if isinstance(nwm_file, list):
            if len(nwm_file) == 0:
                tmpdir = tempfile.TemporaryDirectory()
                logger.info(
                    f"Downloading National Water Model stream network tar file to {tmpdir}"
                )
                try:
                    wget.download(
                        "https://www.nohrsc.noaa.gov/pub/staff/keicher/NWM_live/web/data_tools/NWM_channel_hydrofabric.tar.gz",
                        out=tmpdir.name,
                        bar=wget.bar_adaptive if logger.getEffectiveLevel() < 30 else None,
                    )
                except urllib.error.HTTPError as e:
                    logger.fatal("Could not download NWM_channel_hydrofabric.tar.gz")
                    raise e
                tmpfile = list(pathlib.Path(tmpdir.name).glob("**/*.tar.gz"))[0]
                with tarfile.open(tmpfile, "r:gz") as src:
                    logger.info(
                        f"Extracting National Water Model stream network tar file to {DATADIR}"
                    )

                    src.extractall(DATADIR)
                nwm_file = list(DATADIR.glob('**/*.gdb'))[0]
            elif len(nwm_file) == 1:
                nwm_file = nwm_file[0]
            else:
                raise Exception('Found more than 1 NWM hydrofabric file.')
        self.__nwm_file = nwm_file       


def streamflow_lookup(file, indexes):
    nc = Dataset(file)
    streamflow = nc["streamflow"][:]
    data = []
    # TODO: read scaling factor directly from netcdf file?
    for indxs in indexes:
        # Dataset already considered scale_factor and offset
        data.append(np.sum(streamflow[indxs]))
    return data


class AWSForecatsInventory:
    def __init__(
        self,
        start_date: datetime = None,
        rnday: Union[int, float, timedelta] = timedelta(days=5.0),
        product=None,
        verbose=False,
        fallback=True,
    ):
        """This will download the latest National Water Model data.

        NetCDF files are saved to the system's temporary directory.
        The AWS data goes back 30 days. For requesting hindcast data from
        before we need a different data source
        """
        self.start_date = (
            dates.nearest_cycle()
            if start_date is None
            else dates.nearest_cycle(dates.localize_datetime(start_date))
        )
        # self.start_date = self.start_date.replace(tzinfo=None)
        self.rnday = rnday if isinstance(rnday, timedelta) else timedelta(days=rnday)
        self.product = "medium_range_mem1" if product is None else product
        self.fallback = fallback
        self._files = {
            _: None
            for _ in np.arange(
                self.start_date,
                self.start_date + self.rnday + self.output_interval,
                self.output_interval,
            ).astype(datetime)
        }
        # print(self._files)

        self.data = self.s3.list_objects_v2(
            Bucket=self.bucket,
            Delimiter="/",
            Prefix=f'nwm.{self.start_date.strftime("%Y%m%d")}' f"/{self.product}/",
        )

        for requested_time, _ in self._files.items():
            logger.info(f"Requesting NWM data for time {requested_time}")
            self._files[requested_time] = self.request_data(requested_time)

    def request_data(self, request_time):

        file_metadata = list(
            reversed(
                sorted(
                    [_["Key"] for _ in self.data["Contents"] if "channel" in _["Key"] and 'Contents' in self.data]
                )
            )
        )

        for key in file_metadata[0:240:3]:
            if request_time != self.key2date(key):
                continue
            # print(self.key2date(key))
            # print(request_time)
            filename = pathlib.Path(self.tmpdir.name) / key
            filename.parent.mkdir(parents=True, exist_ok=True)

            with open(filename, "wb") as f:
                logger.info(f"Downloading file {key}, ")
                self.s3.download_fileobj(self.bucket, key, f)
            return filename

    def key2date(self, key):
        base_date_str = f'{key.split("/")[0].split(".")[-1]}'
        timedelta_str = key.split("channel_rt_1.")[-1].split(".")[0].strip("f")
        return (
            datetime.strptime(base_date_str, "%Y%m%d")
            + timedelta(hours=float(timedelta_str))
            # + timedelta(hours=int(key.split('.')[2].strip('tz')))
            # + timedelta(hours=float(timedelta_str))
        )

    def get_nc_pairing_indexes(self, pairings: NWMElementPairings):
        nc_feature_id = Dataset(self.files[0])["feature_id"][:]

        def get_aggregated_features(features):
            aggregated_features = []
            for source_feats in features:
                aggregated_features.extend(list(source_feats))
            in_file = []
            for feature in aggregated_features:
                idx = np.where(nc_feature_id == int(feature))[0]
                in_file.append(idx.item())
            in_file_2 = []
            sidx = 0
            for source_feats in features:
                eidx = sidx + len(source_feats)
                in_file_2.append(in_file[sidx:eidx])
                sidx = eidx
            return in_file_2

        sources = get_aggregated_features(pairings.sources.values())
        sinks = get_aggregated_features(pairings.sinks.values())
        return sources, sinks

    @property
    def bucket(self):
        return "noaa-nwm-pds"

    @property
    def nearest_cycle(self) -> datetime:
        return dates.nearest_cycle(self.start_date)

    @property
    def output_interval(self) -> timedelta:
        return {"medium_range_mem1": timedelta(hours=3)}[self.product]

    @property
    def s3(self):
        try:
            return self._s3
        except AttributeError:
            self._s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
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
            self.output_interval,
        ).astype(datetime)

    @property
    def files(self):
        return sorted(list(pathlib.Path(self.tmpdir.name).glob("**/*.nc")))

    @property
    def requested_product(self):
        return {"medium_range_mem1": "medium_range.channel_rt_1"}[self.product]


class AWSHindcastInventory:

    def __init__(
            self,
            start_date: datetime = None,
            rnday: Union[int, float, timedelta] = timedelta(days=5.),
            product=None,
            verbose=False,
            fallback=True,
    ):
        """This will download the National Water Model retro data.
        A 26-year (January 1993 through December 2018) retrospective 
        simulation using version 2.0 of the NWM. 

        NetCDF files are saved to the system's temporary directory.
        """
        self.start_date = dates.nearest_cycle() if start_date is None \
            else dates.nearest_cycle(dates.localize_datetime(start_date))
        # self.start_date = self.start_date.replace(tzinfo=None)
        self.rnday = rnday if isinstance(rnday, timedelta) \
            else timedelta(days=rnday)
        self.product = 'CHRTOUT_DOMAIN1.comp' if product is None else product
        self.fallback = fallback
        self._files = {_: None for _ in np.arange(
            self.start_date,
            self.start_date + self.rnday + self.output_interval,
            self.output_interval
        ).astype(datetime)}

        paginator = self.s3.get_paginator('list_objects_v2')
        pages = paginator.paginate(
            Bucket=self.bucket,
            Prefix=f'full_physics/{self.start_date.year}'
        )

        self.data = []
        for page in pages:
            for obj in page['Contents']:
                self.data.append(obj)

        self.file_metadata = list(sorted([
            _['Key'] for _ in self.data if 'CHRTOUT_DOMAIN1.comp' in _['Key']
        ]))

        timevector = np.arange(
            datetime(self.start_date.year, 1, 1),
            datetime(self.start_date.year+1, 1, 1),
            np.timedelta64(1, 'h'),
            dtype='datetime64'
        )

        timefile = {pd.to_datetime(str(timevector[i])): self.file_metadata[i] for i in range(len(timevector))}
 
        for requested_time, _ in self._files.items():
            key = timefile.get(requested_time)
            logger.info(f'Requesting NWM data for time {requested_time}')

            self._files[requested_time] = self.request_data(key, requested_time)

    def request_data(self, key, request_time):

        filename = pathlib.Path(self.tmpdir.name) / key
        filename.parent.mkdir(parents=True, exist_ok=True)

        with open(filename, 'wb') as f:
            logger.info(f'Downloading file {key}, ')
            self.s3.download_fileobj(self.bucket, key, f)
        return filename

    def get_nc_pairing_indexes(self, pairings: NWMElementPairings):
        nc_feature_id = Dataset(self.files[0])['feature_id'][:]

        def get_aggregated_features(features):
            aggregated_features = []
            for source_feats in features:
                aggregated_features.extend(list(source_feats))
            in_file = []
            for feature in aggregated_features:
                idx = np.where(nc_feature_id == int(feature))[0]
                in_file.append(idx.item())
            in_file_2 = []
            sidx = 0
            for source_feats in features:
                eidx = sidx + len(source_feats)
                in_file_2.append(in_file[sidx:eidx])
                sidx = eidx
            return in_file_2

        sources = get_aggregated_features(pairings.sources.values())
        sinks = get_aggregated_features(pairings.sinks.values())
        return sources, sinks

    @property
    def bucket(self):
        return 'noaa-nwm-retro-v2.0-pds'

    @property
    def nearest_cycle(self) -> datetime:
        return dates.nearest_cycle(self.start_date)

    @property
    def output_interval(self) -> timedelta:
        return {
            'CHRTOUT_DOMAIN1.comp': timedelta(hours=1)
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
    def files(self):
        return sorted(list(pathlib.Path(self.tmpdir.name).glob('**/*.comp')))


class AWSDataInventory:

    def __new__(
        cls,
        start_date,
        rnday,
        product=None,
        verbose=False,
        fallback=True
    ):
        # AWSHindcastInventory -> January 1993 through December 2018
        if start_date >= dates.localize_datetime(datetime(1993, 1, 1, 0, 0)) \
                and start_date + rnday <= dates.localize_datetime(datetime(2018, 12, 31, 23, 59)):
            return AWSHindcastInventory(start_date, rnday, product, verbose, fallback)

        elif start_date >= dates.nearest_zulu() - timedelta(days=30):
            return AWSForecatsInventory(start_date, rnday, product, verbose, fallback)

        else:
            raise Exception(f'No NWM model data for start_date {start_date} and end_date {start_date+rnday}.')


class NationalWaterModel(SourceSink):

    start_date = dates.StartDate()
    end_date = dates.EndDate()
    # run_days = dates.RunDays()

    def __init__(self, aggregation_radius=None, pairings=None, nwm_file=None):
        self.nwm_file = nwm_file
        self.aggregation_radius = aggregation_radius
        self.pairings = pairings

    def write(
        self,
        output_directory,
        gr3: Gr3,
        start_date: datetime = None,
        end_date: Union[datetime, timedelta] = None,
        overwrite: bool = False,
        nprocs=-1,
        pairings=None,
        product=None,
    ):
        nprocs = -1 if nprocs is None else nprocs
        nprocs = cpu_count() if nprocs == -1 else nprocs
        self.start_date = start_date
        self.end_date = end_date
        pairings = self.pairings if pairings is None else pairings
        self.pairings = NWMElementPairings(gr3, nwm_file=self.nwm_file) if pairings is None else pairings
        self._inventory = AWSDataInventory(
            start_date=self.start_date,
            rnday=self.end_date - self.start_date,
            # product="medium_range_mem1",
        )

        self._timevector = [dates.localize_datetime(d) for d in self.inventory._files]

        src_idxs, snk_idxs = self.inventory.get_nc_pairing_indexes(self.pairings)
        with Pool(processes=nprocs) as pool:
            sources = pool.starmap(
                streamflow_lookup, [(file, src_idxs) for file in self.inventory.files]
            )
            sinks = pool.starmap(
                streamflow_lookup, [(file, snk_idxs) for file in self.inventory.files]
            )
        pool.join()
        self._data = {}
        for i, file in enumerate(self.inventory.files):
            nc = Dataset(file)
            _time = dates.localize_datetime(
                datetime.strptime(nc.model_output_valid_time, "%Y-%m-%d_%H:%M:%S")
            )
            for j, element_id in enumerate(self.pairings.sources.keys()):
                self._data.setdefault(_time, {}).setdefault(element_id, {}).update(
                    {
                        "flow": sources[i][j],
                        "temperature": -9999,
                        "salinity": 0.,
                    }
                )
            for k, element_id in enumerate(self.pairings.sinks.keys()):
                self._data.setdefault(_time, {}).setdefault(element_id, {}).update({"flow": -sinks[i][k]})

        if self.aggregation_radius is not None:
            self.aggregate_by_radius(gr3, self.aggregation_radius)

    @property
    def timevector(self):
        return self._timevector

    @property
    def inventory(self):
        return self._inventory
