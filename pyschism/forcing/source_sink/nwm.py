from abc import ABC, abstractmethod
from collections import defaultdict
from datetime import datetime, timedelta
import json
import logging
from multiprocessing import Pool, cpu_count
import os
import pathlib
import posixpath
import shutil

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
from pyschism.mesh.base import Gr3

from pyschism.forcing.source_sink.base import SourceSink, Sources, Sinks

DATADIR = pathlib.Path(appdirs.user_data_dir("pyschism/nwm"))
DATADIR.mkdir(exist_ok=True, parents=True)

logger = logging.getLogger(__name__)


class NWMElementPairings:
    def __init__(self, hgrid: Gr3, nwm_file=None, workers=-1):

        # TODO: Accelerate using dask: https://blog.dask.org/2017/09/21/accelerating-geopandas-1

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
                    intersections = ring.geometry.intersection(reach.geometry)
                    if isinstance(intersections, MultiPoint):
                        for point in intersections.geoms:
                            data.append({"geometry": point, "reachIndex": i})
                        continue

                    data.append({"geometry": intersections, "reachIndex": i})
                    break

        if len(data) == 0:
            # TODO: change for warning in future.
            raise IOError(
                "No National Water model intersections found on the mesh.")
        intersection = gpd.GeoDataFrame(data, crs=hgrid.crs)
        #TODO: add exporting intersection as an option
        #intersection.to_file('intersections.shp')
        del data

        # 2) Generate element centroid KDTree
        centroids = []
        for element in hgrid.elements.elements.values():
            cent = LinearRing(
                hgrid.nodes.coord[list(
                    map(hgrid.nodes.get_index_by_id, element))]
            ).centroid
            centroids.append((cent.x, cent.y))
        tree = cKDTree(centroids)
        del centroids

        # 3) Match reach/boundary intersection to nearest element centroid
        coords = [
            np.array(inters.geometry.coords) for inters in intersection.itertuples()
        ]
        _, idxs = tree.query(np.vstack(coords), workers=workers)
        del tree

        logger.info(
            "Pairing features to corresponding element took " f"{time()-start}."
        )

        hull = hgrid.hull.multipolygon()

        start = time()
        sources = defaultdict(list)
        sinks = defaultdict(list)
        for row in intersection.itertuples():
            poi = row.geometry
            reach = reaches.iloc[row.reachIndex].geometry
            if not isinstance(reach, LineString):
                reach = ops.linemerge(reach)
            for segment in map(
                LineString, zip(reach.coords[:-1], reach.coords[1:])
            ):
                if segment.intersects(poi.buffer(np.finfo(np.float32).eps)):
                    segment_origin = Point(segment.coords[0])
                    d1 = segment_origin.distance(poi)
                    downstream = segment.interpolate(
                        d1 + np.finfo(np.float32).eps)
                    element = hgrid.elements.gdf.iloc[idxs[row.Index]]
                    if (
                        box(*LineString([poi, downstream]).bounds)
                        .intersection(hull)
                        .intersects(downstream)
                    ):
                        sources[element.id].append(reaches.iloc[row.reachIndex].feature_id)
                    else:
                        sinks[element.id].append(reaches.iloc[row.reachIndex].feature_id)
                    break

        logger.info(
            "Sorting features into sources and sinks took: " f"{time()-start}.")
        self.sources = sources
        self.sinks = sinks

    def make_plot(self):
        # verification plot
        data = []
        egdf = self.hgrid.elements.gdf
        for eid in self.sources.keys():
            eidx = self.hgrid.elements.get_index_by_id(eid)
            data.append({"geometry": egdf.iloc[eidx].geometry})
        src_gdf = gpd.GeoDataFrame(data)
        data = []
        for eid in self.sinks.keys():
            eidx = self._hgrid.elements.get_index_by_id(eid)
            data.append({"geometry": egdf.iloc[eidx].geometry})
        snk_gdf = gpd.GeoDataFrame(data)

        ax = egdf.plot(facecolor="none", edgecolor="black", lw=0.7)
        src_gdf.plot(color="red", ax=ax, alpha=0.5)
        snk_gdf.plot(color="blue", ax=ax, alpha=0.5)
        plt.show()

    def save_json(self, sources=True, sinks=True):

        sources = "sources.json" if sources is True else sources
        if sources:
            logger.info(f"Saving {sources}")
            with open(sources, "w") as fh:
                json.dump(self.sources, fh)

        sinks = "sinks.json" if sinks is True else sinks
        if sinks:
            with open(sinks, "w") as fh:
                logger.info(f"Saving {sinks}")
                json.dump(self.sinks, fh)

    @staticmethod
    def load_json(hgrid, sources=None, sinks=None):
        pairings = NWMElementPairings.__new__(NWMElementPairings)
        logger.info(f"Loading pairing sources: {sources}")
        pairings.sources = json.load(
            open(sources)) if sources is not None else {}
        logger.info(f"Loading pairing sinks: {sinks}")
        pairings.sinks = json.load(open(sinks)) if sinks is not None else {}
        pairings._hgrid = hgrid
        return pairings

    @property
    def sources_gdf(self):
        if not hasattr(self, "_sources_gdf"):
            data = []
            for eid, features in self.sources.items():
                eidx = self.hgrid.elements.get_index_by_id(eid)
                data.append(
                    {
                        "element_id": eid,
                        "geometry": LineString(
                            self.hgrid.elements.gdf.loc[eidx].geometry.exterior.coords
                        ),
                        "features": " ".join(list(map(str, features))),
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
                        "geometry": LineString(
                            self.hgrid.elements.gdf.loc[eidx].geometry.exterior.coords
                        ),
                        "features": " ".join(list(map(str, features))),
                    }
                )
            self._sinks_gdf = gpd.GeoDataFrame(data)
        return self._sinks_gdf

    @property
    def hgrid(self):
        return self._hgrid

    @property
    def gdf(self):
        if not hasattr(self, "_gdf"):
            gdf_coll = []
            for reach_layer in [
                reach_layer
                for reach_layer in fiona.listlayers(self.nwm_file)
                if "reaches" in reach_layer
            ]:
                layer_crs = gpd.read_file(
                    self.nwm_file, rows=1, layer=reach_layer).crs
                bbox = self.hgrid.get_bbox(crs=layer_crs)
                gdf_coll.append(
                    gpd.read_file(
                        self.nwm_file,
                        bbox=(bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymax),
                        layer=reach_layer,
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
        nwm_file = (
            list(DATADIR.glob("**/*hydrofabric*.gdb")
                 ) if nwm_file is None else nwm_file
        )
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
                        bar=wget.bar_adaptive
                        if logger.getEffectiveLevel() < 30
                        else None,
                    )
                except urllib.error.HTTPError as e:
                    logger.fatal(
                        "Could not download NWM_channel_hydrofabric.tar.gz")
                    raise e
                tmpfile = list(pathlib.Path(
                    tmpdir.name).glob("**/*.tar.gz"))[0]
                with tarfile.open(tmpfile, "r:gz") as src:
                    logger.info(
                        f"Extracting National Water Model stream network tar file to {DATADIR}"
                    )

                    def is_within_directory(directory, target):
                        
                        abs_directory = os.path.abspath(directory)
                        abs_target = os.path.abspath(target)
                    
                        prefix = os.path.commonprefix([abs_directory, abs_target])
                        
                        return prefix == abs_directory
                    
                    def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
                    
                        for member in tar.getmembers():
                            member_path = os.path.join(path, member.name)
                            if not is_within_directory(path, member_path):
                                raise Exception("Attempted Path Traversal in Tar File")
                    
                        tar.extractall(path, members, numeric_owner) 
                        
                    
                    safe_extract(src, DATADIR)
                nwm_file = list(DATADIR.glob("**/*.gdb"))[0]
            elif len(nwm_file) == 1:
                nwm_file = nwm_file[0]
            else:
                raise Exception("Found more than 1 NWM hydrofabric file.")
        self.__nwm_file = nwm_file

def get_aggregated_features(nc_feature_id, features):
    aggregated_features = []
    for source_feats in features:
        aggregated_features.extend(list(source_feats))

    in_file=[]
    for feature in aggregated_features:
        idx=np.where(nc_feature_id == int(feature))[0]
        in_file.append(idx.item())

    in_file_2 = []
    sidx = 0
    for source_feats in features:
        eidx = sidx + len(source_feats)
        #in_file_2.append(in_file[sidx:eidx].tolist())
        in_file_2.append(in_file[sidx:eidx])
        sidx = eidx
    return in_file_2

def streamflow_lookup(file, indexes, threshold=-1e-5):
    nc = Dataset(file)
    streamflow = nc["streamflow"][:]
    streamflow[np.where(streamflow < threshold)] = 0.0
    #change masked value to zero
    streamflow[np.where(streamflow.mask)] = 0.0
    data = []
    for indxs in indexes:
        # Note: Dataset already consideres scale factor and offset.
        data.append(np.sum(streamflow[indxs]))
    nc.close()
    return data

class AWSDataInventory(ABC):
    def __new__(
        cls, start_date, rnday, product=None, verbose=False, fallback=True, cache=None
    ):
        # AWSHindcastInventory
        # The latest(as of 12/21/2021) AWSHindcast dataset covers from Feb 1979 through Dec 2020
        if start_date >= dates.localize_datetime(
            datetime(1979, 2, 1, 0, 0)
        ) and start_date + rnday <= dates.localize_datetime(
            datetime(2020, 12, 31, 23, 59)
        ):
            return AWSHindcastInventory.__new__(cls)

        # GOOGLEHindcastInventory -> January 2019 through 30 days earlier than today
        # data before April 20, 2021 (including Apirl 20) is 3-hr interval, after that is hourly
        # Missing data 20211231, 20220101, 20220102, 20220103
        elif start_date >= dates.localize_datetime(
            datetime(2021, 1, 1, 0, 0)
        ) and start_date + rnday < dates.nearest_zulu() - timedelta(days=2):
            return GOOGLEHindcastInventory.__new__(cls)

        elif start_date >= dates.nearest_zulu() - timedelta(days=2):
            return AWSForecastInventory.__new__(cls)

        else:
            raise Exception(
                f"No NWM model data for start_date {start_date} and end_date {start_date+rnday}."
            )

    #@abstractmethod
    #def request_data(self, request_time):
    #    raise NotImplementedError

    #@property
    #@abstractmethod
    #def bucket(self):
    #    raise NotImplementedError

    @property
    def nearest_cycle(self) -> datetime:
        return dates.nearest_cycle(self.start_date)

    #@property
    #def s3(self):
    #    try:
    #        return self._s3
    #    except AttributeError:
    #        self._s3 = boto3.client(
    #            "s3", config=Config(signature_version=UNSIGNED))
    #        return self._s3

    @property
    def tmpdir(self):
        if not hasattr(self, "_tmpdir"):
            self.__tmpdir = tempfile.TemporaryDirectory()
            self._tmpdir = pathlib.Path(self.__tmpdir.name)
        return self._tmpdir

    @property
    def files(self):
        return self._files


class AWSHindcastInventory(AWSDataInventory):
    def __new__(cls):
        return object.__new__(AWSHindcastInventory)

    def __init__(
        self,
        start_date: datetime = None,
        rnday: Union[int, float, timedelta] = timedelta(days=5.0),
        product=None,
        verbose=False,
        fallback=True,
        cache=None,
    ):
        """This will download the National Water Model retro data.
        A 42-year (February 1979 through December 2020) retrospective
        simulation using version 2.1 of the NWM.
        """
        self.product = "CHRTOUT_DOMAIN1.comp" if product is None else product
        self.cache = cache
        self.start_date = (
            dates.nearest_cycle()
            if start_date is None
            else dates.nearest_cycle(dates.localize_datetime(start_date))
        )
        # self.start_date = self.start_date.replace(tzinfo=None)
        self.rnday = rnday if isinstance(
            rnday, timedelta) else timedelta(days=rnday)
        self.fallback = fallback
        self._files = {
            _: None
            for _ in np.arange(
                self.start_date,
                self.start_date + self.rnday + self.output_interval,
                self.output_interval,
            ).astype(datetime)
        }

        end_date = self.start_date + self.rnday

        years = np.arange(self.start_date.year, end_date.year+1)
        
        file_metadata = [] 
        for it, year in enumerate(years):
            paginator = self.s3.get_paginator("list_objects_v2")
            pages = paginator.paginate(
                Bucket=self.bucket, Prefix=f"model_output/{year}"
            )

            self.data = []
            for page in pages:
                for obj in page["Contents"]:
                    self.data.append(obj)

            metadata = sorted([_["Key"] for _ in self.data if "CHRTOUT_DOMAIN1.comp" in _["Key"]])
            [file_metadata.append(i) for i in metadata] 

        timevector = np.arange(
            datetime(self.start_date.year, 1, 1),
            datetime(end_date.year + 1, 1, 1),
            np.timedelta64(1, "h"),
            dtype="datetime64",
        )

        timefile = {
            pd.to_datetime(str(timevector[i])): file_metadata[i]
            for i in range(len(timevector))
        }

        for requested_time in self._files:
            logger.info(f"Requesting NWM data for time {requested_time}")
            self._files[requested_time] = self.request_data(
                timefile.get(requested_time)
            )

    def request_data(self, key):
        filename = self.tmpdir / key
        if filename.is_file() is False:
            cached_file = list(self.tmpdir.glob(f'**/{filename.name}'))
            if len(cached_file) == 1:
                filename = cached_file[0]
                logger.info(f"Using cached file {filename}, ")
            else:
                filename.parent.mkdir(parents=True, exist_ok=True)
                tmpfile = tempfile.NamedTemporaryFile().name
                with open(tmpfile, "wb") as f:
                    logger.info(f"Downloading file {key}, ")
                    self.s3.download_fileobj(self.bucket, key, f)
                shutil.move(tmpfile, filename)
        return filename

    @property
    def bucket(self):
        #return "noaa-nwm-retro-v2.0-pds"
        return "noaa-nwm-retrospective-2-1-pds"

    @property
    def s3(self):
        try:
            return self._s3
        except AttributeError:
            self._s3 = boto3.client(
                "s3", config=Config(signature_version=UNSIGNED))
            return self._s3

    @property
    def output_interval(self) -> timedelta:
        return {"CHRTOUT_DOMAIN1.comp": timedelta(hours=1)}[self.product]

    @property
    def cached_files(self):
        return sorted(list(self.tmpdir.glob("**/*.comp")))

    @property
    def cache(self):
        return self._cache

    @cache.setter
    def cache(self, cache: Union[str, os.PathLike, None, bool]):
        if cache is None or cache is False:
            self._cache = False
        elif cache is True:
            self._cache = pathlib.Path(
                appdirs.user_cache_dir(
                    f"pyschism/nwm/hindcast_data/{self.product}")
            )
            self._cache.mkdir(exist_ok=True, parents=True)
            self._tmpdir = self._cache

        elif isinstance(cache, (str, os.PathLike)):
            self._cache = pathlib.Path(cache)
            self._cache.mkdir(exist_ok=True, parents=True)
            self._tmpdir = self._cache
        else:
            raise TypeError(
                f"Unhandled argument cache={cache} of type {type(cache)}.")

class GOOGLEHindcastInventory(AWSDataInventory):
    def __new__(cls):
        return object.__new__(GOOGLEHindcastInventory)

    def __init__(
        self,
        start_date: datetime = None,
        rnday: Union[int, float, timedelta] = timedelta(days=5.0),
        product='medium_range_mem1',
        verbose=False,
        fallback=True,
        cache=None,
    ):

        self.start_date = dates.nearest_cycle() if start_date is None \
            else dates.nearest_cycle(dates.localize_datetime(start_date))
        self.rnday = rnday if isinstance(rnday, timedelta) \
            else timedelta(days=rnday)
        self.product = product
        self.verbose = verbose
        self.fallback = fallback
        self.cache = cache
        self._files = {_: None for _ in np.arange(
            self.start_date,
            self.start_date + self.rnday + self.output_interval,
            self.output_interval
        ).astype(datetime)}

        timevector=np.arange(datetime(self.start_date.year, 1, 1),
            datetime(self.start_date.year+1,1, 1),
            np.timedelta64(1, 'h'),
            dtype='datetime64')

        for requested_time, _ in self._files.items():

            logger.info(f'Requesting NWM data for time {requested_time}')

            if requested_time.hour == 0:
                requested_time2 = requested_time - timedelta(days=1)
            else:
                requested_time2 = requested_time

            self._files[requested_time] = self.request_data(requested_time, requested_time2)

    def request_data(self, request_time, request_time2):

        fname = self.tmpdir / f'nwm.t00z.{self.product[:12]}.channel_rt_1.' \
            f'{request_time.strftime("%Y%m%d%H")}.conus.nc'
        #logger.info(f'fname is {fname}')

        if fname.is_file():
            cached_file = list(self.tmpdir.glob(f'**/{fname.name}'))
            if len(cached_file) == 1:
                fname = cached_file[0]
                logger.info(f'Use cached file {fname}')
        else:
        
            fname = f'{self.start_date.strftime("%Y%m%d")}/nwm.t00z.' \
                f'{self.product[:12]}.channel_rt_1.{request_time.strftime("%Y%m%d%H")}.conus.nc'

            logger.info(f'Downloading file {request_time}, ')

            it = request_time.strftime("%H")
            if it == '00':
                logger.info(f'Requesting data at 00Z from yesterday!')
                it = str(int(it) + 24)
            it = it.zfill(3)

            url = f'https://storage.googleapis.com/national-water-model/nwm.{request_time2.strftime("%Y%m%d")}' \
                f'/{self.product}/nwm.t00z.{self.product[:12]}.channel_rt_1.f{it}.conus.nc'
            logger.info(f'{url}')
            try:
                wget.download(url, fname)
            except:
                logger.info(f'No data for {request_time}!')
                
        return fname

    @property
    def output_interval(self) -> timedelta:
        return {
            'medium_range_mem1': timedelta(hours=1)
        }[self.product]

    @property
    def cached_files(self):
        return sorted(list(self.tmpdir.glob("**/*.comp")))

    @property
    def cache(self):
        return self._cache

    @cache.setter
    def cache(self, cache: Union[str, os.PathLike, None, bool]):
        if cache is None or cache is False:
            self._cache = False
        elif cache is True:
            self._cache = pathlib.Path(
                appdirs.user_cache_dir(
                    f"pyschism/nwm/hindcast_data/{self.product}")
            )
            self._cache.mkdir(exist_ok=True, parents=True)
            self._tmpdir = self._cache

        elif isinstance(cache, (str, os.PathLike)):
            self._cache = pathlib.Path(cache)
            self._cache.mkdir(exist_ok=True, parents=True)
            self._tmpdir = self._cache
        else:
            raise TypeError(
                f"Unhandled argument cache={cache} of type {type(cache)}.")

class AWSForecastInventory(AWSDataInventory):
    def __new__(cls):
        return object.__new__(AWSForecastInventory)

    def __init__(
        self,
        start_date: datetime = None,
        rnday: Union[int, float, timedelta] = timedelta(days=5.0),
        product=None,
        verbose=False,
        fallback=True,
        cache=None,
    ):
        """This will download the latest National Water Model data.

        NetCDF files are saved to the system's temporary directory.
        The AWS data goes back 30 days. For requesting hindcast data from
        before we need a different data source
        """
        self.product = "medium_range_mem1" if product is None else product
        self.cache = cache
        self.start_date = (
            dates.nearest_cycle()
            if start_date is None
            else dates.nearest_cycle(dates.localize_datetime(start_date))
        )

        yesterday = self.start_date - timedelta(days=1)
       
        #check if previous day's data folder exists
        if not os.path.exists(yesterday.strftime("%Y%m%d")):
            logger.info(f"Downloading NWM data for {yesterday}")
            _ = self.download(yesterday, days=1)
        
        filemaps = self.download(self.start_date, days=rnday.days)

        # self.start_date = self.start_date.replace(tzinfo=None)
        self.rnday = rnday if isinstance(
            rnday, timedelta) else timedelta(days=rnday)
        self.fallback = fallback
        self._files = {
            _: None
            for _ in np.arange(
                self.start_date, # + timedelta(hours=3),
                self.start_date + self.rnday + self.output_interval,
                self.output_interval,
            ).astype(datetime)
        }

        for it, (requested_time, _) in enumerate(self._files.items()):

            logger.info(f"Requesting NWM data for time {requested_time}")
            if it == 0 and requested_time.hour == 0:
                f0 = f'{yesterday.strftime("%Y%m%d")}/nwm.' \
                    + yesterday.strftime("%Y%m%d") \
                    + '/medium_range_mem1/nwm.t00z.medium_range.channel_rt_1.f024.conus.nc'
                logger.info(f'Using data {f0}')
                self.files[requested_time] = f0
                continue

            logger.info(f"Using data {filemaps[requested_time]}")
            self._files[requested_time] = filemaps[requested_time]

    def download(self, nwmdate, days):
        self.data = self.s3.list_objects_v2(
            Bucket=self.bucket,
            Delimiter="/",
            Prefix=f'nwm.{nwmdate.strftime("%Y%m%d")}' f"/{self.product}/",
        )

        file_metadata = list(
                sorted(
                    [
                        _["Key"]
                        for _ in self.data["Contents"]
                        if "channel" in _["Key"] and "Contents" in self.data
                    ]
                )
            )

        filedir = pathlib.Path(nwmdate.strftime("%Y%m%d"))
        filedir.mkdir(exist_ok=True, parents=True)

        filedict = {}
        for key in file_metadata[0:days*24+1]:
            filetime = self.key2date(key)
            filename = nwmdate.strftime("%Y%m%d") + '/' + key
            filesubdir = pathlib.Path(filename)
            filesubdir.parent.mkdir(parents=True, exist_ok=True)

            logger.info(f"Downloading file {key}, ")
            self.s3.download_file(self.bucket, key, filename)
            filedict[filetime] = filename
       
        return filedict

    def key2date(self, key):
        base_date_str = f'{key.split("/")[0].split(".")[-1]}'
        timedelta_str = key.split("channel_rt_1.")[-1].split(".")[0].strip("f")
        return (
            datetime.strptime(base_date_str, "%Y%m%d")
            + timedelta(hours=float(timedelta_str))
            # + timedelta(hours=int(key.split('.')[2].strip('tz')))
            # + timedelta(hours=float(timedelta_str))
        )

    @property
    def bucket(self):
        return "noaa-nwm-pds"
    
    @property
    def s3(self):
        try:
            return self._s3
        except AttributeError:
            self._s3 = boto3.client(
                "s3", config=Config(signature_version=UNSIGNED))
            return self._s3

    @property
    def output_interval(self) -> timedelta:
        return {"medium_range_mem1": timedelta(hours=1)}[self.product]

    @property
    def timevector(self):
        return np.arange(
            self.start_date,
            self.start_date + self.rnday + self.output_interval,
            self.output_interval,
        ).astype(datetime)

    #@property
    #def files(self):
    #    return sorted(list(self.tmpdir.glob("**/*.nc")))

    @property
    def requested_product(self):
        return {"medium_range_mem1": "medium_range.channel_rt_1"}[self.product]

    @property
    def cache(self):
        return self._cache

    @cache.setter
    def cache(self, cache: Union[str, os.PathLike, None, bool]):
        if cache is None or cache is False:
            self._cache = False
        elif cache is True:
            self._cache = pathlib.Path(
                appdirs.user_cache_dir(
                    f"pyschism/nwm/forecast_data/{self.product}")
            )
            self._cache.mkdir(exist_ok=True, parents=True)
            self._tmpdir = self._cache

        elif isinstance(cache, (str, os.PathLike)):
            self._cache = pathlib.Path(cache)
            self._cache.mkdir(exist_ok=True, parents=True)
            self._tmpdir = self._cache
        else:
            raise TypeError(
                f"Unhandled argument cache={cache} of type {type(cache)}.")


class NationalWaterModel(SourceSink):

    start_date = dates.StartDate()
    end_date = dates.EndDate()
    # run_days = dates.RunDays()

    def __init__(
        self, aggregation_radius=None, pairings=None, nwm_file=None, cache=False
    ):
        self.nwm_file = nwm_file
        self.aggregation_radius = aggregation_radius
        self.pairings = pairings
        self.cache = cache
        self._data = None

    def _fetch_data(
        self,
        gr3: Gr3,
        start_date: datetime = None,
        end_date: Union[datetime, timedelta] = None,
        nprocs=-1,
        product=None,
    ):
        nprocs = -1 if nprocs is None else nprocs
        nprocs = cpu_count() if nprocs == -1 else nprocs
        self.start_date = start_date
        self.end_date = end_date
        self.pairings = (
            NWMElementPairings(gr3, nwm_file=self.nwm_file, workers=nprocs)
            if self.pairings is None
            else self.pairings
        )
        self._inventory = AWSDataInventory(
            start_date=self.start_date,
            rnday=self.end_date - self.start_date,
            cache=self.cache,
            # product="medium_range_mem1",
        )

        self._timevector = [
            dates.localize_datetime(d) for d in self.inventory.files.keys()
        ]

        #src_idxs, snk_idxs = self.inventory.get_nc_pairing_indexes(
        #    self.pairings)
        #logger.info(f'Start aggregating NWM timeseries using nprocs={nprocs}')
        start0 = datetime.now()
        #with Pool(processes=nprocs) as pool:
        #    sources = pool.starmap(
        #        streamflow_lookup,
        #        [(file, self.pairings.sources) for file in self.inventory.files.values()],
        #    )
        #    sinks = pool.starmap(
        #        streamflow_lookup,
        #        [(file, self.pairings.sinks) for file in self.inventory.files.values()],
        #    )
        #pool.join()

        sources = []
        sinks = []
        nc_fid0 = Dataset(list(self.inventory.files.values())[0])["feature_id"][:]
        src_idxs = get_aggregated_features(nc_fid0, self.pairings.sources.values())
        snk_idxs = get_aggregated_features(nc_fid0, self.pairings.sinks.values())
        for file in self.inventory.files.values():
            start = datetime.now()
            nc = Dataset(file)
            ncfeatureid=nc['feature_id'][:]
            if not np.all(ncfeatureid == nc_fid0):
                logger.info(f'Indexes of feature_id are changed in  {file}')
                src_idxs=get_aggregated_features(ncfeatureid, self.pairings.sources.values())
                snk_idxs=get_aggregated_features(ncfeatureid, self.pairings.sinks.values())
                nc_fid0 = ncfeatureid

            sources.append(streamflow_lookup(file, src_idxs))
            sinks.append(streamflow_lookup(file, snk_idxs))
            logger.info(f'Processing file {file} took {datetime.now() - start}')
            nc.close()

        source_data = {}
        sink_data = {}
        for i, file in enumerate(self.inventory.files.values()):
            nc = Dataset(file)
            _time = dates.localize_datetime(
                datetime.strptime(nc.model_output_valid_time,
                                  "%Y-%m-%d_%H:%M:%S")
            )
            for j, element_id in enumerate(self.pairings.sources):
                source_data.setdefault(_time, {})[element_id] = {
                    "flow": sources[i][j],
                    "temperature": -9999.0,
                    "salinity": 0.0,
                }

            for k, element_id in enumerate(self.pairings.sinks):
                sink_data.setdefault(_time, {})[element_id] = {
                    "flow": -sinks[i][k],
                }
            nc.close()
        logger.info(f'Timeseries aggregation took {datetime.now() - start0}')
        self._sources = Sources(source_data)
        self._sinks = Sinks(sink_data)
        self._data = {**source_data, **sink_data}

    def write(
        self,
        output_directory,
        gr3: Gr3,
        start_date: datetime = None,
        end_date: Union[datetime, timedelta] = None,
        overwrite: bool = False,
        nprocs=-1,
        product=None,
        msource: Union[str, bool] = True,
        vsource: Union[str, bool] = True,
        vsink: Union[str, bool] = True,
        source_sink: Union[str, bool] = True,
    ):
        if self._data is None:
            self._data = {}
            self._fetch_data(
                gr3,
                start_date=start_date,
                end_date=end_date,
                nprocs=nprocs,
                product=product,
            )

        if self.aggregation_radius is not None:
            self.aggregate_by_radius(gr3, self.aggregation_radius)

        super().write(
            output_directory,
            overwrite=overwrite,
            msource=msource,
            vsource=vsource,
            vsink=vsink,
            source_sink=source_sink,
        )

    @property
    def timevector(self):
        return self._timevector

    @property
    def inventory(self):
        return self._inventory
