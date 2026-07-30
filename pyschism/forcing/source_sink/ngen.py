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
from pyschism.utils.jsonencoder import NpEncoder

from pyschism.forcing.source_sink.base import SourceSink, Sources, Sinks

logger = logging.getLogger(__name__)

class NextGenElementPairings:
    def __init__(self, hgrid: Gr3, hyfab, workers=-1):

        # TODO: Accelerate using dask: https://blog.dask.org/2017/09/21/accelerating-geopandas-1

        self._hyfab_file = hyfab

        logger.info("Computing NextGenElementPairings...")
        self._hgrid = hgrid

        # An STR-Index returns the reaches that are near the boundaries of the
        # mesh. This subsamples the NextGen hydrofabric flowpaths network, but also is not the exact
        # result. This is used to speed-up computations by filtering the input
        # data.
        logger.info("Computing r_index.")
        start = time()
        ngen_r_index = self.gdf.sindex
        logger.info(f"Computing r_index took {time() - start}.")
     
        print(self.gdf)
        # The r-index is used to find intersections between mesh boundary edges
        # and NextGen hydrofabric reaches (approximate results)
        logger.info("Use r_index to filter features.")
        start = time()
        possible_indexes = set()
        for edge in hgrid.hull.edges().itertuples():
            for index in list(ngen_r_index.intersection(edge.geometry.bounds)):
                possible_indexes.add(index)
        possible_matches = self.gdf.iloc[list(possible_indexes)]
        logger.info(f"Filtering features took {time()-start}.")
        del possible_indexes
        del ngen_r_index

        # The hull rings itersections is used to find the exact NextGen reaches
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
                        sources[element.id].append(reaches.iloc[row.reachIndex].id)
                    else:
                        sinks[element.id].append(reaches.iloc[row.reachIndex].id)
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
                json.dump(self.sources, fh, cls=NpEncoder)

        sinks = "sinks.json" if sinks is True else sinks
        if sinks:
            with open(sinks, "w") as fh:
                logger.info(f"Saving {sinks}")
                json.dump(self.sinks, fh, cls=NpEncoder)

    @staticmethod
    def load_json(hgrid, sources=None, sinks=None):
        pairings = NextGenElementPairings.__new__(NextGenElementPairings)
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
#                        "features": " ".join(list(map(str, features))),
    @property
    def hgrid(self):
        return self._hgrid

    @property
    def gdf(self):
        self._gdf = gpd.read_file(self._hyfab_file,layer='flowpaths').to_crs('WGS84')
        return self._gdf

    @property
    def hyfab_file(self):
        return self._hyfab_file

class NextGen(SourceSink):

    def __init__(self, pairings=None):
        self.pairings = pairings

    def Source_nc_writer(self, start, end, dirname=None):
        if dirname is None:
            raise Exception("dirname is required.")
        os.makedirs(dirname, exist_ok=True)

        source_dt = pd.date_range(start=start,end=end,freq='H')[0:-1]

        fout = Dataset(f'{dirname}/source.nc', 'w', format='NETCDF4')
        fout.createDimension('nsources', len(self.pairings.sources_gdf))
        fout.createDimension('nsinks', len(self.pairings.sinks_gdf))
        fout.createDimension('ntracers', 2)
        fout.createDimension('time_msource', len(source_dt))
        fout.createDimension('time_vsource', len(source_dt))
        fout.createDimension('time_vsink', len(source_dt))
        fout.createDimension('one', np.atleast_1d(1))

        fout.createVariable('source_elem', 'i', ('nsources'))
        fout['source_elem'][:] = self.pairings.sources_gdf.element_id

        fout.createVariable('source_ids', 'i', ('nsources'))
        newlist = self.pairings.sources_gdf.features.values
        newlist = [s.replace('wb-','') for s in newlist]
        newlist = [s.split(' ')[0] for s in newlist]
        fout['source_ids'][:] = np.array(newlist,dtype=np.int32)

        fout.createVariable('vsource', 'f', ('time_vsource', 'nsources'))
        fout['vsource'][:,:] = -9999.

        fout.createVariable('msource', 'f', ('time_msource', 'ntracers', 'nsources'))
        fout['msource'][:,:,:] = -9999.

        fout.createVariable('sink_elem', 'i', ('nsinks'))
        fout['sink_elem'][:] = self.pairings.sinks_gdf.element_id

        fout.createVariable('sink_ids', 'i', ('nsinks'))
        newlist = self.pairings.sinks_gdf.features.values
        newlist = [s.replace('wb-','') for s in newlist]
        newlist = [s.split(' ')[0] for s in newlist]
        fout['sink_ids'][:] = np.array(newlist,dtype=np.int32)

        fout.createVariable('vsink', 'f', ('time_vsink', 'nsinks'))
        fout['vsink'][:,:] = -9999.

        fout.createVariable('time_step_vsource', 'f', ('one'))
        fout['time_step_vsource'][:] = 3600.

        fout.createVariable('time_step_msource', 'f', ('one'))
        fout['time_step_msource'][:] = 3600.

        fout.createVariable('time_step_vsink', 'f', ('one'))
        fout['time_step_vsink'][:] = 3600.

        fout.close()

