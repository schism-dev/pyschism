#! /usr/bin/env python

import logging
from multiprocessing import Pool, cpu_count
from time import time

from netCDF4 import Dataset
import numpy as np
from pyschism.mesh import Hgrid
from pyschism.forcing.hydrology.nwm import (
    NWMElementPairings, AWSDataInventory,
    # streamflow_lookup
)


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


if __name__ == '__main__':

    logging.getLogger().setLevel(logging.INFO)

    logging.info('Read hgrid.')
    start = time()
    hgrid = Hgrid.open(
        '/home/jreniel/schism-dev/pyschism-forecast-useast/dev-mesh/static/hgrid.gr3',
        crs='EPSG:4326'
    )
    logging.info(f'Read hgrid took {time() - start}.')

    logging.info('Create NWM element pairings.')
    start = time()
    pairings = NWMElementPairings(hgrid)
    logging.info(f'Create NWM element pairings took {time() - start}.')

    logging.info('Download AWS data.')
    start = time()
    awsdata = AWSDataInventory()
    logging.info(f'Download AWS data took {time() - start}.')

    logging.info('streamflow_lookup parallelized test.')
    start = time()
    with Pool(processes=cpu_count()) as pool:
        res = pool.starmap(
            streamflow_lookup,
            [(file, pairings) for file in awsdata.files])
    pool.join()
    logging.info(f'streamflow_lookup in parallel took {time() - start}.')
