#! /usr/bin/env python
from datetime import datetime
import logging
from multiprocessing import Pool, cpu_count
import os
import pathlib
from time import time

from netCDF4 import Dataset

from pyschism.mesh import Hgrid
from pyschism.forcing.hydrology import (
    Sources,
    Sinks,
    Msource,
    Vsource,
    Vsink,
    SourceSink,
    Hydrology
)
from pyschism.forcing.hydrology.nwm import (
    NWMElementPairings,
    AWSDataInventory,
    streamflow_lookup,
    localize_datetime,
    pivot_time
)


if __name__ == '__main__':

    logging.getLogger().setLevel(logging.INFO)

    logging.info('Read hgrid.')
    start = time()
    hgrid = Hgrid.open(
        os.getenv('NWM_TEST_MESH'),
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

    # Illustrate parallel vs serial call
    nprocs = cpu_count()
    if nprocs > 1:
        logging.info('Running parallalized streamflow_lookup.')
        start = time()
        with Pool(processes=cpu_count()) as pool:
            res = pool.starmap(
                streamflow_lookup,
                [(file, pairings) for file in awsdata.files])
        pool.join()
        logging.info(f'streamflow_lookup in parallel took {time() - start}.')
    else:
        logging.info('streamflow_lookup serial test.')
        start = time()
        res = [streamflow_lookup(file, pairings) for file in awsdata.files]
        logging.info(f'streamflow_lookup in serial took {time() - start}.')

    logging.info('Adding streamflow data as sources and sinks...')
    start = time()
    sources = Sources()
    sinks = Sinks()
    for i, file in enumerate(awsdata.files):
        nc = Dataset(file)
        _time = localize_datetime(datetime.strptime(
            nc.model_output_valid_time,
            "%Y-%m-%d_%H:%M:%S"))
        for j, element_id in enumerate(pairings.sources.keys()):
            sources.add_data(_time, element_id, res[i][0][j], -9999, 0.)
        for k, element_id in enumerate(pairings.sinks.keys()):
            sinks.add_data(_time, element_id, res[i][1][k])
    logging.info(
        f'Adding streamflow data as sources and sinks took {time() - start}.')

    print(str(SourceSink()))
    print(str(Vsource(awsdata.start_date, awsdata.rnday)))
    print(str(Msource(awsdata.start_date, awsdata.rnday)))
    print(str(Vsink(awsdata.start_date, awsdata.rnday)))

    outdir = pathlib.Path('hydro_test')
    outdir.mkdir(exist_ok=True)

    hydro = Hydrology()

    hydro.set_start_date(pivot_time())
    hydro.set_rnday(awsdata.rnday)
    hydro.write(outdir, overwrite=True)
