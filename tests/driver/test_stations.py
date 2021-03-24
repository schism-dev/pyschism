#! /usr/bin/env python
from datetime import timedelta
import logging
import os
from time import time

from pyschism.mesh import Hgrid
from pyschism.stations import Stations


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)

    logging.info('Read hgrid.')
    start = time()
    hgrid = Hgrid.open(
        os.getenv('NWM_TEST_MESH'),
        crs='EPSG:4326'
    )
    logging.info(f'Read hgrid took {time() - start}.')

    stations = Stations.from_file('/home/jreniel/schism-dev/pyschism-forecast-useast/static/station.in', timedelta(minutes=6.0))
    stations.clip(hgrid.hull.multipolygon())
    # breakpoint()
