#! /usr/bin/env python
import pathlib
import logging

from matplotlib.transforms import Bbox

from pyschism.forcing.atmosphere import GlobalForecastSystem
from pyschism.forcing.atmosphere.nws.nws2.sflux import SfluxDataset


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, force=True)
    gfs = GlobalForecastSystem()
    gfs.fetch_data(
        bbox=Bbox([[17.623082, -67.730713], [18.786717, -65.055542]]),
        rnday=0.25,
        prc=False,
        rad=False
    )
    test_gfs = pathlib.Path('test_gfs')
    test_gfs.mkdir(exist_ok=True)
    gfs.write('test_gfs', 1, overwrite=True)
