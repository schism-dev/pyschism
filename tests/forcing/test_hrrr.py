#! /usr/bin/env python
import pathlib
import logging

from matplotlib.transforms import Bbox

from pyschism.forcing.atmosphere.hrrr import HRRR


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, force=True)
    hrrr = HRRR()
    hrrr.fetch_data(
        bbox=Bbox([[17.623082, -67.730713], [18.786717, -65.055542]]),
        rnday=0.25,
        prc=False,
        rad=False
    )
    test_hrrr = pathlib.Path('test_hrrr')
    test_hrrr.mkdir(exist_ok=True)
    hrrr.write('test_hrrr', 1, overwrite=True)
