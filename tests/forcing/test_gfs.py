#! /usr/bin/env python
import unittest

from matplotlib.transforms import Bbox

from pyschism.forcing.atmosphere import GlobalForecastSystem


class GFSTestCase(unittest.TestCase):

    def test_nowcast(self):
        gfs = GlobalForecastSystem()
        gfs.fetch_data(
            bbox=Bbox([[17.623082, -67.730713], [18.786717, -65.055542]])
        )


if __name__ == '__main__':
    unittest.main()
