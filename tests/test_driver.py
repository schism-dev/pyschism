#! /usr/bin/env python
from datetime import timedelta
import pathlib
import tarfile
import tempfile
import unittest
import urllib.request

from pyschism.mesh import Hgrid
from pyschism.driver import ModelConfig
from pyschism.forcing.tides import Tides
from pyschism.forcing.atmosphere import NWS2, GFS, HRRR
from pyschism.forcing.hydrology import NWM


DATA_DIRECTORY = pathlib.Path(__file__).parent.absolute() / 'data'
FORT14 = DATA_DIRECTORY / "NetCDF_Shinnecock_Inlet/fort.14"


class ModelConfigurationTestCase(unittest.TestCase):

    def setUp(self):
        if not FORT14.is_file():
            url = "https://www.dropbox.com/s/1wk91r67cacf132/"
            url += "NetCDF_shinnecock_inlet.tar.bz2?dl=1"
            g = urllib.request.urlopen(url)
            tmpfile = tempfile.NamedTemporaryFile()
            with open(tmpfile.name, 'b+w') as f:
                f.write(g.read())
            with tarfile.open(tmpfile.name, "r:bz2") as tar:
                tar.extractall(DATA_DIRECTORY / "NetCDF_Shinnecock_Inlet")

    def test_basic_config(self):
        hgrid = Hgrid.open(FORT14, crs='epsg:4326')
        tides = Tides()
        tides.use_all()
        nws2 = NWS2(GFS(), HRRR())
        nwm = NWM()
        config = ModelConfig(hgrid, tides=tides, atmosphere=nws2,
                             hydrology=nwm, dramp=timedelta(days=1)
        # coldstart = config.get_coldstart()
        # hotstart = config.get_hotstart()
        driver = config.get_driver()
        driver.run()


if __name__ == '__main__':
    unittest.main()
