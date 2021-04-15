#! /usr/bin/env python
from datetime import timedelta
import logging
import pathlib
import shutil
import tarfile
import tempfile
import unittest
import urllib.request

from pyschism import dates
from pyschism.mesh import Hgrid
from pyschism.driver import ModelConfig
from pyschism.forcing.tides import Tides
from pyschism.forcing.atmosphere import NWS2, GFS, HRRR
from pyschism.forcing.hydrology import NWM


DATA_DIRECTORY = pathlib.Path(__file__).parent.absolute() / 'data'
HGRID = DATA_DIRECTORY / "GulfStreamDevel/hgrid.gr3"

logging.basicConfig(level=logging.INFO, force=True)


class ModelConfigurationTestCase(unittest.TestCase):

    def setUp(self):
        if not HGRID.is_file():
            url = "https://www.dropbox.com/s/mjaxaqeggy721um/"
            url += "Gulf_Stream_develop.tar.gz?dl=1"
            g = urllib.request.urlopen(url)
            tmpfile = tempfile.NamedTemporaryFile()
            with open(tmpfile.name, 'b+w') as f:
                f.write(g.read())
            with tarfile.open(tmpfile.name, "r:gz") as tar:
                tar.extractall(DATA_DIRECTORY / "GulfStreamDevel")

    def test_basic_config_2d(self):

        config = ModelConfig(
            Hgrid.open(HGRID, crs='epsg:4326'),
            tides=Tides(),
            atmosphere=NWS2(
                GFS(),
                # HRRR()
            ),
            hydrology=NWM()
        )

        # create reference dates
        nearest_cycle = dates.nearest_cycle()
        spinup_time = timedelta(days=0.25)

        # create a coldstart object
        coldstart = config.coldstart(
            start_date=nearest_cycle - spinup_time,
            end_date=nearest_cycle,
            timestep=300.,
            dramp=spinup_time,
            dramp_ss=spinup_time,
            drampwind=spinup_time,
            nspool=timedelta(hours=1),
            elev=True,
            dahv=True,
        )

        # optionally run or write the coldstart object
        if shutil.which('pschism_TVD-VL') is not None:
            coldstart.run('/tmp/test/coldstart', overwrite=True)

        else:
            ctmpdir = tempfile.TemporaryDirectory()
            coldstart.write(ctmpdir.name)

        hotstart = config.hotstart(
            coldstart,
            timestep=300.,
            end_date=timedelta(days=2) - timedelta(hours=2),
            nspool=timedelta(hours=1),
            elev=True,
            dahv=True,
        )

        if shutil.which('pschism_TVD-VL') is not None:
            hotstart.run('/tmp/test/hotstart', overwrite=True)
        else:
            htmpdir = tempfile.TemporaryDirectory()
            hotstart.write(htmpdir.name)


if __name__ == '__main__':
    unittest.main()
