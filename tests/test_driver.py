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
FORT14 = DATA_DIRECTORY / "NetCDF_Shinnecock_Inlet/fort.14"

logging.basicConfig(level=logging.INFO, force=True)


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

    def test_basic_config_2d(self):

        import os
        config = ModelConfig(
            Hgrid.open(os.getenv('NWM_TEST_MESH'), crs='epsg:4326'),
            tides=Tides(tidal_database='tpxo'),
            atmosphere=NWS2(
                GFS(),
                HRRR()
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
            # timestep=300.,
            dramp=spinup_time,
            dramp_ss=spinup_time,
            drampwind=spinup_time,
            nspool=timedelta(hours=1),
            elev=True,
            dahv=True,
        )

        # # optionally run or write the coldstart object
        # if shutil.which('pschism_TVD-VL') is not None:
        #     coldstart.run('/tmp/test/coldstart', overwrite=True)

        # else:
        #     tmpdir = tempfile.TemporaryDirectory()
        #     coldstart.write(tmpdir.name)
        coldstart.outdir = '/tmp/test/coldstart'
        hotstart = config.hotstart(
            coldstart,
            timestep=300.,
            end_date=timedelta(days=2) - timedelta(hours=2),
            nspool=timedelta(hours=1),
            elev=True,
            dahv=True,
        )
        if shutil.which('pschism_TVD-VL') is not None:
            # optionally run or write the coldstart object
            hotstart.run('/tmp/test/hotstart', overwrite=True)
        else:
            tmpdir = tempfile.TemporaryDirectory()
            hotstart.write(tmpdir.name)


if __name__ == '__main__':
    unittest.main()
