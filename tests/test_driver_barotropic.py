#! /usr/bin/env python
from datetime import timedelta
import logging
import os
import pathlib
import shutil
import tarfile
import tempfile
import unittest
import urllib.request

from pyschism import dates
from pyschism.driver import ModelConfig
from pyschism.forcing.atmosphere import NWS2, GFS
from pyschism.forcing.hydrology import NWM
from pyschism.forcing.tides import Tides
from pyschism.mesh import Hgrid


logging.basicConfig(level=logging.INFO, force=True)


class ModelConfigurationTestCase(unittest.TestCase):

    def setUp(self):
        hgrid = os.getenv('NWM_TEST_MESH')
        if hgrid is None:
            data_directory = pathlib.Path(__file__).parent.absolute() / 'data'
            hgrid = data_directory / "GulfStreamDevel/hgrid.gr3"
            url = "https://www.dropbox.com/s/mjaxaqeggy721um/"
            url += "Gulf_Stream_develop.tar.gz?dl=1"
            g = urllib.request.urlopen(url)
            tmpfile = tempfile.NamedTemporaryFile()
            with open(tmpfile.name, 'b+w') as f:
                f.write(g.read())
            with tarfile.open(tmpfile.name, "r:gz") as tar:
                tar.extractall(data_directory / "GulfStreamDevel")
        self.hgrid = hgrid

    def test_basic_config(self):
        config = ModelConfig(
            Hgrid.open(self.hgrid, crs='epsg:4326'),
            tides=Tides(),
            atmosphere=NWS2(
                GFS(),
                # HRRR()
            ),
            hydrology=NWM(),
        )

        # create reference dates
        nearest_cycle = dates.nearest_cycle()
        spinup_time = timedelta(days=1)

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
            coldstart.run()

        else:
            ctmpdir = tempfile.TemporaryDirectory()
            coldstart.write(ctmpdir.name)

        hotstart = config.hotstart(
            coldstart,
            end_date=timedelta(days=1),
            timestep=300.,
            nspool=timedelta(hours=1),
            elev=True,
            dahv=True,
        )

        if shutil.which('pschism_TVD-VL') is not None:
            hotstart.run()
        else:
            htmpdir = tempfile.TemporaryDirectory()
            hotstart.write(htmpdir.name)

    def _test_basic_config_passive_ts(self):
        # WIP: ibc=1, ibtp=1
        from pyschism.forcing.baroclinic import RTOFS
        config = ModelConfig(
            Hgrid.open(self.hgrid, crs='epsg:4326'),
            tides=Tides(),
            atmosphere=NWS2(GFS()),
            baroclinic=RTOFS(),
            hydrology=NWM(),
        )

        # create reference dates
        nearest_cycle = dates.nearest_cycle()
        spinup_time = timedelta(days=1)

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
        if shutil.which('_pschism_TVD-VL') is not None:
            coldstart.run('/tmp/test/coldstart', overwrite=True)

        else:
            # ctmpdir = tempfile.TemporaryDirectory()
            coldstart.write(
                # ctmpdir.name
                '/tmp/test/coldstart', overwrite=True
                )

        hotstart = config.hotstart(
            coldstart,
            end_date=timedelta(days=1),
            timestep=300.,
            nspool=timedelta(hours=1),
            elev=True,
            dahv=True,
        )

        if shutil.which('_pschism_TVD-VL') is not None:
            hotstart.run('/tmp/test/hotstart', overwrite=True)
        else:
            # htmpdir = tempfile.TemporaryDirectory()
            hotstart.write(
                # htmpdir.name
                '/tmp/test/hotstart', overwrite=True
                )


if __name__ == '__main__':
    unittest.main()
