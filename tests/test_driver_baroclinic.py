#! /usr/bin/env python
from datetime import timedelta
import logging
# import os
# import pathlib
import shutil
# import tarfile
# import tempfile
import unittest
# import urllib.request

from pyschism import dates
from pyschism.mesh import Hgrid, Vgrid, Fgrid
from pyschism.driver import ModelConfig
from pyschism.forcing.atmosphere import NWS2, GFS, HRRR
from pyschism.forcing.baroclinic import RTOFS, GOFS
from pyschism.forcing.hydrology import NWM
from pyschism.forcing.tides import Tides


logging.basicConfig(level=logging.INFO, force=True)


class ModelConfigurationTestCase(unittest.TestCase):

    def test_basic_config_baroclinic(self):
        config = ModelConfig(
            Hgrid.open('data/baroclinic/hgrid.gr3', crs='epsg:4326'),
            vgrid=Vgrid.open('data/baroclinic/vgrid.in'),
            fgrid=Fgrid.open('data/baroclinic/drag.gr3', crs='epsg:4326'),
            tides=Tides(),
            atmosphere=NWS2(
                GFS(),
                # HRRR()
            ),
            hydrology=NWM(),
            # baroclinic=RTOFS()
            baroclinic=GOFS()
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
