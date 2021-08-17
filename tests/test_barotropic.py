#! /usr/bin/env python
from datetime import timedelta
import logging
import pathlib
import shutil

from pyschism import dates
from pyschism.driver import ModelConfig
from pyschism.forcing.bctides import iettype, ifltype
from pyschism.forcing.nws import NWS2, GFS, HRRR
from pyschism.forcing.source_sink import NWM
from pyschism.mesh import Hgrid


def test_barotropic_forecast():

    test_output_directory = pathlib.Path('/tmp/test')
    if test_output_directory.exists():
        shutil.rmtree(test_output_directory)

    config = ModelConfig(
        Hgrid.open(
            "https://raw.githubusercontent.com/geomesh/test-data/main/NWM/hgrid.ll",
            crs="epsg:4326",
        ),
        iettype=iettype.Iettype3(
            database="tpxo"
            ),
        ifltype=ifltype.Ifltype3(
            database="tpxo"
            ),
        nws=NWS2(
                GFS(),
                HRRR()
            ),
        source_sink=NWM(aggregation_radius=4000.),
    )
    # create reference dates
    nearest_cycle = dates.nearest_cycle()
    spinup_time = timedelta(days=1)
    coldstart = config.coldstart(
        start_date=nearest_cycle - spinup_time,
        end_date=nearest_cycle,
        timestep=300.0,
        dramp=spinup_time,
        dramp_ss=spinup_time,
        drampwind=spinup_time,
        nspool=timedelta(hours=1),
        elev=True,
        dahv=True,
    )

    coldstart.run(test_output_directory / "coldstart")

    hotstart = config.hotstart(
        coldstart,
        end_date=timedelta(days=1) + timedelta(hours=23.),
        timestep=300.0,
        nspool=timedelta(hours=1.),
        elev=True,
        dahv=True,
    )
    hotstart.run(test_output_directory / "hotstart")


if __name__ == "__main__":

    logging.basicConfig(
        format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
        force=True,
    )
    logging.captureWarnings(True)
    logging.getLogger("pyschism").setLevel(logging.DEBUG)
    test_barotropic_forecast()
