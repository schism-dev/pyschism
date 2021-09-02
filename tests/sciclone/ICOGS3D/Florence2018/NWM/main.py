#!/usr/bin/env python
from datetime import datetime
import logging
import pathlib

from pyschism.forcing.source_sink.nwm import NationalWaterModel, NWMElementPairings
from pyschism.mesh import Hgrid

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)

logging.getLogger("pyschism").setLevel(logging.DEBUG)

output_directory = pathlib.Path("outputs")

startdate = datetime(2018, 8, 24)

rnday = 36

hgrid = Hgrid.open('hgrid.ll')

sources_pairings = pathlib.Path("static/sources.json")
sinks_pairings = pathlib.Path("static/sinks.json")

if all([sources_pairings.is_file(), sinks_pairings.is_file()]) is False:
    pairings = NWMElementPairings(hgrid)
    sources_pairings.parent.mkdir(exist_ok=True, parents=True)
    pairings.save_json(sources=sources_pairings, sinks=sinks_pairings)

else:
    pairings = NWMElementPairings.load_json(
        hgrid, sources_pairings, sinks_pairings)

nwm = NationalWaterModel(
    # aggregation_radius=4000,
    pairings=pairings,
    cache='NWM_v2.0'
)
start = datetime.now()
nwm.write(
    output_directory,
    hgrid,
    startdate,
    rnday,
    overwrite=True,
)
print(f'Write out took {datetime.now() - start}')
