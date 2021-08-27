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

output_directory = pathlib.Path("Florence")

startdate = datetime(2018, 8, 24)

rnday = 36
rnday = 2 / 24
import os

hgrid = Hgrid.open(
    # HGRID_PATH,
    os.getenv("NWM_TEST_MESH"),
    crs="epsg:4326",
)

sources_pairings = output_directory / "sources.json"
sinks_pairings = output_directory / "sinks.json"

if all([sources_pairings.is_file(), sinks_pairings.is_file()]) is False:
    pairings = NWMElementPairings(hgrid)
    sources_pairings.parent.mkdir(exist_ok=True, parents=True)
    pairings.save_json(sources=sources_pairings, sinks=sinks_pairings)

else:
    pairings = NWMElementPairings.load_json(hgrid, sources_pairings, sinks_pairings)

nwm = NationalWaterModel(aggregation_radius=4000, pairings=pairings, cache=True)
# print(rnday)
# exit()
nwm.write(
    output_directory,
    hgrid,
    startdate,
    rnday,
    overwrite=True,
)

exit()

# test how much longer it takes depending on the value of rnday
logging.getLogger("pyschism").setLevel(logging.WARNING)


def Fibonacci(n):

    # Check if input is 0 then it will
    # print incorrect input
    if n < 0:
        print("Incorrect input")

    # Check if n is 0
    # then it will return 0
    elif n == 0:
        return 0

    # Check if n is 1,2
    # it will return 1
    elif n == 1 or n == 2:
        return 1

    else:
        return Fibonacci(n - 1) + Fibonacci(n - 2)


for i in range(2, 10):
    start = datetime.now()
    rnday = Fibonacci(i)
    nwm.write(
        output_directory / f"rnday{rnday}",
        hgrid,
        startdate,
        Fibonacci(i),
        overwrite=True,
    )
    print(f"Took {datetime.now() - start} for {rnday}")
