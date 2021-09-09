from datetime import datetime
from time import time

import logging

from pyschism.forcing import source_sink
from pyschism.mesh import Hgrid


logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.captureWarnings(True)

log_level = logging.DEBUG
source_sink.nwm.logger.setLevel(log_level)
source_sink.base.logger.setLevel(log_level)

startdate=datetime(2018, 8, 24)
rnday=40
hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")

t0=time()
nwm=source_sink.nwm.NationalWaterModel()
nwm.write('./', hgrid, startdate, rnday, overwrite=True)
print(f'It took {time()-t0} seconds to generate source/sink')
nwm.pairings.save_json(sources=True, sinks=True)
nwm.pairings.sources_gdf.to_file("sources_run07b_new.shp")
nwm.pairings.sinks_gdf.to_file("sinks_run07b_new.shp")
nwm.pairings.intersection.to_file("intersection_run07b_new.shp")
nwm.pairings.reaches.to_file("reaches_run07b_new.shp")
