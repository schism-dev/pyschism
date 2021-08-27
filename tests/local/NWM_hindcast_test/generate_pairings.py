from pyschism.forcing import source_sink
from pyschism.mesh import Hgrid

import logging

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.captureWarnings(True)

log_level = logging.DEBUG

source_sink.nwm.logger.setLevel(log_level)
source_sink.base.logger.setLevel(log_level)

hgrid = Hgrid.open("../ICOGS3D/Hgrid/hgrid.gr3", crs="epsg:4326")
pairings = source_sink.nwm.NWMElementPairings(hgrid)
pairings.sources_gdf.to_file("sources/sources.shp")
pairings.sinks_gdf.to_file("sinks/sinks.shp")
