from datetime import datetime
from time import time
import pathlib
import logging

from pyschism.forcing.source_sink.ngen import NextGen, NextGenElementPairings
from pyschism.mesh import Hgrid

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.captureWarnings(True)

log_level = logging.DEBUG
logging.getLogger('pyschism').setLevel(log_level)

if __name__ == '__main__':

    startdate = datetime(2017, 12, 1)
    enddate = datetime(2018, 1, 1)
    hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")

    t0 = time()

    #source/sink json files, if not exist, it will call NWMElementPairings to generate.
    sources_pairings = pathlib.Path('./sources.json')
    sinks_pairings = pathlib.Path('./sinks.json')
    output_directory = pathlib.Path('./')

    # Pathway to a NextGen hydrofabric geopackage file to extract
    # the SCHISM mesh sources/sinks with T-Route's flowpaths
    hydrofabric = pathlib.Path('./conus.gpkg')

    # check if source/sink json file exists
    if all([sources_pairings.is_file(), sinks_pairings.is_file()]) is False:
        pairings = NextGenElementPairings(hgrid)
        sources_pairings.parent.mkdir(exist_ok=True, parents=True)
        pairings.save_json(sources=sources_pairings, sinks=sinks_pairings)
    else:
        pairings = NextGenElementPairings.load_json(
            hgrid,
            sources_pairings,
            sinks_pairings)

    #check nc files, if not exist will download
    ngen=NextGen(pairings=pairings)

    # Generate the current Source.nc file required to execute the SCHISM
    # Basic Model Interface (BMI). This file is set up to essentially just initalize
    # SCHISM source/sink element fields and by the external forcing file
    # requirements when a user specifies NWM_BMI=1 within the param.nml file.
    # This function also includes the NextGen T-Route flow path ids as part of
    # the netcdf output that will be utilized as a geospatial field within the
    # SCHISM BMI to connect T-Route flowpaths with SCHISM sources/sinks.
    ngen.Source_nc_write(startdate, enddate, output_directory, overwrite=True)
    print(f'It took {time()-t0} seconds to generate source/sink')
