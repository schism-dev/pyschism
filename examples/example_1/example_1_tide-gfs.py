#! /usr/bin/env python
import pathlib
from datetime import datetime, timedelta, timezone
import logging
from time import time

from pyschism.mesh import Hgrid, Vgrid, Fgrid
from pyschism import ModelDomain, ModelDriver, Stations
from pyschism.forcing import Tides
from pyschism.forcing.atmosphere.nws.nws2 import NWS2
from pyschism.forcing.atmosphere.gfs import GlobalForecastSystem as GFS

#logging.basicConfig(filename='test.log',
#                    level=logging.INFO, 
#                    format='%(asctime)s:%(levelnames'
#                   )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

SIMPLE_SLURM_DRIVER = """#!/bin/bash --login
#SBATCH -D .
#SBATCH -J schism-test
#SBATCH -A nosofs
#SBATCH --mail-type=all
#SBATCH --mail-user=jaime.calzada@noaa.gov
#SBATCH --output=slurm.log
#SBATCH -n 500
#SBATCH --time=02:00:00
#SBATCH --partition=orion

set -e

module load intel/2020 impi/2020 netcdf/4.7.2-parallel

PATH=$HOME/SCHISM/schism/build/bin:$PATH

main() {
  mkdir -p outputs
  time srun pschism_TVD-VL
}

main
"""

# https://eev.ee/blog/2012/05/23/python-faq-descriptors/

PARENT = pathlib.Path(__file__).parent


if __name__ == '__main__':
    # open gr3 file
    logger.info('Reading hgrid file...')
    _tic = time()
    hgrid = Hgrid.open(PARENT / 'hgrid.gr3',crs='EPSG:4326')
    logger.info(f'Reading hgrind file took {time()-_tic}.')

    vgrid = Vgrid()
    fgrid = Fgrid.open(PARENT / 'drag.gr3',crs='EPSG:4326')

    # setup model domain
    domain = ModelDomain(hgrid, vgrid, fgrid)
    logger.info('Model domain setup finished')

    # set tidal boundary conditions
    elevbc = Tides()
    elevbc.use_all()  # activate all forcing constituents
    logger.info('Tidal boundary setup finished')

    # connect the boundary condition to the domain
    domain.add_boundary_condition(elevbc)

    sflux_1 = GFS()
    # sflux_2 = HWRF()

    atmos = NWS2(
             sflux_1,
    #         sflux_2
         )

    domain.set_atmospheric_forcing(atmos)

    #  ------ Param options
    # dt and rnday are required arguments

    # Use int or float for seconds, or timedelta objects for pythonic
    # specifications
    dt = timedelta(seconds=150.)
    rnday = timedelta(days=4.)

    # Use an integer for number of steps or a timedelta to approximate
    # number of steps internally based on timestep
    nspool = timedelta(minutes=30.)

    # The dramp and start_date are optional parameters.
    # start_date is required when using forcings that
    # require some natural dates, (e.g. tides, winds, T&S ...)
    # If using a forcing that requires a start date is invoked, and no start
    # date is provided, the driver will throw an exception.
    # tzinfo is optional, UTC assumed if not provided
    dramp = 0.1 * rnday
    dramp = timedelta(0.0)

    start_date = datetime.utcnow() - dramp

    # Now we add station outputs
#    stations = Stations.from_file(PARENT / 'station.in', timedelta(minutes=6.),
#                                  elev=True, u=True, v=True)

    # init the model driver
    driver = ModelDriver(domain, dt, rnday, dramp=dramp, start_date=start_date,
                         #stations=stations, nspool=nspool, elev=True)
                         nspool=nspool, elev=True)

    # Output requests, as well as most other namelist variables, can be
    # modified post instantiation through their corresponding properties.
    # For example, the following line activates the depth averaged horizontal
    # velocity output request:
    driver.param.schout.dahv = True

    # write files to disk
    outdir = pathlib.Path('staging')
    driver.write(outdir, overwrite=True)

