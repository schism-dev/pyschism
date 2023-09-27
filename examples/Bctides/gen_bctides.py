from time import time
import os
from datetime import datetime, timedelta
import logging
import pathlib

from pyschism.mesh import Hgrid
from pyschism.forcing.bctides import Bctides

logging.basicConfig( 
    format = "[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force = True,
)
logging.getLogger("pyschism").setLevel(logging.DEBUG)

if __name__ == "__main__":

    '''
    Assume files are already located in:
    database='fes2014'
        ~/.local/share/fes2014/eastward_velocity/
        ~/.local/share/fes2014/northward_velocity/  
        ~/.local/share/fes2014/ocean_tide_extrapolated/
    database = 'tpxo'
        ~/.local/share/tpxo/
    '''

    start_date = datetime(2014, 12, 1)
    rnday = 397
    outdir = './'
    hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")

    bctides=Bctides(
        hgrid, 
        flags = [[5, 5, 4, 4], [5, 5, 4, 4], [0, 1, 2, 2]],
        constituents = 'major',
        database = 'fes2014',
        tthconst = 10.0,
        sthconst = 0.0,
    )

    bctides.write(
        outdir, 
        start_date=start_date, 
        rnday=rnday, 
        overwrite=True,
    )
