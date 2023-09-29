from time import time
import os
from datetime import datetime, timedelta
import logging

import numpy as np

from pyschism.mesh import Hgrid
from pyschism.forcing.bctides import Bctides

logging.basicConfig(
    format = "[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
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

    outdir = './'
    hgrid = Hgrid.open('./hgrid.gr3', crs="epsg:4326")
    start_date = datetime(2014, 12, 1)
    rnday = 397
    
    #bctypes for each open boundary, here shows [ocean bnd1, ocean bnd2, river bnd]
    flags = [[5, 5, 4, 4], [5, 5, 4, 4], [0, 1, 1, 2]]
    constituents = 'major'
    database = 'fes2014'
    earth_tidal_potential = True
    #ethconst = [0.5, 0.5, 0.5]
    #vthconst = [-100, -100, -100]
    #tthconst = [10.0, 10.0, 10.0]
    sthconst = [np.nan, np.nan, 0.0]
    tobc = [0.5, 0.5, 1]
    sobc = [0.5, 0.5, 1]
    #relax = [0.5, 0.5]

    bctides=Bctides(
        hgrid = hgrid, 
        flags = flags,
        constituents = constituents, 
        database = database,
        add_earth_tidal = earth_tidal_potential,
        #ethconst = ethconst,
        #vthconst = vthconst,
        #tthconst = tthconst,
        sthconst = sthconst,
        tobc = tobc,
        sobc = sobc,
        #relax = relax,
    )

    bctides.write(
        outdir, 
        start_date=start_date, 
        rnday=rnday, 
        overwrite=True,
    )
