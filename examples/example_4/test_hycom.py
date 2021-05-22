import numpy as np
from datetime import datetime, timedelta
import logging
from time import time

from netCDF4 import Dataset
from matplotlib.transforms import Bbox

from pyschism.mesh import Hgrid
from pyschism.forcing.hycom.hycom import HotStartInventory

logger = logging.getLogger(__name__)

if __name__ ==  '__main__':
    logger.info('Starting!')
    start_date = datetime.now()
    hgrid = Hgrid.open('../ICOGS3D_dep/hgrid.gr3',crs='EPSG:4326')
    hycom = HotStartInventory()
    t0 = time()
    hycom.fetch_data(hgrid, start_date)
    print(f'Took {time()-t0} seconds')
