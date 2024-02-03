from datetime import datetime, timedelta
from time import time
import logging
import pathlib

import numpy as np
import pandas as pd
from matplotlib.transforms import Bbox

from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.nws.nws2.hrrr3 import HRRR
from pyschism.dates import nearest_cycle

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.captureWarnings(True)

log_level = logging.DEBUG
logging.getLogger('pyschism').setLevel(log_level)

if __name__ == "__main__":
    t0 = time()

    #now = datetime.now()
    #last_cycle = np.datetime64(pd.DatetimeIndex([now-timedelta(hours=2)]).floor('6H').values[0], 'h').tolist()
    #print(last_cycle)

    #start = (last_cycle - timedelta(days=1)).replace(hour=0)
    #end = last_cycle + timedelta(days=2)

    #rndelta = end -start
    #rnday = rndelta.total_seconds() / 86400.
    start = datetime(2023, 10, 26)
    rnday = 10
    hgrid = Hgrid.open('./hgrid.gr3', crs='epsg:4326')
    #outdir = path = pathlib.Path('./HRRR_2017_2')

    #pscr - to save downloaded grib2 files and files will be deleted after nc files are created
    pscr='/sciclone/pscr/lcui01/HRRR/'

    hrrr = HRRR(level=2, region='conus', pscr=pscr, bbox=hgrid.bbox)
    hrrr.write(start_date=start, rnday=rnday, air=True, prc=True, rad=True)
    print(f'It took {(time()-t0)/60} mins to process {rnday} days!')
