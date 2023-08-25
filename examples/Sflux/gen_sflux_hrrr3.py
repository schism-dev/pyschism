from datetime import datetime, timedelta
from time import time
import logging

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
    startdate = datetime(2018, 1, 1)
    rnday = 1
    record = 1  #parameter of how many days of time records in each nc file (1-for hindcast, >1 for forecast)
    #startdate = nearest_cycle()
    hgrid = Hgrid.open('./hgrid.gr3', crs='epsg:4326')

    pscr='/sciclone/pscr/lcui01/HRRR/' #tmp dir for saving downloaded grib file
    hrrr = HRRR(start_date=startdate, rnday=rnday, pscr=pscr, record=record, bbox=hgrid.bbox)
    print(f'It took {(time()-t0)/60} mins to process {rnday} days, {record*24} files/day')
