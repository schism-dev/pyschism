from datetime import datetime
from time import time
import pathlib
import logging

from pyschism.forcing.nws.nws2.era5 import ERA5

from pyschism.mesh.hgrid import Hgrid

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.captureWarnings(True)

log_level = logging.DEBUG
logging.getLogger('pyschism').setLevel(log_level)

if __name__ == '__main__':

    startdate=datetime(2017, 4, 26)
    rnday=18

    t0=time()
    hgrid=Hgrid.open('./hgrid.gr3',crs='EPSG:4326')
    bbox = hgrid.get_bbox('EPSG:4326', output_type='bbox')

    er=ERA5()
    outdir = pathlib.Path('./')
    er.write(outdir=outdir, start_date=startdate, rnday=rnday, air=True, rad=True, prc=True, bbox=bbox, overwrite=True)
    print(f'It took {(time()-t0)/60} minutes to generate {rnday} days')
