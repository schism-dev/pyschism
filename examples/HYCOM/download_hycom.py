from datetime import datetime
import logging
from time import time

from matplotlib.transforms import Bbox

from pyschism.forcing.hycom.hycom2schism import DownloadHycom
from pyschism.mesh.hgrid import Hgrid

'''
Download hycom data for Fortran scripts.
Default is to download data for generating initial condition (use hgrid 
    as parameter to cover the whole domain).
Optionally, download data for generating th.nc (bnd=True) and nu.nc (nudge=True) 
    (use bbox as parameter to cut small domain).
'''
logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logger = logging.getLogger('pyschism')
logger.setLevel(logging.INFO)

if __name__ == '__main__':
    date=datetime(2018, 1, 1)
    outdir = './'

    ##example 1 - download data for IC
    #rnday = 1

    #hgrid = Hgrid.open('./hgrid.gr3', crs='epsg:4326')
    #hycom = DownloadHycom(hgrid=hgrid)

    #t0 = time()
    #hycom.fetch_data(date, rnday=rnday, bnd=False, nudge=False, outdir=outdir)
    #print(f'It took {(time()-t0)/60} mins to download')

    #example 2 - download data for bnd
    rnday = 10

    xmin = -55
    xmax = -50
    ymin = 0
    ymax = 53
    bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)

    hycom = DownloadHycom(bbox=bbox)

    t0 = time()
    hycom.fetch_data(date, rnday=rnday, bnd=True, nudge=False, fmt='schism', outdir=outdir)
    print(f'It took {(time()-t0)/60} mins to download')
    #"""

    """example 3 - download data for nudge
    rnday = 20

    xmin = -65
    xmax = -50
    ymin = 0
    ymax = 53
    bbox = Bbox.from_extents(xmin, ymin, xmax, ymax)

    hycom = DownloadHycom(bbox=bbox)

    hycom.fetch_data(date, rnday=rnday, bnd=False, nudge=True, outdir=outdir) 
    print(f'It took {(time()-t0)/60} mins to download') 
    """
