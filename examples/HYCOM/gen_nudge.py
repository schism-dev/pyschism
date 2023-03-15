from datetime import datetime, timedelta
import logging

from pyschism.mesh import Hgrid, Vgrid
from pyschism.forcing.hycom.hycom2schism import Nudge

'''
outputs:
    SAL_nudge.gr3
    SAL_nu.nc
    TEM_nudge.gr3
    TEM_nu.nc
'''

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logger = logging.getLogger('pyschism')
logger.setLevel(logging.INFO)

if __name__ == '__main__':


    hgrid=Hgrid.open('./hgrid.gr3', crs='epsg:4326')
    vgrid='./vgrid.in'
    outdir='./'
    start_date = datetime(2018, 1, 1)
    rnday=30

    #ocean_bnd_ids - segment indices, starting from zero.
    nudge=Nudge(hgrid=hgrid, ocean_bnd_ids=[0, 1, 9])

    #rlmax - max relax distance in m or degree
    #rnu_day - max relax strength in days 
    #restart = True will append to the existing nc file, works when first try doesn't break.
    nudge.fetch_data(outdir, vgrid, start_date, rnday, restart=False, rlmax=1, rnu_day=1)
