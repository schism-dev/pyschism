from datetime import datetime
import logging

from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.hycom.hycom2schism import OpenBoundaryInventory

'''
outputs:
    elev2D.th.nc (elev=True)
    SAL_3D.th.nc (TS=True)
    TEM_3D.th.nc (TS=True)
    uv3D.th.nc   (UV=True)
'''

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logger = logging.getLogger('pyschism')
logger.setLevel(logging.INFO)

if __name__ == "__main__":
    hgrid=Hgrid.open('./hgrid.gr3', crs='epsg:4326')
    vgrid='./vgrid.in'
    outdir='./'
    start_date=datetime(2018, 8, 1)
    rnday=20

    #boundary
    bnd=OpenBoundaryInventory(hgrid, vgrid)

    #ocean_bnd_ids - segment indices, starting from zero.
    bnd.fetch_data(outdir, start_date, rnday, elev2D=True, TS=True, UV=True, ocean_bnd_ids=[0,1,2])
