from pyschism.mesh.hgrid import Hgrid
from pyschism.mesh.gridgr3 import ElevIc

if __name__ == '__main__':

    hgrid = Hgrid.open('./hgrid.gr3', crs='epsg:4326') 

    elevic = ElevIc.default(hgrid=hgrid, offset=0.1)
    elevic.write('elev.ic', overwrite=True)
