from time import time
from datetime import datetime
import pathlib

from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.nws.nws2.hrrr2 import AWSZarrInventory

from pyschism.dates import nearest_cycle

'''
Please use it with caution, some variables may have missing data, see:
https://mesowest.utah.edu/html/hrrr/zarr_documentation/html/zarr_variables.html
'''

if __name__ == "__main__":
    t0 = time()
    startdate = datetime(2023, 10, 3)

    #bbox from hgrid
    #hgrid = Hgrid.open('./hgrid.gr3', crs='epsg:4326')
    #bbox = hgrid.bbox
    #or
    ##bbox from list [xmin, ymin, xmax, ymax]
    bbox = [-98, 7.8, -52.986, 42.09]

    outdir = './'
    
    hrrr = AWSZarrInventory(bbox=bbox)
    hrrr.write(outdir=outdir, start_date=startdate, rnday=2, air=True, prc=True, rad=False)
    print(f'It took {(time() - t0)/60} mins')

