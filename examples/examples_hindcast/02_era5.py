from datetime import datetime
from time import time
import pathlib

from pyschism.forcing.nws.nws2.era5 import ERA5

from pyschism.mesh.hgrid import Hgrid

startdate=datetime(2018,7,1)

t0=time()
hgrid=Hgrid.open('hgrid.gr3',crs='EPSG:4326')
bbox = hgrid.get_bbox('EPSG:4326', output_type='bbox')

er=ERA5()
outdir = pathlib.Path('./ERA5')
er.write(outdir=outdir, start_date=startdate, rnday=10, air=True, rad=True, prc=True, bbox=bbox, overwrite=True)
print(f'It took {(time()-t0)/60} minutes to generate 10 days')
