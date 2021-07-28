from datetime import datetime
from time import time
import pathlib

from pyschism.forcing.atmosphere.era5 import ERA5

from pyschism.mesh.hgrid import Hgrid

startdate=datetime(2018,7,1)

t0=time()
hgrid=Hgrid.open('hgrid.gr3',crs='EPSG:4326')
bbox = hgrid.get_bbox('EPSG:4326', output_type='bbox')

er=ERA5()
outdir = pathlib.Path('./ERA5')
er.gen_sflux(startdate, rnday=10, bbox=bbox, outdir=outdir)
print(f'It took {(time()-t0)/60} minutes to generate 10 days')
