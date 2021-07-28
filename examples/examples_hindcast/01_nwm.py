from datetime import datetime

from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.hydrology.nwmhindcast import NationalWaterModel

startdate=datetime(2018,7,3)
rnday=10

hgrid=Hgrid.open('hgrid.gr3', crs='EPSG:4326')

nwm=NationalWaterModel()
nwm.fetch_data(hgrid, startdate, rnday)
nwm.write('./', overwrite=True)
