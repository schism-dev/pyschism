from datetime import datetime, timedelta
from enum import Enum
import pathlib
import tempfile
from typing import Union
import logging
from time import time

from matplotlib.transforms import Bbox
from netCDF4 import Dataset
import numpy as np
import requests

from pyschism.mesh import Hgrid, Gr3
import geopandas as gpd

logger = logging.getLogger(__name__)

class HYCOM: 
    def __init__(self,hgrid,bbox):
        pass

    def get_data(self, hgrid, start_date: datetime=None,
        rnday: Union[Float, timedelta] = 4,
        ):
  
        bbox = hgrid.get_bbox('EPSG:4326', output_type='bbox')
        #Bbox(x0=-98.0058874, y0=8.5344214, x1=-60.0400016, y1=45.83143077)
        west = bbox.xmin + 360.0
        east = bbox.xmax + 360.0
        north = bbox.ymax
        south = bbox.ymin
      
        loggger.info('Fetching HYCOM data.')

        for i in np.arange(-1,3):
   
            dt = start_date - timedelta(days=i)
            dstr = dt.strftime('%Y-%m-%d-T00:00:00Z')
            var_list = 'surf_el,water_temp,salinity,water_u,water_v'
            url = ('https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0/FMRC/GLBy0.08_930_FMRC_best.ncd'
                +'?var='+var_list
                +'&north='+north+'&south='+south+'&west='+west+'&east='+east
                +'&time='+dstr
                +'&vertCoord=&addLatLon=true&accept=netcdf') 
 
            t0 = time()
            r = requests.get(url)
            with open(fn_out,'wb') as f:
                f.write(r.content) 
        
            logger.info('\nrequests.get took %0.1f seconds' %(time()-t0)) 
 
