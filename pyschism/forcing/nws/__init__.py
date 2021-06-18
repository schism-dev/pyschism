from pyschism.forcing.nws.nws2 import NWS2
from pyschism.forcing.nws.gfs import GlobalForecastSystem
from pyschism.forcing.nws.hrrr import HRRR
from pyschism.forcing.nws.base import NWS


GFS = GlobalForecastSystem


__all__ = [
    'NWS',
    "NWS2",
    'GlobalForecastSystem',
    'GFS',
    'HRRR',
]
