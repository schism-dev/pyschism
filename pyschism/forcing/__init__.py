from pyschism.forcing.bctides.tides import Tides
# from pyschism.forcing.tides.bctides import Bctides
from pyschism.forcing.nws.nws2.gfs import GlobalForecastSystem
from pyschism.forcing.source_sink.nwm import NationalWaterModel
from pyschism.forcing.nws.base import NWS

GFS = GlobalForecastSystem
NWM = NationalWaterModel

__all__ = [
    "Tides",
    "GlobalForecastSystem", 'GFS',
    'NationalWaterModel', 'NWM',
    'NWS',
    # 'Bctides',
]
