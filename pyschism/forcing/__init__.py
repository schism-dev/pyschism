from pyschism.forcing.tides import Tides
# from pyschism.forcing.tides.bctides import Bctides
from pyschism.forcing.atmosphere import GlobalForecastSystem
from pyschism.forcing.hydrology.base import Hydrology
from pyschism.forcing.hydrology.nwm import NationalWaterModel
from pyschism.forcing.atmosphere.nws.nws import NWS

GFS = GlobalForecastSystem
NWM = NationalWaterModel

__all__ = [
    "Tides",
    "GlobalForecastSystem", 'GFS',
    'NationalWaterModel', 'NWM',
    "Hydrology",
    'NWS',
    # 'Bctides',
]
