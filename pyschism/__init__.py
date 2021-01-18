from importlib import util
# import logging

from pyschism.param import Param
from pyschism.stations import Stations
from pyschism.domain import ModelDomain
from pyschism.driver import ModelDriver

__all__ = ["Param", "ModelDomain", 'ModelDriver', 'Stations']


if util.find_spec("colored_traceback") is not None:
    import colored_traceback  # type: ignore[import]
    colored_traceback.add_hook(always=True)

# logging.basicConfig(
#     level=logging.NOTSET,
#     format='[%(asctime)s] %(name)-13s %(levelname)-8s: %(message)s',
#     )
# logging.getLogger('matplotlib').setLevel(logging.WARNING)
# logging.getLogger('fiona').setLevel(logging.CRITICAL)
# logging.getLogger('rasterio').setLevel(logging.WARNING)
