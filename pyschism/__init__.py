from importlib import util

# from pyschism.param.param import Param
from pyschism.stations import Stations
from pyschism.driver import ModelDriver, ModelConfig

__all__ = ['ModelConfig', 'ModelDriver', 'Stations']


if util.find_spec("colored_traceback") is not None:
    import colored_traceback  # type: ignore[import]
    colored_traceback.add_hook(always=True)
