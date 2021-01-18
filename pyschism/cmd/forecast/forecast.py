from argparse import Namespace
from enum import Enum

from pyschism.cmd.forecast.init import ForecastInit
from pyschism.cmd.forecast.update import ForecastUpdate


class Dispatch(Enum):
    INIT = ForecastInit
    UPDATE = ForecastUpdate


class Env(Enum):
    INIT = 'init'
    UPDATE = 'update'


class ForecastCli:

    def __init__(self, args):
        Dispatch[Env(args.action).name].value(args)
        if args.action == "init":
            ForecastUpdate(Namespace(
                project_directory=args.project_directory,
                log_level=args.log_level))
