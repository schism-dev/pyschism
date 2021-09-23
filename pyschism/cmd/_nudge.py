import argparse
import logging
import warnings

from pyproj import CRS

from pyschism.forcing.bctides.nudge import Nudge
from pyschism.forcing.baroclinic import GOFS, RTOFS
from pyschism.mesh import Hgrid, Vgrid

logger = logging.getLogger(__name__)


class NudgeCli:
    def __init__(self, args: argparse.Namespace):
        data_source = {
            "gofs": GOFS,
            "rtofs": RTOFS,
        }[args.baroclinic_database]()
        nudge = Nudge(args.hgrid, rlmax=args.rlmax, rnu_day=args.rnu_day)
        temp = nudge(data_source.temperature, args.vgrid)
        temp.write('TEM_nudge.gr3')
        salt = nudge(data_source.salinity, args.vgrid)
        salt.write('SAL_nudge.gr3')


def add_nudge(subparsers):
    nudge = subparsers.add_parser("nudge")
    nudge.add_argument(
        "hgrid",
        action=HgridAction,
    )
    nudge.add_argument(
        "vgrid",
        action=VgridAction,
    )
    nudge.add_argument(
        "--hgrid-crs",
        action=HgridCrsAction
    )
    nudge.add_argument('--rlmax', type=float, default=1.5)
    nudge.add_argument('--rnu_day', type=float, default=0.25)
    nudge.add_argument(
        '--hycom',
        '--baroclinic-database',
        choices=['rtofs', 'gofs'],
        dest='baroclinic_database',
        default='gofs',
    )
    nudge.add_argument(
        "--log-level",
        choices=[name.lower() for name in logging._nameToLevel],
        default="warning",
    )


class HgridAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hgrid = Hgrid.open(values)
        if len(hgrid.boundaries.open) == 0:
            raise TypeError(f"Hgrid provided {values} contains no open boundaries.")
        setattr(namespace, self.dest, hgrid)


class HgridCrsAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values is not None:
            namespace.hgrid.nodes._crs = CRS.from_user_input(values)
        setattr(namespace, self.dest, values)


class VgridAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, Vgrid.open(values))
