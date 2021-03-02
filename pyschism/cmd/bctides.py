from argparse import Namespace
from datetime import datetime, timedelta
import logging
import pathlib
from time import time

from pyschism.io import Bctides
from pyschism.mesh import Hgrid

_logger = logging.getLogger(__name__)


class BctidesCli:

    def __init__(self, args: Namespace):

        # logging.basicConfig(level=logging._nameToLevel[args.log_level.upper()])
        root_logger = logging.getLogger()
        root_logger.setLevel(logging._nameToLevel[args.log_level.upper()])
        _logger.setLevel(logging._nameToLevel[args.log_level.upper()])

        _logger.info(f'Open hgrid file: {args.hgrid}')
        start = time()
        hgrid = Hgrid.open(args.hgrid, crs=args.hgrid_crs)
        _logger.info(f'Reading hgrid file took {time()-start} seconds.')
        bctides = Bctides(
            hgrid,
            args.start_date,
            args.run_days,
            tidal_database=args.tidal_database,
            velocity=args.include_velocity,
        )

        if args.output_file is not None:
            bctides.write(args.output_file, overwrite=args.overwrite)
        else:
            print(str(bctides))


def add_bctides(subparsers):
    bctides = subparsers.add_parser('bctides')
    bctides.add_argument('hgrid')
    bctides.add_argument(
        'start_date', type=lambda x: datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))
    bctides.add_argument('run_days', type=float)
    bctides.add_argument('--output-file', type=pathlib.Path)
    bctides.add_argument('--include-velocity', action='store_true')
    bctides.add_argument('--tidal-database', '--tidal-db',
                         choices=['hamtide', 'tpxo'],
                         default='hamtide')
    bctides.add_argument('--hgrid-crs')
    bctides.add_argument(
        "--overwrite", action="store_true",
        help="Allow overwrite of output file.")
    bctides.add_argument(
        "--log-level",
        choices=[name.lower() for name in logging._nameToLevel],
        default='warning')
