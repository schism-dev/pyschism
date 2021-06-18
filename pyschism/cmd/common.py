import argparse
import logging
from time import time

from pyschism.mesh import Hgrid, Vgrid

logger = logging.getLogger(__name__)


def add_log_level_to_parser(parser):
    parser.add_argument(
        "--log-level",
        choices=[name.lower() for name in logging._nameToLevel],
        default="warning",
    )


class HgridAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if self.const is None:
            tmp_parser = argparse.ArgumentParser(add_help=False)
            tmp_parser.add_argument('--hgrid-crs')
            tmp_args, _ = tmp_parser.parse_known_args()
            logger.info(f'Opening hgrid from {values}...')
            start = time()
            hgrid = Hgrid.open(values, crs=tmp_args.hgrid_crs)
            logger.info(f'Reading hgrid took {time()-start}...')
            if len(hgrid.boundaries.open) == 0:
                raise TypeError(f"Hgrid provided {values} contains no open boundaries.")
            setattr(namespace, self.dest, hgrid)
        else:
            setattr(namespace, self.dest, self.const)


def add_hgrid_to_parser(parser, const=None):
    parser.add_argument(
        "hgrid",
        action=HgridAction,
        const=const,
        help='Path to the SCHISM hgrid file.'
    )
    parser.add_argument(
        '--hgrid-crs',
        help='Coordinate reference system string of hgrid.'
    )


class VgridAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if self.const is None:  # This is the init
            logger.info(f'Opening vgrid from {values}')
            start = time()
            vgrid = Vgrid.open(values)
            logger.info(f'Opening vgrid took {time()-start}')
            setattr(namespace, self.dest, vgrid)
        else:
            setattr(namespace, self.dest, self.const)


def add_vgrid_to_parser(parser, const=None):
    parser.add_argument(
        "--vgrid",
        action=VgridAction,
        default=Vgrid.default(),
        const=const,
    )
