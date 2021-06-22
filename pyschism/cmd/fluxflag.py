import argparse
import logging

from pyschism.cmd import common
from pyschism.mesh.prop import Fluxflag

logger = logging.getLogger(__name__)


class FluxflagCli:

    def __init__(self, args: argparse.Namespace):

        if args.prop_table is not None:
            fluxflag = Fluxflag.from_prop_table(args.hgrid, args.prop_table)

        if args.plot is True:
            fluxflag.make_plot(show=True)

        if args.output_path is not None:
            fluxflag.write(args.output_path, overwrite=args.overwrite)

    @staticmethod
    def add_subparser_action(subparsers):
        add_prop_options_to_parser(subparsers.add_parser('fluxflag'))


def add_prop_options_to_parser(parser):
    common.add_hgrid_to_parser(parser)
    common.add_log_level_to_parser(parser)
    parser.add_argument('--output-path', '-o')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--plot', action='store_true')

    options = parser.add_mutually_exclusive_group(required=True)
    options.add_argument('--prop-table')
    # options.add_argument('--prop-table')
    # table.add_argument('--regtable', dest='table')
    # files = parser.add_mutually_exclusive_group(required=True)
    # files.add_argument(
    #     '--regfiles',
    #     '-r',
    #     help='Generate default tvdflag.prop file.'
    # )
    # files.add_argument(
    #     '--geometry',
    #     nargs='+',
    # )
    # options.add_argument(
    #     '--regfiles',
    #     '--reg',
    #     '-r',
    #     nargs='+',
    #     help='Uses ACE/gredit reg files.'
    # )
    # options.add_argument(
    #     '--propfile',
    #     '--prop',
    #     '-p',
    #     help='Loads a SCHISM prop file.'
    # )
    # options.add_argument(
    #     '--geometry',
    #     '--geom',
    #     '--gpd',
    #     '--json',
    #     '-g',
    #     nargs='+',
    #     help='Any file or files readable by geopandas. Must contain polygons '
    #          'or multipolygon geometries.'
    # )
