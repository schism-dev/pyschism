import argparse
import logging

from shapely.geometry import MultiPolygon

from pyschism.cmd import common
from pyschism.mesh.prop import Tvdflag, reg2multipoly


logger = logging.getLogger(__name__)


class TvdflagCli:

    def __init__(self, args: argparse.Namespace):

        # these are mutually exclusive:

        if args.default is not None:
            tvdflag = Tvdflag.default(args.hgrid)

        elif args.regfiles is not None:
            polygon_collection = []
            for regfile in args.regfiles:
                logger.info(f'Reading region file {regfile}...')
                for polygon in reg2multipoly(regfile).geoms:
                    polygon_collection.append(polygon)
            logger.info('Generating tvdflag from geometry...')
            tvdflag = Tvdflag.from_geometry(
                args.hgrid,
                MultiPolygon(polygon_collection)
            )

        elif args.geometry is not None:
            raise NotImplementedError('parse geometry')
            tvdflag = Tvdflag.from_geometry(
                args.hgrid,
                MultiPolygon(polygon_collection)
            )

        elif args.propfile is not None:
            element_values = []
            with open(args.propfile) as f:
                for i in range(len(args.hgrid.elements)):
                    element_values.append(int(f.readline().split()[-1]))
            tvdflag = Tvdflag(args.hgrid, element_values)

        if args.plot is True:
            tvdflag.make_plot(show=True)

        if args.output_path is not None:
            tvdflag.write(args.output_path, overwrite=args.overwrite)

    @staticmethod
    def add_subparser_action(subparsers):
        add_prop_options_to_parser(subparsers.add_parser('tvdflag'))


def add_prop_options_to_parser(parser):
    common.add_hgrid_to_parser(parser)
    common.add_log_level_to_parser(parser)
    parser.add_argument('--output-path', '-o')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--plot', action='store_true')
    options = parser.add_mutually_exclusive_group(required=True)
    options.add_argument(
        '--default',
        '-d',
        help='Generate default tvdflag.prop file.'
    )
    options.add_argument(
        '--regfiles',
        '--reg',
        '-r',
        nargs='+',
        help='Uses ACE/gredit reg files.'
    )
    options.add_argument(
        '--propfile',
        '--prop',
        '-p',
        help='Loads a SCHISM prop file.'
    )
    options.add_argument(
        '--geometry',
        '--geom',
        '--gpd',
        '--json',
        '-g',
        nargs='+',
        help='Any file or files readable by geopandas. Must contain polygons '
             'or multipolygon geometries.'
    )
