from argparse import Namespace

from pyschism.mesh.hgrid import Hgrid
from pyschism.mesh.vgrid import Vgrid


class VgridCli:

    def __init__(self, args: Namespace):
        hgrid = Hgrid.open(args.hgrid, crs=args.hgrid_crs)
        vgrid = Vgrid.from_binary(hgrid)
        vgrid.write('/tmp/temp_vgrid.in')

    @staticmethod
    def add_subparser_action(subparsers):
        add_vgrid_options_to_parser(subparsers.add_parser('vgrid'))


def add_vgrid_options_to_parser(parser):
    parser.add_argument('hgrid')
    parser.add_argument('--hgrid-crs')
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Allow overwrite of output file.")
