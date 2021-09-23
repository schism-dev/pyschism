import argparse
import logging

from pyschism.cmd import common


class BootstrapCli:

    def __init__(self, args: argparse.Namespace):
        self.args = args

    @staticmethod
    def add_subparser_action(subparsers):
        add_bootstrap_options_to_parser(subparsers.add_parser('bootstrap'))

    @property
    def args(self):
        return self._args

    @args.setter
    def args(self, args: argparse.Namespace):
        assert isinstance(args, argparse.Namespace)
        self._args = args


def add_bootstrap_options_to_parser(parser):
    parser.add_argument(
        'project_directory',
        help='Directory where project files will be written to.',
    )
    parser.add_argument(
        "hgrid",
        help='Path to the SCHISM hgrid file.'
    )
    parser.add_argument(
        '--hgrid-crs',
        help='Coordinate reference system string of hgrid.'
    )
    parser.add_argument(
        "--vgrid",
    )
    parser.add_argument(
        "--log-level",
        choices=[name.lower() for name in logging._nameToLevel],
        default="warning",
    )
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument(
        '--tides',
        '--tidal-database',
        choices=['hamtide', 'tpxo'],
        default='hamtide',
        dest='tidal_database'
    )
    parser.add_argument(
        '--hycom',
        '--baroclinic-database',
        choices=['rtofs', 'gofs'],
        dest='baroclinic_database',
        default='gofs',
    )
    parser.add_argument(
        '--include-velocity',
        action='store_true'
    )
    parser.add_argument('--Z0', type=float)
    parser.add_argument("--cutoff-depth", type=float, default=50.0)
    common.add_ibctype_to_parser(parser)
