#! /usr/bin/env python
import argparse
from enum import Enum
import pathlib

from pyschism.cmd.forecastd import Forecastd
from pyschism.cmd.viewerd import Viewerd
from pyschism.cmd.plotting import CliOutputPlot


class Dispatch(Enum):
    FORECASTD_DAEMON = Forecastd
    PLOTTING = CliOutputPlot
    VIEWER = Viewerd


class Env(Enum):
    FORECASTD_DAEMON = 'forecastd'
    PLOTTING = 'plot'
    VIEWER = 'viewerd'


def add_forecastd(subparsers):
    forecastd = subparsers.add_parser('forecastd')
    actions = forecastd.add_subparsers(dest='action')
    actions.required = True
    actions.add_parser('start')
    actions.add_parser('stop')
    actions.add_parser('restart')

    add = actions.add_parser('add')
    add.add_argument('hgrid')
    add.add_argument('--vgrid')
    add.add_argument('--fgrid')
    add.add_argument(
        "-c", "--constituents",
        action='append',
        choices=["K1", "O1", "P1", "Q1", "MM", "Mf", "M4", "MN4", "MS4",
                 "2N2", "S1", "all", "major"],
        dest='constituents',
        help="Tidal constituent to be forced in the model. Pass "
             "--use-constituent='all' to use all available constituents "
             "(K1, O1, P1, Q1, MM, Mf, M4, MN4, MS4, 2N2, S1) "
             "or use --constituents='major'. "
             "For a custom list of forcing constituents, pass -c= for each "
             "individual constituent to use (case-insensitive). "
             "Use None for no tidal forcing. Defaults to 'all'."
    )
    add.add_argument("-o", "--output-directory", dest='outdir')
    add.add_argument("--air", "--air-1", dest='air_1')
    add.add_argument("--prc", "--prc-1", dest='prc_1')
    add.add_argument("--rad", "--rad-1", dest='rad_1')
    add.add_argument("--air-2")
    add.add_argument("--prc-2")
    add.add_argument("--rad-2")


def add_viewerd(subparsers):
    viewerd = subparsers.add_parser('viewerd')
    actions = viewerd.add_subparsers(dest='action')
    actions.required = True
    start = actions.add_parser('start')
    start.add_argument('--deploy', action='store_true')
    actions.add_parser('stop')
    actions.add_parser('restart')


def add_autodocd(subparsers):
    autodocd = subparsers.add_parser('autodocd')
    actions = autodocd.add_subparsers(dest='action')
    actions.required = True
    actions.add_parser('start')
    actions.add_parser('stop')
    actions.add_parser('restart')


def add_plot(subparsers):
    plot = subparsers.add_parser('plot')
    # output directory pointer
    plot.add_argument(
        'resource', type=pathlib.Path,
        help='Filename or directory containing SCHISM outputs.')
    # which variable to plot
    plot.add_argument('variable', help='Name of variable to plot.')
    # Do you want to plot a surface or stations?
    output_type = plot.add_subparsers(dest='output_type')
    output_type.required = True
    # surface plot options
    surface_plot = output_type.add_parser('surface', help="Plot SCHISM surface"
                                          " outputs.")
    surface_plot.add_argument("--start")
    # station plot options
    stations = output_type.add_parser('stations', help="Plot SCHISM stations "
                                      "outputs.")
    stations.add_argument('-s', '--station', nargs='*', type=int,
                          help='Station index in station.in to include in '
                          'plot.')


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='mode')
    add_forecastd(subparsers)
    add_plot(subparsers)
    # add_viewerd(subparsers)
    # add_autodocd(subparsers)
    return parser.parse_args()


def main():
    args = parse_args()
    Dispatch[Env(args.mode).name].value(args)


def init():
    if __name__ == '__main__':
        main()


init()
