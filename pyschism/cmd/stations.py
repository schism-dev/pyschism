from argparse import Namespace

from pyschism.outputs.stations import StationsOutput


class StationsCli:

    def __init__(self, args: Namespace):

        if args.action == 'plot':
            StationsOutput(
                    args.outputs_directory,
                    args.stations_file,
            ).plot(args.variable, args.station_index, show=True)

        else:
            raise NotImplementedError(f'Unhandled CLI action: {args.action}.')

    @staticmethod
    def add_subparser_action(subparsers):
        add_stations(subparsers.add_parser('stations'))


def add_stations(parser):
    # define subparser action
    actions = parser.add_subparsers(dest='action')
    # plotting action
    plot = actions.add_parser('plot')
    plot.add_argument('variable')
    plot.add_argument('station_index', type=int)
    plot.add_argument('outputs_directory')
    plot.add_argument('--stations-file')
