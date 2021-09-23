from argparse import Namespace

from pyschism.cmd.fgrid import manning


class FgridCli:

    def __init__(self, args: Namespace):

        if args.action == 'manning':
            manning.ManningsNCli(args)

        else:
            raise NotImplementedError(f'Unhandled CLI action: {args.action}.')

    @staticmethod
    def add_subparser_action(subparsers):
        add_fgrid_options_to_parser(subparsers.add_parser('fgrid'))


def add_fgrid_options_to_parser(parser):
    actions = parser.add_subparsers(dest='action')
    manning.add_manning(actions)
