from argparse import Namespace

from pyschism.cmd.fgrid import manning


class FgridCli:

    def __init__(self, args: Namespace):

        if args.action == 'manning':
            manning.ManningsNCli(args)

        else:
            raise NotImplementedError(f'Unhandled CLI action: {args.action}.')


def fgrid_subparser(subparsers):
    fgrid = subparsers.add_parser('fgrid')
    # define subparser action
    actions = fgrid.add_subparsers(dest='action')
    # creating manning
    manning.add_manning(actions)


add_fgrid = fgrid_subparser
