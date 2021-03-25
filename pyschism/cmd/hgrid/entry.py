from argparse import Namespace


from pyschism.cmd.hgrid.plot import plot_hgrid, add_plot_hgrid


class HgridCli:

    def __init__(self, args: Namespace):

        if args.action == 'plot':
            plot_hgrid(args)

        else:
            raise NotImplementedError(f'Unhandled CLI action: {args.action}.')


def hgrid_subparser(subparsers):

    hgrid = subparsers.add_parser('hgrid')
    actions = hgrid.add_subparsers(dest='action')
    add_plot_hgrid(actions)
