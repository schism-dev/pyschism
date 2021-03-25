from argparse import Namespace


from pyschism.cmd.hgrid.plot import plot_hgrid, add_plot_hgrid
from pyschism.cmd.hgrid.sms2grd import sms2grd, add_sms2grd
# from pyschism.cmd.hgrid.grd2sms import grd2sms, add_grd2sms


class HgridCli:

    def __init__(self, args: Namespace):

        if args.action == 'plot':
            plot_hgrid(args)

        if args.action == 'sms2grd':
            sms2grd(args)

        # if args.action == 'grd2sms':
        #     grd2sms(args)

        else:
            raise NotImplementedError(f'Unhandled CLI action: {args.action}.')


def hgrid_subparser(subparsers):

    hgrid = subparsers.add_parser('hgrid')

    # define subparser action
    actions = hgrid.add_subparsers(dest='action')

    # plotting action
    add_plot_hgrid(actions)
    add_sms2grd(actions)
    # add_grd2sms(actions)
