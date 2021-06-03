#! /usr/bin/env python
import argparse
from datetime import datetime
import logging

from pytz import timezone

from pyschism.cmd.forecast import ForecastCli
from pyschism.cmd.bctides import BctidesCli
from pyschism.cmd.stations import StationsCli
from pyschism.cmd.hgrid import HgridCli
from pyschism.cmd.sms2grd import Sms2grdCli
from pyschism.cmd.grd2sms import Grd2SmsCli
from pyschism.cmd.fgrid.entry import FgridCli
from pyschism.cmd.vgrid import VgridCli


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--log-level",
        choices=['info', 'warning', 'debug'],
        default='info'
    )
    subparsers = parser.add_subparsers(dest='mode')
    ForecastCli.add_subparser_action(subparsers)
    BctidesCli.add_subparser_action(subparsers)
    StationsCli.add_subparser_action(subparsers)
    HgridCli.add_subparser_action(subparsers)
    Sms2grdCli.add_subparser_action(subparsers)
    Grd2SmsCli.add_subparser_action(subparsers)
    FgridCli.add_subparser_action(subparsers)
    VgridCli.add_subparser_action(subparsers)

    args = parser.parse_args()

    logging.basicConfig(
        level={
            'warning': logging.WARNING,
            'info': logging.INFO,
            'debug': logging.DEBUG,
        }[args.log_level],
        format='[%(asctime)s] %(name)s %(levelname)s: %(message)s',
        force=True,
    )

    logging.Formatter.converter = lambda *args: datetime.now(
        tz=timezone('UTC')).timetuple()

    logging.captureWarnings(True)

    if args.mode == 'forecast':
        ForecastCli(args)

    elif args.mode == 'hgrid':
        HgridCli(args)

    elif args.mode == 'fgrid':
        FgridCli(args)

    elif args.mode == 'bctides':
        BctidesCli(args)

    elif args.mode == 'stations':
        StationsCli(args)

    elif args.mode == 'grd2sms':
        Grd2SmsCli(args)

    elif args.mode == 'sms2grd':
        Sms2grdCli(args)

    elif args.mode == 'vgrid':
        VgridCli(args)

    else:
        raise NotImplementedError(f'Unhandled CLI mode: {args.mode}')


if __name__ == '__main__':
    main()
