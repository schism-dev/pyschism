#! /usr/bin/env python
import argparse
from datetime import datetime
import logging

from pytz import timezone

from pyschism.cmd.forecast.forecast import ForecastCli, add_forecast
from pyschism.cmd.bctides import BctidesCli, add_bctides
from pyschism.cmd.stations import StationsCli, add_stations
from pyschism.cmd.hgrid import HgridCli, add_hgrid
from pyschism.cmd.sms2grd import Sms2grdCli, add_sms2grd
# from pyschism.cmd.grd2sms import Grd2smsCli, add_grd2sms
from pyschism.cmd.fgrid.entry import FgridCli, add_fgrid


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--log-level",
        choices=['info', 'warning', 'debug'],
        default='info'
    )
    subparsers = parser.add_subparsers(dest='mode')
    add_forecast(subparsers)
    add_bctides(subparsers)
    add_stations(subparsers)
    add_hgrid(subparsers)
    add_fgrid(subparsers)
    add_sms2grd(subparsers)
    # add_grd2sms(subparsers)
    return parser.parse_args()


def main():
    args = parse_args()

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

    # elif args.mode == 'grd2sms':
    #     Grd2smsCli(args)

    elif args.mode == 'sms2grd':
        Sms2grdCli(args)

    else:
        raise NotImplementedError(f'Unhandled CLI mode: {args.mode}')


if __name__ == '__main__':
    main()
