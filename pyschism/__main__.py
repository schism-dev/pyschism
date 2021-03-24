#! /usr/bin/env python
import argparse
from datetime import datetime
import logging

from pytz import timezone

from pyschism.cmd.forecast.forecast import ForecastCli, add_forecast
from pyschism.cmd.bctides import BctidesCli, add_bctides
from pyschism.cmd.stations import StationsCli, add_stations
from pyschism.cmd.hgrid import HgridCli, add_hgrid


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

    elif args.mode == 'bctides':
        BctidesCli(args)

    elif args.mode == 'stations':
        StationsCli(args)

    elif args.mode == 'hgrid':
        HgridCli(args)

    else:
        raise NotImplementedError(f'Unhandled CLI mode: {args.mode}')


if __name__ == '__main__':
    main()
