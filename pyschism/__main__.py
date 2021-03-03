#! /usr/bin/env python
import argparse

from pyschism.cmd.forecast.forecast import ForecastCli, add_forecast
from pyschism.cmd.bctides import BctidesCli, add_bctides


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='mode')
    # add_forecastd(subparsers)
    add_forecast(subparsers)
    # add_plot(subparsers)
    add_bctides(subparsers)
    # add_viewerd(subparsers)
    # add_autodocd(subparsers)
    return parser.parse_args()


def main():
    args = parse_args()

    if args.mode == 'forecast':
        ForecastCli(args)

    elif args.mode == 'bctides':
        BctidesCli(args)

    else:
        raise NotImplementedError(f'Unhandled CLI mode: {args.mode}')


if __name__ == '__main__':
    main()
