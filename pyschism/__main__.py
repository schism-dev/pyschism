#! /usr/bin/env python
import argparse
from datetime import datetime
import logging

from pytz import timezone

from pyschism.cmd.bctides import BctidesCli
from pyschism.cmd.bootstrap import BootstrapCli
from pyschism.cmd.fgrid.entry import FgridCli
from pyschism.cmd.forecast import ForecastCli
from pyschism.cmd.grd2sms import Grd2SmsCli
from pyschism.cmd.hgrid import HgridCli
from pyschism.cmd.outputs import OutputsCli
from pyschism.cmd.sms2grd import Sms2grdCli
from pyschism.cmd.stations import StationsCli
from pyschism.cmd.tvdflag import TvdflagCli
from pyschism.cmd.vgrid import VgridCli


def init_logger():
    tmp_parser = argparse.ArgumentParser(add_help=False)
    tmp_parser.add_argument(
        "--log-level",
        choices=['info', 'warning', 'debug'],
        default='warning'
    )
    tmp_args, _ = tmp_parser.parse_known_args()
    logging.basicConfig(
        level={
            'warning': logging.WARNING,
            'info': logging.INFO,
            'debug': logging.DEBUG,
        }[tmp_args.log_level],
        format='[%(asctime)s] %(name)s %(levelname)s: %(message)s',
        force=True,
    )

    logging.Formatter.converter = lambda *args: datetime.now(
        tz=timezone('UTC')).timetuple()

    logging.captureWarnings(True)


def main():
    
    init_logger()
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='clitype')
    clitypes = [
        BctidesCli,
        BootstrapCli,
        FgridCli,
        ForecastCli,
        Grd2SmsCli,
        HgridCli,
        OutputsCli,
        Sms2grdCli,
        StationsCli,
        TvdflagCli,
        VgridCli,
    ]

    for clitype in clitypes:
        clitype.add_subparser_action(subparsers)

    args = parser.parse_args()

    for clitype in clitypes:
        if args.clitype == f'{clitype.__name__.lower()[:-3]}':
            clitype(args)
            break


if __name__ == '__main__':
    main()
