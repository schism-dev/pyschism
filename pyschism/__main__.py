#! /usr/bin/env python
import argparse
from datetime import datetime
import inspect
import logging
import sys

from pytz import timezone

from pyschism.cmd import *


def init_logger():
    tmp_parser = argparse.ArgumentParser(add_help=False)
    parser_common.add_log_level_to_parser(tmp_parser)
    tmp_args, _ = tmp_parser.parse_known_args()
    if tmp_args.log_level is not None:
        logging.basicConfig(
            format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
            force=True,
        )

        logging.getLogger("pyschism").setLevel({
                "warning": logging.WARNING,
                "info": logging.INFO,
                "debug": logging.DEBUG,
                "critical": logging.CRITICAL,
                "notset": logging.NOTSET,
            }[tmp_args.log_level])
        logging.Formatter.converter = lambda *args: datetime.now(
            tz=timezone("UTC")
        ).timetuple()

        logging.captureWarnings(True)


def main():

    init_logger()
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="clitype")

    clitypes = [
        clitype
        for cliname, clitype in inspect.getmembers(
            sys.modules[__name__], inspect.isclass
        )
        if "Cli" in cliname
    ]

    for clitype in clitypes:
        clitype.add_subparser_action(subparsers)

    args = parser.parse_args()

    for clitype in clitypes:
        if args.clitype.replace("_", "") == f"{clitype.__name__.lower()[:-3]}":
            clitype(args)
            break


if __name__ == "__main__":
    main()
