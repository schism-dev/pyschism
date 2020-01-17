import sys
import logging
from pyschism.cmd import argument_parser
from pyschism.cmd.basecmd import SchismBaseCommand


class TidalRunCommand(SchismBaseCommand):
    """  Wrapper for SCHISM tidal only run CLI tool.  """


def main():
    parser = argument_parser.get_parser(
        'tidal',
        description="SCHISM tidal only run."
        )
    args = parser.parse_args()
    if len(args.constituents) == 0:
        args.constituents = ['all']
    logging.basicConfig(level=args.log_level)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("paramiko").setLevel(logging.WARNING)
    drv = TidalRunCommand(args)
    exit(drv.run())


# https://medium.com/opsops/how-to-test-if-name-main-1928367290cb
def init():
    if __name__ == "__main__":
        sys.exit(main())


init()
