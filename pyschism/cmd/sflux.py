import argparse
from datetime import timedelta
import pathlib
import sys

import matplotlib.pyplot as plt

from pyschism.cmd import common
from pyschism.cmd.base import CliComponent
from pyschism.enums import Sflux1Types, Sflux2Types
from pyschism.forcing.nws.nws2 import NWS2
from pyschism.mesh import gridgr3

sflux_variables = {
    "air": ["prmsl", "spfh", "stmp", "uwind", "vwind"],
    "prc": ["prate"],
    "rad": ["dlwrf", "dswrf"],
}


class AnimateSubCli:
    def __init__(self, args):
        kwargs = {}

        if args.sflux_1_glob is not None:
            kwargs.setdefault("sflux_1_glob", args.sflux_1_glob)

        if args.sflux_2_glob is not None:
            kwargs.setdefault("sflux_2_glob", args.sflux_2_glob)

        if pathlib.Path(args.sflux_directory).is_dir() is False:

            class NotADirectory(Exception):
                pass

            raise NotADirectory(f"{args.sflux_directory} is not a directory.")

        nws2 = NWS2.read(args.sflux_directory, **kwargs)

        sflux = getattr(nws2, f"sflux_{args.level}")
        for vargroup, vars in sflux_variables.items():
            component = getattr(sflux, f"{vargroup}")
            for variable in vars:
                animations = []
                if getattr(args, variable) is True:
                    animations.append(
                        getattr(component, f"{variable}").animation(
                            save=args.save / f"{variable}.gif",
                            fps=args.fps,
                        )
                    )

            if args.show is True:
                plt.show()


class GenerateSubCli:
    def __init__(self, args):
        if args.file_interval_sflux_1 is not None:
            args.sflux.sflux_1.inventory.file_interval = args.file_interval_sflux_1
        if args.file_interval_sflux_2 is not None:
            args.sflux.sflux_2.inventory.file_interval = args.file_interval_sflux_2
        args.sflux.write(
            args.output_directory,
            start_date=args.start_date,
            end_date=args.end_date,
            bbox=args.hgrid.get_bbox(crs="epsg:4326"),
            overwrite=args.overwrite,
            windrot=args.windrot,
            air=args.air,
            rad=args.rad,
            prc=args.prc,
        )


class SfluxCli(CliComponent):
    def __init__(self, args: argparse.Namespace):
        exec(f"{args.action.capitalize()}SubCli(args)")

    @staticmethod
    def add_subparser_action(subparsers: argparse._SubParsersAction) -> None:
        parser = subparsers.add_parser("sflux")
        subparsers = parser.add_subparsers(dest="action")
        add_animate_to_subparsers(subparsers)
        add_generate_to_subparsers(subparsers)


def add_animate_to_subparsers(subparsers):
    parser = subparsers.add_parser("animate")
    parser.add_argument("sflux_directory")
    parser.add_argument("--level", choices=["1", "2"], default="1")
    parser.add_argument("--sflux-1-glob")
    parser.add_argument("--sflux-2-glob")
    for vargroup, variables in sflux_variables.items():
        arg_group = parser.add_argument_group(f"{vargroup} variables")
        for variable in variables:
            arg_group.add_argument(f"--{variable}", action="store_true", default=False)
    parser.add_argument("--save", type=lambda x: pathlib.Path(x))
    parser.add_argument("--fps", type=float)
    parser.add_argument(
        "--show",
        action="store_true",
        dest="show",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        dest="show",
    )
    parser.set_defaults(show=False if "--save" in " ".join(sys.argv) else True)


def add_generate_to_subparsers(subparsers):
    parser = subparsers.add_parser("generate")
    # add_sflux_to_parser(parser)
    common.add_hgrid_to_parser(parser)
    # add_sflux_to_parser(parser)
    common.add_sflux_to_parser(parser)
    common.add_dates_to_parser(parser)
    parser.add_argument("--output-directory", "-o", required=True, type=pathlib.Path)
    parser.add_argument("--overwrite", action="store_true", default=False)
    air = parser.add_argument_group("air group").add_mutually_exclusive_group()
    air.add_argument("--air", action="store_true", dest="air")
    air.add_argument("--no-air", action="store_false", dest="air")
    prc = parser.add_argument_group("prc group").add_mutually_exclusive_group()
    prc.add_argument("--prc", action="store_true", dest="prc")
    prc.add_argument("--no-prc", action="store_false", dest="prc")
    rad = parser.add_argument_group("rad group").add_mutually_exclusive_group()
    rad.add_argument("--rad", action="store_true", dest="rad")
    rad.add_argument("--no-rad", action="store_false", dest="rad")
    parser.set_defaults(
        air=True,
        prc=True,
        rad=True,
    )
    common.add_log_level_to_parser(parser)


# def add_sflux_to_parser(parser):

#     class SfluxAction(argparse.Action):

#         def __call__(self, parser, namespace, values, option_string=None):
#             if len(values) > 2:
#                 raise ValueError(
#                     "pyschism sflux generate: error: argument sflux_databases: "
#                     "expected at most two arguments."
#                 )
#             sflux_1 = Sflux1Types[values[0].upper()].value(
#                 product="gfs_0p25_1hr" if values[0].lower() == "gfs" else values[0].lower()
#             )
#             if len(values) == 2:
#                 sflux_2 = Sflux2Types[values[1].upper()].value()
#             else:
#                 sflux_2 = None

#             # tmp_parser = argparse.ArgumentParser(add_help=False)
#             # if not bool(set(sys.argv).intersection(['-h', '--help'])):
#             #     common.add_windrot_to_parser(tmp_parser)
#             #     tmp_args = tmp_parser.parse_known_args()[0]
#             #     windrot = gridgr3.Windrot.default(namespace.hgrid) if tmp_args.windrot is None else tmp_args.windrot
#             # else:
#             #     windrot = None
#             setattr(
#                 namespace,
#                 self.dest,
#                 NWS2(
#                     sflux_1,
#                     sflux_2,
#                     windrot=windrot
#                 )
#             )

#     parser.add_argument(
#         'sflux_databases',
#         nargs="+",
#         action=SfluxAction,
#         metavar="sflux_database",
#         help='Names of the two sflux levels. At least one name is required. Valid choices are:'
#              '`gfs`, `hrrr`.',
#     )
