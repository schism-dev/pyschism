import argparse
import pathlib
import sys

import matplotlib.pyplot as plt

from pyschism.cmd.base import CliComponent
from pyschism.forcing.nws.nws2 import NWS2, SfluxDataset

sflux_variables = {
    "air": ["prmsl", "spfh", "stmp", "uwind", "vwind"],
    "prc": ["prate"],
    "rad": ["dlwrf", "dswrf"],
}


class SfluxCli(CliComponent):
    def __init__(self, args: argparse.Namespace):
        kwargs = {}
        if args.sflux_1_glob is not None:
            kwargs.setdefault("sflux_1_glob", args.sflux_1_glob)
        if args.sflux_2_glob is not None:
            kwargs.setdefault("sflux_2_glob", args.sflux_2_glob)
        nws2 = NWS2.read(args.sflux_directory, **kwargs)
        if args.action == "animate":
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

    @staticmethod
    def add_subparser_action(subparsers: argparse._SubParsersAction) -> None:
        parser = subparsers.add_parser("sflux")
        parser.add_argument("sflux_directory")
        parser.add_argument("--level", choices=["1", "2"], default="1")
        parser.add_argument("--sflux-1-glob")
        parser.add_argument("--sflux-2-glob")
        subparsers = parser.add_subparsers(dest="action")
        add_animate_to_subparsers(subparsers)


def add_animate_to_subparsers(subparsers):
    parser = subparsers.add_parser("animate")
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
