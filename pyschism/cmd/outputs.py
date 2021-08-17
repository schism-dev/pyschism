from argparse import Namespace
import sys

from pyschism.cmd import common
from pyschism.outputs.outputs import OutputsCollector


class OutputsCli:

    def __init__(self, args: Namespace):
        outputs = OutputsCollector(args.outputs_directory)
        if args.action == 'plot':
            getattr(outputs, args.variable_name).make_plot(show=True)

        elif args.action == 'animate':
            getattr(outputs, args.variable_name).animation(
                vmin=args.vmin,
                vmax=args.vmax,
                show=args.show,
                fps=args.fps,
                save=args.save,
            )

    @staticmethod
    def add_subparser_action(subparsers):
        add_outputs_options_to_parser(subparsers.add_parser('outputs'))


def add_outputs_options_to_parser(parser):

    parser.add_argument('outputs_directory')
    parser.add_argument('variable_name')
    subparsers = parser.add_subparsers(dest='action')
    add_plot_action(subparsers)
    add_animate_action(subparsers)


def add_plot_action(subparsers):
    parser = subparsers.add_parser('plot')
    parser.add_argument('step')


def add_animate_action(subparsers):
    parser = subparsers.add_parser('animate')
    common.add_vmin_to_parser(parser)
    common.add_vmax_to_parser(parser)
    show = parser.add_mutually_exclusive_group() 
    show.add_argument('--show', dest='show', action='store_true')
    show.add_argument('--no-show', dest='show', action='store_false')
    parser.set_defaults(show=False if "--save" in " ".join(sys.argv) else True)
    parser.add_argument('--fps', type=float)
    parser.add_argument('--save')
    parser.add_argument('--figsize', nargs=2, type=float)
