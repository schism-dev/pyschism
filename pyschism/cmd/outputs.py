from argparse import Namespace

from pyschism.outputs.outputs import OutputsCollector


class OutputsCli:

    def __init__(self, args: Namespace):
        outputs = OutputsCollector(args.outputs_directory)
        getattr(outputs, args.variable_name).make_plot(show=True)

    @staticmethod
    def add_subparser_action(subparsers):
        add_outputs_options_to_parser(subparsers.add_parser('outputs'))


def add_outputs_options_to_parser(parser):
    parser.add_argument('outputs_directory')
    parser.add_argument('variable_name')
    parser.add_argument('step', type=int)
