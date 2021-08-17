from abc import ABC, abstractmethod
import argparse


class CliComponent(ABC):

    @abstractmethod
    def __init__(self, args: argparse.Namespace):
        self.args = args

    @staticmethod
    @abstractmethod
    def add_subparser_action(subparsers: argparse._SubParsersAction) -> None:
        raise NotImplementedError('Method must be implemented by subclass.')
