#! /usr/bin/env python
import unittest
import argparse
from unittest.mock import patch
import sys
from pyschism.cmd import argument_parser


class PlotMeshCmdTestCase(unittest.TestCase):

    def test_add_general_options(self):
        parser = argparse.ArgumentParser()
        with patch.object(sys, 'argv', ["", "NULL"]):
            self.assertIsNone(argument_parser.add_general_options(parser))

    def test_add_general_options_case_2(self):
        parser = argparse.ArgumentParser()
        with patch.object(sys, 'argv', ["", "NULL", "--hostname", "NULL"]):
            self.assertIsNone(argument_parser.add_server_options(parser))

    def test_add_tidal_run_options(self):
        parser = argparse.ArgumentParser()
        cmd = ["", "NULL", "start_date", "end_date", "--spinup-days=0"]
        with patch.object(sys, 'argv', cmd):
            self.assertIsNone(
                argument_parser.add_general_options(parser, 'tidal')
                )

    def test_add_best_track_run_options(self):
        parser = argparse.ArgumentParser()
        cmd = ["", "NULL", "storm_id", "--start-date", "--end-date"]
        with patch.object(sys, 'argv', cmd):
            self.assertIsNone(
                argument_parser.add_general_options(parser, 'best_track')
                )


if __name__ == '__main__':
    unittest.main()
