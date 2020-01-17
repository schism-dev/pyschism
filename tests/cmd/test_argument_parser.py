#! /usr/bin/env python
import unittest
import argparse
from unittest.mock import patch
import sys
from pyschism.cmd import argument_parser


class PlotMeshCmdTestCase(unittest.TestCase):

    def test_add_tidal_run_options(self):
        cmd = ["", "NULL", "start_date", "end_date", "--spinup-days=0"]
        with patch.object(sys, 'argv', cmd):
            self.assertIsInstance(
                argument_parser.get_parser('tidal'),
                argparse.ArgumentParser,
                )

    def test_add_best_track_run_options(self):
        cmd = ["", "NULL", "storm_id", "--start-date", "--end-date"]
        with patch.object(sys, 'argv', cmd):
            self.assertIsInstance(
                argument_parser.get_parser('best_track'),
                argparse.ArgumentParser
                )

    def test_nproc_required(self):
        parser = argparse.ArgumentParser()
        with patch.object(sys, 'argv', ["", "NULL", "--hostname", "NULL"]):
            self.assertIsNone(argument_parser.server(parser))


if __name__ == '__main__':
    unittest.main()
