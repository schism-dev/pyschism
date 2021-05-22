#! /usr/bin/env python
import argparse
# from datetime import datetime
import logging
import os
import pathlib
# import shutil
import tarfile
import tempfile
import unittest
import urllib.request

from pyschism.cmd.bctides import add_bctides, BctidesCli


logging.basicConfig(level=logging.INFO, force=True)


class ModelConfigurationTestCase(unittest.TestCase):

    def setUp(self):
        hgrid2D = os.getenv('NWM_TEST_MESH')
        if hgrid2D is None:
            data_directory = pathlib.Path(__file__).parent.absolute() / 'data'
            hgrid2D = data_directory / "GulfStreamDevel/hgrid.gr3"
            url = "https://www.dropbox.com/s/mjaxaqeggy721um/"
            url += "Gulf_Stream_develop.tar.gz?dl=1"
            g = urllib.request.urlopen(url)
            tmpfile = tempfile.NamedTemporaryFile()
            with open(tmpfile.name, 'b+w') as f:
                f.write(g.read())
            with tarfile.open(tmpfile.name, "r:gz") as tar:
                tar.extractall(data_directory / "GulfStreamDevel")
        self.hgrid2D = hgrid2D

    def _test_bctides_cli_2d(self):
        import tempfile
        tmpdir = tempfile.TemporaryDirectory()
        parser = argparse.ArgumentParser()
        add_bctides(parser.add_subparsers(dest='mode'))
        args = parser.parse_args(
            [
                'bctides',
                f'{self.hgrid2D}',  # hgrid path
                '--run-days=5',  # rndays
                '--hgrid-crs=epsg:4326',
                f'-o={tmpdir.name}',
                '--overwrite'
            ])
        BctidesCli(args)
        # with open(tmpdir.name + '/bctides.in') as f:
        #     for line in f:
        #         print(line, end='')

    def _test_bctides_cli_2d_w_velo(self):
        import tempfile
        tmpdir = tempfile.TemporaryDirectory()
        parser = argparse.ArgumentParser()
        add_bctides(parser.add_subparsers(dest='mode'))
        args = parser.parse_args(
            [
                'bctides',
                f'{self.hgrid2D}',  # hgrid path
                '--run-days=5',  # rndays,
                '--include-velocity',
                '--hgrid-crs=epsg:4326',
                f'-o={tmpdir.name}',
                '--overwrite'
            ])
        BctidesCli(args)
        # with open(tmpdir.name + '/bctides.in') as f:
        #     for line in f:
        #         print(line, end='')

    def test_bctides_cli_2d_w_elev2d(self):
        import tempfile
        tmpdir = tempfile.TemporaryDirectory()
        parser = argparse.ArgumentParser()
        add_bctides(parser.add_subparsers(dest='mode'))
        args = parser.parse_args(
            [
                'bctides',
                f'{self.hgrid2D}',  # hgrid path
                '--run-days=5',  # rndays,
                '--iettype-5',
                # '--ifltype-5',
                '--hgrid-crs=epsg:4326',
                f'-o={tmpdir.name}',
                '--overwrite'
            ])
        BctidesCli(args)

    def _test_bctides_cli_3d(self):
        import tempfile
        tmpdir = tempfile.TemporaryDirectory()
        parser = argparse.ArgumentParser()
        add_bctides(parser.add_subparsers(dest='mode'))
        args = parser.parse_args(
            [
                'bctides',
                f'{self.hgrid}',  # hgrid path
                '--run-days=5',  # rndays
                '--hgrid-crs=epsg:4326',
                f'-o={tmpdir.name}',
                '--overwrite'
            ])
        BctidesCli(args)


if __name__ == '__main__':
    unittest.main()
