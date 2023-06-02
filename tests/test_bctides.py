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
from datetime import datetime, timedelta

from pyschism.cmd.bctides import BctidesCli
from pyschism.driver import ModelConfig
from pyschism.forcing.bctides import iettype, ifltype
from pyschism.mesh import Hgrid


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

    def test_simple_tidal_setup(self):
        # Inputs: Modify As Needed!
        start_date = "2018091000" #in YYYYMMDDHH format
        end_date = "2018091600" #in YYYYMMDDHH format

        tide_source = "tpxo"

        dt = timedelta(seconds=150.)
        nspool = timedelta(minutes=20.)

        start_date = datetime.strptime(start_date, "%Y%m%d%H")
        end_date = datetime.strptime(end_date, "%Y%m%d%H")
        rnday = end_date - start_date

        hgrid = Hgrid.open(
            'https://raw.githubusercontent.com/geomesh/test-data/main/NWM/hgrid.ll',
            crs=4326
        )

        config = ModelConfig(
            hgrid=hgrid,
            iettype=iettype.Iettype3(database=tide_source),
            ifltype=ifltype.Ifltype3(database=tide_source),
        )

        coldstart = config.coldstart(
            start_date=start_date,
            end_date=start_date + rnday,
            timestep=dt,
            nspool=timedelta(hours=1),
            elev=True,
            dahv=True,
        )

        with tempfile.TemporaryDirectory() as dn:
            tmpdir = pathlib.Path(dn)
            coldstart.write(tmpdir, overwrite=True)


if __name__ == '__main__':
    unittest.main()
