#! /usr/bin/env python
import tempfile
import pathlib
from pyschism.mesh.friction import Fgrid
from pyschism.mesh import Mesh
import unittest


class MeshTestCase(unittest.TestCase):

    def setUp(self):
        nodes = {
            '1': ((0., 0.), -5.),
            '2': ((.5, 0.), -4.),
            '3': ((1., 0.), -3.),
            '4': ((1., 1.), -2.),
            '5': ((0., 1.), -1.),
            '6': ((.5, 1.5), 0.),
            '7': ((.33, .33), 1.),
            '8': ((.66, .33), 2.),
            '9': ((.5, .66), 3.),
            '10': ((-1., 1.), 4.),
            '11': ((-1., 0.), 5.),
            }
        elements = {
            '1': ['5', '7', '9'],
            '2': ['1', '2', '7'],
            '3': ['2', '3', '8'],
            '4': ['8', '7', '2'],
            '5': ['3', '4', '8'],
            '6': ['4', '9', '8'],
            '7': ['4', '6', '5'],
            '8': ['5', '10', '11', '1'],
            '9': ['9', '4', '5'],
            '10': ['5', '1', '7']
            }

        boundaries = dict()

        boundaries[None] = {  # "open" boundaries
                0: ['10', '11', '1', '2'],
                1: ['2', '3', '4']
        }

        boundaries[0] = {  # "land" boundaries
            0: ['4', '6'],
            1: ['6',  '5', '10']
        }

        boundaries[1] = {0: ['7', '8', '9', '7']}  # "interior" boundary

        # write hgrid
        tmpdir = tempfile.TemporaryDirectory()
        hgrid = pathlib.Path(tmpdir.name) / 'hgrid.gr3'
        with open(hgrid, 'w') as f:
            f.write('\n')
            f.write(f'{len(elements):d} ')
            f.write(f'{len(nodes):d}\n')
            for id, ((x, y), z) in nodes.items():
                f.write(f"{id} ")
                f.write(f"{x} ")
                f.write(f"{y} ")
                f.write(f"{z}\n")
            for id, geom in elements.items():
                f.write(f"{id} ")
                f.write(f"{len(geom)} ")
                for idx in geom:
                    f.write(f"{idx} ")
                f.write(f"\n")
        vgrid = pathlib.Path(tmpdir.name) / 'vgrid.gr3'
        with open(vgrid, 'w') as f:
            f.write("2 !ivcor\n")
            f.write("2 1 1.e6 !nvrt, kz (# of Z-levels); h_s (transition depth between S and Z)\n")
            f.write("Z levels\n")
            f.write("1  -1.e6\n")
            f.write("S levels\n")
            f.write("40. 1. 1.e-4  !h_c, theta_b, theta_f\n")
            f.write("   1    -1.\n")
            f.write("   2     0.\n")
        fgrid = pathlib.Path(tmpdir.name) / 'fgrid.gr3'
        with open(fgrid, 'w') as f:
            f.write('generic_friction\n')
            f.write(f'{len(elements):d} ')
            f.write(f'{len(nodes):d}\n')
            for id, ((x, y), z) in nodes.items():
                f.write(f"{id} ")
                f.write(f"{x} ")
                f.write(f"{y} ")
                f.write(f"{z}\n")
            for id, geom in elements.items():
                f.write(f"{id} ")
                f.write(f"{len(geom)} ")
                for idx in geom:
                    f.write(f"{idx} ")
                f.write(f"\n")
        self.tmpdir = tmpdir
        self.hgrid = hgrid
        self.vgrid = vgrid
        self.fgrid = fgrid

    def test_open(self):
        self.assertIsInstance(
            Mesh.open(
                self.hgrid,
                self.vgrid,
                self.fgrid
                ),
            Mesh
        )

    def test_make_plot(self):
        m = Mesh.open(self.hgrid)
        m.make_plot()

    def test_default_fgrid(self):
        m = Mesh.open(self.hgrid)
        assert isinstance(m.fgrid, Fgrid)
