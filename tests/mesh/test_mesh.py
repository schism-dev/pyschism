#! /usr/bin/env python
import tempfile
import pathlib
from pyschism.mesh import Mesh, Vgrid
import unittest


class MeshTestCase(unittest.TestCase):

    def setUp(self):
        self.nodes = {
            '1': ((0., 0.), -5.),
            '2': ((.5, 0.), -4.),
            '3': ((0., 1.), -3.),
            '4': ((1., 1.), -2.),
            '5': ((0., 1.), -1.),
            '6': ((.5, 1.5), 0.),
            '7': ((.33, .33), 1.),
            '8': ((.66, .33), 2.),
            '9': ((.5, .66), 3.),
            '10': ((-1., 1.), 4.),
            '11': ((-1., 0.), 5.),
            }
        self.elements = {
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

        self.boundaries = dict()

        self.boundaries[None] = {  # "open" boundaries
                0: ['10', '11', '1', '2'],
                1: ['2', '3', '4']
        }

        self.boundaries[0] = {  # "land" boundaries
            0: ['4', '6'],
            1: ['6',  '5', '10']
        }

        self.boundaries[1] = {0: ['7', '8', '9']}  # "interior" boundary

        self.grd = {
            'nodes': self.nodes,
            'elements': self.elements,
            'boundaries': self.boundaries,
            'description': 'gr3_unittest'
        }

    def test_open(self):
        # write hgrid
        tmpdir = tempfile.TemporaryDirectory()
        hgrid = pathlib.Path(tmpdir.name) / 'hgrid.gr3'
        with open(hgrid, 'w') as f:
            f.write('\n')
            f.write(f'{len(self.elements):d} ')
            f.write(f'{len(self.nodes):d}\n')
            for id, ((x, y), z) in self.nodes.items():
                f.write(f"{id} ")
                f.write(f"{x} ")
                f.write(f"{y} ")
                f.write(f"{z}\n")
            for id, geom in self.elements.items():
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

        self.assertIsInstance(
            Mesh.open(
                hgrid,
                vgrid
                ),
            Mesh
        )
