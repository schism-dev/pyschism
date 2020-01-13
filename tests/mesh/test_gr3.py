#! /usr/bin/env python
import unittest
import pathlib
import tempfile
from pyschism.mesh.gr3 import parse_gr3, write_gr3


class Gr3TestCase(unittest.TestCase):

    def test_gr3(self):
        nodes = {
            0: (0., 0., -99999.),
            1: (.5, 0., -99999.),
            2: (0., 1., -99999.),
            3: (1., 1., -99999.),
            4: (0., 1., -99999.),
            5: (.5, 1.5, -99999.),
            6: (.33, .33, -99999.),
            7: (.66, .33, -99999.),
            8: (.5, .66, -99999.),
            9: (-1., 1, -99999.),
            10: (-1, 0., -99999.),
        }
        elements = {
            0: [4, 6, 8],
            1: [0, 1, 6],
            2: [1, 2, 7],
            3: [7, 6, 1],
            4: [2, 3, 7],
            5: [3, 8, 7],
            6: [3, 5, 4],
            7: [4, 9, 10, 0],
            8: [8, 3, 4],
            9: [4, 0, 6]
        }

        boundaries = dict()
        boundaries[None] = {  # "open" boundaries
            0: [9, 10, 0, 1],
            1: [1, 2, 3]
        }
        boundaries[0] = {  # "land" boundaries
            0: [3, 5],
            1: [5,  4, 9]
        }
        boundaries[1] = {0: [6, 7, 8]}  # "interior" boundary
        grd = {
            'nodes': nodes,
            'elements': elements,
            'boundaries': boundaries,
            'description': 'gr3.py unittest'
        }
        tmpdir = tempfile.TemporaryDirectory()
        write_gr3(grd, pathlib.Path(tmpdir.name) / 'hgrid.gr3')
        parse_gr3(pathlib.Path(tmpdir.name) / 'hgrid.gr3')

    def test_gr3_overwrite(self):
        nodes = {
            0: (0., 0., -99999.),
            1: (.5, 0., -99999.),
            2: (0., 1., -99999.),
            3: (1., 1., -99999.),
            4: (0., 1., -99999.),
            5: (.5, 1.5, -99999.),
            6: (.33, .33, -99999.),
            7: (.66, .33, -99999.),
            8: (.5, .66, -99999.),
            9: (-1., 1, -99999.),
            10: (-1, 0., -99999.),
        }
        elements = {
            0: [4, 6, 8],
            1: [0, 1, 6],
            2: [1, 2, 7],
            3: [7, 6, 1],
            4: [2, 3, 7],
            5: [3, 8, 7],
            6: [3, 5, 4],
            7: [4, 9, 10, 0],
            8: [8, 3, 4],
            9: [4, 0, 6]
        }

        boundaries = dict()
        boundaries[None] = {  # "open" boundaries
            0: [9, 10, 0, 1],
            1: [1, 2, 3]
        }
        boundaries[0] = {  # "land" boundaries
            0: [3, 5],
            1: [5,  4, 9]
        }
        boundaries[1] = {0: [6, 7, 8]}  # "interior" boundary
        grd = {
            'nodes': nodes,
            'elements': elements,
            'boundaries': boundaries,
            'description': 'gr3.py unittest'
        }
        tmpdir = tempfile.TemporaryDirectory()
        write_gr3(grd, pathlib.Path(tmpdir.name) / 'hgrid.gr3')
        self.assertRaises(
            Exception,
            write_gr3,
            grd,
            pathlib.Path(tmpdir.name) / 'hgrid.gr3'
            )


if __name__ == '__main__':
    unittest.main()
