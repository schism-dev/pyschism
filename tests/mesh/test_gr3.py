#! /usr/bin/env python
import unittest
import pathlib
import tempfile
from pyschism.mesh.gr3 import reader, writer


class Gr3TestCase(unittest.TestCase):

    def setUp(self):
        self.nodes = {
            '1': ((0., 0.), -99999.),
            '2': ((.5, 0.), -99999.),
            '3': ((0., 1.), -99999.),
            '4': ((1., 1.), -99999.),
            '5': ((0., 1.), -99999.),
            '6': ((.5, 1.5), -99999.),
            '7': ((.33, .33), -99999.),
            '8': ((.66, .33), -99999.),
            '9': ((.5, .66), -99999.),
            '10': ((-1., 1.), -99999.),
            '11': ((-1., 0.), -99999.),
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
            'description': 'test_writer_unittest'
        }

    def test_write_read(self):
        tmpdir = tempfile.TemporaryDirectory()
        tmpfile = pathlib.Path(tmpdir.name) / 'hgrid.gr3'
        writer(self.grd, pathlib.Path(tmpfile))
        self.assertDictEqual(reader(pathlib.Path(tmpfile)), self.grd)

    def test_overwrite(self):
        tmpdir = tempfile.TemporaryDirectory()
        writer(self.grd, pathlib.Path(tmpdir.name) / 'hgrid.gr3')
        self.assertRaises(
            Exception,
            writer,
            self.grd,
            pathlib.Path(tmpdir.name) / 'hgrid.gr3'
            )


if __name__ == '__main__':
    unittest.main()
