#! /usr/bin/env python
# import tempfile
# import pathlib
from pyschism.mesh.friction import Fgrid
import unittest


class FgridTestCase(unittest.TestCase):

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

        self.fgrid = {
            'nodes': self.nodes,
            'elements': self.elements,
        }

    def test_assert_raises_nchi(self):
        f = Fgrid(**self.fgrid)
        self.assertRaises(Exception, getattr, f, "nchi")
