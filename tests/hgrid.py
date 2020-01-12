#! /usr/bin/env python
from pyschism.mesh import Hgrid
import unittest


class HgridTestCase(unittest.TestCase):

    def test_triangles_only(self):
        x = [0, 1, 0.5]
        y = [0, 0, 0.5]
        triangles = [[0, 1, 2]]
        Hgrid(x, y, triangles)

    def test_quads_only(self):
        x = [0, 1, 1, 0]
        y = [0, 0, 1, 1]
        quads = [[0, 1, 2, 3]]
        Hgrid(x, y, quads)


if __name__ == '__main__':
    unittest.main()
