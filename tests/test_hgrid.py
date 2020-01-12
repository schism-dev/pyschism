#! /usr/bin/env python
import numpy as np
from pyschism.mesh import Hgrid
import unittest


class HgridTestCase(unittest.TestCase):

    def test_triangule_only_mesh(self):
        verts = [(0., 0.),
                 (1., 0.),
                 (0., 1.)]
        values = np.random.rand(len(verts),)
        triangles = [[0, 1, 2]]
        Hgrid(verts, values, triangles=triangles)

    def test_quads_only_mesh(self):
        verts = [(0., 0.),
                 (1., 0.),
                 (1., 1.),
                 (0., 1.)]
        values = np.random.rand(len(verts),)
        quads = [[0, 1, 2, 3]]
        Hgrid(verts, values, quads=quads)

    def test_hybrid_mesh(self):
        verts = [(0., 0.),
                 (1., 0.),
                 (1., 1.),
                 (0., 1.),
                 (0.5, 1.5)]
        values = np.random.rand(len(verts),)
        triangles = [[2, 4, 3]]
        quads = [[0, 1, 2, 3]]
        Hgrid(verts, values, triangles=triangles, quads=quads)


if __name__ == '__main__':
    unittest.main()
