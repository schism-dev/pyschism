#! /usr/bin/env python
import tempfile
import random
from pyschism.mesh import Hgrid
import unittest


class HgridTestCase(unittest.TestCase):

    def setUp(self):
        rand = random.SystemRandom()
        self.values = lambda verts: [
            rand.random() for _ in range(len(verts))]

    def test_triangule_only_mesh(self):
        verts = [(0., 0.),
                 (1., 0.),
                 (0., 1.)]
        triangles = [[0, 1, 2]]
        h = Hgrid(verts, self.values(verts), triangles=triangles)
        self.assertIsInstance(h, Hgrid)

    def test_quads_only_mesh(self):
        verts = [(0., 0.),
                 (1., 0.),
                 (1., 1.),
                 (0., 1.)]
        quads = [[0, 1, 2, 3]]
        h = Hgrid(verts, self.values(verts), quads=quads)
        self.assertIsInstance(h, Hgrid)

    def test_hybrid_mesh(self):
        verts = [(0., 0.),
                 (1., 0.),
                 (1., 1.),
                 (0., 1.),
                 (0.5, 1.5)]
        triangles = [[2, 4, 3]]
        quads = [[0, 1, 2, 3]]
        h = Hgrid(verts, self.values(verts), triangles=triangles, quads=quads)
        self.assertIsInstance(h, Hgrid)

    def test_open_mesh(self):
        verts = [(0., 0.),
                 (1., 0.),
                 (1., 1.),
                 (0., 1.),
                 (0.5, 1.5)]
        triangles = [[2, 4, 3]]
        quads = [[0, 1, 2, 3]]
        values = self.values(verts)
        tmpfile = tempfile.NamedTemporaryFile()
        with open(tmpfile.name, 'w') as f:
            f.write('\n')
            elen = len(triangles) + len(quads)
            f.write(f'{elen:d} ')
            f.write(f'{len(verts):d}\n')
            for i, (x, y) in enumerate(verts):
                f.write(f"{i+1} ")
                f.write(f"{verts[i][0]} ")
                f.write(f"{verts[i][1]} ")
                f.write(f"{values[i]}\n")
            _cnt = 0
            for geom in triangles:
                f.write(f"{_cnt+1} ")
                f.write(f"3 ")
                for tri in geom:
                    f.write(f"{tri} ")
                f.write(f"\n")
                _cnt += 1
            for geom in quads:
                f.write(f"{_cnt+1} ")
                f.write(f"4 ")
                for quad in geom:
                    f.write(f"{quad:d} ")
                f.write(f"\n")
                _cnt += 1
        h = Hgrid.open(tmpfile.name)
        self.assertIsInstance(h, Hgrid)


if __name__ == '__main__':
    unittest.main()
