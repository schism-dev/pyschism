#! /usr/bin/env python
import numpy as np
import tempfile
import pathlib
from pyschism.mesh import Hgrid
import unittest


class HgridTestCase(unittest.TestCase):

    def test_triangule_only_mesh(self):
        nodes = {
            0: (0., 0., np.nan),
            1: (1., 0., np.nan),
            2: (0., 1., np.nan),
        }
        elements = {0: [0, 1, 2]}
        h = Hgrid(nodes, elements)
        self.assertIsInstance(h, Hgrid)

    def test_quads_only_mesh(self):
        nodes = {
            0: (0., 0., np.nan),
            1: (1., 0., np.nan),
            2: (1., 1., np.nan),
            3: (0., 1., np.nan),
        }
        elements = {0: [0, 1, 2, 3]}
        h = Hgrid(nodes, elements)
        self.assertIsInstance(h, Hgrid)

    def test_hybrid_mesh(self):
        nodes = {
            0: (0., 0., np.nan),
            1: (1., 0., np.nan),
            2: (1., 1., np.nan),
            3: (0., 1., np.nan),
            4: (0.5, 1.5, np.nan),
        }
        elements = {
            0: [2, 4, 3],
            1: [3, 0, 1, 2],
        }
        h = Hgrid(nodes, elements)
        self.assertIsInstance(h, Hgrid)

    def test_open_mesh(self):
        nodes = {
            0: (0., 0., np.nan),
            1: (1., 0., np.nan),
            2: (1., 1., np.nan),
            3: (0., 1., np.nan),
            4: (0.5, 1.5, np.nan),
        }
        elements = {
            0: [2, 4, 3],
            1: [0, 1, 2, 3],
        }
        tmpfile = tempfile.NamedTemporaryFile()
        with open(tmpfile.name, 'w') as f:
            f.write('\n')
            f.write(f'{len(elements):d} ')
            f.write(f'{len(nodes):d}\n')
            for id, (x, y, z) in nodes.items():
                f.write(f"{id} ")
                f.write(f"{x} ")
                f.write(f"{y} ")
                f.write(f"{z}\n")
            for id, geom in elements.items():
                f.write(f"{id} ")
                f.write(f"{len(geom)} ")
                for idx in geom:
                    f.write(f"{idx:d} ")
                f.write(f"\n")
        self.assertIsInstance(Hgrid.open(tmpfile.name), Hgrid)

    def test_make_plot(self):
        nodes = {
            0: (0., 0., 0),
            1: (1., 0., 1),
            2: (1., 1., 2),
            3: (0., 1., 3),
            4: (0.5, 1.5, 4),
        }
        elements = {
            0: [2, 4, 3],
            1: [0, 1, 2, 3],
        }
        h = Hgrid(nodes, elements)
        h.make_plot()
        self.assertIsInstance(h, Hgrid)

    def test_make_plot_wet_only(self):
        nodes = {
            0: (0., 0., 0),
            1: (1., 0., -1),
            2: (1., 1., -2),
            3: (0., 1., -3),
            4: (0.5, 1.5, -4),
        }
        elements = {
            0: [2, 4, 3],
            1: [0, 1, 2, 3],
        }
        h = Hgrid(nodes, elements)
        h.make_plot()
        self.assertIsInstance(h, Hgrid)

    def test_dump(self):
        nodes = {
            0: (0., 0., 0),
            1: (1., 0., -1),
            2: (1., 1., -2),
            3: (0., 1., -3),
            4: (0.5, 1.5, -4),
        }
        elements = {
            0: [2, 4, 3],
            1: [0, 1, 2, 3],
        }
        h = Hgrid(nodes, elements)
        tmpdir = tempfile.TemporaryDirectory()
        h.dump(pathlib.Path(tmpdir.name) / 'test_hgrid.gr3')
        self.assertIsInstance(h, Hgrid)


if __name__ == '__main__':
    unittest.main()
