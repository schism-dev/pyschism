#! /usr/bin/env python
import tempfile
import pathlib
from unittest.mock import patch
from pyschism.mesh import Hgrid
from pyschism.mesh.friction import Fgrid
import unittest


class HgridTestCase(unittest.TestCase):

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

    def test_triangles_only(self):
        self.assertIsInstance(
            Hgrid(
                self.nodes,
                {id: geom for geom in self.elements.values() if len(geom) == 3}
                ),
            Hgrid
            )

    def test_quads_only(self):
        self.assertIsInstance(
            Hgrid(
                self.nodes,
                {id: geom for geom in self.elements.values() if len(geom) == 4}
                ),
            Hgrid
            )

    def test_hybrid(self):
        self.assertIsInstance(Hgrid(self.nodes, self.elements), Hgrid)

    def test_open(self):
        tmpfile = tempfile.NamedTemporaryFile()
        with open(tmpfile.name, 'w') as f:
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
        self.assertIsInstance(Hgrid.open(tmpfile.name), Hgrid)

    @patch('matplotlib.pyplot.show')
    def test_make_plot(self, mock):
        h = Hgrid(self.nodes, self.elements)
        h.make_plot(
            show=True,
            extent=[0, 1, 0, 1],
            title='test',
            cbar_label='elevation [m]',
            vmax=0.
            )
        self.assertIsInstance(h, Hgrid)

    def test_plot_boundary(self):
        h = Hgrid(**self.grd)
        h.plot_boundary(None, 0)
        self.assertIsInstance(h, Hgrid)

    def test_make_plot_wet_only(self):
        nodes = {
            0: ((0., 0.), 0.),
            1: ((1., 0.), -1.),
            2: ((1., 1.), -2.),
            3: ((0., 1.), -3.),
            4: ((0.5, 1.5), -4.),
        }
        elements = {
            0: [2, 4, 3],
            1: [0, 1, 2, 3],
        }
        h = Hgrid(nodes, elements)
        h.make_plot()
        self.assertIsInstance(h, Hgrid)

    def test_write(self):
        h = Hgrid(self.nodes, self.elements)
        tmpdir = tempfile.TemporaryDirectory()
        h.write(pathlib.Path(tmpdir.name) / 'test_hgrid.gr3')
        self.assertIsInstance(h, Hgrid)

    def test_set_friction(self):
        h = Hgrid(self.nodes, self.elements)
        fric = h.set_friction(0.025, 'manning')
        self.assertIsInstance(fric, Fgrid)

    def test_add_custom_boundary_custom(self):
        h = Hgrid(self.nodes, self.elements)
        h.add_boundary_type('ibtype')
        indexes = [('2', '7'), ('3', '8'), ('4', '9')]
        props = {'flow': [1, 2, 3]}
        h.set_boundary_data('ibtype', 0, indexes, properties=props)

    def test_add_boundary_custom_raise(self):
        h = Hgrid(self.nodes, self.elements)
        h.add_boundary_type('ibtype')
        indexes = [('2', '7'), ('3', '10000'), ('4', '9')]
        props = {'flow': [1, 2, 3]}
        self.assertRaises(
            AssertionError,
            h.set_boundary_data,
            'ibtype',
            0,
            indexes,
            **props
        )

    def test_add_boundary(self):
        h = Hgrid(self.nodes, self.elements)
        h.add_boundary_type('ibtype')
        indexes = ['1']
        h.set_boundary_data('ibtype', 0, indexes)

    def test_add_boundary_raise(self):
        h = Hgrid(self.nodes, self.elements)
        h.add_boundary_type('ibtype')
        indexes = ['10000']
        self.assertRaises(
            AssertionError,
            h.set_boundary_data,
            'ibtype',
            0,
            indexes,
        )

    def test_plot_boundaries(self):
        h = Hgrid(self.nodes, self.elements, self.boundaries)
        h.plot_boundaries()


if __name__ == '__main__':
    unittest.main()
