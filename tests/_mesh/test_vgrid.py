#! /usr/bin/env python
import tempfile
import pathlib
import unittest
from pyschism.mesh import Vgrid
from pyschism.mesh.hgrid import Hgrid
from pyschism.mesh.vgrid import LSC2
import numpy as np

class VgridTestCase(unittest.TestCase):
    def setUp(self):
        nodes = {
            '1': ((0., 0.), 1.5),
            '2': ((.5, 0.), 2.5),
            '3': ((1., 0.), 3.5),
            '4': ((1., 1.), 0.),
            '5': ((0., 1.), 1.),
            '6': ((.5, 1.5), -9.),
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

        hsm=[1,2,3,6]
        nv=[2,4,5,5]
        theta_b=0
        theta_f=2.5
        h_c=2

        self.hgrid=Hgrid.open(hgrid,crs=4326)
        self.hsm=hsm
        self.nv=nv
        self.theta_b=theta_b
        self.theta_f=theta_f
        self.h_c=h_c

    def test_init(self):
        v = Vgrid()
        v.boilerplate_2D
        self.assertIsInstance(v, Vgrid)

    def test_LSC2(self):
        lsc2_obj = LSC2(self.hsm,self.nv,self.h_c,self.theta_b,self.theta_f)
        self.assertIsInstance(lsc2_obj, LSC2)

    def test_calc_m_grid(self):
        lsc2_obj = LSC2(self.hsm,self.nv,self.h_c,self.theta_b,self.theta_f)
        lsc2_obj.calc_m_grid()
        self.assertIsInstance(lsc2_obj.m_grid.shape, (4,5))

    def test_make_m_plot(self):
        lsc2_obj = LSC2(self.hsm,self.nv,self.h_c,self.theta_b,self.theta_f)
        lsc2_obj.calc_m_grid()
        lsc2_obj.make_m_plot()

    def test_calc_lsc2_att(self):
        gr3=self.hgrid
        lsc2_obj = LSC2(self.hsm,self.nv,self.h_c,self.theta_b,self.theta_f)
        lsc2_obj.calc_m_grid()
        lsc2_obj.calc_lsc2_att(gr3)
        self.assertIsInstance(lsc2_obj._znd.shape,(11,5))
        self.assertIsInstance(lsc2_obj._snd.shape,(11,5))
        self.assertIsInstance(lsc2_obj._nlayer,np.array([4,5,5,2,2,2,2,4,5,5,5]))

    def test_lsc2_write(self):
        tmpdir = tempfile.TemporaryDirectory()
        gr3=self.hgrid
        lsc2_obj = LSC2(self.hsm,self.nv,self.h_c,self.theta_b,self.theta_f)
        lsc2_obj.calc_m_grid()
        lsc2_obj.calc_lsc2_att(gr3)
        lsc2_obj.write(pathlib.Path(tmpdir.name) / 'vgrid_lsc2.in')

    def test_write(self):
        tmpdir = tempfile.TemporaryDirectory()
        v = Vgrid()
        v.write(pathlib.Path(tmpdir.name) / 'vgrid.in')
        self.assertIsInstance(v, Vgrid)

    def test_write_overwrite_raises(self):
        tmpfile = tempfile.NamedTemporaryFile()
        v = Vgrid()
        self.assertRaises(
            Exception,
            v.write,
            pathlib.Path(tmpfile.name)
        )
