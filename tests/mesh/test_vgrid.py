#! /usr/bin/env python
import tempfile
import pathlib
import unittest
from pyschism.mesh import Vgrid


class VgridTestCase(unittest.TestCase):

    def test_init(self):
        v = Vgrid()
        v.boilerplate_2D
        self.assertIsInstance(v, Vgrid)

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
