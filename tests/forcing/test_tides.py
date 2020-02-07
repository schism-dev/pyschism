#! /usr/bin/env python
import unittest
from datetime import datetime, timedelta
from pyschism.forcing import Tides


class TidesTestCase(unittest.TestCase):

    def setUp(self):
        self.tf = Tides()
        self.tf.start_date = datetime(2017, 9, 23)
        self.tf.end_date = datetime(2017, 9, 25)
        self.tf.spinup_time = timedelta(days=7)
        self.tf.use_constituent('M2')

    def test_use_ll(self):
        self.tf.use_all()

    def test___call__(self):
        self.tf('M2')

    def test___iter__(self):
        for c, v in self.tf:
            pass

    def test___len__(self):
        len(self.tf)

    def test_use_major(self):
        self.tf.use_major()

    def test_drop_constituent(self):
        self.tf.drop_constituent('M2')

    def test_get_active_constituents(self):
        self.tf.get_active_constituents()

    def test_get_nodal_factor(self):
        for f in self.tf.orbital_frequencies:
            self.tf.get_nodal_factor(f)

    def test_get_greenwich_factor(self):
        for f in self.tf.orbital_frequencies:
            self.tf.get_greenwich_factor(f)

    def test_get_nodal_factor_raise(self):
        self.assertRaises(TypeError, self.tf.get_nodal_factor, '')

    def test_get_greenwich_factor_raise(self):
        self.assertRaises(TypeError, self.tf.get_greenwich_factor, '')

    def test_reset_dates(self):
        new_date = self.tf.start_date
        del(self.tf.start_date)
        self.tf.end_date = new_date

    def test_use_all(self):
        self.tf.use_all()

    def test_deafault_spiunp_time(self):
        del(self.tf.spinup_time)
        self.tf.spinup_time


if __name__ == '__main__':
    unittest.main()
