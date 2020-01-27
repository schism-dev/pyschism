#! /usr/bin/env python
import unittest
from unittest.mock import patch
import sys
import pathlib
from pyschism.cmd import tidal_run


class TidalRunCmdTestCase(unittest.TestCase):

    def _test_tidal_run_mesh(self):
        parent = pathlib.Path(__name__).parent
        hgrid = parent.resolve() / "examples/example_1/hgrid.gr3"
        cmd = ["tidal_run"]
        cmd += [str(hgrid.resolve())]
        cmd += ["2017-9-18T12:00"]
        cmd += ["2017-9-23T12:00"]
        cmd += ["--spinup-days=2"]
        print(" ".join(cmd))
        with patch.object(sys, 'argv', cmd):
            self.assertEqual(tidal_run.main(), 0)

    def test_init(self):
        with patch.object(tidal_run, "main", return_value=0):
            with patch.object(tidal_run, "__name__", "__main__"):
                with patch.object(tidal_run.sys, 'exit') as mock_exit:
                    tidal_run.init()
                    assert mock_exit.call_args[0][0] == 0


if __name__ == '__main__':
    unittest.main()
