#! /usr/bin/env python
import unittest
from unittest.mock import patch
import sys
import tempfile
import warnings
import pathlib
from pyschism.cmd import plot_mesh


class PlotMeshCmdTestCase(unittest.TestCase):

    @patch('matplotlib.pyplot.show')
    def test_plot_mesh_command(self, mock):
        tmpdir = tempfile.TemporaryDirectory()
        outdir = pathlib.Path(tmpdir.name)
        nodes = {
            0: (0., 0., -1.),
            1: (1., 0., 0.),
            2: (1., 1., 0.),
            3: (0., 1., 0.),
            4: (0.5, 1.5, 0.),
        }
        elements = {
            0: [2, 4, 3],
            1: [0, 1, 2, 3],
        }
        with open(outdir / 'hgrid.gr3', 'w') as f:
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
        cmd = ["plot_mesh"]
        cmd += [str(outdir / 'hgrid.gr3')]
        cmd += ["--plot-elements"]
        cmd += ["--plot-boundaries"]
        cmd += [f"--save-path={str(outdir / 'hgrid.png')}"]
        with patch.object(sys, 'argv', cmd):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.assertEqual(plot_mesh.main(), 0)

    def test_init(self):
        with patch.object(plot_mesh, "main", return_value=0):
            with patch.object(plot_mesh, "__name__", "__main__"):
                with patch.object(plot_mesh.sys, 'exit') as mock_exit:
                    plot_mesh.init()
                    assert mock_exit.call_args[0][0] == 0


if __name__ == '__main__':
    unittest.main()
