import os
import pathlib
import subprocess
import tempfile
import shutil
from typing import Union


class Vgrid:

    def __init__(self):
        """Represents a SCHISM vertical grid."""
        self._vgrid = self._get_2D_string()

    @classmethod
    def from_binary(cls, outdir: Union[str, os.PathLike], hgrid, *args, **kwargs):
        _tmpdir = tempfile.TemporaryDirectory()
        tmpdir = pathlib.Path(_tmpdir.name)
        hgrid.write(tmpdir / 'hgrid.gr3')
        subprocess.check_call(['gen_vqs'], cwd=tmpdir)
        outdir = pathlib.Path(outdir)
        print(f'write vgrid to dir {outdir}')
        #outdir = pathlib.Path('./')
        shutil.copy2(tmpdir / 'vgrid.in', outdir / 'vgrid.in')
        
        obj = cls()
        #obj._vgrid = open(tmpdir / 'vgrid.in').read()
        return obj

    @staticmethod
    def open(path):
        raise NotImplementedError('Vgrid.open()')

    def __str__(self):
        return self._vgrid

    def write(self, path, overwrite=False):
        if path.is_file() and not overwrite:
            raise Exception(
                'File exists, pass overwrite=True to allow overwrite.')

        with open(path, 'w') as f:
            f.write(str(self))

    def is_2D(self):
        return self._vgrid == self._get_2D_string()

    def is_3D(self):
        return not self.is_2D()

    def _get_2D_string(self):
        return """2 !ivcor
2 1 1.e6 !nvrt, kz (# of Z-levels); h_s (transition depth between S and Z)
Z levels
1  -1.e6
S levels
40. 1. 1.e-4  !h_c, theta_b, theta_f
   1    -1.
   2    0."""
