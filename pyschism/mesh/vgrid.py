import os
import pathlib
import subprocess
import tempfile
import shutil
from typing import Union

import numpy as np

from pyschism.mesh.hgrid import Hgrid


class Vgrid:

    def __init__(self):
        """Represents a SCHISM vertical grid."""
        self._vgrid = self._get_2D_string()
        pass

    @classmethod
    def from_binary(cls, outdir: Union[str, os.PathLike], hgrid, *args, **kwargs):
        _tmpdir = tempfile.TemporaryDirectory()
        tmpdir = pathlib.Path(_tmpdir.name)
        hgrid = Hgrid.open(hgrid, crs='EPSG:4326')
        hgrid.write(tmpdir / 'hgrid.gr3')
        subprocess.check_call(['gen_vqs'], cwd=tmpdir)
        outdir = pathlib.Path(outdir)
        print(f'write vgrid to dir {outdir}')
        if not outdir.is_file():
            shutil.copy2(tmpdir / 'vgrid.in', outdir)
        else:
            shutil.copy2(tmpdir / 'vgrid.in', outdir / 'vgrid.in')
        
        obj = cls()
        #with open(path) as f:
        #    obj._vgrid = f.read()
        return obj

    @classmethod
    def open(cls, path):
        path = pathlib.Path(path)
        if path.name != 'vgrid.in':
            raise TypeError('Not a valid vgrid.in file')
        obj = cls()
        with open(path) as f:
            obj._vgrid = f.read()
        return obj

    '''
    read_vgrid is based on:
    https://github.com/wzhengui/pylibs/blob/master/Utility/schism_file.py
    '''
    @staticmethod
    def read_vgrid(fname):
        fid=open(fname, 'r')
        lines=fid.readlines()
        fid.close()

        ivcor=int(lines[0].strip().split()[0])
        nvrt=int(lines[1].strip().split()[0])

        if ivcor == 1:
            lines=lines[2:]
            kbp=np.array([int(i.split()[1])-1 for i in lines])
            NP=len(kbp)
            #print(NP)
            sigma=-np.ones([NP,nvrt])
            for i, line in enumerate(lines):
                sigma[i,kbp[i]:]=np.array(line.strip().split()[2:]).astype('float')
        elif ivcor == 2:
            kz, h_s=lines[1].strip().split()[1:3]
            kz=int(kz)
            h_s=float(kz)
            
            #read z grid
            ztot=[]
            irec=2
            for i in np.arange(kz):
                irec=irec+1
                ztot.append(lines[irec].strip().split()[1])
            ztot=np.array(ztot).astype('float')

            #read s grid
            sigma=[]
            irec=irec+2
            nsigma=nvrt-kz+1
            h_c, theta_b, theta_f=np.array(lines[irec].strip().split()[:3]).astype('float')
            for i in np.arange(nsigma):
                irec=irec+1
                sigma.append(lines[irec].strip().split()[1])
            sigma=np.array(sigma).astype('float')

        return sigma

    def __str__(self):
        return self._vgrid

    def write(self, path, overwrite=False):
        if path.is_file() and not overwrite:
            raise Exception(
                'File exists, pass overwrite=True to allow overwrite.')

        with open(path, 'w') as f:
            f.write(str(self))

    def is2D(self):
        return self._vgrid == self._get_2D_string()

    def is3D(self):
        return not self.is2D()

    def _get_2D_string(self):
        return """2 !ivcor
2 1 1.e6 !nvrt, kz (# of Z-levels); h_s (transition depth between S and Z)
Z levels
1  -1.e6
S levels
40. 1. 1.e-4  !h_c, theta_b, theta_f
   1    -1.
   2    0."""
