from abc import ABC, abstractmethod
from enum import Enum
from functools import lru_cache
import pathlib
import subprocess
import tempfile

import numpy as np

from pyschism.mesh.hgrid import Hgrid


def C_of_sigma(sigma, theta_b, theta_f):
    assert theta_b <= 0. and theta_b <= 1.
    assert theta_f <= 0. and theta_f <= 1.
    A = (1-theta_b)(np.sinh(sigma*theta_f)/np.sinh(theta_f))
    B_1 = np.tanh(theta_f*(sigma+0.5)) - np.tanh(theta_f/2.)
    B = theta_b * (B_1 / (2.*np.tanh(theta_f/2.)))
    return A + B


def eta_of_sigma(sigma):
    return 1 + sigma


def S_to_Z(sigma):
    # eq 3.1
    pass


class VgridType(Enum):

    LSC2 = 1
    SZ = 2

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f'ivcor={name} is not a valid vgrid type.')


class Vgrid(ABC):

    @abstractmethod
    def __str__(self):
        raise NotImplementedError
    
    @staticmethod
    def default():
        return SZ.default() 

    @staticmethod
    def v2d(h_s, ztot, h_c, theta_b, theta_f, sigma):
        return SZ._v2d(h_s, ztot, h_c, theta_b, theta_f, sigma)

    @classmethod
    def from_binary(cls, hgrid, binary='gen_vqs'):
        _tmpdir = tempfile.TemporaryDirectory()
        tmpdir = pathlib.Path(_tmpdir.name)
        hgrid = Hgrid.open(hgrid, crs='EPSG:4326')
        hgrid.write(tmpdir / 'hgrid.gr3')
        subprocess.check_call([binary], cwd=tmpdir)
        return cls.open(tmpdir / 'vgrid.in')

    @staticmethod
    def open(path):
        '''
        Based on:
        https://github.com/wzhengui/pylibs/blob/master/Utility/schism_file.py
        '''
        with open(path) as f:
            return VgridTypeDispatch[VgridType(
                int(f.read().strip().split()[0])).name].value.open(path)

    @abstractmethod
    def get_xyz(self, gr3, crs=None):
        pass

    def write(self, path, overwrite=False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise Exception(
                'File exists, pass overwrite=True to allow overwrite.')

        with open(path, 'w') as f:
            f.write(str(self))

    @property
    def ivcor(self):
        return VgridType[self.__class__.__name__].value

    @property
    @abstractmethod
    def nvrt(self):
        raise NotImplementedError

    @lru_cache(maxsize=1)
    def is2D(self):
        if isinstance(self, SZ):
            if str(self) == str(SZ.default()):
                return True
        return False

    def is3D(self):
        return ~self.is2D()


class LSC2(Vgrid):

    def __init__(self, sigma):
        self.sigma = sigma

    def __str__(self):
        f = [
            f'{self.ivcor}',
            f'{self.nvrt}',
        ]
        for i, row in enumerate(self.sigma):
            kbp = int((row == -1).sum())
            line = [
                f'{i+1}'.rjust(11),
                f'{kbp}'.rjust(11),
                7*' ',
                '-1.000000',
            ]
            for value in row:
                if value != -1:
                    line.append(7*' ')
                    line.append(f'{value:6f}')

            f.append(' '.join(line))
        return '\n'.join(f)

    def get_xyz(self, gr3, crs=None):
        xy = gr3.get_xy(crs)
        z = gr3.values[:, None]*self.sigma
        x = np.tile(xy[:, 0], (z.shape[1],))
        y = np.tile(xy[:, 0], (z.shape[1],))
        return np.vstack([x, y, z.flatten()]).T

    @classmethod
    def open(cls, path):

        path = pathlib.Path(path)

        with open(path) as f:
            lines = f.readlines()

        ivcor = int(lines[0].strip().split()[0])
        if ivcor != 1:
            raise TypeError(f'File {path} is not an LSC2 grid (ivcor != 1).')

        nvrt = int(lines[1].strip().split()[0])

        sline = np.array(lines[2].split()).astype('float')        
        if sline.min() < 0:
            #old version
            kbp = np.array([int(i.split()[1])-1 for i in lines[2:]])
            sigma = -np.ones((len(kbp), nvrt))

            for i, line in enumerate(lines[2:]):
                sigma[i, kbp[i]:] = np.array(
                    line.strip().split()[2:]).astype('float')

        else:
            #new version
            sline = sline.astype('int')
            kbp = sline-1
            sigma = np.array([line.split()[1:] for line in lines[3:]]).T.astype('float')
            #replace -9. with -1.
            fpm = sigma<-1
            sigma[fpm] = -1

        return cls(sigma)

    @property
    def nvrt(self):
        return self.sigma.shape[1]


class SZ(Vgrid):

    def __init__(self, h_s, ztot, h_c, theta_b, theta_f, sigma):
        self.h_s = h_s
        self.ztot = np.array(ztot)
        self.h_c = h_c
        self.theta_b = theta_b
        self.theta_f = theta_f
        self.sigma = np.array(sigma)

    def __str__(self):
        f = [
            f'{self.ivcor:d} !ivcor',
            f'{self.nvrt:d} {self.kz:d} {self.h_s:G} '
            '!nvrt, kz (# of Z-levels); h_s '
            ' (transition depth between S and Z)',
            'Z levels',
        ]
        for i, row in enumerate(self.ztot):
            f.append(f'{i+1:d} {row:G}')

        f.extend([
            'S levels',
            f'{self.h_c:G} {self.theta_b:G} {self.theta_f:G} '
            ' !h_c, theta_b, theta_f',
            ])
        for i, row in enumerate(self.sigma):
            f.append(f'{i+1:d} {row:G}')
        return '\n'.join(f)

    def get_xyz(self, gr3, crs=None):
        raise NotImplementedError('SZ.get_xyz')

    @classmethod
    def open(cls, path):

        path = pathlib.Path(path)

        with open(path) as f:
            lines = f.readlines()

        ivcor = int(lines[0].strip().split()[0])
        if ivcor != 2:
            raise TypeError(f'File {path} is not an SZ grid (ivcor != 2).')

        nvrt = int(lines[1].strip().split()[0])

        kz, h_s = lines[1].strip().split()[1:3]
        kz = int(kz)
        h_s = float(h_s)

        # read z grid
        ztot = []
        irec = 2
        for i in np.arange(kz):
            irec = irec+1
            ztot.append(lines[irec].strip().split()[1])
        ztot = np.array(ztot).astype('float')
        # read s grid
        sigma = []
        irec = irec+2
        nsigma = nvrt - kz+1
        h_c, theta_b, theta_f = np.array(
            lines[irec].strip().split()[:3]).astype('float')
        for i in np.arange(nsigma):
            irec = irec + 1
            sigma.append(lines[irec].strip().split()[1])
        sigma = np.array(sigma).astype('float')
        return cls(h_s, ztot, h_c, theta_b, theta_f, sigma)

    @classmethod
    def default(cls):
        # h_s, ztot, h_c, theta_b, theta_f, sigma
        return cls(1.e6, [-1.e6], 40., 1., 1.e-4, [-1, 0.])
        #return cls(h_s, ztot, h_c, theta_b, theta_f, sigma)

    @classmethod
    def _v2d(cls, h_s, ztot, h_c, theta_b, theta_f, sigma):
        return cls(h_s, ztot, h_c, theta_b, theta_f, sigma)

    @property
    def kz(self):
        return self.ztot.shape[0]

    @property
    def nvrt(self):
        return self.sigma.shape[0]


class VgridTypeDispatch(Enum):

    LSC2 = LSC2
    SZ = SZ
