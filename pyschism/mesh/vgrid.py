from abc import ABC, abstractmethod
from enum import Enum
from functools import lru_cache
import pathlib
import subprocess
import tempfile

import numpy as np

from pyschism.mesh.hgrid import Hgrid

from matplotlib.pyplot import *


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

    def __init__(self, hsm, nv, h_c, theta_b, theta_f):
        self.hsm = np.array(hsm)
        self.nv = np.array(nv)
        self.h_c = h_c
        self.theta_b = theta_b
        self.theta_f = theta_f
        self.m_grid = None
        self._znd = None
        self._snd = None
        self._nlayer = None

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
        if type(gr3) == Hgrid:
            gr3 = gr3
        else:
            gr3=Hgrid.open(gr3)
        xy = gr3.get_xy(crs)
        z = gr3.values[:, None]*self.sigma
        x = np.tile(xy[:, 0], (z.shape[1],))
        y = np.tile(xy[:, 0], (z.shape[1],))
        return np.vstack([x, y, z.flatten()]).T

    def calc_m_grid(self):
        '''
        create master grid
        Adapted from:
        https://github.com/wzhengui/pylibs/blob/master/pyScripts/gen_vqs.py
        '''
        if self.m_grid:
            pass
        else:
            z_mas=np.ones([self.nhm,self.nvrt])*np.nan; eta=0.0
            for m, [hsmi,nvi] in enumerate(zip(self.hsm,self.nv)):
                #strethcing funciton
                hc=min(hsmi,self.h_c)
                for k in np.arange(nvi):
                    sigma= k/(1-nvi)  #zi=-sigma #original sigma coordiante
                    #compute zcoordinate
                    cs=(1-self.theta_b)*np.sinh(self.theta_f*sigma)/np.sinh(self.theta_f)+\
                        self.theta_b*(np.tanh(self.theta_f*(sigma+0.5))-\
                                np.tanh(self.theta_f*0.5))/2/np.tanh(self.theta_f*0.5)
                    z_mas[m,k]=eta*(1+sigma)+hc*sigma+(hsmi-hc)*cs

                #normalize z_mas
                z_mas[m]=-(z_mas[m]-z_mas[m,0])*hsmi/(z_mas[m,nvi-1]-z_mas[m,0])
            s_mas=np.array([z_mas[i]/self.hsm[i] for i in np.arange(self.nhm)])

            self.m_grid = z_mas

    def make_m_plot(self):
        '''
        plot master grid
        Adapted from:
        https://github.com/wzhengui/pylibs/blob/master/pyScripts/gen_vqs.py
        '''        
        #check master grid
        for i in np.arange(self.nhm-1):
            if min(self.m_grid[i,:self.nv[i]]-self.m_grid[i+1,:self.nv[i]])<0: \
                print('check: master grid layer={}, hsm={}, nv={}'.\
                      format(i+1,self.hsm[i+1],self.nv[i+1]))

        #plot master grid
        figure(figsize=[10,5])
        for i in np.arange(self.nhm): plot(i*np.ones(self.nvrt),\
                                           self.m_grid[i],'k-',lw=0.3)
        for k in np.arange(self.nvrt): plot(np.arange(self.nhm),\
                                            self.m_grid.T[k],'k-',lw=0.3)
        setp(gca(),xlim=[-0.5,self.nhm-0.5],ylim=[-self.hsm[-1],0.5])
        gcf().tight_layout()

    def calc_lsc2_att(self, gr3, crs=None):
        '''
        master grid to lsc2:
        compute vertical layers at nodes
        gr3 is either file or Hgrid Object
        Adapted from:
        https://github.com/wzhengui/pylibs/blob/master/pyScripts/gen_vqs.py
        '''
        if type(gr3) == Hgrid:
            gd = gr3
        else:
            gd=Hgrid.open(gr3)
        dp = gd.values*-1
        fpz=dp<self.hsm[0]
        dp[fpz]=self.hsm[0]
        
        #find hsm index for all points
        rat=np.ones(len(gd.nodes.id))*np.nan
        nlayer=np.zeros(len(gd.nodes.id)).astype('int')
        ind1=np.zeros(len(gd.nodes.id)).astype('int')
        ind2=np.zeros(len(gd.nodes.id)).astype('int')
        for m, hsmi in enumerate(self.hsm):
            if m==0:
                fp=dp<=self.hsm[m]
                ind1[fp]=0; ind2[fp]=0
                rat[fp]=0; nlayer[fp]=self.nv[0]
            else:
                fp=(dp>self.hsm[m-1])*(dp<=self.hsm[m])
                ind1[fp]=m-1
                ind2[fp]=m
                rat[fp]=(dp[fp]-self.hsm[m-1])/(self.hsm[m]-self.hsm[m-1])
                nlayer[fp]=self.nv[m]

        #Find the last non NaN node and fills the NaN values with it
        last_non_nan = (~np.isnan(self.m_grid)).cumsum(1).argmax(1)
        z_mas=np.array([np.nan_to_num(z_mas_arr,nan=z_mas_arr[last_non_nan[i]])\
                        for i, z_mas_arr in enumerate(self.m_grid)])
        znd=z_mas[ind1]*(1-rat[:,None])+z_mas[ind2]*rat[:,None]; #z coordinate
        for i in np.arange(len(gd.nodes.id)):
            znd[i,nlayer[i]-1]=-dp[i]
            znd[i,nlayer[i]:]=np.nan
        snd=znd/dp[:,None]; #sigma coordinate

        #check vgrid
        for i in np.arange(len(gd.nodes.id)):
            for k in np.arange(self.nvrt-1):
                if znd[i,k]<=znd[i,k+1]:
                    raise TypeError(f'wrong vertical layers')

        self._znd = znd
        self._snd = snd
        self._nlayer = nlayer


    def write(self, path, overwrite=False):
        '''
        write mg2lsc2 into vgrid.in
        '''
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise Exception(
                'File exists, pass overwrite=True to allow overwrite.')

        with open(path, 'w') as fid:
            fid.write('           1 !average # of layers={:0.2f}\n          {} !nvrt\n'.format\
                    (np.mean(self._nlayer),self.nvrt))
            bli=[]#bottom level index
            for i in np.arange(len(self._nlayer)):
                nlayeri=self._nlayer[i]; si=np.flipud(self._snd[i,:nlayeri])
                bli.append(self.nvrt-nlayeri+1)
                fstr=f"         {self.nvrt-nlayeri+1:2}"
                fid.write(fstr)
            for i in range(self.nvrt):
                fid.write(f'\n         {i+1}')
                for n,bl in enumerate(bli):
                    si=np.flipud(self._snd[n])
                    if bl <= i+1:
                        fid.write(f"      {si[i]:.6f}")
                    else:
                        fid.write(f"      {-9.:.6f}")
            fid.close()


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
        return self.nv[-1]

    @property
    def nhm(self):
        return self.hsm.shape[0]


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
