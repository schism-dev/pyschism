from abc import abstractmethod
import logging

from numba import prange, jit
import numpy as np

from pyschism.mesh.gridgr3 import Gr3Field


logger = logging.getLogger(__name__)


class Nudge(Gr3Field):

    def __init__(self, bctides, data_source, rlmax=1.5, rnu_day=0.25):

        @jit(nopython=True, parallel=True)
        def compute_nudge(lon, lat, opbd, out):

            nnode = lon.shape[0]

            rnu_max = 1./rnu_day/86400.

            for idn in prange(nnode):
                if idn in opbd:
                    rnu = rnu_max
                    distmin = 0.
                else:
                    distmin = np.finfo(np.float64).max
                    for j in opbd:
                        rl2 = np.sqrt(
                            np.square(lon[idn] - lon[j-1])
                            + np.square(lat[idn] - lat[j-1])
                            )
                        if rl2 < distmin:
                            distmin = rl2
                rnu = 0.
                if distmin <= rlmax:
                    rnu = (1-distmin/rlmax)*rnu_max
                out[idn] = rnu
        out = np.zeros(bctides.hgrid.values.shape)
        xy = bctides.hgrid.get_xy(crs='epsg:4326')
        opbd = []
        for boundary in bctides.gdf.itertuples():
            forcing = getattr(boundary, self.bctype)
            if forcing.nudge is True:
                opbd.extend(list(boundary.indexes))
        opbd = np.array(opbd)
        if len(opbd) > 0:
            logger.info(f'Begin compute_nudge for {self.name}.')
            from time import time
            start = time()
            compute_nudge(xy[:, 0], xy[:, 1], opbd, out)
            logger.info(f'compute_nudge took {time()-start} seconds.')
        super().__init__(
            nodes={i: (coord, out[i]) for i, coord in enumerate(bctides.hgrid.coords)},
            elements=bctides.hgrid.elements.elements,
            description=f"{rlmax}, {rnu_day}",
            crs=bctides.hgrid.crs,
        )

    @property
    @abstractmethod
    def name(self):
        ''''''

    @property
    @abstractmethod
    def bctype(self):
        ''''''


class TEM_Nudge(Nudge):

    @property
    def name(self):
        return 'temperature'

    @property
    def bctype(self):
        return 'itetype'


class SAL_Nudge(Nudge):

    @property
    def name(self):
        return 'salinity'

    @property
    def bctype(self):
        return 'isatype'


TempNudge = TEM_Nudge
SaltNudge = SAL_Nudge
