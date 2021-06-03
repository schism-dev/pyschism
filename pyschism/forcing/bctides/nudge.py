import logging

from numba import prange, jit
import numpy as np

from pyschism.mesh.gridgr3 import Gr3Field


logger = logging.getLogger(__name__)


class Nudge(Gr3Field):

    def __init__(self, hgrid, rlmax=1.5, rnu_day=0.25):

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

        opbd = []
        for row in hgrid.boundaries.open.itertuples():
            opbd.extend(row.indexes)
        opbd = np.array(opbd)
        out = np.zeros(hgrid.values.shape)
        xy = hgrid.get_xy(crs='epsg:4326')
        logger.info('Begin compute_nudge')
        from time import time
        start = time()
        compute_nudge(xy[:, 0], xy[:, 1], opbd, out)
        logger.info(f'compute_nudge took {time()-start} seconds.')
        super().__init__(
            nodes={i: (coord, out[i]) for i, coord in enumerate(hgrid.coords)},
            elements=hgrid.elements.elements,
            description=f"{rlmax}, {rnu_day}",
            crs=hgrid.crs,
        )

    def __call__(self, data_source, vgrid):
        if data_source.name == 'temperature':
            return TempNudge(self, data_source, vgrid)

    @classmethod
    def default(cls, hgrid):
        return cls(hgrid)


class TempNudge(Nudge):

    def __init__(self, nudge, data_source, vgrid):
        super().__init__(
            nudge.nodes.to_dict(),
            nudge.elements.elements,
            nudge.description,
            nudge.crs
        )
        self.vgrid = vgrid
