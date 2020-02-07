import numpy as np
from functools import lru_cache
from collections.abc import Iterable
from pyschism.forcing.bctypes import BoundaryCondition
from pyschism.mesh.hgrid import Hgrid
from pyschism.mesh.vgrid import Vgrid
from pyschism.mesh.friction import (
    Fgrid,
    ManningsN,
    DragCoefficient,
    RoughnessLength
)


class Mesh:
    """
    Class representing a SCHISM computational domain.
    This class combines the horizontal grid (hgrid), vertical grid (vgrid) and
    friction/drag grids (fgrid).
    """

    def __init__(
        self,
        hgrid,
        vgrid=None,
        fgrid=None,
    ):
        self._hgrid = hgrid
        self._vgrid = vgrid
        self._fgrid = fgrid

    @classmethod
    def open(cls, hgrid, vgrid=None, fgrid=None, crs=None):
        if vgrid is not None:
            vgrid = Vgrid.open(vgrid)

        if fgrid is not None:
            fgrid = Fgrid.open(fgrid, crs)

        return cls(
            Hgrid.open(hgrid, crs),
            vgrid,
            fgrid
            )

    def set_friction(self, ftype, value):

        # certify ftype
        ftypes = {
            'manning': ManningsN,
            'drag': DragCoefficient,
            'rough': RoughnessLength
        }
        msg = f"ftype argument must be one of {ftypes.keys()}"
        assert ftype.lower() in ftypes, msg

        # certify value
        msg = "value argument must be an instance of type "
        msg += f"{int}, {float} or an iterable ."
        assert isinstance(value, (Iterable, int, float, Fgrid)), msg

        if isinstance(value, (int, float)):
            if ftype == 'manning':
                self._fgrid = ftypes[ftype].constant(self.hgrid, value)

        return self.fgrid

    def set_boundary_forcing(self, forcing, id=None):
        if id is None:
            for i in range(len(self.open_boundaries)):
                self.set_boundary_forcing(forcing, i)
        else:
            assert isinstance(forcing, BoundaryCondition)
            self._open_boundaries[id][forcing.vartype] = forcing

    def make_plot(self, **kwargs):
        if self.vgrid.is3D():
            msg = "Plotting not yet supported for 3D meshes."
            raise NotImplementedError(msg)
        elif self.vgrid.is2D():
            self.hgrid.make_plot(**kwargs)

    @property
    def hgrid(self):
        return self._hgrid

    @property
    def vgrid(self):
        return self._vgrid

    @property
    def fgrid(self):
        return self._fgrid

    @property
    def crs(self):
        return self.hgrid.crs

    @property
    def ics(self):
        if self.hgrid.crs is None:
            msg = "Can't determine ics parameter. No projection information "
            msg += "has been provided for the hgrid."
            raise Exception(msg)
        if self.crs.is_geographic:
            return 2
        else:
            return 1

    @property
    def slam0(self):
        return np.median(self.hgrid.get_x("EPSG:4326"))

    @property
    def sfea0(self):
        return np.median(self.hgrid.get_y("EPSG:4326"))

    @property
    def open_boundaries(self):
        open_boundaries = self._open_boundaries.copy()
        for id in self._open_boundaries:
            indexes = list(map(
                    self.hgrid.get_node_index,
                    self.hgrid.boundaries[None][id]['indexes']))
            open_boundaries[id].update({'indexes': indexes})
        return open_boundaries

    @property
    @lru_cache
    def _open_boundaries(self):
        _open_boundaries = self.hgrid.boundaries[None].copy()
        for id in _open_boundaries:
            _open_boundaries[id]['eta'] = {}
            _open_boundaries[id]['uv'] = {}
            _open_boundaries[id]['stt'] = {}
        return _open_boundaries

    @property
    def _hgrid(self):
        return self.__hgrid

    @property
    def _vgrid(self):
        return self.__vgrid

    @property
    def _fgrid(self):
        return self.__fgrid

    @_hgrid.setter
    def _hgrid(self, hgrid):
        assert isinstance(hgrid, Hgrid)
        self.__hgrid = hgrid

    @_vgrid.setter
    def _vgrid(self, vgrid):
        if vgrid is not None:
            assert isinstance(vgrid, Vgrid)
        else:
            vgrid = Vgrid()
        self.__vgrid = vgrid

    @_fgrid.setter
    def _fgrid(self, fgrid):
        if fgrid is None:
            fgrid = ManningsN.constant(self.hgrid, 0.)
        assert isinstance(fgrid, Fgrid)
        self.__fgrid = fgrid
