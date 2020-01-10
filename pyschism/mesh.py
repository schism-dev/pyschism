from pyschism.hgrid import Hgrid
from pyschism.vgrid import Vgrid

__all__ = [
    "Hgrid",
    "Vgrid"
]


class SchismMesh(Hgrid):
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
        return cls(Hgrid.open(hgrid))

    @property
    def hgrid(self):
        return self._hgrid

    @property
    def vgrid(self):
        return self._vgrid

    @property
    def _hgrid(self):
        return self.__hgrid

    @property
    def _vgrid(self):
        return self.__vgrid

    @_hgrid.setter
    def _hgrid(self, hgrid):
        assert isinstance(hgrid, Hgrid)
        self.__hgrid = hgrid

    @_vgrid.setter
    def _vgrid(self, vgrid):
        assert isinstance(vgrid, vgrid)
        self.__vgrid = vgrid
