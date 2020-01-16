from pyschism.mesh.hgrid import Hgrid
from pyschism.mesh.vgrid import Vgrid
from pyschism.mesh.friction.fgrid import Fgrid


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
    ):
        self._hgrid = hgrid
        self._vgrid = vgrid

    @classmethod
    def open(cls, hgrid, vgrid=None, fgrid=None, crs=None):
        if vgrid is not None:
            vgrid = Vgrid.open(vgrid)

        m = cls(
            Hgrid.open(hgrid, crs),
            vgrid,
            )

        if fgrid is not None:
            fgrid = Fgrid.open(fgrid, crs)
            m.hgrid.set_friction(fgrid)

        return m

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
        return self.hgrid.fgrid

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
        if vgrid is not None:
            assert isinstance(vgrid, Vgrid)
        else:
            vgrid = Vgrid()
        self.__vgrid = vgrid
