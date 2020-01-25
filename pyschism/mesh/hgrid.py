from collections.abc import Iterable
from pyschism.mesh import gr3
from pyschism.mesh import Gmesh
from pyschism.mesh.friction import (
    Fgrid,
    ManningsN,
    DragCoefficient,
    RoughnessLength
)


class Hgrid(Gmesh):
    """
    Class that represents the unstructured planar mesh used by SCHISM.
    """

    def __init__(
        self,
        nodes,
        elements,
        boundaries=None,
        crs=None,
        description=None,
    ):
        grd = {
            'nodes': nodes,
            'elements': elements,
            'description': description,
            'boundaries': boundaries,
        }
        super().__init__(**gr3.to_gmesh(grd), crs=crs)

    @staticmethod
    def open(path, crs=None):
        return Hgrid(**gr3.reader(path), crs=crs)

    def set_friction(self, value, ftype='manning'):

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

        if isinstance(value, Fgrid):
            self._fgrid = value

        elif isinstance(value, (int, float)):
            if ftype == 'manning':
                self._fgrid = ftypes[ftype].constant(self, value)

        return self.fgrid

    @property
    def fgrid(self):
        try:
            return self.__fgrid
        except AttributeError:
            self._fgrid = ManningsN.constant(self, 0.025)
            return self.__fgrid

    @property
    def _fgrid(self):
        return self.__fgrid

    @_fgrid.setter
    def _fgrid(self, fgrid):
        assert isinstance(fgrid, Fgrid)
        self.__fgrid = fgrid
