from pyschism.mesh import gr3
from pyschism.mesh.gmesh import Gmesh


class Fgrid(Gmesh):
    """
    Base class for all friction types (e.g. manning.gr3, drag.gr3, etc...)
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
        return Fgrid(**gr3.reader(path), crs=crs)