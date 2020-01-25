from pyschism.mesh import grd
from pyschism.mesh.gmesh import Gmesh


class Fgrid(Gmesh):
    """
    Base class for all friction types (e.g. manning.grd, drag.grd, etc...)
    """
    def __init__(
        self,
        nodes,
        elements,
        boundaries=None,
        crs=None,
        description=None,
    ):
        super().__init__(**grd.to_gmesh({
            'nodes': nodes,
            'elements': elements,
            'description': description,
            'boundaries': boundaries}),
        crs=crs)

    @staticmethod
    def open(path, crs=None):
        return Fgrid(**grd.reader(path), crs=crs)