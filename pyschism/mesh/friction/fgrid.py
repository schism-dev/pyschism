from pyschism.mesh import grd
from pyschism.mesh.base import EuclideanMesh2D


class Fgrid(EuclideanMesh2D):
    """
    Base class for all friction types (e.g. manning.grd, drag.grd, etc...)
    """
    def __init__(
        self,
        nodes,
        elements,
        crs=None,
        description=None,
    ):
        msh = grd.euclidean_mesh({
            'nodes': nodes,
            'elements': elements,
            'description': description,
            'crs': crs})
        msh.update({"crs": crs})
        super().__init__(**msh)

    @staticmethod
    def open(path, crs=None):
        msh = grd.reader(path)
        msh.update({"crs": crs})
        return Fgrid(**msh)

    @property
    def nchi(self):
        msg = "Child class must implement nchi parameter"
        raise NotImplementedError(msg)
