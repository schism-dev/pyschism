import os
from typing import Union

from pyproj import CRS

from pyschism.mesh import grd
from pyschism.mesh.base import EuclideanMesh2D


class Fgrid(EuclideanMesh2D):
    """
    Base class for all friction types (e.g. manning.grd, drag.grd, etc...)
    """
    def __init__(self, path: os.PathLike, crs: Union[str, CRS] = None):
        msh = grd.reader(path)
        msh.update({"crs": crs})
        super().__init__(**grd.euclidean_mesh(**msh))

    @property
    def nchi(self):
        msg = "Child class must implement nchi parameter"
        raise NotImplementedError(msg)
