from pyschism.mesh.gmesh import Gmesh


class Fgrid(Gmesh):
    """
    Base class for all friction types (e.g. manning.gr3, drag.gr3, etc...)
    """

    def __init__(
        self,
        nodes,
        elements,
        values=None,
        crs=None,
        description=None,
    ):
        super().__init__(*self._gr3_to_mesh(nodes, elements), crs, description)
