from pyschism.mesh.friction.fgrid import Fgrid


class ManningsN(Fgrid):
    """  Class for representing Manning's n values.  """

    @classmethod
    def constant(cls, mesh, value):
        return cls(
            mesh.nodes,
            mesh.elements,
            crs=mesh.crs,
            description=mesh.description + "_mannings"
            )
