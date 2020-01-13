from pyschism.mesh.friction.fgrid import Fgrid


class ManningsN(Fgrid):
    """  Class for representing Manning's n values.  """

    @classmethod
    def constant(cls, mesh, value):
        return cls(
            mesh._coords,
            mesh._triangles,
            mesh._quads,
            values=len(mesh.coords)*[value],
            crs=mesh.crs,
            description=mesh.description + "_mannings"
            )
