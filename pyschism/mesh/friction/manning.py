from pyschism.mesh.friction.fgrid import Fgrid


class ManningsN(Fgrid):
    """  Class for representing Manning's n values.  """

    @classmethod
    def constant(cls, mesh, value):
        nodes = dict()
        for id, (x, y, _) in mesh.nodes.items():
            nodes[id] = (x, y, value)
        description = mesh.description + "_mannings"
        return cls(nodes, mesh.elements, crs=mesh.crs, description=description)
