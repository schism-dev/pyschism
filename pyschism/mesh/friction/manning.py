from pyschism.mesh.friction.fgrid import Fgrid


class ManningsN(Fgrid):
    """  Class for representing Manning's n values.  """

    @classmethod
    def constant(cls, mesh, value):
        return cls(
            {id: (coords, float(value))
             for (id, coords) in mesh._coords.items()},
            mesh._elements,
            crs=mesh.crs,
            description=mesh.description + " mannings"
            )

    @property
    def nchi(self):
        return -1

    @property
    def hmin_man(self):
        try:
            return self.__hmin_man
        except AttributeError:
            return 1.

    @hmin_man.setter
    def hmin_man(self, hmin_man):
        self.__hmin_man = float(hmin_man)

    @property
    def fname(self):
        return 'manning.gr3'
