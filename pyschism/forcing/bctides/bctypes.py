from abc import abstractmethod

from pyschism.forcing.base import ModelForcing


class BoundaryForcing(ModelForcing):

    @abstractmethod
    def get_boundary_string(self, boundary) -> str:
        pass


Bctype = BoundaryForcing
