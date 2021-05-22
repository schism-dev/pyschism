from abc import ABC, abstractmethod


class Bctype(ABC):

    @abstractmethod
    def get_boundary_string(self, boundary) -> str:
        pass
