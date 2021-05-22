from abc import abstractmethod
from typing import List

from pyschism.forcing.bctides.bctypes import Bctype


class Itetype(Bctype):

    @property
    @abstractmethod
    def itetype(self) -> int:
        pass


class UniformTimeHistoryTemperature(Itetype):

    # def __init__(self, time_history: Dict[str, List[float]], tobc: float = 1.):
    #     self.time_history = time_history
    #     super().__init__(tobc)

    @property
    def itetype(self) -> int:
        return 1


class ConstantTemperature(Itetype):

    # def __init__(self, value, tobc: float = 1.):
    #     self.value = value
    #     super().__init__(tobc)

    @property
    def itetype(self) -> int:
        return 2


class TemperatureInitialConditions(Itetype):

    # def __init__(self, tobc: float = 1., tth):
    #     raise NotImplementedError(f'{self.__class__.__name__}')

    @property
    def itetype(self):
        return 3


class SpatiallyVaryingTimeHistoryTemperature(Itetype):

    def __init__(self, data_source, tobc: List[float] = None):
        self.data_source = data_source
        self.tobc = tobc

    def get_boundary_string(self, hgrid, boundary):
        if self.tobc is None:
            return '1.'

    @property
    def itetype(self):
        return 4


Itetype0 = Itetype
Itetype1 = UniformTimeHistoryTemperature
Itetype2 = ConstantTemperature
Itetype3 = TemperatureInitialConditions
Itetype4 = SpatiallyVaryingTimeHistoryTemperature
