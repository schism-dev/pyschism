from abc import abstractmethod
from typing import List

from pyschism.forcing.bctides.bctypes import Bctype


class Isatype(Bctype):

    @property
    @abstractmethod
    def isatype(self) -> int:
        pass


class UniformTimeHistorySalinity(Isatype):

    def __init__(self, time_history):
        raise NotImplementedError(f'{self.__class__.__name__}')
        self.time_history = time_history

    @property
    def isatype(self) -> int:
        return 1


class ConstantSalinity(Isatype):

    def __init__(self, value):
        raise NotImplementedError(f'{self.__class__.__name__}')
        self.value = value

    @property
    def isatype(self) -> int:
        return 2


class SalinityInitialConditions(Isatype):

    def __init__(self):
        raise NotImplementedError(f'{self.__class__.__name__}')

    @property
    def isatype(self):
        return 3


class SpatiallyVaryingTimeHistorySalinity(Isatype):

    def __init__(self, data_source, tobc: List[float] = None):
        self.data_source = data_source
        self.tobc = tobc

    def get_boundary_string(self, hgrid, boundary):
        if self.tobc is None:
            return '1.'

    @property
    def isatype(self):
        return 4


Isatype1 = UniformTimeHistorySalinity
Isatype2 = ConstantSalinity
Isatype3 = SalinityInitialConditions
Isatype4 = SpatiallyVaryingTimeHistorySalinity
