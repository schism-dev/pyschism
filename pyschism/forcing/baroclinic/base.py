from abc import ABC, abstractmethod
from datetime import datetime, timedelta
from typing import Union, Dict

from netCDF4 import Dataset

from pyschism import dates
from pyschism.mesh.base import Gr3


class BaroclinicComponent(ABC):

    @abstractmethod
    def interpolate(self, gr3: Gr3, date: datetime, level: int = 0, step: int = 0):
        """Used to generate *.ic files"""

    @property
    @abstractmethod
    def output_interval(self):
        """Dataset output frequency."""


class Temperature(BaroclinicComponent):
    pass


class Salinity(BaroclinicComponent):
    pass


class BaroclinicForcing:

    def __init__(self, resource):
        self.resource

    @property
    def resource(self):
        return self._resource

    @resource.setter
    def resource(self, resource):
        self._resource = resource

    @property
    def temperature(self) -> Temperature:
        return self._temperature

    @property
    def salinity(self) -> Salinity:
        return self._salinity
