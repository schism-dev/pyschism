from abc import ABC, abstractmethod
# from datetime import datetime, timedelta
# from typing import Union, Dict

# from netCDF4 import Dataset

# from pyschism import dates
# from pyschism.mesh.base import Gr3


class BaroclinicComponent(ABC):

    @property
    @abstractmethod
    def output_interval(self):
        """Dataset output frequency."""


class BaroclinicForcing:

    def __init__(self, resource):
        self.resource

    @property
    def resource(self):
        return self._resource

    @resource.setter
    def resource(self, resource):
        self._resource = resource
