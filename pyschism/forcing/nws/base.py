from pyschism.forcing.base import ModelForcing
from abc import abstractmethod
# from datetime import datetime


class NWS(ModelForcing):

    def __str__(self):
        """Returns string used in param.nml"""
        return f"{self.dtype.value}"

    def __call__(self, model_driver):
        self._start_date = model_driver.param.opt.start_date
        self._rnday = model_driver.param.core.rnday
        model_driver.param.opt.nws = self

    @abstractmethod
    def write(self, path, overwrite=False):
        """Provides a method for writting SCHISM atmospherics files to disk.

        Since the output is different for each NWS type, the derived class
        must implement this method.
        """

    @property
    @abstractmethod
    def dtype(self):
        """Returns the NWSType of this object.

        Derived classes must return the corresponding instance of
        :enum:`pyschism.enums.NWSType`
        """
