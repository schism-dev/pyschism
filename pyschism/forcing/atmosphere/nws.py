from abc import ABC, abstractmethod


class NWS(ABC):

    def __str__(self):
        """Returns string used in param.nml"""
        return f"{self.dtype.value}"

    @abstractmethod
    def write(self):
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
