import cf  # type: ignore[import]
import pathlib
from typing import Union

from pyschism.enums import NWSType
from pyschism.forcing.atmosphere.nws import NWS


class NWS2(NWS):

    def __init__(self, field_list: cf.FieldList):
        """Loads a field_list to use as NWS2 input.

        This class is used by :class:`pyschism.forcing.Atmosphere`

        Args:
            field_list: File system path of directory to load.
        """
        if not isinstance(field_list, cf.FieldList):
            raise TypeError(f'field_list must be of type {cf.FieldList}')
        # TODO: Sanity check on field_list object.
        self.field_list = field_list

    @property
    def dtype(self):
        """Returns the datatype of the object"""
        return NWSType(2)


class Sflux(NWS2):

    def __init__(self, level_1, air_1=None, prc_1=None, rad_1=None,
                 level_2=None, air_2=None, prc_2=None, rad_2=None):
        """Creates symlinks to the atmospheric wind files.

        It will not

        """
        self._level_1 = level_1
        self._air_1 = air_1
        self._prc_1 = prc_1
        self._rad_1 = rad_1
        self._level_2 = level_2
        self._air_2 = air_2
        self._prc_2 = prc_2
        self._rad_2 = rad_2
        # TODO: Run sanity check here
        super().__init__(cf.read())

    @staticmethod
    def load(path: Union[str, pathlib.Path]):
        pass


class SfluxServerFiles(NWS2):

    def __init__(self, level_1, air_1=None, prc_1=None, rad_1=None,
                 level_2=None, air_2=None, prc_2=None, rad_2=None):
        """Creates symlinks to the atmospheric wind files.

        It will not

        """
        # set level_1

        self._level_1 = level_1
        self._air_1 = air_1
        self._prc_1 = prc_1
        self._rad_1 = rad_1
        self._level_2 = level_2
        self._air_2 = air_2
        self._prc_2 = prc_2
        self._rad_2 = rad_2
        # Send an empty cf.FileList()
        super().__init__(cf.FileList())
