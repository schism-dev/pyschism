import os
import pathlib
from typing import Union

import f90nml  # type: ignore[import]

from pyschism.enums import NWSType
from pyschism.forcing.atmosphere.nws.nws import NWS
from pyschism.forcing.atmosphere.nws.nws2.sflux import SfluxDataset
from pyschism.mesh import Gr3, gridgr3


SFLUX_DEFAULTS = f90nml.read(
    pathlib.Path(__file__).parent / 'sflux_inputs.txt')


class NWS2(NWS):

    def __init__(
            self,
            sflux_1: SfluxDataset,
            sflux_2: SfluxDataset = None,
            windrot: gridgr3.Windrot = None
    ):
        """Loads SfluxDataset to use as NWS2 input. """

        self.sflux_1 = sflux_1
        self.sflux_2 = sflux_2
        self.windrot = windrot

    def fetch_data(self, start_date=None, rnday=None, bbox=None, **kwargs):
        if hasattr(self.sflux_1, 'fetch_data'):
            self.sflux_1.fetch_data(start_date, rnday, bbox=bbox, **kwargs)
        if self.sflux_2 is not None:
            if hasattr(self.sflux_2, 'fetch_data'):
                self.sflux_2.fetch_data(
                    start_date,
                    rnday,
                    bbox=bbox,
                    **kwargs
                    )

    def __str__(self):
        data = []
        data = '\n'.join(data)
        return f'&sflux_inputs\n{data}/\n'

    def write(self, path: Union[str, os.PathLike], overwrite: bool = False,
              windrot: bool = True):
        # write sflux namelist
        path = pathlib.Path(path)
        if path.name != 'sflux':
            path /= 'sflux'
        path.mkdir(exist_ok=True)
        with open(path / 'sflux_inputs.txt', 'w') as f:
            f.write(str(self))
        # write sflux data
        self.sflux_1.write(
            path, 1, overwrite,
            # start_date=self.start_date,
            # rnday=self._rnday
        )
        if self.sflux_2 is not None:
            self.sflux_2.write(
                path, 2, overwrite,
                # start_date=self.start_date,
                # rnday=self._rnday
                )
        # # write windrot data
        if windrot is not False and self.windrot is not None:
            windrot = 'windrot_geo2proj.gr3' if windrot is True else windrot
            self.windrot.write(path.parent / 'windrot_geo2proj.gr3', overwrite)

    @property
    def dtype(self) -> NWSType:
        """Returns the datatype of the object"""
        return NWSType(2)

    @property
    def sflux_1(self) -> SfluxDataset:
        return self._sflux_1

    @sflux_1.setter
    def sflux_1(self, sflux_1):
        if not isinstance(sflux_1, SfluxDataset):
            raise TypeError(f'Argument sflux_1 must be of type {SfluxDataset},'
                            f' not type {type(sflux_1)}')
        self._sflux_1 = sflux_1

    @property
    def sflux_2(self) -> SfluxDataset:
        return self._sflux_2

    @sflux_2.setter
    def sflux_2(self, sflux_2):
        if sflux_2 is not None:
            if not isinstance(sflux_2, SfluxDataset):
                raise TypeError(
                    f'Argument sflux_2 must be of type {SfluxDataset}, not '
                    f'type {type(sflux_2)}.')
        self._sflux_2 = sflux_2

    @property
    def windrot(self):
        return self._windrot

    @windrot.setter
    def windrot(self, windrot: Union[gridgr3.Windrot, None]):
        if not isinstance(windrot, (gridgr3.Windrot, type(None))):
            raise TypeError(
                f'Argument windrot must be of type {gridgr3.Windrot} or None, '
                f'not type {type(windrot)}.')
        self._windrot = windrot

    @property
    def timevector(self):
        if self.sflux_2 is None:
            return self.sflux_1.timevector
        rnday_1 = self.sflux_1.timevector[-1] - self.sflux_1.timevector[0]
        rnday_2 = self.sflux_2.timevector[-1] - self.sflux_2.timevector[0]
        if rnday_2 <= rnday_1:
            return self.sflux_2.timevector
        else:
            return self.sflux_1.timevector
