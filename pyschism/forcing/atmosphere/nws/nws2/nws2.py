import os
import pathlib
from typing import Union

import f90nml  # type: ignore[import]

from pyschism.enums import NWSType
from pyschism.mesh.base import Gr3
from pyschism.forcing.atmosphere.nws.nws import NWS
from pyschism.forcing.atmosphere.nws.nws2.sflux import SfluxDataset


SFLUX_DEFAULTS = f90nml.read(
    pathlib.Path(__file__).parent / 'sflux_inputs.txt')


class Sflux1Descriptor:

    def __set__(self, obj, sflux_1: SfluxDataset):
        if not isinstance(sflux_1, SfluxDataset):
            raise TypeError(f'Argument sflux_1 must be of type {SfluxDataset},'
                            f' not type {type(sflux_1)}')
        obj.__dict__['sflux_1'] = sflux_1

    def __get__(self, obj, val):
        return obj.__dict__['sflux_1']


class Sflux2Descriptor:

    def __set__(self, obj, sflux_2: Union[SfluxDataset, None]):
        if sflux_2 is not None:
            if not isinstance(sflux_2, SfluxDataset):
                raise TypeError(f'Argument sflux_2 must be of type {SfluxDataset},'
                                f' not type {type(sflux_2)}.')
            obj.__dict__['sflux_2'] = sflux_2

    def __get__(self, obj, val):
        return obj.__dict__.get('sflux_2')


class WindRotDescriptor:

    def __set__(self, obj, grd: Gr3):
        if not isinstance(grd, Gr3):
            raise TypeError(f'Argument grd must be of type {Gr3}, '
                            f'not type {type(grd)}.')
        data = grd.to_dict()
        data.pop("boundaries", None)
        data["nodes"] = {id: (coord, 0.) for id, (coord, _)
                         in data["nodes"].items()}
        obj.__dict__['windrot'] = Gr3(**data)

    def __get__(self, obj, val):
        return obj.__dict__.get('windrot')


class NWS2(NWS):

    _sflux_1 = Sflux1Descriptor()
    _sflux_2 = Sflux2Descriptor()
    _windrot = WindRotDescriptor()

    def __init__(self, sflux_1: SfluxDataset,
                 sflux_2: SfluxDataset = None):
        """Loads SfluxDataset to use as NWS2 input. """

        # for key, item in SFLUX_DEFAULTS['sflux_inputs'].items():
        #     print(key, item)
        #     setattr(self, key, item)

        self._sflux_1 = sflux_1
        self._sflux_2 = sflux_2

    def __call__(self, model_driver):
        super().__call__(model_driver)
        if hasattr(self.sflux_1, 'fetch_data'):
            self.sflux_1.fetch_data(
                self._start_date, self._rnday,
                bbox=model_driver.model_domain.hgrid.get_bbox('EPSG:4326'))
        if self.sflux_2 is not None:
            if hasattr(self.sflux_2, 'fetch_data'):
                self.sflux_2.fetch_data(
                    self._start_date, self._rnday,
                    bbox=model_driver.model_domain.hgrid.get_bbox('EPSG:4326'))
        self._windrot = model_driver.model_domain.hgrid

    def __str__(self):
        data = []
        data = '\n'.join(data)
        return f'&sflux_inputs\n{data}/\n'

    def write(self, path: Union[str, os.PathLike], overwrite: bool = False,
              wind_rot: bool = True):
        # write sflux namelist
        path = pathlib.Path(path)
        if path.name != 'sflux':
            path /= 'sflux'
        path.mkdir(exist_ok=True)
        with open(path / 'sflux_inputs.txt', 'w') as f:
            f.write(str(self))
        # write sflux data
        self.sflux_1.write(
            path, 1, overwrite, start_date=self._start_date, rnday=self._rnday)
        if self.sflux_2 is not None:
            self.sflux_2.write(
                path, 2, overwrite, start_date=self._start_date,
                rnday=self._rnday)
        # # write windrot data
        if self._windrot is not None:
            if wind_rot is True:
                self._windrot.write(path.parent / 'windrot_geo2proj.gr3',
                                    overwrite)

    @property
    def dtype(self) -> NWSType:
        """Returns the datatype of the object"""
        return NWSType(2)

    @property
    def sflux_1(self) -> SfluxDataset:
        return self._sflux_1

    @property
    def sflux_2(self) -> SfluxDataset:
        return self._sflux_2
