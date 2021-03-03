from datetime import timedelta
import logging
from typing import Union

from pyschism.enums import (
    IofHydroVariables,
    IofDvdVariables,
    IofWwmVariables,
    IofGenVariables,
    IofAgeVariables,
    IofSedVariables,
    IofEcoVariables,
    IofIcmVariables,
    IofCosVariables,
    IofFibVariables,
    IofSed2dVariables,
    IofMarshVariables,
    IofIceVariables,
    IofAnaVariables,
    SchoutType
)


_logger = logging.getLogger(__name__)


class OutputVariableDescriptor:

    def __init__(self, iof_type, name, index):
        self._iof_type = iof_type
        self._name = name
        self._index = index

    def __get__(self, obj, val):
        return bool(getattr(obj, f'_{self._iof_type}')[self._index])

    def __set__(self, obj, val: bool):
        if not isinstance(val, bool):
            raise TypeError(f'Argument to {self._name} must be boolean, not '
                            f'type {type(val)}.')
        iof = getattr(obj, f'_{self._iof_type}')
        iof[self._index] = int(val)


class SurfaceOutputVars:

    def __init__(self):
        self._surface_output_vars = {
            'iof_hydro': [(var.value, i) for i, var
                          in enumerate(IofHydroVariables)],
            'iof_dvd': [(var.value, i) for i, var
                        in enumerate(IofDvdVariables)],
            'iof_wwm': [(var.value, i) for i, var
                        in enumerate(IofWwmVariables)],
            'iof_gen': [(var.value, i) for i, var
                        in enumerate(IofGenVariables)],
            'iof_age': [(var.value, i) for i, var
                        in enumerate(IofAgeVariables)],
            'iof_sed': [(var.value, i) for i, var
                        in enumerate(IofSedVariables)],
            'iof_eco': [(var.value, i) for i, var
                        in enumerate(IofEcoVariables)],
            'iof_icm': [(var.value, i) for i, var
                        in enumerate(IofIcmVariables)],
            'iof_cos': [(var.value, i) for i, var
                        in enumerate(IofCosVariables)],
            'iof_fib': [(var.value, i) for i, var
                        in enumerate(IofFibVariables)],
            'iof_sed2d': [(var.value, i) for i, var
                          in enumerate(IofSed2dVariables)],
            'iof_marsh': [(var.value, i) for i, var
                          in enumerate(IofMarshVariables)],
            'iof_ice': [(var.value, i) for i, var
                        in enumerate(IofIceVariables)],
            'iof_ana': [(var.value, i) for i, var
                        in enumerate(IofAnaVariables)],
        }

    def __get__(self, obj, val):
        return self._surface_output_vars


class SchoutMeta(type):

    surface_output_vars = SurfaceOutputVars()

    def __new__(meta, name, bases, attrs):
        for iof_type in meta.surface_output_vars.keys():
            attrs[f'_{iof_type}'] = len(SchoutType[iof_type].value)*[0]

        for iof_type, vardata in meta.surface_output_vars.items():
            for name, index in vardata:
                attrs[name] = OutputVariableDescriptor(iof_type, name, index)
        attrs['surface_output_vars'] = meta.surface_output_vars
        return type(name, bases, attrs)


class Nhot:

    def __set__(self, obj, nhot: int):
        if nhot not in [0, 1]:
            raise ValueError(f"nhot must be 0 or 1, not {nhot}")
        obj.__dict__['nhot'] = nhot

    def __get__(self, obj, val):
        return obj.__dict__.get('nhot')


class NhotWrite:

    def __set__(self, obj, nhot_write: int):
        obj.__dict__['nhot_write'] = nhot_write
        obj.__dict__['nhot'] = 1

    def __get__(self, obj, val):
        return obj.__dict__.get('nhot_write')


class IoutSta:

    def __set__(self, obj, iout_sta: int):
        if iout_sta not in [0, 1]:
            raise ValueError(f"iout_sta must be 0 or 1, not {iout_sta}")
        obj.__dict__['iout_sta'] = iout_sta

    def __get__(self, obj, val):
        return obj.__dict__.get('iout_sta')


class NspoolSta:

    def __set__(self, obj, nspool_sta: Union[int, timedelta]):
        obj.__dict__['nspool_sta'] = nspool_sta
        obj.__dict__['iout_sta'] = 1

    def __get__(self, obj, val):
        return obj.__dict__.get('nspool_sta')


class SCHOUT(metaclass=SchoutMeta):
    """ Provides error checking implementation for SCHOUT group """
    _iout_sta = IoutSta()
    _nhot = Nhot()
    nhot_write = NhotWrite()
    nspool_sta = NspoolSta()

    def __init__(self, **outputs):
        _logger.info('Initializing SCHOUT.')
        for key, val in outputs.items():
            setattr(self, key, val)

    def __iter__(self):
        for outvar in self._surface_output_vars:
            yield outvar, getattr(self, outvar)

    def __str__(self):
        schout = ["&SCHOUT"]
        if self.nhot_write is not None:
            schout.append(f"  nhot={self._nhot}")
            schout.append(f"  nhot_write={self.nhot_write}")
        if self.nspool_sta is not None:
            schout.append(f"  iout_sta={self._iout_sta}")
            schout.append(f"  nspool_sta={self.nspool_sta}")
        for var in dir(self):
            if var.startswith('_iof'):
                for i, state in enumerate(getattr(self, var)):
                    if state == 1:
                        schout.append(f'  {var[1:]}({i+1})={state}')
        schout.append('/')
        return '\n'.join(schout)
