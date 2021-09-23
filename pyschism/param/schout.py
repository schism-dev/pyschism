from datetime import timedelta
import logging
from typing import Union

from pyschism.enums import (
    # IofWetdryVariables,
    # IofZcorVariables,
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


class SchoutMeta(type):

    surface_output_vars = SurfaceOutputVars()

    def __new__(meta, name, bases, attrs):
        for iof_type in meta.surface_output_vars.keys():
            attrs[f'_{iof_type}'] = len(SchoutType[iof_type].value)*[0]

        for iof_type, vardata in meta.surface_output_vars.items():
            for name, index in vardata:
                attrs[name] = OutputVariableDescriptor(iof_type, name, index)
        output_vars = []
        for iof_, outputs in meta.surface_output_vars.items():
            for name, id in outputs:
                output_vars.append(name)
        attrs['surface_output_vars'] = output_vars
        return type(name, bases, attrs)


class SCHOUT(
        metaclass=SchoutMeta
):
    """ Provides error checking implementation for SCHOUT group """

    def __init__(
            self,
            nhot_write: int = None,
            nspool_sta: int = None,
            **outputs
    ):
        """
        nhot_write:
            - if -1 will write last timestep (default)
            - if None it will be disabled.
            - if int interpreted as iteration
            - if timedelta it will be rounded to the nearest iteration
        """

        self.nhot_write = nhot_write
        self.nspool_sta = nspool_sta

        for key, val in outputs.items():
            setattr(self, key, val)

    def __iter__(self):
        for outvar in self._surface_output_vars:
            yield outvar, getattr(self, outvar)

    def __str__(self):
        schout = ["&SCHOUT"]
        if self.nhot_write is not None:
            schout.append(f"  nhot={self.nhot}")
            schout.append(f"  nhot_write={self.nhot_write}")
        if self.nspool_sta is not None:
            nspool_sta = self.nspool_sta
            if isinstance(nspool_sta, timedelta):
                nspool_sta = int(round(nspool_sta.total_seconds() / self.dt))
            if isinstance(nspool_sta, float):
                nspool_sta = int(
                    round(timedelta(hours=nspool_sta) / self.dt))
            if isinstance(nspool_sta, (int, float)):
                if nspool_sta <= 0:
                    raise ValueError("nspool_sta must be positive.")
            schout.append(f"  iout_sta={self.iout_sta}")
            schout.append(f"  nspool_sta={nspool_sta}")
        for var in dir(self):
            if var.startswith('iof'):
                for i, state in enumerate(getattr(self, var)):
                    if state == 1:
                        schout.append(f'  {var[1:]}({i+1})={state}')
        schout.append('/')
        return '\n'.join(schout)

    def to_dict(self):
        data = {}
        if self.nhot_write is not None:
            data['nhot'] = self.nhot
            data['nhot_write'] = self.nhot_write
        if self.nspool_sta is not None:
            nspool_sta = self.nspool_sta
            if isinstance(nspool_sta, timedelta):
                nspool_sta = int(round(nspool_sta.total_seconds() / self.dt))
            if isinstance(nspool_sta, float):
                nspool_sta = int(
                    round(timedelta(hours=nspool_sta) / self.dt))
            if isinstance(nspool_sta, (int, float)):
                if nspool_sta <= 0:
                    raise ValueError("nspool_sta must be positive.")
            data['iout_sta'] = self.iout_sta
            data['nspool_sta'] = nspool_sta
        for var in dir(self):
            if var.startswith('iof'):
                _var = var[1:]
                data[_var] = len(getattr(self, var)) * [0]
                for i, state in enumerate(getattr(self, var)):
                    if state == 1:
                        data[_var][i] = state
        return data

    @property
    def nhot_write(self):
        return self._nhot_write

    @nhot_write.setter
    def nhot_write(self, nhot_write: Union[int, None]):
        if nhot_write is not None:
            if not isinstance(nhot_write, int):
                raise TypeError(
                    f'Argument nhot_write must be of type {int} or None, not '
                    f'type {type(nhot_write)}.')
        self._nhot_write = nhot_write

    @property
    def nhot(self) -> Union[int, None]:
        if not hasattr(self, '_nhot') and self.nhot_write is not None:
            return 1
        else:
            return self._nhot

    @nhot.setter
    def nhot(self, nhot: Union[int, None]):
        if nhot not in [0, 1]:
            raise ValueError('Argument nhot must be 0, 1.')
        self._nhot = nhot

    @nhot.deleter
    def nhot(self):
        if hasattr(self, '_nhot'):
            del self._nhot

    @property
    def nspool_sta(self):
        return self._nspool_sta

    @nspool_sta.setter
    def nspool_sta(self, nspool_sta: Union[int, None]):
        if nspool_sta is not None:
            if not isinstance(nspool_sta, int):
                raise TypeError(
                    f'Argument nspool_sta must be of type {int} or None, not '
                    f'type {type(nspool_sta)}.')
        self._nspool_sta = nspool_sta

    @property
    def iout_sta(self) -> Union[int, None]:
        if not hasattr(self, '_iout_sta') and self.nspool_sta is not None:
            return 1
        else:
            return self._iout_sta

    @iout_sta.setter
    def iout_sta(self, iout_sta: Union[int, None]):
        if iout_sta not in [0, 1]:
            raise ValueError('Argument iout_sta must be 0, 1.')
        self._iout_sta = iout_sta

    @iout_sta.deleter
    def iout_sta(self):
        if hasattr(self, '_iout_sta'):
            del self._iout_sta
