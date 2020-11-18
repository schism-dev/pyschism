from datetime import datetime, timedelta
import pathlib
from typing import Union

from pyschism.enums import Stratification
from pyschism.domain import ModelDomain
from pyschism.param.core import CORE
from pyschism.param.opt import OPT
from pyschism.param.schout import SCHOUT


class OptDescriptor:

    def __set__(self, obj, opt: OPT):
        # friction parameters
        opt.nchi = obj.model_domain.fgrid
        # set coordinate system
        opt.ics = obj.model_domain.ics
        # set coriolis
        opt.ncor = obj.model_domain.ncor
        # set atmospheric forcing
        if obj.model_domain.nws is not None:
            opt.nws = obj.model_domain.nws
        # TODO: Set the remaining options:
        # msc2
        # mdc2
        # ntracer_gen
        # ntracer_age
        # sed_class
        # eco_class
        self._opt = opt

    def __get__(self, obj, val):
        return self._opt


class NhotWriteDescriptor:

    def __set__(self, obj, nhot_write: Union[int, bool, timedelta, None]):

        if not isinstance(nhot_write, (int, bool, timedelta, type(None))):
            raise TypeError(f"Argument nhot_write must be of type {int}, "
                            f"{bool}, {timedelta}, or None.")

        if nhot_write is True:
            nhot_write = int(round(obj.core.rnday / obj.core.dt))

        elif isinstance(nhot_write, timedelta):
            nhot_write = int(round(nhot_write / obj.core.dt))

        if nhot_write is not None:
            if nhot_write % obj.core.ihfskip != 0:
                raise ValueError("nhot_write must be a multiple of ihfskip")
            obj.schout.nhot_write = nhot_write

    def __get__(self, obj, val):
        return obj.schout.nhot_write


class Param:

    _opt = OptDescriptor()
    _nhot_write = NhotWriteDescriptor()

    def __init__(
            self,
            model_domain: ModelDomain,
            dt: Union[int, float, timedelta],
            rnday: Union[int, float, timedelta],
            dramp: Union[int, float, timedelta] = None,
            start_date: datetime = None,
            ibc: Union[Stratification, int, str] = Stratification.BAROTROPIC,
            drampbc: Union[int, float, timedelta] = None,
            nspool: Union[int, float, timedelta] = None,
            ihfskip: Union[int, timedelta] = None,
            nhot_write: Union[int, timedelta, bool] = None,
            **surface_outputs):
        self._model_domain = model_domain
        self._core = CORE(ibc, rnday, dt, nspool, ihfskip)
        self._opt = OPT(dramp, drampbc, start_date)
        self._schout = SCHOUT(**surface_outputs)
        self._nhot_write = nhot_write

    def __str__(self):
        return f"{str(self.core)}\n\n{str(self.opt)}\n\n{str(self.schout)}\n"

    def write(self, path, overwrite=False, use_template=False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise IOError(f"File {path} exists and overwrite=False")
        with open(path, 'w') as f:
            f.write(str(self))

    @property
    def core(self):
        return self._core

    @property
    def opt(self):
        return self._opt

    @property
    def schout(self):
        return self._schout

    @property
    def model_domain(self):
        return self._model_domain
