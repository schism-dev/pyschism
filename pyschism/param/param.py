from datetime import timedelta
import logging
import pathlib
from typing import Union

import f90nml

# from pyschism.domain import ModelDomain
from pyschism.enums import Stratification
from pyschism.param.core import CORE
from pyschism.param.opt import OPT
from pyschism.param.schout import SCHOUT

# from pyschism.stations import Stations

# PARAM_TEMPLATE = pathlib.Path(__file__).parent / 'param.nml.template'
from pyschism.param.schism_init import GitParamTemplate

logger = logging.getLogger(__name__)


class Param:
    def __init__(
        self,
        dt: Union[int, float, timedelta] = None,
        rnday: Union[int, float, timedelta] = None,
        ibc: Union[Stratification, int, str] = Stratification.BAROTROPIC,
        nspool: Union[int, float, timedelta] = None,
        ihfskip: Union[int, timedelta] = None,
        template=None,
    ):

        self.core = CORE()
        self.core.ibc = ibc
        self.core.rnday = rnday
        self.core.dt = dt
        self.core.nspool = nspool
        self.core.ihfskip = ihfskip
        self.opt = OPT()
        self.schout = SCHOUT()
        self.template = template

    def __str__(self):
        return f"{str(self.core)}\n\n{str(self.opt)}\n\n{str(self.schout)}\n"

    def write(self, path, overwrite=False, use_template=False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise IOError(f"File {path} exists and overwrite=False")
        if use_template:
            PARAM_TEMPLATE = GitParamTemplate().path if use_template is True else use_template
            f90nml.patch(PARAM_TEMPLATE, self.to_dict(), path)
        else:
            with open(path, "w") as f:
                f.write(str(self))

    def to_dict(self):
        return {
            "CORE": self.core.to_dict(),
            "OPT": self.opt.to_dict(),
            "SCHOUT": self.schout.to_dict(),
        }

    # @property
    # def core(self):
    #     return self._core

    # @property
    # def opt(self):
    #     return self._opt

    # @property
    # def schout(self):
    #     return self._schout

    # @property
    # def model_domain(self):
    #     return self.__model_domain

    # @property
    # def nhot_write(self):
    #     return self.schout.nhot_write

    # @property
    # def stations(self):
    #     return self.__stations

    # @property
    # def _model_domain(self):
    #     return self.__model_domain

    # @_model_domain.setter
    # def _model_domain(self, model_domain):
    #     assert isinstance(model_domain, ModelDomain), \
    #         f"Argument model_domain must be of type {ModelDomain}, " \
    #         f"not {type(model_domain)}."
    #     self.__model_domain = model_domain

    # @property
    # def _opt(self):
    #     return self.__opt

    # @_opt.setter
    # def _opt(self, opt: OPT):
    #     # friction parameters
    #     opt.nchi = self.model_domain.fgrid
    #     # set coordinate system
    #     opt.ics = self.model_domain.ics
    #     # set coriolis
    #     opt.ncor = self.model_domain.ncor
    #     # set atmospheric forcing
    #     if self.model_domain.nws is not None:
    #         opt.nws = self.model_domain.nws
    #     # TODO: Set the remaining options:
    #     # msc2
    #     # mdc2
    #     # ntracer_gen
    #     # ntracer_age
    #     # sed_class
    #     # eco_class
    #     self.__opt = opt

    # @property
    # def _nhot_write(self):
    #     return self.__nhot_write

    # @_nhot_write.setter
    # def _nhot_write(self, nhot_write: Union[int, timedelta, bool, None]):

    #     if not isinstance(nhot_write, (int, bool, timedelta, type(None))):
    #         raise TypeError(f"Argument nhot_write must be of type {int}, "
    #                         f"{bool}, {timedelta}, or None.")

    #     if nhot_write is True:
    #         nhot_write = int(round(self.core.rnday / self.core.dt))

    #     elif isinstance(nhot_write, timedelta):
    #         nhot_write = int(round(nhot_write / self.core.dt))

    #     if nhot_write is not None:
    #         if nhot_write % self.core.ihfskip != 0:
    #             raise ValueError("nhot_write must be a multiple of ihfskip")
    #         self.schout.nhot_write = nhot_write

    # @property
    # def _stations(self):
    #     return self.__stations

    # @_stations.setter
    # def _stations(self, stations: Union[Stations, None]):
    #     assert isinstance(stations, (Stations, type(None))), \
    #         f"Argument stations must be of type {Stations} or None, " \
    #         f"not type {type(stations)}."
    #     if isinstance(stations, Stations):
    #         stations.transform_to(self.model_domain.hgrid.crs)
    #         stations.clip(
    #             self.model_domain.hgrid.get_bbox()
    #         )
    #         self.schout.iout_sta = 1
    #         self.schout.nspool_sta = stations.nspool_sta
    #     self.__stations = stations
