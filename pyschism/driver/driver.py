from datetime import datetime, timedelta
from enum import Enum
import os
import pathlib
import subprocess
from typing import Union, List

import numpy as np

from pyschism import dates
from pyschism.enums import Stratification
from pyschism.forcing import Tides, NWS, Hydrology
from pyschism.mesh import Hgrid, Vgrid, Fgrid, ManningsN, gridgr3
from pyschism.param import Param
from pyschism.server.base import ServerConfig
from pyschism.stations import Stations


def raise_type_error(argname, obj, cls):
    raise TypeError(
        f'Argument {argname} must be of type {cls}, not '
        f'type {type(obj)}.')


class ModelForcings:

    def fetch_data(self, hgrid: Hgrid, start_date, rnday):

        if not isinstance(hgrid, Hgrid):
            raise_type_error('hgrid', hgrid, Hgrid)

        if self.tides is not None:
            self.tides.fetch_data(hgrid, start_date, rnday)

        if self.atmosphere is not None:
            self.atmosphere.fetch_data(hgrid, start_date, rnday)

        if self.hydrology is not None:
            for forcing in self.hydrology:
                forcing.fetch_data(hgrid, start_date, rnday)

        if self.baroclinic is not None:
            self.baroclinic.fetch_data(hgrid, start_date, rnday)

        if self.waves is not None:
            self.waves.fetch_data(hgrid, start_date, rnday)

    @property
    def tides(self):
        return self._tides

    @tides.setter
    def tides(self, tides: Union[Tides, None]):
        if tides is not None:
            if not isinstance(tides, Tides):
                raise_type_error('tides', tides, Tides)
        self._tides = tides

    @property
    def atmosphere(self):
        return self._atmosphere

    @atmosphere.setter
    def atmosphere(self, atmosphere: Union[NWS, None]):
        if atmosphere is not None:
            if not isinstance(atmosphere, NWS):
                raise_type_error('atmosphere', atmosphere, NWS)
        self._atmosphere = atmosphere

    @property
    def hydrology(self):
        return self._hydrology

    @hydrology.setter
    def hydrology(self, hydrology: Union[Hydrology, List[Hydrology], None]):
        if hydrology is not None:
            if not isinstance(hydrology, Hydrology):
                hydrology = list(hydrology)
            for forcing in hydrology:
                raise_type_error('hydrology', forcing, Hydrology)
        self._hydrology = hydrology

    @property
    def baroclinic(self):
        return self._baroclinic

    @baroclinic.setter
    def baroclinic(self, baroclinic: None):
        if baroclinic is not None:
            raise NotImplementedError(
                'baroclinic forcing not yet implemented.')
        self._baroclinic = baroclinic

    @property
    def waves(self):
        return self._waves

    @waves.setter
    def waves(self, waves: None):
        if waves is not None:
            raise NotImplementedError(
                'waves forcing not yet implemented.')
        self._waves = waves


class Gr3FieldTypes(Enum):
    ALBEDO = gridgr3.Albedo
    DIFFMIN = gridgr3.Diffmin
    DIFFMAX = gridgr3.Diffmax
    WATERTYPE = gridgr3.Watertype
    FLUXFLAG = gridgr3.Fluxflag
    TVDFLAG = gridgr3.Tvdflag
    WINDROT = gridgr3.Windrot

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f'{name} is not a valid {gridgr3.Gr3Field} type.')


class Gridgr3Descriptor:

    def __init__(self, gridgr3_type: Gr3FieldTypes):
        self.name = gridgr3_type.name
        self.type = gridgr3_type.value
        self.gr3field = None

    def __set__(self, obj, val: Union[gridgr3.Gr3Field, None]):
        if not isinstance(val, self.type) and val is not None:
            raise ValueError(
                f'Argument {self.name.lower()} must be of type {self.type} '
                f'not type {type(val)}.')
        self.gr3field = val

    def __get__(self, obj, val):
        return self.gr3field


class ModelConfigMeta(type):

    def __new__(meta, name, bases, attrs):
        attrs['forcings'] = ModelForcings()
        attrs['start_date'] = dates.StartDate()
        attrs['end_date'] = dates.EndDate()
        attrs['spinup_time'] = dates.SpinupTime()
        for gr3field_type in Gr3FieldTypes:
            name = gr3field_type.name.lower()
            attrs[name] = Gridgr3Descriptor(gr3field_type)
        return type(name, bases, attrs)


class ModelConfig(metaclass=ModelConfigMeta):
    """Class representing a SCHISM model configuration.

    This class combines the horizontal grid (hgrid), vertical grid (vgrid)
    and friction/drag grids (fgrid). Additionally, this class holds
    information about forcings.

    Args:
        hgrid: :class:`pyschism.mesh.Hgrid` instance.
        vgrid: :class:`pyschism.mesh.Vgrid` instance.
        fgrid: :class:`pyschism.mesh.Fgrid` derived instance.
    """

    def __init__(
            self,
            hgrid: Hgrid,
            vgrid: Vgrid = None,
            fgrid: Fgrid = None,
            albedo: gridgr3.Albedo = None,
            diffmin: gridgr3.Diffmin = None,
            diffmax: gridgr3.Diffmax = None,
            watertype: gridgr3.Watertype = None,
            fluxflag: gridgr3.Fluxflag = None,
            tvdflag: gridgr3.Tvdflag = None,
            windrot: gridgr3.Windrot = None,
            tides: Tides = None,
            atmosphere: NWS = None,
            hydrology: Union[Hydrology, List[Hydrology]] = None,
            baroclinic=None,
            waves=None,
    ):
        self.hgrid = hgrid
        self.vgrid = vgrid
        self.fgrid = fgrid
        self.albedo = albedo
        self.diffmin = diffmin
        self.diffmax = diffmax
        self.watertype = watertype
        self.fluxflag = fluxflag
        self.tvdflag = tvdflag
        self.windrot = windrot
        self.forcings.tides = tides
        self.forcings.atmosphere = atmosphere
        self.forcings.hydrology = hydrology
        self.forcings.baroclinic = baroclinic
        self.forcings.waves = waves

    def coldstart(
            self,
            timestep: Union[float, timedelta] = 150.,
            start_date: datetime = None,
            rnday: Union[datetime, timedelta] = None,
            dramp: Union[float, timedelta] = None,
            drampbc: Union[float, timedelta] = None,
            nspool: Union[int, timedelta] = None,
            ihfskip: int = None,
            stations: Stations = None,
            server_config: ServerConfig = None,
            **surface_outputs
    ) -> 'ModelDriver':
        return ModelDriver(
            self,
            timestep=timestep,
            start_date=start_date,
            rnday=rnday,
            dramp=dramp,
            drampbc=drampbc,
            server_config=server_config,
            **surface_outputs,
        )

    def hotstart(
            self,
            hotstart: Union[str, os.PathLike, 'CombineHotstart'],
            timestep: Union[float, timedelta] = 150.,
            rnday: Union[datetime, timedelta] = None,
            nspool: Union[int, timedelta] = None,
            ihfskip: int = None,
            stations: Stations = None,
            **surface_outputs
    ) -> 'ModelDriver':
        return ModelDriver(
            self,
            timestep=timestep,
            start_date=start_date,
            rnday=rnday,
            combine_hotstart=hotstart,
            nspool=nspool,
            ihfskip=ihfskip,
            stations=stations,
        )

    @property
    def hgrid(self):
        return self._hgrid

    @hgrid.setter
    def hgrid(self, hgrid: Hgrid):
        if not isinstance(hgrid, Hgrid):
            raise_type_error('hgrid', hgrid, Hgrid)
        self._hgrid = hgrid

    @property
    def vgrid(self):
        return self._vgrid

    @vgrid.setter
    def vgrid(self, vgrid: Union[Vgrid, None]):
        if vgrid is None:
            vgrid = Vgrid()
        if not isinstance(vgrid, Vgrid):
            raise_type_error('vgrid', vgrid, Vgrid)
        self._vgrid = vgrid

    @property
    def fgrid(self):
        return self._fgrid

    @fgrid.setter
    def fgrid(self, fgrid: Union[Fgrid, None]):
        if fgrid is None:
            if self.vgrid.is2D():
                fgrid = ManningsN.linear_with_depth(self.hgrid)
            else:
                raise NotImplementedError('Must provide fgrid for 3D runs.')
        if not isinstance(fgrid, Fgrid):
            raise_type_error('fgrid', fgrid, Fgrid)
        if self.vgrid.is2D() is True and not isinstance(fgrid, ManningsN):
            raise TypeError(
                f'2D model must use {ManningsN} but got {type(fgrid)}.')
        self._fgrid = fgrid

    @property
    def timestep(self):
        return self._timestep

    @timestep.setter
    def timestep(self, timestep: Union[float, timedelta, None]):
        if timestep is None:
            self._timestep = timedelta(seconds=150.)
        elif not isinstance(timestep, timedelta):
            self._timestep = timedelta(seconds=float(timestep))
        return self._timestep


class ModelDriver:

    def __init__(
            self,
            config: ModelConfig,
            dt: Union[float, timedelta],
            rnday: Union[float, timedelta],
            ihfskip: int = None,
            dramp: Union[float, timedelta] = None,
            start_date: datetime = None,
            drampbc: Union[float, timedelta] = None,
            stations: Stations = None,
            nspool: Union[int, timedelta] = None,
            # nhot_write: Union[int, timedelta, bool] = None,
            server_config: ServerConfig = None,
            combine_hotstart: Union[str, os.PathLike,
                                    'CombineHotstart'] = None,
            # cutoff_depth: float = 50.,
            **surface_outputs
    ):
        self.config = config
        self.param = Param()

        # set core parameters
        self.param.core.dt = dt
        self.param.core.rnday = rnday
        self.param.core.nspool = nspool
        self.param.core.ihfskip = ihfskip
        if config.vgrid.is2D():
            self.param.core.ibc = Stratification.BAROTROPIC
        else:
            self.param.core.ibc = Stratification.BAROCLINIC
            raise NotImplementedError('Must set ibtp')
        # TODO: must also set msc2/mdc2 here.

        # set opt
        self.param.opt.start_date = start_date
        self.param.opt.ics = 2 if self.config.hgrid.crs.is_geographic() else 1
        self.param.opt.dramp = dramp
        self.param.opt.drampbc = drampbc

        self.combine_hotstart = combine_hotstart

    def run(self):
        pass

    @property
    def config(self) -> ModelConfig:
        return self._config

    @config.setter
    def config(self, config: ModelConfig):
        if not isinstance(config, ModelConfig):
            raise_type_error('config', config, ModelConfig)
        self._config = config

    @property
    def combine_hotstart(self):
        return self._combine_hotstart

    @combine_hotstart.setter
    def combine_hotstart(
            self,
            combine_hotstart: Union[str, os.PathLike, 'CombineHotstart', None]
    ):

        if combine_hotstart is None:
            pass

        elif not isinstance(combine_hotstart, 'CombineHotstart'):
            combine_hotstart = CombineHotstartBinary(combine_hotstart)
            self.param.opt.ihot = 1

        self._combine_hotstart = combine_hotstart

    def write(self):
        pass
        # """Writes to disk the full set of input files necessary to run SCHISM.
        # """
        # outdir = pathlib.Path(output_directory)

        # if not (outdir / 'outputs').exists():
        #     (outdir / 'outputs').mkdir(parents=True)

        # if hgrid:
        #     hgrid = 'hgrid.gr3' if hgrid is True else hgrid
        #     self.model_domain.hgrid.write(outdir / hgrid, overwrite)
        #     if self.model_domain.ics == 2:
        #         original_dir = os.getcwd()
        #         os.chdir(outdir)  # pushd
        #         try:
        #             os.remove('hgrid.ll')
        #         except OSError:
        #             pass
        #         os.symlink(hgrid, 'hgrid.ll')
        #         os.chdir(original_dir)  # popd

        # if self.model_domain.albedo is not None:
        #     # self.model_domain.albedo = Albedo.constant(self.model_domain.hgrid, 0.15)
        #     self.model_domain.albedo.write(outdir / 'albedo.gr3', overwrite)

        # if diffmax:
        #     self.diffmax = Diffmax.constant(self.model_domain.hgrid, 1.0)
        #     self.diffmax.write(outdir / 'diffmax.gr3', overwrite)

        # if diffmin:
        #     self.diffmin = Diffmax.constant(self.model_domain.hgrid, 1.0e-6)
        #     self.diffmin.write(outdir / 'diffmin.gr3', overwrite)

        # if watertype:
        #     self.watertype = Diffmax.constant(self.model_domain.hgrid, 1.0)
        #     self.watertype.write(outdir / 'watertype.gr3', overwrite)

        # if fluxflag:
        #     self.fluxflag = Fluxflag.constant(self.model_domain.hgrid, -1)
        #     with open(outdir / 'fluxflag.prop', 'w+') as fid:
        #         fid.writelines(self.fluxflag)

        # if tvdflag:
        #     # Hard-wire the polygon at this point.
        #     coords = [(-75.340506, 40.843483), (-75.533474, 40.436019), (-75.796036, 39.535807),
        #               (-75.672664, 39.339972), (-75.305709,
        #                                         39.460000), (-75.103251, 39.636884),
        #               (-74.692008, 39.744277), (-74.391485,
        #                                         40.009603), (-74.359851, 40.252818),
        #               (-74.514858, 40.745565), (-74.834362,
        #                                         40.957194), (-75.210807, 40.935083),
        #               (-75.283565, 40.925607)]
        #     poly = Polygon(coords)
        #     self.tvdflag = Tvdflag.define_by_region(
        #         hgrid=self.model_domain.hgrid, region=poly, value=1)
        #     with open(outdir / 'tvd.prop', 'w+') as fid:
        #         fid.writelines(self.tvdflag)

        # if rtofs:
        #     self.start_date = nearest_cycle_date()
        #     self.hotstart = HotStartInventory()
        #     self.hotstart.fetch_data(
        #         outdir, self.model_domain.hgrid, self.start_date)
        #     self.obnd = OpenBoundaryInventory()

        # if vgrid:
        #     vgrid = 'vgrid.in' if vgrid is True else vgrid
        #     self.model_domain.vgrid = Vgrid.from_binary(
        #         self.model_domain.hgrid)
        #     # TODO: Unhardcode indexes
        #     self.obnd.fetch_data(outdir, self.start_date, rnday=3,
        #                          idx_min=2687, idx_max=2714, jdx_min=1181,
        #                          jdx_max=1634)

        # if fgrid:
        #     fgrid = f'{self.model_domain.fgrid.fname}' if fgrid is True \
        #             else fgrid
        #     self.model_domain.fgrid.write(outdir / fgrid, overwrite)
        # if param:
        #     param = 'param.nml' if param is True else param
        #     self.param.write(outdir / param, overwrite,
        #                      use_template=use_param_template)
        # if bctides:
        #     bctides = 'bctides.in' if bctides is True else bctides
        #     self.bctides.write(outdir / bctides, overwrite)
        # if nws:
        #     if self.nws is not None:
        #         if isinstance(self.nws, NWS2):
        #             self.nws.write(outdir / 'sflux', overwrite,
        #                            wind_rot=wind_rot)
        #         else:
        #             self.nws.write(outdir, overwrite)
        # if stations:
        #     if self.stations is not None:
        #         stations = 'station.in' if stations is True else stations
        #         self.stations.write(outdir / stations, overwrite)

        # for hydrology in self._hydrology:
        #     hydrology.write(outdir, overwrite)

        # if self.hotstart_file is not None:
        #     hotstart_lnk = outdir / 'hotstart.nc'
        #     if overwrite is True:
        #         if hotstart_lnk.exists():
        #             hotstart_lnk.unlink()
        #     os.symlink(os.path.relpath(
        #         self.hotstart_file, outdir), hotstart_lnk)

        # self.makefile.write(outdir / 'Makefile', overwrite)


class CombineHotstart:
    pass


class CombineHotstartBinary(CombineHotstart):

    def __init__(self, path: Union[str, os.PathLike], iteration=None):
        path = pathlib.Path(path)
        if iteration is None:
            combine_hotstart = path.glob('hotstart_[0-9][0-9][0-9][0-9]_*.nc')
            increments = set([file.name.split('_')[-1].split('.nc')[0]
                              for file in combine_hotstart])
            iteration = np.max([int(increment) for increment in increments])
        subprocess.check_call(
            ["combine_hotstart7", '-i', f'{iteration}'], cwd=path)
        self.path = path / f"hotstart_it={iteration}.nc"


# class ModelDriver:
#     """Main driver class used to generate SCHISM input files

#     This __init__ calls initialization routines for all SCHISM data inputs
#     and prepares them to be dumped to files by the write() method.

#     Arguments:
#         domain: :class:`pyschism.domain.ModelDomain`
#         # param: :class:`pyschism.param.param.Param`
#     """

#     def __init__(
#             self,
#             config: ModelConfig,
#             timestep: Union[float, timedelta],
#             start_date: datetime = None,
#             end_date: Union[datetime, timedelta] = None,
#             spinup_time: Union[datetime, timedelta] = None,
#             server_config: ServerConfig = None,
#     ):

#         self.config = config

#         # if not isinstance(config, ModelConfig):
#         #     raise

#         # set param instance
#         # self.param = param

#     # @classmethod
#     # def coldstart(cls, config: ModelConfig):
#     #     obj = cls(config)
#         # # init param.nml
#         # self._param = Param(domain, dt, rnday, dramp, start_date, ibc,
#         #                     drampbc, nspool, ihfskip, nhot_write, stations,
#         #                     **surface_outputs)

#         # # init bctides.in
#         # self._bctides = Bctides(domain, self.param,
#         #                         cutoff_depth=cutoff_depth)

#         # # init atmospheric data.
#         # self._nws = domain.nws

#         # # init hydrology data
#         # self._hydrology = domain.hydrology

#         # # init hotstart file
#         # self._combine_hotstart = combine_hotstart

#         # if domain.albedo is not None:
#         #     self.param.opt.albedo = 1

#         # # do same here

#         # # init Makefile drivers
#         # self._makefile = MakefileDriver(server_config=server_config)

#     # @classmethod
#     # def hotstart(cls, previous_directory):
#     #     pass

#     def write(
#             self,
#             output_directory,
#             overwrite=False,
#             hgrid=True,
#             vgrid=True,
#             fgrid=True,
#             param=True,
#             bctides=True,
#             nws=True,
#             wind_rot=True,
#             stations=True,
#             use_param_template=True,
#             albedo=True,
#             diffmax=True,
#             diffmin=True,
#             watertype=True,
#             fluxflag=True,
#             tvdflag=True,
#             rtofs=True,
#     ):
#         """Writes to disk the full set of input files necessary to run SCHISM.
#         """
#         outdir = pathlib.Path(output_directory)

#         if not (outdir / 'outputs').exists():
#             (outdir / 'outputs').mkdir(parents=True)

#         if hgrid:
#             hgrid = 'hgrid.gr3' if hgrid is True else hgrid
#             self.domain.hgrid.write(outdir / hgrid, overwrite)
#             if self.domain.ics == 2:
#                 original_dir = os.getcwd()
#                 os.chdir(outdir)  # pushd
#                 try:
#                     os.remove('hgrid.ll')
#                 except OSError:
#                     pass
#                 os.symlink(hgrid, 'hgrid.ll')
#                 os.chdir(original_dir)  # popd

#         if self.domain.albedo is not None:
#             # self.domain.albedo = Albedo.constant(self.domain.hgrid, 0.15)
#             self.domain.albedo.write(outdir / 'albedo.gr3', overwrite)

#         if diffmax:
#             self.diffmax = Diffmax.constant(self.domain.hgrid, 1.0)
#             self.diffmax.write(outdir / 'diffmax.gr3', overwrite)

#         if diffmin:
#             self.diffmin = Diffmax.constant(self.domain.hgrid, 1.0e-6)
#             self.diffmin.write(outdir / 'diffmin.gr3', overwrite)

#         if watertype:
#             self.watertype = Diffmax.constant(self.domain.hgrid, 1.0)
#             self.watertype.write(outdir / 'watertype.gr3', overwrite)

#         if fluxflag:
#             self.fluxflag = Fluxflag.constant(self.domain.hgrid, -1)
#             with open(outdir / 'fluxflag.prop', 'w+') as fid:
#                 fid.writelines(self.fluxflag)

#         if tvdflag:
#             # Hard-wire the polygon at this point.
#             coords = [(-75.340506, 40.843483), (-75.533474, 40.436019), (-75.796036, 39.535807),
#                       (-75.672664, 39.339972), (-75.305709,
#                                                 39.460000), (-75.103251, 39.636884),
#                       (-74.692008, 39.744277), (-74.391485,
#                                                 40.009603), (-74.359851, 40.252818),
#                       (-74.514858, 40.745565), (-74.834362,
#                                                 40.957194), (-75.210807, 40.935083),
#                       (-75.283565, 40.925607)]
#             poly = Polygon(coords)
#             self.tvdflag = Tvdflag.define_by_region(
#                 hgrid=self.domain.hgrid, region=poly, value=1)
#             with open(outdir / 'tvd.prop', 'w+') as fid:
#                 fid.writelines(self.tvdflag)

#         if rtofs:
#             self.start_date = nearest_cycle_date()
#             self.hotstart = HotStartInventory()
#             self.hotstart.fetch_data(
#                 outdir, self.domain.hgrid, self.start_date)
#             self.obnd = OpenBoundaryInventory()
#             self.obnd.fetch_data(outdir, self.start_date, rnday=3,
#                                  idx_min=2687, idx_max=2714, jdx_min=1181, jdx_max=1634)

#         if vgrid:
#             vgrid = 'vgrid.in' if vgrid is True else vgrid
#             self.domain.vgrid = Vgrid.from_binary(
#                 self.domain.hgrid)

#             #self.domain.vgrid.write(outdir / vgrid, overwrite)
# # lcui
#         if fgrid:
#             fgrid = f'{self.domain.fgrid.fname}' if fgrid is True \
#                     else fgrid
#             self.domain.fgrid.write(outdir / fgrid, overwrite)
#         if param:
#             param = 'param.nml' if param is True else param
#             self.param.write(outdir / param, overwrite,
#                              use_template=use_param_template)
#         if bctides:
#             bctides = 'bctides.in' if bctides is True else bctides
#             self.bctides.write(outdir / bctides, overwrite)
#         if nws:
#             if self.nws is not None:
#                 if isinstance(self.nws, NWS2):
#                     self.nws.write(outdir / 'sflux', overwrite,
#                                    wind_rot=wind_rot)
#                 else:
#                     self.nws.write(outdir, overwrite)
#         if stations:
#             if self.stations is not None:
#                 stations = 'station.in' if stations is True else stations
#                 self.stations.write(outdir / stations, overwrite)

#         for hydrology in self._hydrology:
#             hydrology.write(outdir, overwrite)

#         if self.hotstart_file is not None:
#             hotstart_lnk = outdir / 'hotstart.nc'
#             if overwrite is True:
#                 if hotstart_lnk.exists():
#                     hotstart_lnk.unlink()
#             os.symlink(
#                 os.path.relpath(self.hotstart_file, outdir),
#                 hotstart_lnk
#             )

#         self.makefile.write(outdir / 'Makefile', overwrite)

#     @property
#     def config(self):
#         return self._config

#     @config.setter
#     def config(self, config: ModelConfig):
#         raise_type_error('config', config, ModelConfig)
#         self._config = config

#     @property
#     def param(self):
#         return self._param

#     @param.setter
#     def param(self, param: Union[Param, None]):
#         if not isinstance(param, Param):
#             raise_type_error('param', param, Param)
#         self._param = param


# --- references:

    # @property
    # def param(self):
    #     return self._param

    # @property
    # def domain(self):
    #     return self.param.domain

    # @property
    # def stations(self):
    #     return self.param.stations

    # @property
    # def bctides(self):
    #     return self._bctides

    # @property
    # def nws(self):
    #     return self._nws

    # @property
    # def hydrology(self):
    #     return self._hydrology

    # @property
    # def hotstart_file(self):
    #     if self._combine_hotstart is not None:
    #         return self._combine_hotstart.path

    # @property
    # def makefile(self):
    #     return self._makefile

    # @property
    # def _nws(self):
    #     return self.__nws

    # @_nws.setter
    # def _nws(self, nws):
    #     if nws is not None:
    #         nws(self)
    #     self.__nws = nws

    # @property
    # def _combine_hotstart(self):
    #     return self.__combine_hotstart

    # @_combine_hotstart.setter
    # def _combine_hotstart(self, combine_hotstart):
    #     if combine_hotstart is not None:
    #         combine_hotstart = CombineHotstartBinary(combine_hotstart)
    #         self.param.opt.ihot = 1
    #     self.__combine_hotstart = combine_hotstart

    # @property
    # def _hydrology(self):
    #     return self.__hydrology

    # @_hydrology.setter
    # def _hydrology(self, hydrology):
    #     for forcing in hydrology:
    #         if callable(forcing):
    #             forcing(self)
    #     self.__hydrology = hydrology


# from datetime import datetime, timedelta
# import os
# import pathlib
# import subprocess
# from typing import Union

# from netCDF4 import Dataset
# import numpy as np
# from shapely.geometry import Polygon, MultiPolygon, Point, MultiPoint

# from pyschism.domain import ModelDomain
# from pyschism.driver.makefile import MakefileDriver
# from pyschism.enums import Stratification
# from pyschism.forcing.tides.bctides import Bctides
# from pyschism.forcing.atmosphere import NWS2
# from pyschism.forcing.hycom.hycom import HotStartInventory, OpenBoundaryInventory
# from pyschism.param import Param
# from pyschism.server import ServerConfig
# from pyschism.stations import Stations
# from pyschism.mesh.gridgr3 import Albedo, Diffmax, Diffmin, Watertype
# from pyschism.mesh.prop import Fluxflag, Tvdflag
# from pyschism.mesh.vgrid import Vgrid
# from pyschism.dates import pivot_time, localize_datetime, nearest_cycle_date


# from functools import lru_cache
# import os
# from typing import Union, List

# import numpy as np  # type: ignore[import]
# from pyproj import CRS  # type: ignore[import]

# from pyschism.forcing.tides.bctypes import BoundaryCondition
# from pyschism.forcing.tides.tides import Tides
# from pyschism.forcing.atmosphere.nws import NWS
# from pyschism.forcing.hydrology.base import Hydrology
# # from pyschism.forcing.atmosphere.nws.nws2 import NWS2
# from pyschism.mesh import Hgrid, Vgrid, Fgrid
# from pyschism.enums import Coriolis


# class ModelConfig:

#     def __init__(
#             self,
#             hgrid: Hgrid,
#             vgrid: Vgrid = None,
#             fgrid: Fgrid = None,
#             param: Param = None,
#     ):
    # """Class representing a SCHISM computational domain.

    # This class combines the horizontal grid (hgrid), vertical grid (vgrid)
    # and friction/drag grids (fgrid). Additionally, this class holds
    # information about forcings.
    # Args:
    #     hgrid: :class:`pyschism.mesh.Hgrid` instance.
    #     vgrid: :class:`pyschism.mesh.Vgrid` instance.
    #     fgrid: :class:`pyschism.mesh.Fgrid` derived instance.
    # """
    # self._hgrid = hgrid
    # self._vgrid = vgrid
    # self._fgrid = fgrid
    # self._open_boundaries = OpenBoundaries(hgrid)
    # self._ncor = Coriolis.AUTO
    # self._nws: Union[NWS, None] = None
    # self._hydrology: List[Hydrology] = []

    # @staticmethod
    # def open(hgrid: Union[str, os.PathLike], fgrid: Union[str, os.PathLike],
    #          vgrid: os.PathLike = None, hgrid_crs: Union[str, CRS] = None,
    #          fgrid_crs: Union[str, CRS] = None):
    #     """Open files from disk"""
    #     return ModelDomain(
    #         Hgrid.open(hgrid, hgrid_crs),
    #         Vgrid.open(vgrid) if vgrid is not None else Vgrid(),
    #         Fgrid.open(fgrid, fgrid_crs))

    # def add_boundary_condition(self, forcing: BoundaryCondition, id=None):
    #     if id is None:
    #         for i in range(len(self.open_boundaries)):
    #             self.add_boundary_condition(forcing, i)
    #     else:
    #         if not isinstance(forcing, BoundaryCondition):
    #             raise TypeError("Argument must be of type "
    #                             f"{BoundaryCondition} but got type "
    #                             f"{type(forcing)}")
    #         self.open_boundaries[id]['forcing'] = forcing

    # def set_atmospheric_forcing(self, atmospheric_forcing: NWS):
    #     self._nws = atmospheric_forcing

    # def set_coriolis(self, ncor: Coriolis):
    #     self._ncor = ncor

    # def add_hydrology(self, hydrology: Hydrology):
    #     assert isinstance(hydrology, Hydrology), \
    #         f"Argument hydrology must be of type {Hydrology}, " \
    #         f"not type {type(hydrology)}."
    #     # if self.vgrid.is_2D:
    #     #     self.elev_ic = ElevIc()
    #     self._hydrology.append(hydrology)

# class ElecIc(Gr3Field):

#     def __init__(self, offset=-0.1):

#         mask = np.logical_and(
#                 self.values > 0.,
#                 nc['eta2'][:] < hgrid.values
#             )
#         idxs = np.where(mask)
#         nc['eta2'][idxs] = hgrid.values[idxs] + offset

    # @lru_cache(maxsize=1)
    # def get_active_potential_constituents(self):
    #     # PySCHISM allows the user to input the tidal potentials individually
    #     # for each boundary, however, SCHISM supports only a global
    #     # specification. Here, we collect all the activated tidal potentials
    #     # on each boundary and activate them all globally
    #     const = dict()
    #     for id in self.open_boundaries():
    #         forcing = self.open_boundaries[id]['forcing']
    #         if isinstance(forcing, Tides):
    #             for active in forcing.get_active_potential_constituents():
    #                 const[active] = True
    #     return tuple(const.keys())

    # @lru_cache(maxsize=1)
    # def get_active_forcing_constituents(self):
    #     # PySCHISM allows the user to input the tidal forcings individually
    #     # for each boundary, however, SCHISM supports only a global
    #     # specification. Here, we collect all the activated tidal forcings
    #     # on each boundary and activate them all globally
    #     const = dict()
    #     for id in self.open_boundaries():
    #         forcing = self.open_boundaries[id]['forcing']
    #         if isinstance(forcing, Tides):
    #             for active in forcing.get_active_forcing_constituents():
    #                 const[active] = True
    #     return tuple(const.keys())

    # def make_plot(self, **kwargs):
    #     if self.vgrid.is3D():
    #         raise NotImplementedError(
    #             "Plotting not yet supported for 3D meshes.")
    #     elif self.vgrid.is2D():
    #         self.hgrid.make_plot(**kwargs)

    # @property
    # def hgrid(self):
    #     return self._hgrid

    # @property
    # def vgrid(self):
    #     return self._vgrid

    # @property
    # def fgrid(self):
    #     return self._fgrid

    # @property
    # def nws(self):
    #     return self._nws

    # @property
    # def hydrology(self):
    #     return self._hydrology

    # @property
    # def ics(self):
    #     if self.hgrid.crs is None:
    #         return None
    #     elif self.hgrid.crs.is_geographic:
    #         return 2
    #     else:
    #         return 1

    # @property
    # def ncor(self):
    #     return self._ncor

    # @property
    # def open_boundaries(self):
    #     return self._open_boundaries

    # @property
    # def bctides(self):
    #     return self._bctides

    # @property
    # def _hgrid(self):
    #     return self.__hgrid

    # @_hgrid.setter
    # def _hgrid(self, hgrid: Hgrid):
    #     assert isinstance(hgrid, Hgrid), \
    #         f"Argument hgrid must be of type {Hgrid}, not type {type(hgrid)}."

    #     self.__hgrid = hgrid

    # @property
    # def _ncor(self):
    #     return self.__ncor

    # @_ncor.setter
    # def _ncor(self, ncor: Coriolis):
    #     if not isinstance(ncor, Coriolis):
    #         raise TypeError(f"ncor must be of type {Coriolis}, not type "
    #                         f"{type(ncor)}.")
    #     if ncor == Coriolis.AUTO:
    #         if self.ics == 1:
    #             self.sfea0 = np.median(self.hgrid.get_y("EPSG:4326"))
    #         elif self.ics == 2:
    #             pass  # nothing to do for ics=2
    #         else:
    #             raise ValueError(
    #                 f'Unknown hgrid.ics parameter {self.ics}')

    #     elif ncor == Coriolis.CORICOEFF:
    #         self.coricoef = 0.

    #     elif ncor == Coriolis.RLATITUDE:
    #         self.rlatitude = 46.

    #     else:
    #         raise NotImplementedError(
    #             f"Unknown value for Coriolis enum type {ncor}.")
    #     self.__ncor = ncor


# class CombineHotstartBinary:

#     def __init__(self, path: Union[str, os.PathLike], iteration=None):
#         path = pathlib.Path(path)
#         if iteration is None:
#             combine_hotstart = path.glob('hotstart_[0-9][0-9][0-9][0-9]_*.nc')
#             increments = set([file.name.split('_')[-1].split('.nc')[0]
#                               for file in combine_hotstart])
#             iteration = np.max([int(increment) for increment in increments])
#         subprocess.check_call(
#             ["combine_hotstart7", '-i', f'{iteration}'], cwd=path)
#         self.path = path / f"hotstart_it={iteration}.nc"

#     def add_elev_ic(self, hgrid, offset=-0.1):
#         nc = Dataset(self.path, 'r+')
#         mask = np.logical_and(
#             hgrid.values > 0.,
#             nc['eta2'][:] < hgrid.values
#         )
#         idxs = np.where(mask)
#         nc['eta2'][idxs] = hgrid.values[idxs] + offset
#         nc.close()
#


# class OpenBoundaries:

#     def __init__(self, hgrid: Hgrid):
#         open_boundaries = {}
#         for bnd in hgrid.boundaries.ocean().itertuples():
#             open_boundaries[bnd.id] = {
#                 'indexes': bnd.indexes, 'forcing': None}
#         self._hgrid = hgrid
#         self._open_boundaries = open_boundaries

#     def __call__(self):
#         return self._open_boundaries

#     def __len__(self):
#         return len(self._hgrid.boundaries.ocean())

#     def __getitem__(self, id):
#         return self._open_boundaries[id]

#     def __iter__(self):
#         for id, data in self._open_boundaries.items():
#             yield id, data

# class HgridDescriptor:

#     def __set__(self, obj, hgrid: Hgrid):
#         if not isinstance(hgrid, Hgrid):
#             raise TypeError(f'Argument hgrid must be of type {Hgrid}, '
#                             f'not {type(hgrid)}')
#         obj.__dict__['hgrid'] = hgrid

#     def __get__(self, obj, val):
#         return obj.__dict__['hgrid']


# class VgridDescriptor:

#     def __set__(self, obj, vgrid: Union[Vgrid, None]):
#         if not isinstance(vgrid, (Vgrid, type(None))):
#             raise TypeError(f'Argument vgrid must be of type {Vgrid} or None, '
#                             f'not {type(vgrid)}')
#         if vgrid is None:
#             vgrid = Vgrid()
#         obj.__dict__['vgrid'] = vgrid

#     def __get__(self, obj, val):
#         return obj.__dict__['vgrid']


# class FgridDescriptor:

#     def __set__(self, obj, fgrid: Fgrid):
#         if not isinstance(fgrid, Fgrid):
#             raise TypeError(f'Argument fgrid must be of type {Fgrid} or None, '
#                             f'not {type(fgrid)}')
#         obj.__dict__['fgrid'] = fgrid

#     def __get__(self, obj, val):
#         return obj.__dict__['fgrid']


# class OpenBoundariesDescriptor:

#     def __get__(self, obj, val):
#         open_boudaries = obj.__dict__.get('open_boudaries')
#         if open_boudaries is None:
#             open_boudaries = {}
#             for bnd in obj.hgrid.boundaries.ocean.itertuples():
#                 open_boudaries[bnd.id] = {
#                     'indexes': bnd.indexes, 'forcing': None}
#             obj.__dict__['open_boudaries'] = open_boudaries
#         return open_boudaries


# class NwsDescriptor:

#     def __set__(self, obj, nws: NWS):
#         if not isinstance(nws, NWS):
#             raise TypeError(f"Argument nws must be of type {NWS}, not "
#                             f"type {type(nws)}.")
#         obj.__dict__['nws'] = nws

#     def __get__(self, obj, val):
#         return obj.__dict__.get('nws')


# class NcorDescriptor:

#     def __set__(self, obj, ncor: Coriolis):
#         if not isinstance(ncor, Coriolis):
#             raise TypeError(f"ncor must be of type {Coriolis}, not type "
#                             f"{type(ncor)}.")
#         if ncor == Coriolis.AUTO:
#             if obj.hgrid.ics == 1:
#                 obj.sfea0 = np.median(obj.hgrid.get_y("EPSG:4326"))
#             elif obj.hgrid.ics == 2:
#                 pass  # nothing to do for ics=2
#             else:
#                 raise ValueError(
#                     f'Unknown hgrid.ics parameter {obj.hgrid.ics}')

#         elif ncor == Coriolis.CORICOEFF:
#             obj.coricoef = 0.

#         elif ncor == Coriolis.RLATITUDE:
#             obj.rlatitude = 46.

#         else:
#             raise NotImplementedError(
#                 f"Unknown value for Coriolis enum type {ncor}.")
#         obj.__dict__['ncor'] = ncor

#     def __get__(self, obj, val):
#         return obj.__dict__.get('ncor', Coriolis.AUTO)


# class Sfea0Descriptor:

#     def __set__(self, obj, sfea0: float):
#         obj.__dict__['sfea0'] = sfea0

#     def __get__(self, obj, val):
#         return obj.__dict__.get('sfea0')


# class CoricoeffDescriptor:

#     def __set__(self, obj, coricoeff: float):
#         obj.__dict__['coricoeff'] = coricoeff

#     def __get__(self, obj, val):
#         return obj.__dict__.get('coricoeff', 0.)


# class RlatitudeDescriptor:

#     def __set__(self, obj, rlatitude: float):
#         obj.__dict__['rlatitude'] = rlatitude

#     def __get__(self, obj, val):
#         return obj.__dict__.get('rlatitude', 46.)

    # _hgrid = HgridDescriptor()
    # _vgrid = VgridDescriptor()
    # _fgrid = FgridDescriptor()
    # _open_boundaries = OpenBoundariesDescriptor()
    # _nws = NwsDescriptor()
    # _ncor = NcorDescriptor()
    # _ics = IcsDescriptor()
    # sfea0 = Sfea0Descriptor()
    # coricoef = CoricoeffDescriptor()
    # rlatitude = RlatitudeDescriptor()
