from datetime import datetime, timedelta
import os
import pathlib
import subprocess
from typing import Union

from netCDF4 import Dataset
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point, MultiPoint

from pyschism.domain import ModelDomain
from pyschism.driver.makefile import MakefileDriver
from pyschism.enums import Stratification
from pyschism.forcing.tides.bctides import Bctides
from pyschism.forcing.atmosphere import NWS2
from pyschism.forcing.hycom.hycom import HotStartInventory, OpenBoundaryInventory
from pyschism.param import Param
from pyschism.server import ServerConfig
from pyschism.stations import Stations
<<<<<<< HEAD
from pyschism.mesh.gridgr3 import Albedo, Diffmax, Diffmin, Watertype

=======
from pyschism.mesh.gridgr3 import Albedo,Diffmax,Diffmin,Watertype
from pyschism.mesh.prop import Fluxflag, Tvdflag
from pyschism.mesh.vgrid import Vgrid
from pyschism.dates import pivot_time, localize_datetime, nearest_cycle_date
>>>>>>> origin/ICOGS3D

class CombineHotstartBinary:

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

    def add_elev_ic(self, hgrid, offset=-0.1):
        nc = Dataset(self.path, 'r+')
        mask = np.logical_and(
                hgrid.values > 0.,
                nc['eta2'][:] < hgrid.values
            )
        idxs = np.where(mask)
        nc['eta2'][idxs] = hgrid.values[idxs] + offset
        nc.close()


class ModelDriver:

    def __init__(
            self,
            model_domain: ModelDomain,
            dt: Union[float, timedelta],
            rnday: Union[float, timedelta],
            ihfskip: int = None,
            dramp: Union[float, timedelta] = None,
            start_date: datetime = None,
            ibc: Union[Stratification, int, str] = Stratification.BAROTROPIC,
            drampbc: Union[float, timedelta] = None,
            stations: Stations = None,
            nspool: Union[int, timedelta] = None,
            nhot_write: Union[int, timedelta, bool] = None,
            server_config: ServerConfig = None,
            combine_hotstart: Union[str, os.PathLike] = None,
            cutoff_depth: float = 50.,
            **surface_outputs
    ):
        """Main driver class used to generate SCHISM input files

        This __init__ calls initialization routines for all SCHISM data inputs
        and prepares them to be dumped to files by the write() method.

        Arguments:
            mesh: :class:`pyschism.domain.ModelDomain`
            param: :class:`pyschism.param.param.Param`
        """
        # init param.nml
        self._param = Param(model_domain, dt, rnday, dramp, start_date, ibc,
                            drampbc, nspool, ihfskip, nhot_write, stations,
                            **surface_outputs)

        # init bctides.in
        self._bctides = Bctides(model_domain, self.param,
                                cutoff_depth=cutoff_depth)

        # init atmospheric data.
        self._nws = model_domain.nws

        # init hydrology data
        self._hydrology = model_domain.hydrology

        # init hotstart file
        self._combine_hotstart = combine_hotstart

        if model_domain.albedo is not None:
            self.param.opt.albedo = 1

        # do same here

        # init Makefile drivers
        self._makefile = MakefileDriver(server_config=server_config)

    def write(
            self,
            output_directory,
            overwrite=False,
            hgrid=True,
            vgrid=True,
            fgrid=True,
            param=True,
            bctides=True,
            nws=True,
            wind_rot=True,
            stations=True,
<<<<<<< Updated upstream
            use_param_template=True,
            albedo = True,
            diffmax = True,
            diffmin = True,
            watertype = True,
            fluxflag = True,
            tvdflag = True,
            rtofs = True,
=======
            use_param_template=False,
            albedo=True,
            diffmax=True,
            diffmin=True,
            watertype=True,
>>>>>>> Stashed changes
    ):
        """Writes to disk the full set of input files necessary to run SCHISM.
        """
        outdir = pathlib.Path(output_directory)

        if not (outdir / 'outputs').exists():
            (outdir / 'outputs').mkdir(parents=True)

        if hgrid:
            hgrid = 'hgrid.gr3' if hgrid is True else hgrid
            self.model_domain.hgrid.write(outdir / hgrid, overwrite)
            if self.model_domain.ics == 2:
                original_dir = os.getcwd()
                os.chdir(outdir)  # pushd
                try:
                    os.remove('hgrid.ll')
                except OSError:
                    pass
                os.symlink(hgrid, 'hgrid.ll')
                os.chdir(original_dir)  # popd

<<<<<<< Updated upstream
        if self.model_domain.albedo is not None:
            # self.model_domain.albedo = Albedo.constant(self.model_domain.hgrid, 0.15)
            self.model_domain.albedo.write(outdir / 'albedo.gr3', overwrite)

        if diffmax:
<<<<<<< HEAD
            self.diffmax = Diffmax.constant(self.model_domain.hgrid, 1.0e-6)
            self.diffmax.write(outdir / 'diffmax.gr3', overwrite)

        if diffmin:
            self.diffmin = Diffmax.constant(self.model_domain.hgrid, 1.0)
            self.diffmin.write(outdir / 'diffmin.gr3', overwrite)

        if watertype:
            self.watertype = Diffmax.constant(self.model_domain.hgrid, 1.0)
            self.watertype.write(outdir / 'watertype.gr3', overwrite)
=======
            self.diffmax = Diffmax.constant(self.model_domain.hgrid, 1.0)
            self.diffmax.write(outdir / 'diffmax.gr3',overwrite)

        if diffmin:
            self.diffmin = Diffmax.constant(self.model_domain.hgrid, 1.0e-6)
            self.diffmin.write(outdir / 'diffmin.gr3',overwrite)

        if watertype:
            self.watertype = Diffmax.constant(self.model_domain.hgrid, 1.0)
            self.watertype.write(outdir / 'watertype.gr3',overwrite)

        if fluxflag:
            self.fluxflag = Fluxflag.constant(self.model_domain.hgrid, -1)
            with open(outdir / 'fluxflag.prop', 'w+') as fid:
                fid.writelines(self.fluxflag)

        if tvdflag:
            #Hard-wire the polygon at this point. 
            coords = [(-75.340506, 40.843483), (-75.533474, 40.436019), (-75.796036, 39.535807), \
                (-75.672664, 39.339972), (-75.305709, 39.460000), (-75.103251, 39.636884), \
                (-74.692008, 39.744277), (-74.391485, 40.009603), (-74.359851, 40.252818), \
                (-74.514858, 40.745565), (-74.834362, 40.957194), (-75.210807, 40.935083),
                (-75.283565, 40.925607) ]
            poly = Polygon(coords)
            self.tvdflag = Tvdflag.define_by_region(hgrid=self.model_domain.hgrid, region=poly, value=1)
            with open(outdir / 'tvd.prop', 'w+') as fid:
                fid.writelines(self.tvdflag)

        if rtofs:
            self.start_date = nearest_cycle_date()
            self.hotstart = HotStartInventory()
            self.hotstart.fetch_data(outdir, self.model_domain.hgrid, self.start_date)
            self.obnd = OpenBoundaryInventory()
            self.obnd.fetch_data(outdir, self.start_date, rnday=3, \
                idx_min=2687, idx_max=2714, jdx_min=1181, jdx_max=1634)
>>>>>>> origin/ICOGS3D

        if vgrid:
            vgrid = 'vgrid.in' if vgrid is True else vgrid
            self.model_domain.vgrid = Vgrid.from_binary(self.model_domain.hgrid)
             
            #self.model_domain.vgrid.write(outdir / vgrid, overwrite)
#lcui
=======
        if vgrid:
            vgrid = 'vgrid.in' if vgrid is True else vgrid
            self.model_domain.vgrid.write(outdir / vgrid, overwrite)

>>>>>>> Stashed changes
        if fgrid:
            fgrid = f'{self.model_domain.fgrid.fname}' if fgrid is True \
                    else fgrid
            self.model_domain.fgrid.write(outdir / fgrid, overwrite)

        if param:
            param = 'param.nml' if param is True else param
            self.param.write(outdir / param, overwrite,
                             use_template=use_param_template)
        if bctides:
            bctides = 'bctides.in' if bctides is True else bctides
            self.bctides.write(outdir / bctides, overwrite)

        if nws:
            if self.nws is not None:
                if isinstance(self.nws, NWS2):
                    self.nws.write(outdir / 'sflux', overwrite,
                                   wind_rot=wind_rot)
                else:
                    self.nws.write(outdir, overwrite)
        if stations:
            if self.stations is not None:
                stations = 'station.in' if stations is True else stations
                self.stations.write(outdir / stations, overwrite)

        for hydrology in self._hydrology:
            hydrology.write(outdir, overwrite)

        if self.hotstart_file is not None:
            hotstart_lnk = outdir / 'hotstart.nc'
            if overwrite is True:
                if hotstart_lnk.exists():
                    hotstart_lnk.unlink()
            os.symlink(os.path.relpath(
                self.hotstart_file, outdir), hotstart_lnk)

        # TODO: Fix.
        if self.model_domain.albedo is not None:
            # self.model_domain.albedo = Albedo.constant(self.model_domain.hgrid, 0.15)
            self.model_domain.albedo.write(outdir / 'albedo.gr3', overwrite)

        if diffmax:
            self.diffmax = Diffmax.constant(self.model_domain.hgrid, 1.0e-6)
            self.diffmax.write(outdir / 'diffmax.gr3', overwrite)

        if diffmin:
            self.diffmin = Diffmax.constant(self.model_domain.hgrid, 1.0)
            self.diffmin.write(outdir / 'diffmin.gr3', overwrite)

        if watertype:
            self.watertype = Diffmax.constant(self.model_domain.hgrid, 1.0)
            self.watertype.write(outdir / 'watertype.gr3', overwrite)

        self.makefile.write(outdir / 'Makefile', overwrite)

    @property
    def param(self):
        return self._param

    @property
    def model_domain(self):
        return self.param.model_domain

    @property
    def stations(self):
        return self.param.stations

    @property
    def bctides(self):
        return self._bctides

    @property
    def nws(self):
        return self._nws

    @property
    def hydrology(self):
        return self._hydrology

    @property
    def hotstart_file(self):
        if self._combine_hotstart is not None:
            return self._combine_hotstart.path

    @property
    def makefile(self):
        return self._makefile

    @property
    def _nws(self):
        return self.__nws

    @_nws.setter
    def _nws(self, nws):
        if nws is not None:
            nws(self)
        self.__nws = nws

    @property
    def _combine_hotstart(self):
        return self.__combine_hotstart

    @_combine_hotstart.setter
    def _combine_hotstart(self, combine_hotstart):
        if combine_hotstart is not None:
            combine_hotstart = CombineHotstartBinary(combine_hotstart)
            self.param.opt.ihot = 1
        self.__combine_hotstart = combine_hotstart

    @property
    def _hydrology(self):
        return self.__hydrology

    @_hydrology.setter
    def _hydrology(self, hydrology):
        for forcing in hydrology:
            if callable(forcing):
                forcing(self)
        self.__hydrology = hydrology
