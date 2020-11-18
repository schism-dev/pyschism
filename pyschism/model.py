from datetime import datetime, timedelta
import os
import pathlib
import subprocess
from typing import Union

import numpy as np  # type: ignore[import]

from pyschism.domain import ModelDomain
from pyschism.bash import BashDriver
from pyschism.stations import Stations
from pyschism.param import Param
from pyschism.server.base import ServerConfig
from pyschism.enums import Stratification
from pyschism.forcing.tides.bctides import Bctides
from pyschism.forcing.atmosphere import NWS2
# from pyschism.forcing.atmosphere.nws.nws2.dataset import SfluxDownloadableDataset


# SIMPLISTIC_BASH_DRIVER = r"""#!/usr/bin/bash
# set -e
# touch outputs/mirror.out outputs/fatal.error
# tail -f outputs/mirror.out -f outputs/fatal.error &
# PID=$!
# mpiexec -n $(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}' ) pschism_TVD-VL
# kill -9 $PID
# """


class BctidesDescriptor:

    def __get__(self, obj, val):
        obj.__dict__['bctides'] = Bctides(obj.model_domain, obj.param)
        return obj.__dict__['bctides']


class ModelDomainDescriptor:

    def __set__(self, obj, model_domain: ModelDomain):
        if not isinstance(model_domain, ModelDomain):
            raise TypeError("Argument model_domain must be of type "
                            f"{ModelDomain}, not {type(model_domain)}.")
        obj.__dict__['model_domain'] = model_domain

    def __get__(self, obj, val):
        return obj.__dict__['model_domain']


class NwsDescriptor:

    def __get__(self, obj, val):
        nws = obj.__dict__.get('nws')
        if nws is None:
            nws = obj.model_domain.nws
            if nws is not None:
                nws(obj)
            obj.__dict__['nws'] = nws
        return nws


class StationsDescriptor:
    def __set__(self, obj, stations: Union[Stations, None]):
        if not isinstance(stations, (Stations, type(None))):
            raise TypeError('stations argument must be of type '
                            f'{Stations}')
        if stations is not None:
            stations.clip(
                obj.model_domain.hgrid.get_multipolygon(
                    stations.crs))
            if len(stations.get_active_vars()) > 0 and \
                    len(stations.stations) > 0:
                nspool_sta = stations.nspool_sta
                if isinstance(nspool_sta, timedelta):
                    nspool_sta = int(round(nspool_sta / obj.param.core.dt))
                obj.param.schout._nspool_sta = nspool_sta
            obj.__dict__['stations'] = stations

    def __get__(self, obj, val):
        return obj.__dict__.get('stations')


class ServerConfigDescriptor:

    def __init__(self):
        self.server_config = None

    def __set__(self, obj, server_config: Union[ServerConfig, None]):
        if not isinstance(server_config, (ServerConfig, type(None))):
            raise TypeError("Argument server_config must be of type "
                            f"{ServerConfig} or None, not type "
                            f"{type(server_config)}.")
        self.server_config = server_config

    def __get__(self, obj, val):
        return self.server_config


class DriverFileDescriptor:

    def __init__(self):
        self._driver_file = None

    def __get__(self, obj, val):
        if self._driver_file is None:
            self._driver_file = BashDriver(server_config=obj._server_config)

        return self._driver_file


class HotstartFile:

    def __set__(self, obj, hotstart_file: Union[str, os.PathLike]):
        obj.__dict__['hotstart_file'] = hotstart_file
        obj.param.opt.ihot = 1

    def __get__(self, obj, val):
        return obj.__dict__.get('hotstart_file')


class CombineHotstartDescriptor:

    def __set__(self, obj, outputs: Union[str, os.PathLike, None]):
        combine_hotstart = obj.__dict__.get('combine_hotstart')
        if outputs is not None:
            combine_hotstart = pathlib.Path(outputs).glob(
                'hotstart_[0-9][0-9][0-9][0-9]_*.nc')
            increments = set([file.name.split('_')[-1].split('.nc')[0]
                              for file in combine_hotstart])
            latest = np.max([int(increment) for increment in increments])
            subprocess.check_call(
                ["combine_hotstart7", '-i', f'{latest}'],
                cwd=outputs)
            obj._hotstart_file = pathlib.Path(outputs) / \
                f"hotstart_it={latest}.nc"
            obj.__dict__['combine_hotstart'] = combine_hotstart

    def __get__(self, obj, val):
        raise obj.__dict__.get('combine_hotstart')


class ModelDriver:

    _model_domain = ModelDomainDescriptor()
    _bctides = BctidesDescriptor()
    _nws = NwsDescriptor()
    _stations = StationsDescriptor()
    _server_config = ServerConfigDescriptor()
    _driver_file = DriverFileDescriptor()
    _combine_hotstart = CombineHotstartDescriptor()
    _hotstart_file = HotstartFile()

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
            **surface_outputs):
        """Main driver class used to generate SCHISM input files

        Arguments:
            mesh: :class:`pyschism.domain.ModelDomain`
            param: :class:`pyschism.param.param.Param`
        """

        self._model_domain = model_domain
        self._param = Param(model_domain, dt, rnday, dramp, start_date, ibc,
                            drampbc, nspool, ihfskip, nhot_write,
                            **surface_outputs)
        self._stations = stations
        self._server_config = server_config
        self._combine_hotstart = combine_hotstart

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
            driver_file=True,
            job_monitor_file=True,
    ):
        """Writes to disk the full set of input files necessary to run SCHISM.
        """
        outdir = pathlib.Path(output_directory)
        (outdir / 'outputs').mkdir(parents=True, exist_ok=True)
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
        if vgrid:
            vgrid = 'vgrid.in' if vgrid is True else vgrid
            self._model_domain.vgrid.write(outdir / vgrid, overwrite)
        if fgrid:
            fgrid = f'{self._model_domain.fgrid.fname}' if fgrid is True \
                    else fgrid
            self._model_domain.fgrid.write(outdir / fgrid, overwrite)
        if param:
            param = 'param.nml' if param is True else param
            self._param.write(outdir / param, overwrite)
        if bctides:
            bctides = 'bctides.in' if bctides is True else bctides
            self._bctides.write(outdir / bctides, overwrite)
        if nws:
            if self._nws is not None:
                if isinstance(self._nws, NWS2):
                    self._nws.write(outdir / 'sflux', overwrite,
                                    wind_rot=wind_rot)
                else:
                    self._nws.write(outdir, overwrite)
        if stations:
            if self._stations is not None:
                stations = 'station.in' if stations is True else stations
                self._stations.write(outdir / stations, overwrite)
        if self._hotstart_file is not None:
            hotstart_lnk = outdir / 'hotstart.nc'
            if overwrite is True:
                if hotstart_lnk.exists():
                    hotstart_lnk.unlink()
            os.symlink(os.path.relpath(
                self._hotstart_file, outdir), hotstart_lnk)
        if driver_file:
            driver_file = 'driver.sh' if driver_file is True else driver_file
            self._driver_file.write(outdir / driver_file, overwrite)
        if job_monitor_file:
            pass

    @property
    def param(self):
        return self._param

    @property
    def model_domain(self):
        return self._model_domain
