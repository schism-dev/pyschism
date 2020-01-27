import tempfile
from datetime import datetime, timedelta
import pathlib
from pyschism.mesh import Mesh
from pyschism.driver.param import Param
from pyschism.forcing import (
    TidalForcing,
    WindForcing,
)


class SchismRun:

    def __init__(
        self,
        mesh,
        start_date,
        end_date,
        spinup_time,
        tidal_forcing=None,
        wind_forcing=None,
        wave_forcing=None,
        tracers=None,
        mode='barotropic',
        transport=False,
        netcdf=True,
    ):
        self._mesh = mesh
        self._start_date = start_date
        self._end_date = end_date
        self._spinup_time = spinup_time
        self._tidal_forcing = tidal_forcing
        self._wind_forcing = wind_forcing
        self._wave_forcing = wave_forcing
        self._tracers = tracers
        self._mode = mode
        self._transport = transport
        self._netcdf = netcdf

    def run(
        self,
        outdir=None,
        nproc=-1,
        overwrite=False,
        server_config=None,
    ):

        if outdir is None:
            self._outdir_tmpdir = tempfile.TemporaryDirectory()
            outdir = pathlib.Path(self._outdir_tmpdir.name)
        else:
            outdir = pathlib.Path(outdir)
            if outdir.exists() and not overwrite:
                msg = f"{outdir} exists and overwrite is not enabled."
                raise IOError(msg)

        if not outdir.exists():
            outdir.mkdir(parents=True)

        # local run
        if server_config is None:
            self._run_local(
                nproc=nproc,
                outdir=outdir,
                overwrite=overwrite
            )

        # server run
        else:
            server_config.run(
                driver=self,
                outdir=outdir,
                overwrite=overwrite,
            )

        self._load_outdir(outdir)

        return self._output_collection

    def write(
        self,
        output_directory,
        overwrite=False,
        hgrid='hgrid.gr3',
        vgrid='vgrid.in',
        fgrid='fgrid.gr3',
        param='param.in',
        bctides='bctides.in',
    ):
        outdir = pathlib.Path(output_directory).absolute()
        if not outdir.exists():
            outdir.mkdir(parents=True)
        if hgrid:
            self.mesh.hgrid.write(outdir / hgrid, overwrite)
        if vgrid:
            self.mesh.vgrid.write(outdir / vgrid, overwrite)
        if fgrid:
            self.mesh.fgrid.write(outdir / fgrid, overwrite)
        if param:
            self.param.write(outdir / param, overwrite)
        if bctides:
            self.bctides.write(outdir / bctides, overwrite)

    @property
    def mesh(self):
        return self._mesh

    @property
    def tidal_forcing(self):
        return self._tidal_forcing

    @property
    def wind_forcing(self):
        return self._wind_forcing

    @property
    def mode(self):
        return self._mode

    @property
    def transport(self):
        return self._transport

    @property
    def param(self):
        return Param(
            self.start_date,
            self.end_date,
            self.spinup_time,
            )

    @property
    def start_date(self):
        return self._start_date

    @property
    def end_date(self):
        return self._end_date

    @property
    def spinup_time(self):
        return self._spinup_time

    def _run_local(self, nproc, outdir, overwrite):
        self.write(outdir, overwrite)

    def _run_coldstart(self, nproc, wdir):
        self._stage_files('coldstart', nproc, wdir)

    @property
    def _mesh(self):
        return self.__mesh

    @property
    def _start_date(self):
        return self.__start_date

    @property
    def _end_date(self):
        return self.__end_date

    @property
    def _spinup_time(self):
        return self.__spinup_time

    @property
    def _tidal_forcing(self):
        return self.__tidal_forcing

    @property
    def _wind_forcing(self):
        return self.__wind_forcing

    @property
    def _mode(self):
        return self.__mode

    @property
    def _transport(self):
        return self.__transport

    @_mesh.setter
    def _mesh(self, mesh):
        msg = f"mesh argument must be of type {Mesh}."
        assert isinstance(mesh, Mesh), msg
        self.__mesh = mesh

    @_start_date.setter
    def _start_date(self, start_date):
        msg = f"start_date must be a {datetime} instance."
        assert isinstance(start_date, datetime), msg
        self.__start_date = start_date

    @_end_date.setter
    def _end_date(self, end_date):
        msg = f"end_date must be a {datetime} instance."
        assert isinstance(end_date, datetime), msg
        self.__end_date = end_date

    @_spinup_time.setter
    def _spinup_time(self, spinup_time):
        msg = f"spinup_time must be a {timedelta} instance."
        assert isinstance(spinup_time, timedelta), msg
        self.__spinup_time = spinup_time

    @_tidal_forcing.setter
    def _tidal_forcing(self, tidal_forcing):
        if tidal_forcing is not None:
            assert isinstance(tidal_forcing, TidalForcing)
            tidal_forcing.start_date = self.start_date
            tidal_forcing.end_date = self.end_date
            tidal_forcing.spinup_time = self.spinup_time
        self.__tidal_forcing = tidal_forcing

    @_wind_forcing.setter
    def _wind_forcing(self, wind_forcing):
        if wind_forcing is not None:
            assert isinstance(wind_forcing, WindForcing)
            wind_forcing.start_date = self.start_date
            wind_forcing.end_date = self.end_date
            wind_forcing.spinup_time = self.spinup_time
        self.__wind_forcing = wind_forcing

    @_mode.setter
    def _mode(self, mode):
        msg = "Run mode choices are 'barotropic' or 'baroclinic'."
        assert mode.lower() in ['barotropic', 'baroclinic'], msg
        self.__mode = mode

    @_transport.setter
    def _transport(self, transport):
        assert isinstance(transport, bool)
        self.__transport = transport
