import tempfile
import pathlib
from pyschism.forcing import (
    TidalForcing,
    WindForcing,
)
from pyschism.mesh import Mesh


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
        mode='barotropic',
        use_transport=False,
        netcdf=True,
    ):
        self._mesh = mesh
        self._tidal_forcing = tidal_forcing
        self._wind_forcing = wind_forcing
        self._wave_forcing = wave_forcing
        self._start_date = start_date
        self._end_date = end_date
        self._spinup_time = spinup_time
        self._mode = mode
        self._use_transport = use_transport
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

    def dump(
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
            self.mesh.hgrid.dump(outdir / hgrid, overwrite)
        if vgrid:
            self.mesh.vgrid.dump(outdir / vgrid, overwrite)
        if fgrid:
            self.mesh.fgrid.dump(outdir / fgrid, overwrite)
        if param:
            self.param.dump(outdir / param, overwrite)
        if bctides:
            self.bctides.dump(outdir / bctides, overwrite)

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
    def use_transport(self):
        return self._use_transport

    def _run_local(self, nproc, outdir, overwrite):
        self.dump(outdir, overwrite)

    def _run_coldstart(self, nproc, wdir):
        self._stage_files('coldstart', nproc, wdir)

    @property
    def _mesh(self):
        return self.__mesh

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
    def _use_transport(self):
        return self.__use_transport

    @_mesh.setter
    def _mesh(self, mesh):
        msg = f"mesh argument must be of type {Mesh}."
        assert isinstance(mesh, Mesh), msg
        self.__mesh = mesh

    @_tidal_forcing.setter
    def _tidal_forcing(self, tidal_forcing):
        if tidal_forcing is not None:
            assert isinstance(tidal_forcing, TidalForcing)
        self.__tidal_forcing = tidal_forcing

    @_wind_forcing.setter
    def _wind_forcing(self, wind_forcing):
        if wind_forcing is not None:
            assert isinstance(wind_forcing, WindForcing)
        self.__wind_forcing = wind_forcing

    @_mode.setter
    def _mode(self, mode):
        msg = "Run mode choices are 'barotropic' or 'baroclinic'."
        assert mode.lower() in ['barotropic', 'baroclinic'], msg
        self.__mode = mode

    @_use_transport.setter
    def _use_transport(self, use_transport):
        assert isinstance(use_transport, bool)
        self.__use_transport = use_transport
