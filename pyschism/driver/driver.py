import tempfile
import pathlib
from functools import lru_cache
from pyschism.mesh import Mesh
from pyschism.driver.param import Param
from pyschism.forcing.bctides import Bctides


class SchismRun:

    def __init__(
        self,
        mesh,
        param,
    ):
        self._mesh = mesh
        self._param = param

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
    def param(self):
        return self._param

    @property
    def mesh(self):
        return self._mesh

    @property
    def start_date(self):
        return self.param.start_date

    @property
    def end_date(self):
        return self.param.end_date

    @property
    def spinup_time(self):
        return self.param.spinup_time

    @property
    @lru_cache(maxsize=None)
    def bctides(self):
        return Bctides(self.mesh, self.start_date, self.end_date, self.spinup_time)

    def _run_local(self, nproc, outdir, overwrite):
        self.write(outdir, overwrite)

    def _run_coldstart(self, nproc, wdir):
        self._stage_files('coldstart', nproc, wdir)

    @property
    def _mesh(self):
        return self.__mesh

    @property
    def _param(self):
        return self.__param

    @_mesh.setter
    def _mesh(self, mesh):
        msg = f"mesh argument must be of type {Mesh}."
        assert isinstance(mesh, Mesh), msg
        self.__mesh = mesh

    @_param.setter
    def _param(self, param):
        msg = f"param argument must be of type {Param}."
        assert isinstance(param, Param), msg
        # set friction
        param.opt.nchi = self.mesh.fgrid.nchi
        if param.opt.nchi == 1:
            param.opt.hmin_man = self.mesh.fgrid.hmin_man
        if param.opt.nchi == 1:
            param.opt.dbz_min = self.mesh.fgrid.dbz_min
            param.opt.dbz_decay = self.mesh.fgrid.dbz_decay
        # set coordinate system
        if self.mesh.crs.is_geographic:
            param.opt.ics = self.mesh.ics
            param.opt.slam0 = self.mesh.slam0
            param.opt.sfea0 = self.mesh.sfea0
        else:
            param.opt.ics = self.mesh.ics
        # set coriolis
        if self.mesh.ics == 2:
            param.opt.ncor = 1
        else:
            param.opt.ncor = 1
        self.__param = param

    # def run_local(self):

    #     def config():
    #         src = pathlib.Path(self.args.config_file).resolve()
    #         dst = self.output_directory / src.name
    #         if dst.is_file() and not self.overwrite:
    #             msg = "Destination file exists and overwrite=False"
    #             raise Exception(msg)
    #         shutil.copy2(src, self.output_directory)

    #     def param():
    #         self.param.write(self.staging_directory)

    #     def mesh():
    #         self.mesh.write(
    #             hgrid=self.output_path / 'hgrid.gr3',
    #             vgrid=self.output_path / 'vgrid.in',
    #             fgrid=self.output_path / 'fgrid.gr3',
    #             overwrite=self.overwrite)

    #     def bctides():
    #         pass

    #     def pschism_exec(self):
    #         raise NotImplementedError

    #     config()
    #     param()
    #     mesh()
    #     bctides()
    #     pschism_exec()
    #     return 0

    # def run_server(self):
    #     raise NotImplementedError
