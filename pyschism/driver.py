import os
import pathlib
from datetime import timedelta


from .mesh import Mesh
from .param import Param
from .forcing.bctides import Bctides

SIMPLISTIC_BASH_DRIVER = """#!/usr/bin/bash
set -e
touch outputs/mirror.out outputs/fatal.error
tail -f outputs/mirror.out -f outputs/fatal.error &
mpiexec -n $(($(nproc)/2)) pschism_TVD-VL


"""


class SchismRun:

    def __init__(self, mesh: Mesh, param: Param):

        # assert Mesh
        assert isinstance(mesh, Mesh), f"mesh argument must be of type {Mesh}."
        self.__mesh = mesh

        # assert Param
        assert isinstance(param, Param), f'param is not an istance of {Param}'
        self.__param = param

        # --------------------------------
        # Set initial parameters for opt |
        # --------------------------------
        # set friction
        self.param.opt.nchi = self.mesh.fgrid.nchi
        if self.param.opt.nchi == -1:
            self.param.opt.hmin_man = self.mesh.fgrid.hmin_man
        if self.param.opt.nchi == 1:
            self.param.opt.dbz_min = self.mesh.fgrid.dbz_min
            self.param.opt.dbz_decay = self.mesh.fgrid.dbz_decay

        # set coordinate system
        self.param.opt.ics = self.mesh.ics

        if self.param.opt.ncor is None:
            # set coriolis parameters (if not explicitly set by user)
            self.param.opt.set_ncor(1, sfea0=self.mesh.sfea0)

        # ------------------------------------
        # Set initial parameters for bctides |
        # ------------------------------------
        # check if start_date was given in case tidal forcings are requested.
        afc = self.mesh.get_active_forcing_constituents()
        if len(afc) > 0 and self.param.opt.start_date is None:
            raise Exception('start_date argument is required for simulating '
                            'tidal forcing.')
        start_date = self.param.opt.start_date
        end_date = start_date + timedelta(days=self.param.core.rnday)
        self.__bctides = Bctides(
            self.mesh,
            start_date=start_date,
            end_date=end_date
        )

    def write(
        self,
        output_directory,
        overwrite=False,
        hgrid=True,
        vgrid=True,
        fgrid=True,
        param=True,
        bctides=True,
        driver_file=True,
        job_monitor_file=True
    ):
        outdir = pathlib.Path(output_directory)
        (outdir / 'outputs').mkdir(parents=True, exist_ok=True)
        if hgrid:
            hgrid = 'hgrid.gr3' if hgrid is True else hgrid
            self.mesh.hgrid.write(outdir / hgrid, overwrite)
            if self.mesh.ics == 2:
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
            self.mesh.vgrid.write(outdir / vgrid, overwrite)
        if fgrid:
            fgrid = f'{self.mesh.fgrid.fname}' if fgrid is True else fgrid
            self.mesh.fgrid.write(outdir / fgrid, overwrite)
        if param:
            param = 'param.nml' if param is True else param
            self.param.write(outdir / param, overwrite)
        if bctides:
            bctides = 'bctides.in' if bctides is True else bctides
            self.bctides.write(outdir / bctides, overwrite)
        if driver_file:
            driver_file = 'driver.sh' if driver_file is True else driver_file
            with open(outdir / driver_file, 'w') as f:
                f.write(SIMPLISTIC_BASH_DRIVER)
        if job_monitor_file:
            pass

    @property
    def param(self):
        return self.__param

    @property
    def mesh(self):
        return self.__mesh

    @property
    def bctides(self):
        return self.__bctides

    @property
    def end_date(self):
        return self.__end_date
        return self.start_date + timedelta(days=self.core.rnday)
