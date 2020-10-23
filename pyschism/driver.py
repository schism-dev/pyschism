import os
import pathlib
from datetime import timedelta


from .mesh import Mesh
from .param import Param
from .forcing.tides.bctides import Bctides

SIMPLISTIC_BASH_DRIVER = r"""#!/usr/bin/bash
set -e
touch outputs/mirror.out outputs/fatal.error
tail -f outputs/mirror.out -f outputs/fatal.error &
PID=$!
mpiexec -n $(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}' ) pschism_TVD-VL
kill -9 $PID
"""


class SchismRun:

    def __init__(self, mesh: Mesh, param: Param):
        """Main driver class used to generate SCHISM input files

        Arguments:
            mesh: :class:`pyschism.Mesh`
            param: :class:`pyschism.Param`
        """

        # assert Mesh
        assert isinstance(mesh, Mesh), f"mesh argument must be of type {Mesh}."
        self._mesh = mesh

        # assert Param
        assert isinstance(param, Param), f'param is not an istance of {Param}'
        self._param = param

        # ------------------------------------------------
        # Set initial Mesh-dependent parameters for core |
        # ------------------------------------------------
        # TODO: Must set (when applicable)
        # msc2
        # mdc2
        # ntracer_gen
        # ntracer_age
        # sed_class
        # eco_class

        # --------------------------------
        # Set initial parameters for opt |
        # --------------------------------
        # set friction
        self._param.opt.nchi = self._mesh.fgrid.nchi
        if self._param.opt.nchi == -1:
            self._param.opt.hmin_man = self._mesh.fgrid.hmin_man
        if self._param.opt.nchi == 1:
            self._param.opt.dbz_min = self._mesh.fgrid.dbz_min
            self._param.opt.dbz_decay = self._mesh.fgrid.dbz_decay

        # set coordinate system
        self._param.opt.ics = self._mesh.ics

        if self._param.opt.ncor is None:
            # set coriolis parameters (if not explicitly set by user)
            self._param.opt.set_ncor(1, sfea0=self._mesh.sfea0)

        # ------------------------------------
        # Set initial parameters for bctides |
        # ------------------------------------
        self._bctides = Bctides(self._mesh, self._param)

        # -----------------------------------
        # Set parameters for station outputs|
        # -----------------------------------
        if self._param.stations is not None:
            self._param.stations.clip(self._mesh.hgrid.get_multipolygon(
                self._param.stations.crs))
            if len(self._param.stations.get_active_vars()) > 0 and \
                    len(self._param.stations.stations) > 0:
                self._param.schout.iout_sta = 1
                nspool_sta = self._param.stations.nspool_sta
                if isinstance(nspool_sta, timedelta):
                    nspool_sta = int(round(nspool_sta / self._param.core.dt))
                self._param.schout.nspool_sta = nspool_sta

    def write(
        self,
        output_directory,
        overwrite=False,
        hgrid=True,
        vgrid=True,
        fgrid=True,
        param=True,
        bctides=True,
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
            self._mesh.hgrid.write(outdir / hgrid, overwrite)
            if self._mesh.ics == 2:
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
            self._mesh.vgrid.write(outdir / vgrid, overwrite)
        if fgrid:
            fgrid = f'{self._mesh.fgrid.fname}' if fgrid is True else fgrid
            self._mesh.fgrid.write(outdir / fgrid, overwrite)
        if param:
            param = 'param.nml' if param is True else param
            self._param.write(outdir / param, overwrite)
        if bctides:
            bctides = 'bctides.in' if bctides is True else bctides
            self._bctides.write(outdir / bctides, overwrite)
        if stations:
            stations = 'station.in' if stations is True else stations
            self._param.stations.write(outdir / stations, overwrite)
        if driver_file:
            driver_file = 'driver.sh' if driver_file is True else driver_file
            with open(outdir / driver_file, 'w') as f:
                f.write(SIMPLISTIC_BASH_DRIVER)
        if job_monitor_file:
            pass
