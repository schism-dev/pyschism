#! /usr/bin/env python
import pathlib
from datetime import datetime, timezone, timedelta
from pyschism import Param, Mesh, forcing, driver


def main():
    # setup mesh
    mesh = Mesh.open(
        pathlib.Path('hgrid.gr3'),  # hgrid
        vgrid=None,  # optional parameter
        crs="EPSG:4326"
        )

    # set mesh friction
    mesh.set_friction('manning', 0.025)

    # set boundary conditions
    elevbc = forcing.Tides()
    elevbc.use_all()
    mesh.set_boundary_forcing(elevbc)

    # setup model config
    param = Param()

    # set model run timings
    param.start_date = datetime(
        2017, 9, 18, 12, 00,
        tzinfo=timezone(timedelta(hours=-4)))
    param.run_time = timedelta(days=5)
    param.spinup_time = timedelta(days=7)

    # set model output options
    param.nspool = timedelta(minutes=30.)
    param.elev = True
    param.dahv = True

    # init the model driver
    model = driver.SchismRun(mesh, param)

    # execute model
    model.run(nproc=-1, workdir='staging')

    print(model.outputs)


if __name__ == '__main__':
    main()
