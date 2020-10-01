#! /usr/bin/env python
import pathlib
from datetime import datetime, timezone, timedelta
from pyschism import Param, Mesh, forcing, SchismRun


def main():

    # setup mesh
    mesh = Mesh.open('hgrid.gr3')

    # set mesh friction
    mesh.set_friction('manning', 0.025)

    # set tidal boundary conditions
    elevbc = forcing.Tides()
    # activate all forcing constituents
    elevbc.use_all()

    # connect the boundary condition to the mesh
    mesh.set_boundary_forcing(elevbc)

    #  ------ Param options
    # rnday and nspool are required arguments

    # Use a float or int to specify run days or use a timedelta for a
    # pythonic specification.
    rnday = timedelta(seconds=864000)

    # Use an integer for number of steps or a timedelta to approximate
    # number of steps internally based on timestep
    nspool = timedelta(minutes=30.)

    # The start date is an optional parameter, except when using forcings that
    # require some natural dates, (e.g. tides, winds, T&S ...)
    # If using a forcing that requires a start date is invoked, and no start
    # date is provided, the driver will throw an ecxeption.
    # tzinfo is optional, UTC assumed if not provided
    start_date = datetime(2017, 9, 18, 12, 00,
                          tzinfo=timezone(timedelta(hours=-4)))

    # Output requests can be activated as part of the optional arguments.
    # dt is another optional argument (timestep). dt is omitted in this example
    # and thus takes the default value of 150. seconds internally
    param = Param(rnday, nspool, start_date=start_date, elev=True)

    # Output requests, as well as other namelist variables, can be modified
    # post-hoc through their corresponding properties.
    # For example, the following line activates the depth averaged horizontal
    # velocity output request:
    param.schout.dahv = True

    # init the model driver
    driver = SchismRun(mesh, param)

    # write files to disk
    driver.write('staging', overwrite=True)


if __name__ == '__main__':
    main()
