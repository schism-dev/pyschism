from time import time
import os
import argparse
from datetime import datetime, timedelta
import logging
import json

import numpy as np

from pyschism.mesh import Hgrid
from pyschism.forcing.bctides import Bctides

logging.basicConfig(
    format = "[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.getLogger("pyschism").setLevel(logging.DEBUG)

def list_of_strings(arg):
    return arg.split(',')

if __name__ == "__main__":

    '''
    Assume files are already located in:
    database='fes2014'
        ~/.local/share/fes2014/eastward_velocity/
        ~/.local/share/fes2014/northward_velocity/  
        ~/.local/share/fes2014/ocean_tide_extrapolated/
    database = 'tpxo'
        ~/.local/share/tpxo/
    '''
    parser = argparse.ArgumentParser(description="Create bctides.in for SCHISM with command-line arguments! e.g. python test_bctides.py hgrid.ll 2014-12-01 397 '[[5,5,4,4],[5,5,4,4],[0,1,1,2]]' major fes2014")
    
    #Add arguments
    parser.add_argument('hgrid', type=str, help='hgrid.ll (lon/lat) file')
    parser.add_argument('start_date', type=datetime.fromisoformat, help='model startdate')
    parser.add_argument('rnday', type=float, help='model rnday')
    parser.add_argument('bctypes', type=str, help="JSON format for Flags for each open boundary, '[[5,5,4,4],[5,5,4,4],[0,1,1,2]]'")
    parser.add_argument('constituents', type=list_of_strings, help="Choose tidal constituents to be included, major, minor, or list of constituents ['K1', 'O1', 'M2']")
    parser.add_argument('database', type=str, help='Tidal datbase: tpxo or fes2014')

    #Parse the command-line arguments
    args = parser.parse_args()
    hgrid_filename = args.hgrid
    start_date = args.start_date
    rnday = args.rnday
    bctypes = args.bctypes
    constituents = args.constituents
    database = args.database

    # Parse the JSON string into a Python data structure
    try:
        flags = json.loads(bctypes)
        print("Parsed bctype list:", flags)
    except json.JSONDecodeError:
        raise TypeError("Invalid JSON format for bctype list.")

    #earth tidal potential
    add_earth_tidal_potential = input("Would you like to add earth tidal potential? Y/N: ")
    if add_earth_tidal_potential.lower() == "y":
        earth_tidal_potential = True
    else:
        earth_tidal_potential = False

    #Check if constant values needed
    ethconst = []
    vthconst = []
    tthconst = []
    sthconst = []
    tobc = []
    sobc = []
    relax = []

    for ibnd, flag in enumerate(flags):
        iettype, ifltype, itetype, isatype = [i for i in flag]
        if iettype == 2:
            val = input(f"Elevation value at boundary {ibnd+1}: ")
            ethconst.append(float(val))
        else:
            ethconst.append(np.nan)

        if ifltype == 2:
            val = input(f"Discharge value (negative for inflow) at boundary {ibnd+1}: ")
            vthconst.append(float(val))
        elif ifltype == -4:
            val = input(f"Relaxation constants (between 0 and 1) for inflow at boundary {ibnd+1}: ")
            relax.append(float(val))
            val = input(f"Relaxation constants (between 0 and 1) for outflow at boundary {ibnd+1}: ")
            relax.append(float(val))
            
        else:
            vthconst.append(np.nan)

        if itetype == 2:
            val = input(f"Temperature value at boundary {ibnd+1}: ")
            tthconst.append(float(val))
            val = input(f"Nuding factor (between 0 and 1) for temperature at boundary {ibnd+1}: ")
            tobc.append(float(val))
        elif itetype == 1 or itetype == 3 or itetype == 4:
            tthconst.append(np.nan)
            val = input(f"Nuding factor (between 0 and 1) for temperature at boundary {ibnd+1}: ")
            tobc.append(float(val))
        else:
            tthconst.append(np.nan)
            tobc.append(np.nan)

        if isatype == 2:
            val = input(f"Salinity value at boundary {ibnd+1}: ")
            sthconst.append(float(val))
            val = input(f"Nuding factor (between 0 and 1) for salinity at boundary {ibnd+1}: ")
            sobc.append(float(val))
        elif isatype == 1 or isatype == 3 or isatype == 4:
            sthconst.append(np.nan)
            val = input(f"Nuding factor (between 0 and 1) for salinity at boundary {ibnd+1}: ")
            sobc.append(float(val))
        else:
            sthconst.append(np.nan)
            sobc.append(np.nan)

    outdir = './'
    hgrid = Hgrid.open(hgrid_filename, crs="epsg:4326")

    bctides=Bctides(
        hgrid = hgrid, 
        flags = flags,
        constituents = constituents, 
        database = database,
        add_earth_tidal = earth_tidal_potential,
        ethconst = ethconst,
        vthconst = vthconst,
        tthconst = tthconst,
        sthconst = sthconst,
        tobc = tobc,
        sobc = sobc,
        relax = relax,
    )

    bctides.write(
        outdir, 
        start_date=start_date, 
        rnday=rnday, 
        overwrite=True,
    )
