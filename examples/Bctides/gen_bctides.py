from time import time
import os
import argparse
from datetime import datetime, timedelta
import logging
import pathlib
import json

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
    parser = argparse.ArgumentParser(description="Create bctides.in for SCHISM with command-line arguments! e.g. python test_bctides.py hgrid.gr3 2014-12-01 397 '[[5,5,4,4],[5,5,4,4],[0,1,1,2]]' major fes2014")
    
    #Add arguments
    parser.add_argument('hgrid', type=str, help='hgrid (lon/lat) file')
    parser.add_argument('start_date', type=datetime.fromisoformat, help='model startdate')
    parser.add_argument('rnday', type=float, help='model rnday')
    parser.add_argument('bctypes', type=str, help="JSON format for Flags for each open boundary, '[5,5,4,4],[5,5,4,4],[0,1,1,2]'")
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

    #start_date = datetime(2014, 12, 1)
    #rnday = 397
    outdir = './'
    hgrid = Hgrid.open(hgrid_filename, crs="epsg:4326")

    bctides=Bctides(
        hgrid = hgrid, 
        #flags = [[5, 5, 4, 4], [5, 5, 4, 4], [0, 1, 2, 2]],
        flags = flags,
        constituents = constituents, 
        database = database,
        tthconst = 10.0,
        sthconst = 0.0,
    )

    bctides.write(
        outdir, 
        start_date=start_date, 
        rnday=rnday, 
        overwrite=True,
    )
