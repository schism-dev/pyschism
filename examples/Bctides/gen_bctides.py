from time import time
import os
from datetime import datetime, timedelta
import logging
import pathlib
import f90nml

from pyschism.mesh import Hgrid
from pyschism.forcing.bctides import Bctides, iettype, ifltype, isatype, itetype
from pyschism.forcing.source_sink.nwm import NationalWaterModel, NWMElementPairings
from pyschism.forcing.nws import NWS2, GFS, HRRR

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
    #setup logging
    logging.basicConfig(
        format = "[%(asctime)s] %(name)s %(levelname)s: %(message)s",
        force=True,
    )
    logging.getLogger("pyschism").setLevel(logging.DEBUG)

    start_date = datetime(2017, 12, 1)
    rnday = 10
    outdir = './'
    hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")

    #elevation
    #constituents: major [M2, S2, N2, K2, K1, O1, P1, Q1]
    iet3 = iettype.Iettype3(constituents='major', database='fes2014')
    iet4 = iettype.Iettype4()
    iet5 = iettype.Iettype5(iettype3=iet3, iettype4=iet4)

    #velocity
    ifl3 = ifltype.Ifltype3(constituents='major', database='fes2014')
    ifl4 = ifltype.Ifltype4()
    ifl5 = ifltype.Ifltype5(ifltype3=ifl3, ifltype4=ifl4)

    #salinity 
    isa4 = isatype.Isatype4() #time-varying salinity
    isa2 = isatype.Isatype2(0, 1) # constant salinity (=0)

    #temperature
    ite4 = itetype.Itetype4() #time-varying temperature
    ite2 = itetype.Itetype2(0, 1) #constant temperature (=0)

    #example1 - only one ocean open boundary and type is [5, 5, 4, 4]
    #bctides = Bctides(hgrid, iettype=iet5, ifltype=ifl5, isatype=isa4, itetype=ite4)

    #example2 - two ocean open boundaries, one river boundary and type is [[5, 5, 4, 4], [5, 5, 4, 4], [0, 0, 3, 3]]
    bctides=Bctides(
        hgrid, 
        iettype={'1': iet5, '2': iet5}, 
        ifltype={'1': ifl5, '2': ifl5}, 
        isatype={'1': isa4, '2':isa4, '3': isa2}, 
        itetype={'1': ite4, '2':ite4, '3': ite2}
    )

    bctides.write(
        outdir, 
        start_date=start_date, 
        end_date=rnday, 
        bctides=True, 
        overwrite=True,
    )
