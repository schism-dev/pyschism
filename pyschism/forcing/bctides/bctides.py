from datetime import datetime, timedelta
from functools import cached_property, lru_cache 
import pathlib
import logging

from pyschism.forcing.bctides import Tides

logger = logging.getLogger(__name__)
class Bctides:

    def __init__(
        self,
        hgrid,
        flags,
        constituents = 'major',
        database = 'tpxo',
        cutoff_depth: float = 50.0,
        ethconst = None,
        vthconst = None,
        tthconst = None,
        sthconst = None,
        tobc = None,
        sobc = None,
    ):

        self.hgrid = hgrid
        self.cutoff_depth = cutoff_depth
        self.tides = Tides(constituents=constituents, tidal_database=database)
        self.flags = flags
        self.ethconst = ethconst 
        self.vthconst = vthconst 
        self.tthconst = tthconst 
        self.sthconst = sthconst 
        self.tobc = tobc if tobc is not None else 1.0
        self.sobc = sobc if sobc is not None else 1.0

    def __str__(self):
 
        #first line in the bctides.in is a note, not used in the schism code
        f = [
            f"!{str(self.start_date)} UTC",
        ]

        #get earth tidal potential and frequency
        f.append(f"{self.ntip} {self.cutoff_depth} !number of earth tidal potential, cut-off depth for applying tidal potential")
        if self.ntip > 0:
            for constituent in self.tides.get_active_potential_constituents(): 
                forcing = self.tides(self.start_date, self.rnday, constituent)
                f.append(
                    " ".join([
                        f"{constituent}\n",
                        f"{forcing[0]:G}",
                        f"{forcing[1]:G}",
                        f"{forcing[2]:G}",
                        f"{forcing[3]:G}",
                        f"{forcing[4]:G}",
                    ])
                )

        #get tidal boundary 
        f.append(f"{self.nbfr:d} !nbfr")
        if self.nbfr > 0:
            for constituent in self.tides.get_active_forcing_constituents():
                forcing = self.tides(self.start_date, self.rnday, constituent)
                f.append(
                    " ". join([
                        f"{constituent}\n", 
                        f"{forcing[2]:G}",
                        f"{forcing[3]:G}",
                        f"{forcing[4]:G}",

                    ])
                )
 
        #get amplitude and phase for each open boundary         
        f.append(f"{len(self.gdf)} !nope")
        for ibnd, (boundary, flag) in enumerate(zip(self.gdf.itertuples(), self.flags)):
            logger.info(f"Processing boundary {ibnd}:")
            #number of nodes and flags
            line = [
                f"{len(boundary.indexes)}",
                *[str(digit) for digit in flag],
                f"!open bnd {ibnd+1}",
            ]
            f.append(" ".join(line))

            #It only accounts for elevation, velocity, temperature, salinity
            #TODO: add information for tracers
            if len(flag) > 4:
                raise NotImplementedError(f"Tracer module is not implemented yet!") 
            iettype, ifltype, itetype, isatype = [i for i in flag]

            #elevation boundary
            logger.info(f"Elevation type: {iettype}")
            if iettype == 1:
                logger.warning(f'time history of elevation is read in from elev.th (ASCII)!')
            elif iettype == 2:
                logger.info("This boundary is forced by a constant elevation!")
                if self.ethconst is None:
                    raise IOError("You are choosing type 2 for elevation, please specify a value!")
                else:
                    f.append("{str(self.ethconst)} !T")
            elif iettype == 4:
                logger.warning('time history of elevation is read in from elev2D.th.nc (netcdf)')
            elif iettype == 3 or iettype == 5:
                if iettype == 5:
                    logger.warning(f'Combination of 3 and 4, time history of elevation is read in from elev2D.th.nc!')
                for constituent in self.tides.get_active_forcing_constituents():
                    f.append(f"{constituent}")
                    vertices = self.hgrid.get_xy(crs=self.hgrid.crs)[boundary.indexes, :]
                    amp, phase = self.tides.get_elevation(constituent, vertices)
                    for i in range(len(boundary.indexes)):
                        f.append(f"{amp[i]: .6f} {phase[i]: .6f}")
            elif iettype == 0:
                logger.warning(f"elevations are not specified for this boundary (in this case the discharge must be specified)")
            else:
                raise IOError(f"Invalid type {iettype} for elevation!")

            #velocity
            logger.info(f"Velocity type: {ifltype}")
            if ifltype == 0:
                logger.info("Velocity is not sepcified, not input needed!")
            elif ifltype == 1:
                logger.warning("time history of discharge is read in from flux.th (ASCII)!")
            elif ifltype == 2:
                logger.info("This boundary is forced by a constant discharge!")
                if self.vthconst is None:
                    raise IOError("You are choosing type 2 for velocity, please specify dischage value!")
                else:
                    f.append("{str(self.vthconst)}")
            elif ifltype == 3 or ifltype == 5:
                if ifltype == 5:
                    logger.warning(f'Combination of 3 and 4, time history of velocity is read in from uv.3D.th.nc!')
                for constituent in self.tides.get_active_forcing_constituents():
                    f.append(f"{constituent}")
                    vertices = self.hgrid.get_xy(crs=self.hgrid.crs)[boundary.indexes, :]
                    uamp, uphase, vamp, vphase = self.tides.get_velocity(constituent, vertices)
                    for i in range(len(boundary.indexes)):
                        f.append(f"{uamp[i]: .6f} {uphase[i]: .6f} {vamp[i]: .6f} {vphase[i]: 6f}")
            elif ifltype == 4 or -4:
                logger.warning("time history of velocity (not discharge!) is read in from uv3D.th.nc (netcdf)")
                if ifltype == -4:
                    if self.rel1 is None or self.rel2 is None:
                        raise IOError(f"Please specify relaxation constants for inflow and outflow (between 0 and 1 with 1 being strongest nudging)")
                    else:
                        logger.info(f"relax value 1: {self.rel1}, relax value 2: {self.rel2}")
            elif ifltype == -1: 
                raise NotImplementedError(f"Velocity type {ifltype} not implemented yet!")
                #logger.info(f"Flather type radiation b.c. (iettype must be 0 in this case)!")
                #f.append(['vn_mean'])
                
                #TODO: add mean normal velocity at the node (at all levels)[
            else:
                raise IOError(f"Invalid type {ifltype} for veloicty!")

            #temperature
            logger.info(f"Temperature type: {itetype}")
            if itetype == 0:
                logger.warning("Temperature is not sepcified, not input needed!")
            elif itetype == 1:
                logger.warning("time history of temperature will be read in from TEM_1.th!")
                logger.info(f"Using nudging factor: {self.tobc}!")
                f.append(f"{self.tobc} !nudging factor for T")
            elif itetype == 2:
                logger.info("This boundary is forced by a constant temperature!")
                if self.tthconst is None:
                    raise IOError("You are choosing type 2 for temperature, please specify temperature value!")
                else:
                    f.append(f"{str(self.tthconst)} !T")

                logger.info(f"Using nudging factor: {self.tobc}!")
                f.append(f"{self.tobc} !nudging factor for T")
            elif itetype == 3:
                logger.info("Using initial temperature profile for inflow")
                logger.info(f"Using nudging factor: {self.tobc}!")
                f.append(f"{self.tobc} !nudging factor for T")
            elif itetype == 4:
                logger.warning("time history of temperature is read in from TEM_3D.th.nc (netcdf)!")
                logger.info(f"Using nudging factor: {self.tobc}!")
                f.append(f"{self.tobc} !nudging factor for T")
            else:
                raise IOError(f"Invalid type {itetype} for salinity!")        

            #salinity
            logger.info(f"Salinity type: {isatype}")
            if isatype == 0:
                logger.info("Salinity is not sepcified, not input needed!")
            elif isatype == 1:
                logger.warning("time history of salinity will be read in from SAL_1.th!")
                logger.info(f"Using nudging factor: {self.sobc}!")
                f.append(f"{self.sobc} !nudging factor for S")
            elif isatype == 2:
                logger.info("This boundary is forced by a constant salinity!")
                if self.sthconst is None:
                    raise IOError("You are choosing type 2 for salinity, please specify salinity value!")
                else:
                    f.append(f"{str(self.sthconst)} !S")

                logger.info(f"Using nudging factor: {self.sobc}!")
                f.append(f"{self.sobc} !nudging factor for S")
            elif isatype == 3:
                logger.info("Using initial salinity profile for inflow")
                logger.info(f"Using nudging factor: {self.sobc}!")
                f.append(f"{self.sobc} !nudging factor for S")
            elif isatype == 4:
                logger.warning("time history of salinity is read in from SAL_3D.th.nc (netcdf)!")
                logger.info(f"Using nudging factor: {self.sobc}!")
                f.append(f"{self.sobc} !nudging factor for S")
            else:
                raise IOError(f"Invalid type {isatype} for salinity!") 
                

        return "\n".join(f)


    def write(
        self,
        output_directory,
        start_date: datetime = None,
        rnday = None,
        constituents = ['K1', 'O1'],
        overwrite: bool = True, 
    ):

        if start_date is not None:
            self.start_date = start_date
        else:
            raise IOError("Please specify the start_date!")

        if rnday is not None:
            self.rnday = rnday
        else:
            raise IOError("Please specify the number of days!")

        if constituents is not None:
            self.constituents = constituents
        else:
            raise IOError("Please specify tidal constituents!")

        output_directory = pathlib.Path(output_directory)
        output_directory.mkdir(exist_ok=overwrite, parents=True) 
        bctides = output_directory / "bctides.in"
        if bctides.exists() and not overwrite:
            raise IOError("path exists and overwrite is False") 
        with open(bctides, "w") as f:
            f.write(str(self))

    @cached_property
    def gdf(self):
       return self.hgrid.boundaries.open.copy()
        
    @property
    def ntip(self):
        return len(self.tides.get_active_potential_constituents())

    @property
    def nbfr(self):
        return len(self.tides.get_active_potential_constituents())

