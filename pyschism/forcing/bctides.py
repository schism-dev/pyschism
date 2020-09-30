from abc import ABC, abstractmethod
from datetime import datetime, timedelta
from enum import Enum
from functools import lru_cache

from . import bctypes
from .tides import Tides
from ..mesh.hgrid import Hgrid
from ..mesh.mesh import Mesh


class NullWritter:

    def __init__(self, *argv):
        pass

    def __str__(self):
        return ""


class NotImplementedWritter:
    
    def __init__(self, *argv):
        raise NotImplementedError('NotImplementedWritter')

class TidalVariableWritter(ABC):
    
    def __init__(self, boundary, bctides):
        self.forcing = boundary.forcing
        self.indexes = boundary.indexes
        self.active_constituents = bctides.active_forcing_constituents
        self.mesh = bctides.mesh

    @abstractmethod
    def __str__(self):
        raise NotImplementedError

class TidalElevationWritter(TidalVariableWritter):

    def __str__(self):
        f = ""
        for constituent in self.active_constituents:
            f += f'{constituent}\n'
            vertices = self.mesh.hgrid.get_xy(
                    crs='EPSG:4326')[self.indexes, :]
            amp, phase = self.forcing.get_elevation(constituent, vertices)
            for i in range(len(vertices)):
                f += f'{amp[i]:.8e} {phase[i]:.8e}\n'
        return f

class TidalVelocityWritter(TidalVariableWritter):

    def __str__(self):
        f = ''
        for constituent in self.active_constituents:
            f += f'{constituent}\n'
            vertices = self.mesh.hgrid.get_xy(
                    crs='EPSG:4326')[self.indexes, :]
            uamp, uphase, vamp, vphase = self.forcing.get_velocity(
                    constituent, vertices)
            for i in range(len(vertices)):
                f += f'{uamp[i]:.8e} {uphase[i]:.8e} ' \
                     f'{vamp[i]:.8e} {vphase[i]:.8e}\n'
        return f


class iettypeWritter(Enum):
    NONE = NullWritter
    TIME_VARYING = NotImplementedWritter
    CONSTANT = NotImplementedWritter
    TIDAL = TidalElevationWritter
    SPACE_TIME_VARYING = NotImplementedWritter
    TIDAL_AND_SPACE_TIME_VARYING = NotImplementedWritter

class ifltypeWritter(Enum):
    NONE = NullWritter
    TIME_VARYING = NotImplementedWritter
    CONSTANT = NotImplementedWritter
    TIDAL = TidalVelocityWritter
    SPACE_TIME_VARYING = NotImplementedWritter
    TIDAL_AND_SPACE_TIME_VARYING = NotImplementedWritter

class itetypeWritter(Enum):
    NONE = NullWritter
    TIME_VARYING = NotImplementedWritter
    CONSTANT = NotImplementedWritter
    INITIAL_PROFILE_FOR_INFLOW = NotImplementedWritter
    INPUT_3D = NotImplementedWritter

class isatypeWritter(Enum):
    NONE = NullWritter
    TIME_VARYING = NotImplementedWritter
    CONSTANT = NotImplementedWritter
    INITIAL_PROFILE_FOR_INFLOW = NotImplementedWritter
    INPUT_3D = NotImplementedWritter

class itrtypeWritter(Enum):
    NONE = NullWritter
    TIME_VARYING = NotImplementedWritter
    CONSTANT = NotImplementedWritter
    INITIAL_PROFILE_FOR_INFLOW = NotImplementedWritter
    INPUT_3D = NotImplementedWritter

class Bctides:

    def __init__(self, mesh: Mesh, start_date: datetime, end_date: datetime,
                 spinup_time: timedelta, cutoff_depth: float = 50.):
        self.mesh = mesh
        self.start_date = start_date
        self.end_date = end_date
        self.spinup_time = spinup_time
        self.cutoff_depth = cutoff_depth

    def write(self, path, overwrite=False):
        with open(path, 'w') as f:
            f.write(self.bctides)

    @property
    def bctides(self):
        f = f"{self.forcing_start_date}\n" \
            f"{self.ntip} {self.cutoff_depth}\n"
        if self.ntip > 0:
            for constituent in self.active_potential_constituents:
                forcing = self.tidal_forcing(constituent)
                f += f'{constituent} \n' \
                     f'{forcing[0]:G} ' \
                     f"{forcing[1]:G} " \
                     f'{forcing[2]:G} ' \
                     f'{forcing[3]:G} ' \
                     f'{forcing[4]:G}\n'
        f += f'{self.tidal_forcing.nbfr:d}\n'
        for constituent in self.active_forcing_constituents:
            forcing = self.tidal_forcing(constituent)
            f += f'{constituent} \n' \
                 f"{forcing[2]:G} " \
                 f'{forcing[3]:G} ' \
                 f'{forcing[4]:G}\n'
        f += f"{len(self.mesh.open_boundaries)}\n"  # nope
        for id in self.mesh.open_boundaries:
            boundary = self.mesh.open_boundaries[id]
            f += f"{len(boundary.indexes)} " \
                 f'{boundary.forcing.bctype}\n' \
                 f'{iettypeWritter[boundary.forcing.iettype.name].value(boundary, self)}' \
                 f'{ifltypeWritter[boundary.forcing.ifltype.name].value(boundary, self)}' \
                 f'{itetypeWritter[boundary.forcing.itetype.name].value(boundary, self)}' \
                 f'{isatypeWritter[boundary.forcing.isatype.name].value(boundary, self)}' \
                 f'{itrtypeWritter[boundary.forcing.itrtype.name].value(boundary, self)}'
        return f

    @property  # type: ignore[misc]
    @lru_cache(maxsize=None)
    def forcing_start_date(self):
        return self.start_date - self.spinup_time

    @property
    @lru_cache(maxsize=None)
    def active_potential_constituents(self):
        # PySCHISM allows the user to input the tidal potentials individually
        # for each boundary, however, SCHISM supports only a global
        # specification. Here, we collect all the activated tidal potentials
        # on each boundary and activate them all globally
        const = dict()
        for id in self.mesh.open_boundaries:
            forcing = self.mesh.open_boundaries[id].forcing
            if isinstance(forcing, Tides):
                for active in forcing.get_active_potential_constituents():
                    const[active] = True
        return tuple(const.keys())

    @property
    @lru_cache(maxsize=None)
    def active_forcing_constituents(self):
        # PySCHISM allows the user to input the tidal forcings individually
        # for each boundary, however, SCHISM supports only a global
        # specification. Here, we collect all the activated tidal forcings
        # on each boundary and activate them all globally
        const = dict()
        for id in self.mesh.open_boundaries:
            forcing = self.mesh.open_boundaries[id].forcing
            if isinstance(forcing, Tides):
                for active in forcing.get_active_forcing_constituents():
                    const[active] = True
        return tuple(const.keys())

    @property  # type: ignore[misc]
    @lru_cache(maxsize=None)
    def ntip(self):
        return len(self.active_potential_constituents)

    @property
    @lru_cache(maxsize=None)
    def tidal_forcing(self):
        tides = Tides()
        tides.start_date = self.start_date
        tides.end_date = self.end_date
        tides.spinup_time = self.spinup_time
        for const in tides.constituents:
            tides.use_constituent(
                    const,
                    potential=True if const in 
                              self.active_potential_constituents else False,
                    forcing=True if const in
                            self.active_forcing_constituents else False
                    )
        return tides

