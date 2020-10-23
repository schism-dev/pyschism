from abc import ABC, abstractmethod
from datetime import timedelta, timezone
from enum import Enum

from .tides import Tides
from ...param.param import Param
from ...mesh import Mesh


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

    def __init__(self, mesh: Mesh, param: Param, cutoff_depth: float = 50.):

        # check if start_date was given in case tidal forcings are requested.
        # Note: This is done twice so that this class can be used independently
        # from Param to just write bctides files
        afc = mesh.get_active_forcing_constituents()
        if len(afc) > 0 and param.opt.start_date is None:
            raise Exception('start_date argument is required for simulating '
                            'tidal forcing.')

        self.__mesh = mesh
        self.__param = param
        self.cutoff_depth = cutoff_depth

        # PySCHISM allows the user to input the tidal potentials and forcings
        # individually at each boundary, however, SCHISM supports only a global
        # specification. Here, we collect all the activated tidal potentials
        # on each boundary and activate them all globally

        # set active tidal potential constituents
        const = dict()
        for id in self.mesh.open_boundaries:
            forcing = self.mesh.open_boundaries[id].forcing
            if isinstance(forcing, Tides):
                for active in forcing.get_active_potential_constituents():
                    const[active] = True
        self.__active_potential_constituents = tuple(const.keys())

        # set active tidal forcing constituents
        for id in self.mesh.open_boundaries:
            forcing = self.mesh.open_boundaries[id].forcing
            if isinstance(forcing, Tides):
                for active in forcing.get_active_forcing_constituents():
                    const[active] = True
        self.__active_forcing_constituents = tuple(const.keys())

        # init the main tidal forcing object
        tides = Tides()
        for const in tides.constituents:
            tides.use_constituent(
                    const,
                    potential=True if const in
                    self.active_potential_constituents else False,
                    forcing=True if const in
                    self.active_forcing_constituents else False
                    )
        self.__tidal_forcing = tides

    def __str__(self):
        f = f"{self.start_date}\n" \
            f"{self.ntip} {self.cutoff_depth}\n"
        if self.ntip > 0:
            for constituent in self.active_potential_constituents:
                forcing = self.tidal_forcing(
                    self.start_date, self.rnday, constituent)
                f += f'{constituent} \n' \
                     f'{forcing[0]:G} ' \
                     f"{forcing[1]:G} " \
                     f'{forcing[2]:G} ' \
                     f'{forcing[3]:G} ' \
                     f'{forcing[4]:G}\n'
        f += f'{self.tidal_forcing.nbfr:d}\n'
        for constituent in self.active_forcing_constituents:
            forcing = self.tidal_forcing(
                self.start_date, self.rnday, constituent)
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

    def write(self, path, overwrite=False):
        with open(path, 'w') as f:
            f.write(str(self))

    @property
    def mesh(self):
        return self.__mesh

    @property
    def start_date(self):
        return self.__param.opt.start_date

    @property
    def rnday(self):
        return self.__param.core.rnday

    @property
    def active_potential_constituents(self):
        return self.__active_potential_constituents

    @property
    def active_forcing_constituents(self):
        return self.__active_forcing_constituents

    @property
    def ntip(self):
        return len(self.active_potential_constituents)

    @property
    def tidal_forcing(self):
        return self.__tidal_forcing

    @property
    def start_date_utc(self):
        if self.start_date.tzinfo is not None and \
                self.start_date.tzinfo.utcoffset(self.start_date) is not None:
            return self.__start_date.astimezone(timezone(timedelta(0)))
        else:
            return self.__start_date
