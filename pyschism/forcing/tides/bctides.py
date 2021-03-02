from abc import ABC, abstractmethod
from datetime import timedelta, timezone
from enum import Enum
from functools import lru_cache
import logging

from pyschism.domain import ModelDomain
from pyschism.forcing.tides.tides import Tides
from pyschism.param.param import Param

_logger = logging.getLogger(__name__)

class NullWritter:

    def __init__(self, *argv):
        pass

    def __str__(self):
        return ""


class NotImplementedWritter:

    def __init__(self, *argv):
        raise NotImplementedError('Writter for variable is not implemented.')


class TidalVariableWritter(ABC):

    def __init__(self, boundary, bctides):
        self.forcing = boundary['forcing']
        self.indexes = boundary['indexes']
        self.active_constituents = bctides.get_active_forcing_constituents()
        self._model_domain = bctides._model_domain

    @abstractmethod
    def __str__(self):
        raise NotImplementedError(f'str({self.__name__})')


class TidalElevationWritter(TidalVariableWritter):

    def __str__(self):
        f = ""
        for constituent in self.active_constituents:
            f += f'{constituent}\n'
            vertices = self._model_domain.hgrid.get_xy(
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
            vertices = self._model_domain.hgrid.get_xy(
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

    def __init__(self, model_domain: ModelDomain, param: Param,
                 cutoff_depth: float = 50.):
        """Provides an interface to write bctides.in to file. """
        _logger.info('Initializing Bctides.')
        # check if start_date was given in case tidal forcings are requested.
        # Note: This is done twice so that this class can be used independently
        # from Param to just write bctides files
        afc = model_domain.get_active_forcing_constituents()
        if len(afc) > 0 and param.opt.start_date is None:
            raise Exception('start_date argument is required for simulating '
                            'tidal forcing.')

        self._model_domain = model_domain
        self._param = param
        self._cutoff_depth = cutoff_depth

        # init the main tidal forcing object
        tides = Tides()
        for const in tides.all_constituents:
            tides.use_constituent(
                const,
                potential=True if const in
                self.get_active_potential_constituents() else False,
                forcing=True if const in
                self.get_active_forcing_constituents() else False
            )
        self.__tidal_forcing = tides

    def __str__(self):
        f = f"{self.start_date}\n" \
            f"{self.ntip} {self._cutoff_depth}\n"
        if self.ntip > 0:
            for constituent in self.get_active_potential_constituents():
                forcing = self.tidal_forcing(
                    self.start_date, self.rnday, constituent)
                f += f'{constituent} \n' \
                     f'{forcing[0]:G} ' \
                     f"{forcing[1]:G} " \
                     f'{forcing[2]:G} ' \
                     f'{forcing[3]:G} ' \
                     f'{forcing[4]:G}\n'
        f += f'{self.tidal_forcing.nbfr:d}\n'
        for constituent in self.get_active_forcing_constituents():
            forcing = self.tidal_forcing(
                self.start_date, self.rnday, constituent)
            f += f'{constituent} \n' \
                 f"{forcing[2]:G} " \
                 f'{forcing[3]:G} ' \
                 f'{forcing[4]:G}\n'
        f += f"{len(self._model_domain.open_boundaries)}\n"  # nope
        for id, data in self._model_domain.open_boundaries:
            f += f"{len(data['indexes'])} " \
                 f'{str(data["forcing"])}\n' \
                 f'{iettypeWritter[data["forcing"].iettype.name].value(data, self)}' \
                 f'{ifltypeWritter[data["forcing"].ifltype.name].value(data, self)}' \
                 f'{itetypeWritter[data["forcing"].itetype.name].value(data, self)}' \
                 f'{isatypeWritter[data["forcing"].isatype.name].value(data, self)}' \
                 f'{itrtypeWritter[data["forcing"].itrtype.name].value(data, self)}'
        return f

    @lru_cache(maxsize=1)
    def get_active_potential_constituents(self):
        # PySCHISM allows the user to input the tidal potentials and forcings
        # individually at each boundary, however, SCHISM supports only a global
        # specification. Here, we collect all the activated tidal potentials
        # on each boundary and activate them all globally
        # set active tidal potential constituents
        const = dict()
        for id, data in self._model_domain.open_boundaries:
            forcing = data['forcing']
            if isinstance(forcing, Tides):
                for active in forcing.get_active_potential_constituents():
                    const[active] = True
        return tuple(const.keys())

    @lru_cache(maxsize=1)
    def get_active_forcing_constituents(self):
        # set active tidal forcing constituents
        const = dict()
        for id, data in self._model_domain.open_boundaries:
            forcing = data['forcing']
            if isinstance(forcing, Tides):
                for active in forcing.get_active_forcing_constituents():
                    const[active] = True
        return tuple(const.keys())

    def write(self, path, overwrite=False):
        with open(path, 'w') as f:
            f.write(str(self))

    @property
    def start_date(self):
        return self._param.opt.start_date

    @property
    def rnday(self):
        return self._param.core.rnday

    @property
    def ntip(self):
        return len(self.get_active_potential_constituents())

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
