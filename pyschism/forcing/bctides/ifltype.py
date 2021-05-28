from abc import abstractmethod

from pyschism.forcing.bctides.bctypes import Bctype
from pyschism.forcing.tides import Tides


class Ifltype(Bctype):

    @property
    @abstractmethod
    def ifltype(self) -> int:
        '''Returns integer representig SCHISM ifltype code for bctides.in'''


class UniformTimeHistoryVelocity(Ifltype):

    def __init__(self, time_history):
        raise NotImplementedError(f'{self.__class__.__name__}')
        self.time_history = time_history

    @property
    def ifltype(self):
        return 1


class ConstantVelocity(Ifltype):

    def __init__(self, value):
        raise NotImplementedError(f'{self.__class__.__name__}')
        self.value = value

    @property
    def ifltype(self) -> int:
        return 2


def get_boundary_string(self, hgrid, boundary):
    f = []
    for constituent in \
            self.tides.get_active_forcing_constituents():
        f.append(f'{constituent}')
        vertices = hgrid.get_xy(crs='EPSG:4326')[boundary.indexes, :]
        uamp, uphase, vamp, vphase = self.tides.get_velocity(
            constituent, vertices)
        for i in range(len(vertices)):
            f.append(f'{uamp[i]:.8e} {uphase[i]:.8e} '
                     f'{vamp[i]:.8e} {vphase[i]:.8e}')
    return '\n'.join(f)


class TidalVelocity(Ifltype):

    def __init__(self, tides: Tides):
        if not isinstance(tides, Tides):
            raise TypeError(
                f'Argument tides must be an isinstance of {Tides} not type'
                f'{type(tides)}.')
        self.tides = tides

    def get_boundary_string(self, hgrid, boundary):
        return get_boundary_string(self, hgrid, boundary)

    @property
    def ifltype(self):
        return 3


class SpatiallyVaryingTimeHistoryVelocity(Ifltype):

    def __init__(self, data_source):
        self.data_source = data_source

    def write(
            self,
            path,
            hgrid,
            vgrid,
            start_date,
            run_days,
            overwrite: bool = False
    ):
        self.data_source.velocity.write(
            path,
            hgrid,
            vgrid,
            start_date,
            run_days,
            overwrite
        )

    @property
    def ifltype(self):
        return 4


class TidalAndSpatiallyVaryingVelocityTimeHistory(Ifltype):

    def __init__(self, tides, data_source):
        if not isinstance(tides, Tides):
            raise TypeError(
                f'Argument tides must be an isinstance of {Tides} not type'
                f'{type(tides)}.')
        self.tides = tides
        self.data_source = data_source

    def write(
            self,
            path,
            hgrid,
            vgrid,
            start_date,
            run_days,
            overwrite: bool = False
    ):
        self.data_source.velocity.write(
            path,
            hgrid,
            vgrid,
            start_date,
            run_days,
            overwrite
        )

    @property
    def ifltype(self):
        return 5

    def get_boundary_string(self, hgrid, boundary):
        return get_boundary_string(self, hgrid, boundary)


class ZeroVelocity(Ifltype):

    def __init__(self):
        raise NotImplementedError(f'{self.__class__.__name__}')

    @property
    def ifltype(self):
        return -1


class InflowVelocity(Ifltype):

    def __init__(self):
        raise NotImplementedError(f'{self.__class__.__name__}')

    @property
    def ifltype(self):
        return -4


class OutflowVelocity(Ifltype):

    def __init__(self):
        raise NotImplementedError(f'{self.__class__.__name__}')

    @property
    def ifltype(self):
        return -5


Ifltype1 = UniformTimeHistoryVelocity
Ifltype2 = ConstantVelocity
Ifltype3 = TidalVelocity
Ifltype4 = SpatiallyVaryingTimeHistoryVelocity
Ifltype5 = TidalAndSpatiallyVaryingVelocityTimeHistory
Ifltype_1 = ZeroVelocity
Flather = ZeroVelocity
