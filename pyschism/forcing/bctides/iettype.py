from abc import abstractmethod
from typing import List

from pyschism.forcing.bctides.bctypes import Bctype
from pyschism.forcing.tides import Tides


class Iettype(Bctype):

    @property
    @abstractmethod
    def iettype(self) -> int:
        pass


class UniformTimeHistoryElevation(Iettype):

    def __init__(self, time_history: List[List[float]], boundaries):
        self.time_history = time_history
        self.boundaries = boundaries

    def get_boundary_string(self, hgrid, boundary):
        return ""

    @property
    def iettype(self) -> int:
        return 1


class ConstantElevation(Iettype):

    def __init__(self, values: List[float]):
        self.values = List[float]

    def get_boundary_string(self, hgrid, boundary):
        return f'{self.values[boundary.index]:G}'

    @property
    def iettype(self) -> int:
        return 2


def get_boundary_string(self, hgrid, boundary):
    f = []
    for constituent in self.tides.get_active_forcing_constituents():
        f.append(f'{constituent}')
        vertices = hgrid.get_xy(crs='EPSG:4326')[boundary.indexes, :]
        amp, phase = self.tides.get_elevation(
            constituent, vertices)
        for i in range(len(amp)):
            f.append(f'{amp[i]:.8e} {phase[i]:.8e}')
    return '\n'.join(f)


class TidalElevation(Iettype):

    def __init__(self, tides: Tides):
        self.tides = tides

    def get_boundary_string(self, hgrid, boundary):
        return get_boundary_string(self, hgrid, boundary)

    @property
    def iettype(self):
        return 3


class SpatiallyVaryingTimeHistoryElevation(Iettype):

    def __init__(self, data_source):
        self.data_source = data_source
        # self.start_date = start_date
        # self.run_days = run_days

    def get_boundary_string(self, hgrid, boundary):
        return ''

    def write(
            self,
            path,
            hgrid,
            start_date,
            run_days,
            overwrite: bool = False
    ):
        self.data_source.elevation.write(
            path,
            hgrid,
            start_date,
            run_days,
            overwrite
        )

    @property
    def iettype(self):
        return 4


class TidalAndSpatiallyVaryingElevationTimeHistory(Iettype):

    def __init__(self, tides, data_source):
        self.tides = tides
        self.data_source = data_source

    def get_boundary_string(self, hgrid, boundary):
        return get_boundary_string(self, hgrid, boundary)

    def write(
            self,
            path,
            hgrid,
            start_date,
            run_days,
            overwrite: bool = False
    ):
        self.data_source.elevation.write(
            path,
            hgrid,
            start_date,
            run_days,
            overwrite
        )

    @property
    def iettype(self):
        return 5


class ZeroElevation(Iettype):

    def __init__(self):
        raise NotImplementedError(f'{self.__class__.__name__}')

    def get_boundary_string(self):
        return ''

    @property
    def iettype(self):
        return -1


Iettype1 = UniformTimeHistoryElevation
Iettype2 = ConstantElevation
Iettype3 = TidalElevation
Iettype4 = SpatiallyVaryingTimeHistoryElevation
Iettype5 = TidalAndSpatiallyVaryingElevationTimeHistory
Iettype_1 = ZeroElevation
