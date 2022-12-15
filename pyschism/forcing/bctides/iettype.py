from abc import abstractmethod
from enum import Enum
from typing import List

from pyschism.forcing.bctides.bctypes import Bctype
from pyschism.forcing.bctides.tides import Tides
from pyschism.forcing import hycom


class Iettype(Bctype):

    @property
    @abstractmethod
    def iettype(self) -> int:
        pass

    @property
    def forcing_digit(self):
        return self.iettype


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
        return f"{self.values[boundary.index]:G}"

    @property
    def iettype(self) -> int:
        return 2


class TidalElevation(Iettype):
    def __init__(self, constituents='all', database="hamtide"):
        self.tides = Tides(tidal_database=database, constituents=constituents)
        
    def add_constituent(
            self,
            name: str,
            amplitude,
            angular_frequency,
            phase=0.,
            **kwargs,
    ):
        self.tides.add_constituent(
            name=name,
            angular_frequency=angular_frequency,
            elevation_amplitude=amplitude,
            elevation_phase=phase,
            velocity_amplitude=self.tides._amplitudes['velocity'].get(name, 0.),
            velocity_phase=self.tides._phases['velocity'].get(name, 0.),
            **kwargs
        )

    def get_boundary_string(self, hgrid, boundary, global_constituents: List[str] = None):
        f = []
        crs = 'epsg:4326' if hgrid.crs is not None else None
        if global_constituents is None:
            for constituent in self.tides.get_active_forcing_constituents():
                f.append(f"{constituent}")
                vertices = hgrid.get_xy(crs=crs)[boundary.indexes, :]
                amp, phase = self.tides.get_elevation(constituent, vertices)
                for i in range(len(boundary.indexes)):
                    f.append(f"{amp[i]:.8e} {phase[i]:.8e}")
        else:
            required_constituent = self.tides.get_active_forcing_constituents()
            for constituent in global_constituents:
                f.append(f"{constituent}")
                if constituent not in required_constituent:
                    for i in range(len(boundary.indexes)):
                        f.append(f"{0:.8e} {0:.8e}")
                else:
                    vertices = hgrid.get_xy(crs=crs)[boundary.indexes, :]
                    amp, phase = self.tides.get_elevation(constituent, vertices)
                    for i in range(len(boundary.indexes)):
                        f.append(f"{amp[i]:.8e} {phase[i]:.8e}")

        return "\n".join(f)

    @property
    def iettype(self):
        return 3


class HarmonicElevation(TidalElevation):

    def __init__(self):
        self.tides = Tides(constituents=None)

class SpatiallyVaryingTimeHistoryElevation(Iettype):
    class BaroclinicDatabases(Enum):
        #RTOFS = hycom.RTOFS
        GOFS = hycom.GOFS

    def __init__(self, data_source="gofs"):

        if isinstance(data_source, str):
            data_source = self.BaroclinicDatabases[data_source.upper()].value()

        if not isinstance(data_source, hycom.Hycom):
            raise TypeError(
                "Argument data_source must be of type str or type "
                f"{type(hycom.Hycom)}, not type {type(data_source)}."
            )

        self.data_component: hycom.HycomComponent = data_source.elevation

    def get_boundary_string(self, hgrid, boundary, global_constituents: List[str] = None):
        return ""

    def write(self, path, hgrid, start_date, run_days, overwrite: bool = False):
        self.data_component.write(path, hgrid, start_date, run_days, overwrite)

    @property
    def iettype(self):
        return 4


class TidalAndSpatiallyVaryingElevationTimeHistory(Iettype):
    def __init__(
        self, iettype3: TidalElevation, iettype4: SpatiallyVaryingTimeHistoryElevation
    ):
        self.iettype3 = iettype3
        self.iettype4 = iettype4

    def get_boundary_string(self, hgrid, boundary, global_constituents: List[str] = None):
        return self.iettype3.get_boundary_string(hgrid, boundary, global_constituents)

    def write(self, path, hgrid, start_date, run_days, overwrite: bool = False):
        self.data_component.write(path, hgrid, start_date, run_days, overwrite)

    @property
    def iettype(self):
        return 5

    @property
    def tides(self):
        return self.iettype3.tides

    @property
    def data_component(self):
        return self.iettype4.data_component


class ZeroElevation(Iettype):
    def __init__(self):
        raise NotImplementedError(f"{self.__class__.__name__}")

    def get_boundary_string(self):
        return ""

    @property
    def iettype(self):
        return -1


class Iettype1(UniformTimeHistoryElevation):
    pass


class Iettype2(ConstantElevation):
    pass


class Iettype3(TidalElevation):
    pass


class Iettype4(SpatiallyVaryingTimeHistoryElevation):
    pass


class Iettype5(TidalAndSpatiallyVaryingElevationTimeHistory):
    pass


class Iettype_1(ZeroElevation):
    pass
