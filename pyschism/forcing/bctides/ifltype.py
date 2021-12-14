from abc import abstractmethod
from enum import Enum

from pyschism.forcing.bctides.bctypes import Bctype
from pyschism.forcing.bctides.tides import Tides
from pyschism.forcing import hycom


class Ifltype(Bctype):
    @property
    @abstractmethod
    def ifltype(self) -> int:
        """Returns integer representig SCHISM ifltype code for bctides.in"""


class UniformTimeHistoryVelocity(Ifltype):
    def __init__(self, time_history):
        raise NotImplementedError(f"{self.__class__.__name__}")
        self.time_history = time_history

    @property
    def ifltype(self):
        return 1


class ConstantVelocity(Ifltype):
    def __init__(self, value):
        raise NotImplementedError(f"{self.__class__.__name__}")
        self.value = value

    @property
    def ifltype(self) -> int:
        return 2


class TidalVelocity(Ifltype):
    def __init__(self, constituents="all", database="hamtide"):
        self.tides = Tides(tidal_database=database, constituents=constituents)

    def get_boundary_string(self, hgrid, boundary, global_constituents=None):
        f = []
        if global_constituents is None:
            for constituent in self.tides.get_active_forcing_constituents():
                f.append(f"{constituent}")
                vertices = hgrid.get_xy(crs="EPSG:4326")[boundary.indexes, :]
                uamp, uphase, vamp, vphase = self.tides.get_velocity(constituent, vertices)
                for i in range(len(vertices)):
                    f.append(
                        f"{uamp[i]:.8e} {uphase[i]:.8e} {vamp[i]:.8e} {vphase[i]:.8e}"
                    )
        else:
            required_constituent = self.tides.get_active_forcing_constituents()
            for constituent in global_constituents:
                f.append(f"{constituent}")
                if constituent not in required_constituent:
                    for i in range(len(boundary.indexes)):
                        f.append(f"{0:.8e} {0:.8e} {0:.8e} {0:.8e}")
                else:
                    vertices = hgrid.get_xy(crs="EPSG:4326")[boundary.indexes, :]
                    uamp, uphase, vamp, vphase = self.tides.get_velocity(constituent, vertices)
                    for i in range(len(boundary.indexes)):
                        f.append(f"{uamp[i]:.8e} {uphase[i]:.8e} {vamp[i]:.8e} {vphase[i]:.8e}")

        return "\n".join(f)

    @property
    def ifltype(self):
        return 3


class SpatiallyVaryingTimeHistoryVelocity(Ifltype):
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

        self.data_component: hycom.HycomComponent = data_source.velocity

    def write(self, path, hgrid, vgrid, start_date, run_days, overwrite: bool = False):
        self.data_component.write(path, hgrid, vgrid, start_date, run_days, overwrite)

    def get_boundary_string(self, *args, **kwargs):
        return ""

    @property
    def ifltype(self):
        return 4


class TidalAndSpatiallyVaryingVelocityTimeHistory(Ifltype):
    def __init__(
        self,
        ifltype3: TidalVelocity,
        ifltype4: SpatiallyVaryingTimeHistoryVelocity,
    ):
        self.ifltype3 = ifltype3
        self.ifltype4 = ifltype4

    def write(self, path, hgrid, vgrid, start_date, run_days, overwrite: bool = False):
        self.data_component.write(path, hgrid, vgrid, start_date, run_days, overwrite)

    @property
    def ifltype(self):
        return 5

    def get_boundary_string(self, *args, **kwargs):
        return self.ifltype3.get_boundary_string(*args, **kwargs)

    @property
    def tides(self):
        return self.ifltype3.tides

    @property
    def data_component(self):
        return self.ifltype4.data_component


class ZeroVelocity(Ifltype):
    def __init__(self):
        raise NotImplementedError(f"{self.__class__.__name__}")

    @property
    def ifltype(self):
        return -1


class InflowVelocity(Ifltype):
    def __init__(self):
        raise NotImplementedError(f"{self.__class__.__name__}")

    @property
    def ifltype(self):
        return -4


class OutflowVelocity(Ifltype):
    def __init__(self):
        raise NotImplementedError(f"{self.__class__.__name__}")

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
