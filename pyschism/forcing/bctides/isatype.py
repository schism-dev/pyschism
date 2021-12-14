from abc import abstractmethod
from enum import Enum

from pyschism.forcing.bctides.bctypes import Bctype
from pyschism.forcing import hycom


class Isatype(Bctype):
    @property
    @abstractmethod
    def isatype(self) -> int:
        pass


class UniformTimeHistorySalinity(Isatype):
    def __init__(self, time_history):
        raise NotImplementedError(f"{self.__class__.__name__}")
        self.time_history = time_history

    @property
    def isatype(self) -> int:
        return 1


class ConstantSalinity(Isatype):
    def __init__(self, value):
        raise NotImplementedError(f"{self.__class__.__name__}")
        self.value = value

    @property
    def isatype(self) -> int:
        return 2


class SalinityInitialConditions(Isatype):
    def __init__(self):
        raise NotImplementedError(f"{self.__class__.__name__}")

    @property
    def isatype(self):
        return 3


class SpatiallyVaryingTimeHistorySalinity(Isatype):
    class BaroclinicDatabases(Enum):
        #RTOFS = hycom.RTOFS
        GOFS = hycom.GOFS

    def __init__(
        self,
        data_source="gofs",
        nudge: bool = True,
        rlmax=1.5,
        rnu_day=0.25,
    ):
        if isinstance(data_source, str):
            data_source = self.BaroclinicDatabases[data_source.upper()].value()
        if not isinstance(data_source, hycom.Hycom):
            raise TypeError(
                "Argument data_source must be of type str or type "
                f"{type(hycom.Hycom)}, not type {type(data_source)}."
            )
        self.data_source = data_source
        self.data_component = data_source.salinity
        self.nudge = nudge
        self.rlmax = rlmax
        self.rnu_day = rnu_day

    def get_boundary_string(self, *args, **kwargs):
        return "1."

    @property
    def isatype(self):
        return 4


Isatype1 = UniformTimeHistorySalinity
Isatype2 = ConstantSalinity
Isatype3 = SalinityInitialConditions
Isatype4 = SpatiallyVaryingTimeHistorySalinity
