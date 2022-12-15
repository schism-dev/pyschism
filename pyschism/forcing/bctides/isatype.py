from abc import abstractmethod
from enum import Enum

from pyschism.forcing.bctides.bctypes import Bctype
from pyschism.forcing import hycom


class Isatype(Bctype):

    @property
    @abstractmethod
    def isatype(self) -> int:
        pass

    @property
    def forcing_digit(self):
        return self.isatype


class UniformTimeHistorySalinity(Isatype):
    def __init__(self, time_history):
        raise NotImplementedError(f"{self.__class__.__name__}")
        self.time_history = time_history

    @property
    def isatype(self) -> int:
        return 1


class ConstantSalinity(Isatype):
    def __init__(self, value: float, nudging_factor: float):
        self.value = value
        if not (nudging_factor >= 0) and (nudging_factor <= 1):
            raise ValueError(f'Argument `nudging_factor` must be get 0 and let 1, but got {nudging_factor}.')
        self.nudging_factor = nudging_factor
        super().__init__()


    @property
    def isatype(self) -> int:
        return 2

    def get_boundary_string(self, *args, **kwargs) -> str:
        boundary_string = [
            f'{self.value:0.6f}',
            f'{self.nudging_factor:0.6f}',
        ]
        return '\n'.join(boundary_string)


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


class Isatype1(UniformTimeHistorySalinity):
    pass


class Isatype2(ConstantSalinity):
    pass


class Isatype3(SalinityInitialConditions):
    pass


class Isatype4(SpatiallyVaryingTimeHistorySalinity):
    pass
