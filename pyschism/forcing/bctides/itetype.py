from abc import abstractmethod
from enum import Enum

# from typing import Union

from pyschism.forcing.bctides.bctypes import Bctype
from pyschism.forcing import hycom

# from pyschism.forcing.bctides.nudge import TempNudge
# from pyschism.mesh.base import Gr3


class Itetype(Bctype):
    @property
    @abstractmethod
    def itetype(self) -> int:
        pass

    @property
    def forcing_digit(self):
        return self.itetype


class UniformTimeHistoryTemperature(Itetype):

    # def __init__(self, time_history: Dict[str, List[float]], tobc: float = 1.):
    #     self.time_history = time_history
    #     super().__init__(tobc)

    @property
    def itetype(self) -> int:
        return 1


class ConstantTemperature(Itetype):

    def __init__(self, value: float, nudging_factor: float):
        self.value = value
        if not (nudging_factor >= 0) and (nudging_factor <= 1):
            raise ValueError(
                    'Argument `nudging_factor` must be >= 0 and <= 1,'
                    f'but got {nudging_factor}.')
        self.nudging_factor = nudging_factor
        super().__init__()

    @property
    def itetype(self) -> int:
        return 2

    def get_boundary_string(self, *args, **kwargs) -> str:
        boundary_string = [
            f'{self.value:0.6f}',
            f'{self.nudging_factor:0.6f}',
        ]
        return '\n'.join(boundary_string)


class TemperatureInitialConditions(Itetype):

    # def __init__(self, tobc: float = 1., tth):
    #     raise NotImplementedError(f'{self.__class__.__name__}')

    @property
    def itetype(self):
        return 3


class SpatiallyVaryingTimeHistoryTemperature(Itetype):
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
        self.data_component = data_source.temperature
        self.nudge = bool(nudge)
        self.rlmax = rlmax
        self.rnu_day = rnu_day

    def get_boundary_string(self, *args, **kwargs):
        return "1."

    # def write(self, *args, **kwargs):
    #     if self.nudge is not None:
    #         self.nudge.write(*args, **kwargs)

    @property
    def itetype(self):
        return 4


class Itetype0(Itetype):
    pass


class Itetype1(UniformTimeHistoryTemperature):
    pass


class Itetype2(ConstantTemperature):
    pass


class Itetype3(TemperatureInitialConditions):
    pass


class Itetype4(SpatiallyVaryingTimeHistoryTemperature):
    pass
