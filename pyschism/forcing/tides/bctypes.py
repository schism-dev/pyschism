from enum import Enum


class BcType(Enum):
    ELEVATION = 'iettype'
    FLOW = 'ifltype'
    TEMPERATURE = 'itetype'
    SALINITY = 'isatype'
    TRACER = 'itrtype'


class InitialElevationType(Enum):
    NONE = 0
    TIME_VARYING = 1
    CONSTANT = 2
    TIDAL = 3
    SPACE_TIME_VARYING = 4
    TIDAL_AND_SPACE_TIME_VARYING = 5


class InitialFlowType(Enum):
    NONE = 0
    TIME_VARYING = 1
    CONSTANT = 2
    TIDAL = 3
    SPACE_TIME_VARYING = 4
    SPACE_TIME_VARYING_3D = -4
    TIDAL_AND_SPACE_TIME_VARYING_3D = 5
    FLANTHER = -1


class InitialTemperatureType(Enum):
    NONE = 0
    TIME_VARYING = 1
    CONSTANT = 2
    INITIAL_PROFILE_FOR_INFLOW = 3
    INPUT_3D = 4


class InitialSalinityType(Enum):
    NONE = 0
    TIME_VARYING = 1
    CONSTANT = 2
    INITIAL_PROFILE_FOR_INFLOW = 3
    INPUT_3D = 4


class InitialTracerType(Enum):
    NONE = 0
    TIME_VARYING = 1
    CONSTANT = 2
    INITIAL_PROFILE_FOR_INFLOW = 3
    INPUT_3D = 4


class InitialConditionType(Enum):
    ELEVATION = InitialElevationType
    FLOW = InitialFlowType
    TEMPERATURE = InitialTemperatureType
    SALINITY = InitialSalinityType
    TRACER = InitialTracerType


class BoundaryCondition:

    def __init__(
            self,
            iettype: InitialElevationType = InitialElevationType.NONE,
            ifltype: InitialFlowType = InitialFlowType.NONE,
            itetype: InitialTemperatureType = InitialTemperatureType.NONE,
            isatype: InitialSalinityType = InitialSalinityType.NONE,
            itrtype: InitialTracerType = InitialTracerType.NONE
    ):
        self.iettype = iettype
        self.ifltype = ifltype
        self.itetype = itetype
        self.isatype = isatype
        self.itrtype = itrtype

    def __str__(self):
        return f"{self.iettype.value} " \
               f"{self.ifltype.value} " \
               f"{self.itetype.value} " \
               f"{self.isatype.value} " \
               f"{self.itrtype.value}"
