from pyschism.forcing.source_sink.base import (
    # Hydrology,
    # Sources,
    # Sinks,
    Msource,
    Vsource,
    Vsink,
    SourceSink,
)
from pyschism.forcing.source_sink.nwm import NationalWaterModel
from pyschism.forcing.source_sink.ngen import NextGen


NWM = NationalWaterModel
NGen = NextGen

__all__ = [
    # "Hydrology",
    "NationalWaterModel",
    "NextGen",
    # "Sources",
    # "Sinks",
    "Msource",
    "Vsource",
    "Vsink",
    "SourceSink",
    'NWM',
    'NGen',    
]
