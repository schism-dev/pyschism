import abc


class BoundaryCondition(metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def bctype(self):
        """
        <int>
        1: time history file
        2: constant
        3: initial condition
        4: time history 3D
        5: 3 + 4
        -1: uv2D
        -4,-5: UV3D
        See page 62 of SCHISM v5.7 manual for details
        """

    @property
    @abc.abstractproperty
    def vartype(self):
        """
        <str>: <'elev' | 'vel' | 'stt'>
        """


class ElevBc(BoundaryCondition):

    @property
    def vartype(self):
        return 'elev'


class VelBc(BoundaryCondition):

    @property
    def vartype(self):
        return 'vel'


class STT(BoundaryCondition):

    @property
    def vartype(self):
        return 'stt'
