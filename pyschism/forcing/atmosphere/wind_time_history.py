from .nws import NWS, NWSType


class WindTimeHistory(NWS):

    def __init__(self):
        super().__init__(NWSType.TIME_HISTORY)
        raise NotImplementedError
