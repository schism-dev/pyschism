

class Fgrid:
    """
    Base class for all friction types (e.g. manning.gr3, drag.gr3, etc...)
    """

    def dump(self):
        raise NotImplementedError
