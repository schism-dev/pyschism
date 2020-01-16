

class Param:

    def write(self, path, overwrite=False):
        raise NotImplementedError(self.nml)

    @property
    def nml(self):
        return {
            "CORE": self.core,
            "OPT": self.opt,
            "SCHOUT": self.schout
        }
