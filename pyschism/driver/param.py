import pathlib
import f90nml


class Param:

    def __init__(self, start_date, end_date, spinup_time):
        self._start_date = start_date
        self._end_date = end_date
        self._spinup_time = spinup_time

    def write(self, path, overwrite=False):
        raise NotImplementedError(self.nml)

    @property
    def param(self):
        param = {}
        self._add_core(param)
        self._add_opt_params(param)
        self._add_opt_params(param)
        param = {
            "CORE"
        }
        return {
            "CORE": self.params['core'],
            "OPT": self.opt,
            "SCHOUT": self.schout
        }

    @property
    def core(self):
        return {
            "ipre": self.ipre,
            "ibtp": self.ibtp,
            "rnday": self.rnday,
            "dt": self.dt,
            "msc2": self.msc2,
            "mdc2": self.mdc2,
            "ntracer_gen": self.ntracer_gen,
            "ntracer_age": self.ntracer_age,
            "sed_class": self.sed_class,
            "eco_class": self.eco_class,
            "nspool": self.nspool,
            "ihfskip": self.ihfskip
        }

    @property
    def opt(self):
        return {
            "ipre2": self.ipre2,
            "start_year": self.start_year,
            "start_month": self.start_month,
            "start_day": self.start_day,
            "start_hour": self.start_hour,
            "utc_start": self.utc_start,
            "ics": self.ics,
            "ihot": self.ihot,
        }

    @property
    def defaults(self):
        return f90nml.read(pathlib.Path(__name__) / 'params.nml')
