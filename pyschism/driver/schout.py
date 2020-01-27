

class SCHOUT:
    """ SCHOUT group """

    def __init__(self, nml):
        self.nml = nml

    @property
    def nhot(self):
        return self.nml['schout']['nhot']

    @property
    def nhot_write(self):
        return self.nml['schout']['nhot_write']

    @property
    def iout_sta(self):
        return self.nml['schout']['iout_sta']

    @property
    def nspool_sta(self):
        return self.nml['schout']['nspool_sta']

    @property
    def iof_hydro(self):
        return self.nml['schout']['iof_hydro']

    @property
    def iof_wwm(self):
        return self.nml['schout']['iof_wwm']

    @property
    def iof_gen(self):
        return self.nml['schout']['iof_gen']

    @property
    def iof_age(self):
        return self.nml['schout']['iof_age']

    @property
    def iof_sed(self):
        return self.nml['schout']['iof_sed']

    @property
    def iof_eco(self):
        return self.nml['schout']['iof_eco']

    @property
    def iof_icm(self):
        return self.nml['schout']['iof_icm']

    @property
    def iof_cos(self):
        return self.nml['schout']['iof_cos']

    @property
    def iof_fib(self):
        return self.nml['schout']['iof_fib']

    @property
    def iof_sed2d(self):
        return self.nml['schout']['iof_sed2d']

    @property
    def iof_marsh(self):
        return self.nml['schout']['iof_marsh']

    @property
    def iof_ice(self):
        return self.nml['schout']['iof_ice']

    @property
    def iof_ana(self):
        return self.nml['schout']['iof_ana']

    @nhot.setter
    def nhot(self, nhot):
        self.nml['schout']['nhot'] = nhot

    @nhot_write.setter
    def nhot_write(self, nhot_write):
        self.nml['schout']['nhot_write'] = nhot_write

    @iout_sta.setter
    def iout_sta(self, iout_sta):
        self.nml['schout']['iout_sta'] = iout_sta

    @nspool_sta.setter
    def nspool_sta(self, nspool_sta):
        self.nml['schout']['nspool_sta'] = nspool_sta

    @iof_hydro.setter
    def iof_hydro(self, iof_hydro):
        self.nml['schout']['iof_hydro'] = iof_hydro

    @iof_wwm.setter
    def iof_wwm(self, iof_wwm):
        self.nml['schout']['iof_wwm'] = iof_wwm

    @iof_gen.setter
    def iof_gen(self, iof_gen):
        self.nml['schout']['iof_gen'] = iof_gen

    @iof_age.setter
    def iof_age(self, iof_age):
        self.nml['schout']['iof_age'] = iof_age

    @iof_sed.setter
    def iof_sed(self, iof_sed):
        self.nml['schout']['iof_sed'] = iof_sed

    @iof_eco.setter
    def iof_eco(self, iof_eco):
        self.nml['schout']['iof_eco'] = iof_eco

    @iof_icm.setter
    def iof_icm(self, iof_icm):
        self.nml['schout']['iof_icm'] = iof_icm

    @iof_cos.setter
    def iof_cos(self, iof_cos):
        self.nml['schout']['iof_cos'] = iof_cos

    @iof_fib.setter
    def iof_fib(self, iof_fib):
        self.nml['schout']['iof_fib'] = iof_fib

    @iof_sed2d.setter
    def iof_sed2d(self, iof_sed2d):
        self.nml['schout']['iof_sed2d'] = iof_sed2d

    @iof_marsh.setter
    def iof_marsh(self, iof_marsh):
        self.nml['schout']['iof_marsh'] = iof_marsh

    @iof_ice.setter
    def iof_ice(self, iof_ice):
        self.nml['schout']['iof_ice'] = iof_ice

    @iof_ana.setter
    def iof_ana(self, iof_ana):
        self.nml['schout']['iof_ana'] = iof_ana
