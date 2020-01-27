class CORE:
    """ CORE group """

    def __init__(self, nml):
        self.nml = nml

    @property
    def ipre(self):
        return self.nml['core']['ipre']

    @property
    def ibc(self):
        return self.nml['core']['ibc']

    @property
    def ibtp(self):
        return self.nml['core']['ibtp']

    @property
    def rnday(self):
        return self.nml['core']['rnday']

    @property
    def dt(self):
        return self.nml['core']['dt']

    @property
    def msc2(self):
        return self.nml['core']['msc2']

    @property
    def mdc2(self):
        return self.nml['core']['mdc2']

    @property
    def ntracer_gen(self):
        return self.nml['core']['ntracer_gen']

    @property
    def ntracer_age(self):
        return self.nml['core']['ntracer_age']

    @property
    def sed_class(self):
        return self.nml['core']['sed_class']

    @property
    def eco_class(self):
        return self.nml['core']['eco_class']

    @property
    def nspool(self):
        return self.nml['core']['nspool']

    @property
    def ihfskip(self):
        return self.nml['core']['ihfskip']

    @ipre.setter
    def ipre(self, ipre):
        assert ipre in [0, 1]
        self.nml['core']['ipre'] = ipre

    @ibc.setter
    def ibc(self, ibc):
        assert ibc in [0, 1]
        self.nml['core']['ibc'] = ibc

    @ibtp.setter
    def ibtp(self, ibtp):
        self.nml['core']['ibtp'] = ibtp

    @rnday.setter
    def rnday(self, rnday):
        self.nml['core']['rnday'] = rnday

    @dt.setter
    def dt(self, dt):
        self.nml['core']['dt'] = dt

    @msc2.setter
    def msc2(self, msc2):
        self.nml['core']['msc2'] = msc2

    @mdc2.setter
    def mdc2(self, mdc2):
        self.nml['core']['mdc2'] = mdc2

    @ntracer_gen.setter
    def ntracer_gen(self, ntracer_gen):
        self.nml['core']['ntracer_gen'] = ntracer_gen

    @ntracer_age.setter
    def ntracer_age(self, ntracer_age):
        self.nml['core']['ntracer_age'] = ntracer_age

    @sed_class.setter
    def sed_class(self, sed_class):
        self.nml['core']['sed_class'] = sed_class

    @eco_class.setter
    def eco_class(self, eco_class):
        self.nml['core']['eco_class'] = eco_class

    @nspool.setter
    def nspool(self, nspool):
        self.nml['core']['nspool'] = nspool

    @ihfskip.setter
    def ihfskip(self, ihfskip):
        self.nml['core']['ihfskip'] = ihfskip
