class CORE:
    """ CORE group """

    def __init__(self, nml):
        self._nml = nml

    @property
    def ipre(self):
        return self._nml['core']['ipre']

    @property
    def ibc(self):
        return self._nml['core']['ibc']

    @property
    def ibtp(self):
        return self._nml['core']['ibtp']

    @property
    def rnday(self):
        return self._nml['core']['rnday']

    @property
    def dt(self):
        return self._nml['core']['dt']

    @property
    def msc2(self):
        return self._nml['core']['msc2']

    @property
    def mdc2(self):
        return self._nml['core']['mdc2']

    @property
    def ntracer_gen(self):
        return self._nml['core']['ntracer_gen']

    @property
    def ntracer_age(self):
        return self._nml['core']['ntracer_age']

    @property
    def sed_class(self):
        return self._nml['core']['sed_class']

    @property
    def eco_class(self):
        return self._nml['core']['eco_class']

    @property
    def nspool(self):
        return self._nml['core']['nspool']

    @property
    def ihfskip(self):
        return self._nml['core']['ihfskip']

    @ipre.setter
    def ipre(self, ipre):
        assert ipre in [0, 1]
        self._nml['core']['ipre'] = ipre

    @ibc.setter
    def ibc(self, ibc):
        assert ibc in [0, 1]
        self._nml['core']['ibc'] = ibc

    @ibtp.setter
    def ibtp(self, ibtp):
        self._nml['core']['ibtp'] = ibtp

    @rnday.setter
    def rnday(self, rnday):
        self._nml['core']['rnday'] = rnday

    @dt.setter
    def dt(self, dt):
        assert dt > 0., "dt (timestep) must a positive float or integer."
        self._nml['core']['dt'] = dt

    @msc2.setter
    def msc2(self, msc2):
        self._nml['core']['msc2'] = msc2

    @mdc2.setter
    def mdc2(self, mdc2):
        self._nml['core']['mdc2'] = mdc2

    @ntracer_gen.setter
    def ntracer_gen(self, ntracer_gen):
        self._nml['core']['ntracer_gen'] = ntracer_gen

    @ntracer_age.setter
    def ntracer_age(self, ntracer_age):
        self._nml['core']['ntracer_age'] = ntracer_age

    @sed_class.setter
    def sed_class(self, sed_class):
        self._nml['core']['sed_class'] = sed_class

    @eco_class.setter
    def eco_class(self, eco_class):
        self._nml['core']['eco_class'] = eco_class

    @nspool.setter
    def nspool(self, nspool):
        assert nspool >= 0, "nspool must be zero or positive integer."
        self._nml['core']['nspool'] = nspool

    @ihfskip.setter
    def ihfskip(self, ihfskip):
        assert ihfskip >= 0, "ihfskip must be zero or positive integer."
        self._nml['core']['ihfskip'] = ihfskip
