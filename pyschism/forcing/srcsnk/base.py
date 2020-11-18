
class Elements:

    def __get__(self, obj, val):
        pass


class MsourceDescriptor:

    def __get__(self, obj, val):
        pass


class VsourceDescriptor:

    def __get__(self, obj, val):
        pass


class VsinkDescriptor:

    def __get__(self, obj, val):
        pass


class SourceSink:

    elements = Elements()
    msource = MsourceDescriptor()
    vsource = VsourceDescriptor()
    vsink = VsinkDescriptor()

    def __init__(self, msource, vsource, vsink):
        self.msource = msource
        self.vsource = vsource
        self.vsink = vsink
