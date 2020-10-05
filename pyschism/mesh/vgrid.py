import pathlib


class Vgrid:

    def __str__(self):
        return """2 !ivcor
2 1 1.e6 !nvrt, kz (# of Z-levels); h_s (transition depth between S and Z)
Z levels
1  -1.e6
S levels
40. 1. 1.e-4  !h_c, theta_b, theta_f
   1    -1.
   2    0."""

    @classmethod
    def open(cls, path):
        return cls()

    def write(self, path, overwrite=False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            msg = 'File exists, pass overwrite=True to allow overwrite.'
            raise Exception(msg)

        with open(path, 'w') as f:
            f.write(str(self))

    def is_2D(self):
        return True

    def is_3D(self):
        return False
