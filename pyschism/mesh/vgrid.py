import os
import pathlib

from typing import Union


class Vgrid:

    def __init__(
        self,
        vgrid: Union[str, os.PathLike] = None
    ):
        """Represents a SCHISM vertical grid.

        Args:
            vgrid (optional): A path to a file or None. If vgrid is None, it
            is assummed that the user wants to use 2D.

        WARNING: This class only support 2D at the moment and it will ignore
        any inputs. This is so it can be used as a placeholder that outputs at
        least a 2D grid for minimalistic model configuration.
        """
        pass

    @staticmethod
    def open(path):
        raise NotImplementedError('Vgrid.open()')

    def __str__(self):
        return """2 !ivcor
2 1 1.e6 !nvrt, kz (# of Z-levels); h_s (transition depth between S and Z)
Z levels
1  -1.e6
S levels
40. 1. 1.e-4  !h_c, theta_b, theta_f
   1    -1.
   2    0."""

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
