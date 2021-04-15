from abc import ABC, abstractmethod
from datetime import datetime, timedelta
import os
import pathlib
import shutil
import subprocess
from typing import Union

import f90nml

from pyschism import dates


class Hotstart(ABC):

    @property
    @abstractmethod
    def time(self) -> datetime:
        """Datetime object for the hotstart snapshot"""

    @property
    @abstractmethod
    def iteration(self) -> int:
        """Iteration number for hotstart snapshot"""

    @property
    @abstractmethod
    def path(self) -> pathlib.Path:
        """Path for output file."""

    @property
    def status(self):
        return self.path.exists()

    @staticmethod
    def combine(path):
        return CombineHotstartBinary(path)

    def move(self, dest, overwrite=False):
        self._path = shutil.move(self.path, dest)

    def symlink(self, dest, overwrite=False):
        hotstart_lnk = dest / 'hotstart.nc'
        if overwrite is True:
            if hotstart_lnk.exists():
                hotstart_lnk.unlink()
        os.symlink(os.path.relpath(
            self.hotstart.path, dest), hotstart_lnk)


class CombineHotstartBinary(Hotstart):

    def __init__(self, path: Union[str, os.PathLike], iteration=None,
                 binary='combine_hotstart7'):
        path = pathlib.Path(path)
        self.binary = binary
        if iteration is not None:
            self.iteration = iteration
        else:
            param = path / 'param.out.nml'
            if param.exists():
                param = f90nml.read(param)
            else:
                param = f90nml.read(path.parent / 'param.nml')
            self.iteration = int(
                timedelta(days=param['core']['rnday']) /
                timedelta(seconds=param['core']['dt']))
        self._path = pathlib.Path(path) / f"hotstart_it={self.iteration}.nc"

    def run(self):
        subprocess.check_call([
            f'{self.binary}', '-i', f'{self.iteration}'], cwd=self.path.parent)

    @property
    def time(self):
        param = f90nml.read(self.path.parent.parent / 'param.nml')
        start_hour = float(param['opt']['start_hour'])
        hours = int(start_hour)
        minutes = int((start_hour*60) % 60)
        seconds = int((start_hour*3600) % 60)
        return dates.localize_datetime(
            datetime(
                int(param['opt']['start_year']),
                int(param['opt']['start_month']),
                int(param['opt']['start_day']),
                hours, minutes, seconds
            )
            + timedelta(days=param['core']['rnday'])
            - timedelta(days=param['opt']['utc_start']))

    @property
    def iteration(self):
        return self._iteration

    @iteration.setter
    def iteration(self, iteration: int):
        self._iteration = int(iteration)

    @property
    def path(self):
        return self._path

    @property
    def binary(self) -> pathlib.Path:
        return self._binary

    @binary.setter
    def binary(self, binary):
        self._binary = pathlib.Path(binary)
