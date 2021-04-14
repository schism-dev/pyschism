from abc import ABC, abstractmethod
from datetime import datetime, timedelta
import os
import pathlib
import shutil
import subprocess
from typing import Union

import f90nml
import numpy as np

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

    def __init__(self, path: Union[str, os.PathLike], iteration=None):
        path = pathlib.Path(path)
        if iteration is None:
            combine_hotstart = path.glob('hotstart_[0-9][0-9][0-9][0-9]_*.nc')
            increments = set([file.name.split('_')[-1].split('.nc')[0]
                              for file in combine_hotstart])
            iteration = np.max(
                [int(increment) for increment in increments])
        subprocess.check_call(
            ["combine_hotstart7", '-i', f'{iteration}'], cwd=path)
        self._path = path / f"hotstart_it={iteration}.nc"
        self._iteration = iteration

    @property
    def time(self):
        # return Param.read(self.path.parent / 'param.nml').start_date
        param = f90nml.read(self.path.parent / 'param.out.nml')
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

    @property
    def path(self):
        return self._path
