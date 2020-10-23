from datetime import datetime, timezone, timedelta
import pathlib
from typing import Union, List
import warnings

import f90nml
import numpy as np


class StationsOutput:

    def __init__(self, outputs: Union[str, pathlib.Path],
                 stations_file: Union[str, pathlib.Path] = None,
                 ):
        outputs = pathlib.Path(outputs)
        if not outputs.is_dir():
            raise Exception(f'{str(outputs)} is not a directory.')
        self.__station_id: List = []
        if stations_file is None:
            stations_file = outputs.resolve() / '../station.in'
            if not stations_file.is_file():
                stations_file = None
            else:
                with open(stations_file) as f:
                    f.readline()
                    for station in range(int(f.readline())):
                        line = f.readline()
                        if '!' in line:
                            self.__station_id.append(line.split('!')[-1])
                        else:
                            self.__station_id.append(None)
        self._manifest = list(outputs.glob('staout_*'))

        # find the start_date. Let's first check if it's on the param.nml file
        param_file = outputs.resolve() / '../param.nml'
        if param_file.is_file():
            param = f90nml.read(param_file)
            start_year = param['opt'].get('start_year')
            start_month = param['opt'].get('start_month')
            start_day = param['opt'].get('start_day')
            start_hour = param['opt'].get('start_hour')
            utc_start = param['opt'].get('utc_start')
            if None in [start_year, start_month, start_day, start_hour,
                        utc_start]:
                warnings.warn('Could not determine start date automatically.')
                self.__start_date = None
            else:
                self.__start_date = datetime(start_year, start_month,
                                             start_day, int(start_hour), 0,
                                             tzinfo=timezone(timedelta(
                                                hours=-utc_start)))
        self.__rndays = None
        for file in self._manifest:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                data = np.loadtxt(file)
                if len(w) != 0:
                    pass
                else:
                    self.__rndays = timedelta(seconds=data[-1, 0])
                    break

        if self.rndays is None:
            raise Exception("Ouptut directory doesn't contain any station "
                            "output data.")

    def set_start_date(self, start_date: datetime):
        self.__start_date = start_date

    def get_station_id(self, index):
        return self.__station_id[index].strip()

    @property
    def start_date(self):
        return self.__start_date

    @property
    def rndays(self):
        return self.__rndays
