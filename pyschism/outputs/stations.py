from datetime import datetime, timezone, timedelta
import pathlib
from typing import Union, List
import warnings

import f90nml
import matplotlib.pyplot as plt
import numpy as np

from pyschism.utils.coops import CoopsDataCollector
from pyschism.enums import StationOutputIndex, StationOutputVariables


class StationsOutput:

    def __init__(self, outputs: Union[str, pathlib.Path],
                 stations_file: Union[str, pathlib.Path] = None,
                 ):
        outputs = pathlib.Path(outputs)
        if not outputs.is_dir():
            raise Exception(f'{str(outputs)} is not a directory.')
        self._station_id: List = []
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
                            self._station_id.append(line.split('!')[-1])
                        else:
                            self._station_id.append(None)
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
                self._start_date = None
            else:
                self._start_date = datetime(start_year, start_month,
                                            start_day, int(start_hour), 0,
                                            tzinfo=timezone(
                                                timedelta(hours=-utc_start)))
        self._rndays = None
        for file in self._manifest:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                data = np.loadtxt(file)
                if len(w) != 0:
                    pass
                else:
                    self._rndays = timedelta(seconds=data[-1, 0])
                    break

        if self.rndays is None:
            raise Exception("Ouptut directory doesn't contain any station "
                            "output data.")

    def set_start_date(self, start_date: datetime):
        self._start_date = start_date

    def get_station_id(self, index):
        return self._station_id[index].strip()

    def plot(self, variable, station_index=None, show=False):
        filenames = [path.name for path in self._manifest]
        var_index = StationOutputIndex[
            StationOutputVariables(variable).name].value + 1

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            data = np.loadtxt(
                self._manifest[filenames.index(f'staout_{var_index}')])
            if len(w) != 0:
                raise AttributeError(f'Empty record for variable {variable}.')

        if self.start_date is not None:
            start_date = self.start_date
            dates = [start_date + timedelta(seconds=t) for t in data[:, 0]]
        else:
            dates = data[:, 0]

        coops = CoopsDataCollector()
        if station_index is None:
            for i in range(data.shape[1] - 2):
                plt.figure(i)
                if self._obs:
                    station_id = self.get_station_id(i)
                    if station_id is not None:
                        obs = coops.fetch(
                            station_id, variable, self.start_date,
                            self.rndays)
                        plt.plot(obs['datetime'], obs['values'])
                        plt.title(obs['name'])
                        # fig, ax = plt.subplots(figsize=(8, 3))
                        # import matplotlib.dates as mdates
                        # plt.plot_date(obs['datetime'], obs['values'], ls='solid', lw=1.5, aa=True,marker='None',color='g')
                        # ax = plt.gca()
                        # dayFmt = mdates.DateFormatter('%a-%b-%d')
                        # hrFmt = mdates.DateFormatter('%H:00')
                        # ax.xaxis.set_major_formatter(dayFmt)
                        # ax.xaxis.set_major_locator(mdates.DayLocator())
                        # ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[12]))
                        # ax.xaxis.set_minor_formatter(hrFmt)
                        # plt.grid(b=True, which='both', linestyle='dotted')
                        # plt.show()
                plt.plot(dates, data[:, i+1])

        else:
            plt.plot(data[:, 0], data[:, station_index])

        if show is True:
            plt.show()
        # import pdb; pdb.set_trace()
        # with open(self.manifest[filenames.index(f'staout_{index}')]) as f:
        #     for line in f:
        #         print(len(line.split()))
        #         exit()

    @property
    def start_date(self):
        return self._start_date

    @property
    def rndays(self):
        return self._rndays