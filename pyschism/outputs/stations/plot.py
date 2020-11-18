from datetime import timezone, timedelta
import warnings

import matplotlib.pyplot as plt
import numpy as np

from pyschism.outputs import StationsOutput
from pyschism.utils.coops import CoopsDataCollector
from pyschism.enums import StationOutputIndex, StationOutputVariables


class PlotOutputStations:

    def __init__(self, stations: StationsOutput, obs=True):
        self.__stations = stations
        self._obs = obs
        if self._obs:
            self._coops = CoopsDataCollector()

    def plot(self, variable, station_index=None):
        filenames = [path.name for path in self.stations._manifest]
        var_index = StationOutputIndex[
            StationOutputVariables(variable).name].value + 1

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            data = np.loadtxt(
                self.stations._manifest[filenames.index(
                    f'staout_{var_index}')])
            if len(w) != 0:
                raise AttributeError(f'Empty record for variable {variable}.')

        if self.stations.start_date is not None:
            start_date = self.stations.start_date
            dates = [start_date + timedelta(seconds=t) for t in data[:, 0]]
        else:
            dates = data[:, 0]

        if station_index is None:
            for i in range(data.shape[1] - 2):
                plt.figure(i)
                if self._obs:
                    station_id = self.stations.get_station_id(i)
                    if station_id is not None:
                        obs = self._coops.fetch(
                            station_id, variable, self.stations.start_date,
                            self.stations.rndays)
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
        plt.show()
        # import pdb; pdb.set_trace()
        # with open(self.manifest[filenames.index(f'staout_{index}')]) as f:
        #     for line in f:
        #         print(len(line.split()))
        #         exit()

    @property
    def stations(self):
        return self.__stations
