from datetime import datetime, timedelta
import pickle
import copy
from numbers import Integral, Real

from pandas import read_csv, to_datetime, to_numeric
import pandas as pd
from matplotlib import pyplot
import numpy

def running_mean1(X, N):
    Y = X * numpy.nan
    if (X.ndim > 1):  # multi-dimensional
        for i, x in enumerate(X.transpose()):
            Y[:, i] = numpy.convolve(x, numpy.ones((N,))/N)[(N-1):]
    else:  # 1-d array
        Y = numpy.convolve(X, numpy.ones((N,))/N)[(N-1):]

    return Y

def running_mean(X, N):
    cumsum = numpy.cumsum(numpy.insert(X, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)


class TimeHistory:


    """Class for manipulating *.th"""

    def __init__(self, file_name=None, start_time_str="2000-01-01 00:00:00", mask_val=None, sec_per_time_unit=1.0, data_array=None, columns=None):

        if mask_val is None:
            mask_val = -9999999
        self.mask_val = mask_val

        """initialize """
        self.source_file = file_name
        self.start_time_str = start_time_str
        self.n_time = -99999999
        self.n_station = -99999999
        self.time = None
        self.data = None
        self.delta_t = -99999999.9
        self.datetime = None
        self.sec_per_time_unit = sec_per_time_unit

        if file_name is None:
            self.df = pd.DataFrame(data_array)
        else:
            self.df = read_csv(file_name, delim_whitespace=True,  index_col=False, header=None)
        if columns is not None:
            self.df.columns = columns

        self.df.rename(columns={0: 'datetime'}, inplace=True)
        self.df_propagate()  # populate vars from dataframe

    def df_propagate(self):
        [self.n_time, self.n_station] = list(self.df.shape)
        self.n_station -= 1  # first col is time
        if isinstance(self.df['datetime'][0], (Real, Integral)):
            self.time = self.df['datetime'].to_numpy()
            self.datetime = [datetime.strptime(self.start_time_str, '%Y-%m-%d %H:%M:%S')
                             + timedelta(seconds=x*self.sec_per_time_unit) for x in self.time]
            self.df['datetime'] = self.datetime
        elif isinstance(self.df['datetime'][0], datetime):
            time_delta = self.df['datetime'] - self.df['datetime'][0]
            self.time = numpy.array([x.total_seconds() for x in time_delta])
            self.datetime = self.df['datetime']
        else:
            raise Exception('unknown datatype for the first column, neither datetime or float')


        # mask invalid values
        self.data = self.df.iloc[:, 1:].to_numpy(dtype=float)
        self.data[abs(self.data-float(self.mask_val)) < 1e-5] = numpy.nan
        print("Number of times: " + str(self.n_time))
        print("Number of stations: " + str(self.n_station))

        self.delta_t = self.time[1] - self.time[0]
        if self.delta_t < 50:  # probably in days
            self.time *= 86400.0
            self.delta_t *= 86400.0
            print("Time unit converted from days to seconds")
        else:
            print("Time unit: seconds")

    def get_running_mean(self, station_idx, window):
        """ sort time series by time average"""
        if station_idx is None:
            valid_idx = range(0, self.n_station)
        else:
            valid_idx = station_idx
        print(valid_idx)

        self.data = self.df.iloc[:, 1:].to_numpy(dtype=float)
        self.data[abs(self.data-float(self.mask_val)) < 1e-5] = numpy.nan
        data_rm = running_mean(self.data[:, valid_idx], window)

        return [data_rm]

    def get_time_average(self, station_idx, start_time_str=None, end_time_str=None):
        if station_idx == []:
            valid_idx = range(0, self.n_station)
        else:
            valid_idx = station_idx
        print(valid_idx)

        if start_time_str is None:
            idx1 = 0
        else:
            idx1 = self.get_snapshot(start_time_str)

        if end_time_str is None:
            idx2 = self.n_time
        else:
            idx2 = self.get_snapshot(end_time_str)

        self.data = self.df.iloc[:, 1:].to_numpy(dtype=float)
        self.data[abs(self.data-float(self.mask_val)) < 1e-5] = numpy.nan
        data_mean = numpy.mean(self.data[idx1:idx2, valid_idx], axis=0)

        return data_mean

    def sort_time_average(self, station_idx):
        """ sort time series by time average"""
        if station_idx == []:
            valid_idx = range(0, self.n_station)
        else:
            valid_idx = station_idx
        print(valid_idx)

        self.data = self.df.iloc[:, 1:].to_numpy(dtype=float)
        self.data[abs(self.data-float(self.mask_val)) < 1e-5] = numpy.nan
        data_mean = numpy.mean(self.data[:, valid_idx], axis=0)
        id_sort = numpy.argsort(data_mean)
        print("Preview on sorted: ")
        print(data_mean[id_sort[-1:-11:-1]])

        return [id_sort, data_mean]

    def get_snapshot(self, this_time_str=None, this_time_s_1970=None):
        """ get a snapshot closest to the input time"""
        if this_time_str is None:
            raise Exception("needs this_time_str to get snapshot")
        this_time_sec = timedelta(datetime.strftime(this_time_str, "%Y-%m-%d %H:%M:%S")
                                  - self.datetime[0]).seconds

        t_idx = numpy.argmin(numpy.abs(self.time - this_time_sec))
        return(t_idx)

    def get_max_station(self):
        """ get max among all stations """
        self.data = self.df.iloc[:, 1:].to_numpy(dtype=float)
        self.data[abs(self.data-float(self.mask_val)) < 1e-5] = numpy.nan
        data_sum = numpy.sum(self.data, axis=0)
        max_idx = numpy.argmax(data_sum)
        print("Station with the largest time-average value: " + str(max_idx))

        return max_idx

    def simple_plot(self, idx, **args):
        self.data = self.df.iloc[:, 1:].to_numpy(dtype=float)
        self.data[abs(self.data-float(self.mask_val)) < 1e-5] = numpy.nan
        pyplot.plot(self.datetime, self.data[:, idx], **args)

    def time_unit_conversion(self, unit_conv):
        """convert time unit"""
        """Col 0 is time: unit_conv from sec to day is 1/86400"""
        self.time = self.time * unit_conv

    def data_unit_conversion(self, unit_conv):
        """convert data unit"""
        self.data = self.df.iloc[:, 1:].to_numpy(dtype=float)
        self.data[abs(self.data-float(self.mask_val)) < 1e-5] = numpy.nan
        self.data = self.data * unit_conv

    def export_subset(self, station_idx=None, time_idx=None, i_reset_time=False, subset_filename=None):
        import os
        """extract a subset from the original *.th"""
        """by station_idx and time_idx"""
        """[] idx means no operation"""
        """if i_reset_time == True, then the subset *.th starts from 0 time"""

        if subset_filename is None:
            subset_filename = os.path.abspath(self.source_file) + '.subset'

        if station_idx is None:
            station_idx = numpy.array(range(self.n_station))
        else:
            station_idx = numpy.array(station_idx)

        if time_idx is None:
            time_idx = numpy.array(range(self.n_time))
        else:
            time_idx = numpy.array(time_idx)

        self.data = self.df.iloc[:, 1:].to_numpy(dtype=float)
        self.data[abs(self.data-float(self.mask_val)) < 1e-5] = numpy.nan

        sub_data = copy.deepcopy(self.data)
        sub_time = copy.deepcopy(self.time)
        if (len(station_idx) > 0):
            sub_data = sub_data[:, station_idx]
        if (len(time_idx) > 0):
            sub_data = sub_data[time_idx, :]
            sub_time = self.time[time_idx]

        if i_reset_time:
            new_start_time_str = datetime.strftime(datetime.strptime(self.start_time_str, "%Y-%m-%d %H:%M:%S") +timedelta(seconds=sub_time[0]*self.sec_per_time_unit), "%Y-%m-%d %H:%M:%S")
            sub_time = sub_time - sub_time[0]
        else:
            new_start_time_str = self.start_time_str

        with open(subset_filename, 'w') as fout:
            for i, _ in enumerate(sub_time):
                fout.write(str(sub_time[i]) + " " +
                           ' '.join(map(str, sub_data[i, :])) +
                           "\n")
        
        return TimeHistory(file_name=subset_filename, start_time_str=new_start_time_str, mask_val=self.mask_val, sec_per_time_unit=self.sec_per_time_unit,
                           data_array=None, columns=None)

    def writer(self, out_file_name):
        self.df_propagate()
        """ write self """
        with open(out_file_name, 'w') as fout:
            for i, _ in enumerate(self.time):
                fout.write(str(self.time[i]) + " " +
                           ' '.join(map(str, self.data[i, :])) +
                           "\n")

    @classmethod
    def timestamp_microsecond(cls, epoch, utc_time):
        """ convert date to microseconds """
        time_diff = utc_time - epoch
        assert time_diff.resolution == timedelta(microseconds=1)
        return (time_diff.days * 86400 + time_diff.seconds) * 10**6 + time_diff.microseconds
