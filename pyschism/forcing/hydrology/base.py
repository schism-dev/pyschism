import logging
import os
import pathlib
from typing import Union

import pandas as pd


_logger = logging.getLogger(__name__)


class TimeHistoryFile:

    def __init__(self, df, start_date, rnday, filename):
        self.df = df
        self.filename = filename
        self.start_date = start_date
        self.rnday = rnday

    def write(self, path: Union[str, os.PathLike], overwrite: bool = False):
        elements = self.df['element_id'].unique()
        times = self.df['time'].unique()
        with open(pathlib.Path(path) / self.filename, 'w') as f:
            for time in times:
                relative_time = (time - self.start_date).total_seconds()
                if relative_time < 0:
                    continue
                line = [f"{relative_time:G}"]
                for element_id in elements:
                    pd_series = self.df["flow"][
                        (self.df["time"] == time) &
                        (self.df["element_id"] == element_id)]
                    line.append(f'{pd_series.values[0]:.4e}')
                f.write(' '.join(line))
                f.write('\n')


class Vsource(TimeHistoryFile):

    def __init__(self, df, start_date, rnday, filename='vsource.th'):
        super().__init__(df, start_date, rnday, filename)


class Vsink(TimeHistoryFile):

    def __init__(self, df, start_date, rnday, filename='vsink.th'):
        super().__init__(df, start_date, rnday, filename)


class Msource(TimeHistoryFile):

    def __init__(self, df, start_date, rnday, filename='msource.th'):
        super().__init__(df, start_date, rnday, filename)

    def write(self, path: Union[str, os.PathLike], overwrite: bool = False):
        elements = self.df['element_id'].unique()
        times = self.df['time'].unique()
        with open(pathlib.Path(path) / self.filename, 'w') as f:
            for time in times:
                relative_time = (time - self.start_date).total_seconds()
                if relative_time < 0:
                    continue
                line = [f"{relative_time:G}"]
                for element_id in elements:
                    pd_series = self.df["temperature"][
                        (self.df["time"] == time) &
                        (self.df["element_id"] == element_id)]
                    line.append(f'{pd_series.values[0]:.4e}')
                for element_id in elements:
                    pd_series = self.df["salinity"][
                        (self.df["time"] == time) &
                        (self.df["element_id"] == element_id)]
                    line.append(f'{pd_series.values[0]:.4e}')
                f.write(' '.join(line))
                f.write('\n')


class SourcesDataFrame:
    '''This descriptor holds the shared dataframe for the API components.'''

    def __get__(self, obj, val):
        if not hasattr(self, '_df'):
            self._df = pd.DataFrame(
                columns=[
                    'element_id',
                    'time',
                    'flow',
                    'temperature',
                    'salinity'
                ])
        return self._df

    def __set__(self, obj, df: pd.DataFrame):
        self._df = df


class Sources:
    '''Public interface used for adding sources to the model.'''

    _df = SourcesDataFrame()

    def __len__(self):
        return len(self.df)

    def add_data(self, time, element_id, flow, temperature, salinity):

        # TODO: sources/sinks.add_data is slow, it should be vectorized
        # main reason is because it's tallying multiple datapoints that
        # corresponds to the same element, so it has to check the dataframe
        # for each datapoint added.

        # TODO: Consider additional inputs besides TS?
        # Check if datapoint already exists (auto-tally)
        tally = ((self.df['time'] == time) &
                 (self.df['element_id'] == element_id))
        if tally.any():
            # previous_flow = self.df.loc[tally, 'flow']
            self.df.loc[tally, 'flow'] += flow
            # TODO: What happens if we have two different flows that both are
            # assigned to the same element? Example: 100 m^3/s @ 1 psu then
            # another flow on the same element of 1 m^3/s @ 100 psu. How do we
            # combine these on a single element? Flow is just simple summation,
            # but TS is different.
            if (self.df.loc[tally, 'salinity'] != salinity).any():
                raise NotImplementedError('Two different values of salinity '
                                          'for same time/element')

            if (self.df.loc[tally, 'temperature'] != temperature).any():
                raise NotImplementedError('Two different values of temperature'
                                          ' for same time/element')

        else:
            self._df = self.df.append(
                pd.Series(
                    [element_id, time, flow, temperature, salinity],
                    index=self.df.columns,
                ),
                ignore_index=True,
                # TODO: in_place=True, ?
            )

    @property
    def df(self):
        return self._df


class SinksDataFrame:
    '''This descriptor holds the shared dataframe for the API components.'''

    def __get__(self, obj, val):
        if not hasattr(self, '_df'):
            self._df = pd.DataFrame(
                columns=[
                    'element_id',
                    'time',
                    'flow'
                ])
        return self._df

    def __set__(self, obj, df: pd.DataFrame):
        self._df = df


class Sinks:

    _df = SinksDataFrame()

    def __len__(self):
        return len(self.df)

    def add_data(self, time, element_id, flow):
        tally = ((self.df['time'] == time) &
                 (self.df['element_id'] == element_id))
        if tally.any():
            self.df.loc[tally, 'flow'] += flow
        else:
            self._df = self.df.append(
                pd.Series(
                    [element_id, time, flow],
                    index=self.df.columns,
                ),
                ignore_index=True,
            )

    @property
    def df(self):
        return self._df


class SourceSink:

    def __init__(
            self,
            sources: Sources,
            sinks: Sinks,
            filename='source_sink.in'
    ):
        self._sources = sources
        self._sinks = sinks
        self.filename = filename

    def write(self, path: Union[str, os.PathLike], overwrite: bool = False):
        sources = self.sources.df['element_id'].unique()
        sinks = self.sinks.df['element_id'].unique()
        data = []
        with open(pathlib.Path(path) / self.filename, 'w') as f:
            data.append(f'{len(sources)}')
            for element_id in sources:
                data.append(f'{element_id}')
            data.append('\n')
            data.append(f'{len(sinks)}')
            for element_id in sinks:
                data.append(f'{element_id}')
            f.write('\n'.join(data))

    @property
    def sources(self):
        return self._sources

    @property
    def sinks(self):
        return self._sinks


class Hydrology:

    def __init__(self):
        self._sources = Sources()
        self._sinks = Sinks()

    def __call__(self, model_driver):
        self._model_driver = model_driver

    def write(
            self,
            path: Union[str, os.PathLike],
            overwrite: bool = False,
            msource: Union[str, bool] = True,
            vsource: Union[str, bool] = True,
            vsink: Union[str, bool] = True,
            source_sink: Union[str, bool] = True,
    ):

        start_date = self.model_driver.param.opt.start_date
        rnday = self.model_driver.param.core.rnday

        # write source sink
        if source_sink is True:
            fname = 'source_sink.in'
        elif isinstance(source_sink, str):
            fname = source_sink
        if source_sink is not False:
            SourceSink(self.sources, self.sinks, fname).write(path, overwrite)

        if vsource is True:
            fname = 'vsource.th'
        elif isinstance(vsource, str):
            fname = vsource
        if vsource is not False:
            Vsource(
                self.sources.df, start_date, rnday, fname).write(
                path, overwrite)

        if msource is True:
            fname = 'msource.th'
        elif isinstance(msource, str):
            fname = msource
        if msource is not False:
            Msource(
                self.sources.df, start_date, rnday, fname).write(
                path, overwrite)

        if vsink is True:
            fname = 'vsink.th'
        elif isinstance(vsink, str):
            fname = vsink
        if vsink is not False:
            Vsink(
                self.sinks.df, start_date, rnday, fname).write(
                path, overwrite)

    @property
    def sources(self):
        return self._sources

    @property
    def sinks(self):
        return self._sinks

    @property
    def model_driver(self):
        return self._model_driver
