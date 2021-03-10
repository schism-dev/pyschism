from datetime import datetime
import logging
import os
import pathlib
from typing import Union

import numpy as np
import pytz

from pyschism.dates import localize_datetime


_logger = logging.getLogger(__name__)


class Sources:
    '''Public interface used for adding sources to the model.'''

    _data = {}

    def __len__(self):
        return len(self._data)

    def add_data(
            self,
            time: datetime,
            element_id: str,
            flow: float,
            temperature: float = np.nan,
            salinity: float = np.nan,
    ):

        time = localize_datetime(time).astimezone(pytz.utc)
        data_for_element = self._data.get(time, {}).get('element_id', {})

        # TODO: What happens if we have two different flows that both are
        # assigned to the same element? Example: 100 m^3/s @ 1 psu then
        # another flow on the same element of 1 m^3/s @ 100 psu. How do we
        # combine these on a single element? Flow is just simple summation,

        _tmp = data_for_element.get('temperature', np.nan)
        if not np.isnan(_tmp):
            if _tmp != temperature:
                raise NotImplementedError(
                    'Two different values of temperature for same '
                    'time/element.')

        _sal = data_for_element.get('salinity', np.nan)
        if not np.isnan(_sal):
            if _sal != salinity:
                raise NotImplementedError(
                    'Two different values of salinity for same time/element.')

        self._data.setdefault(time, {}).setdefault(element_id, {}).update(
            {
                'flow': np.nansum(
                    [data_for_element.get('flow', np.nan), flow]),
                'temperature': temperature,
                'salinity': salinity
            })

    @property
    def elements(self):
        unique_elements = set()
        for elements in self._data.values():
            for element_id in elements.keys():
                unique_elements.add(element_id)
        return list(map(str, sorted(list(map(int, unique_elements)))))

    @property
    def timevector(self):
        return list(sorted(self._data.keys()))


class Sinks:

    _data = {}

    def __len__(self):
        return len(self._data)

    def add_data(
            self,
            time: datetime,
            element_id: str,
            flow: float,
    ):

        time = localize_datetime(time).astimezone(pytz.utc)
        data_for_element = self._data.get(time, {}).get('element_id', {})
        if flow > 0:
            raise ValueError('Argument flow for must be smaller than 0.')
        self._data.setdefault(time, {}).setdefault(element_id, {}).update(
            {
                'flow': np.nansum(
                    [data_for_element.get('flow', np.nan), flow]),
            })

    @property
    def elements(self):
        unique_elements = set()
        for elements in self._data.values():
            for element_id in elements.keys():
                unique_elements.add(element_id)
        return list(map(str, sorted(list(map(int, unique_elements)))))

    @property
    def timevector(self):
        return list(sorted(self._data.keys()))


class TimeHistoryFile:

    def __init__(self, start_date, rnday, filename):

        self.filename = filename
        self.start_date = start_date
        self.rnday = rnday

    def __str__(self):
        data = []
        for time in self._source_sink.timevector:
            relative_time = (time - self.start_date).total_seconds()
            if relative_time < 0:
                continue
            line = [f"{relative_time:G}"]
            for element_id in self._source_sink.elements:
                line.append(
                    f'{self._source_sink._data[time][element_id]["flow"]:.4e}')
            data.append(' '.join(line))
        return '\n'.join(data)

    def write(self, path: Union[str, os.PathLike], overwrite: bool = False):
        path = pathlib.Path(path)
        if path.exists() and overwrite is not True:
            raise IOError('File exists and overwrite is not True.')
        open(path / self.filename, 'w').write(str(self))


class Vsource(TimeHistoryFile):

    _sources = Sources()

    def __init__(self, start_date, rnday, filename='vsource.th'):
        super().__init__(start_date, rnday, filename)
        self._source_sink = self._sources


class Vsink(TimeHistoryFile):

    _sinks = Sinks()

    def __init__(self, start_date, rnday, filename='vsink.th'):
        super().__init__(start_date, rnday, filename)
        self._source_sink = self._sinks


class Msource(TimeHistoryFile):

    _sources = Sources()

    def __init__(self, start_date, rnday, filename='msource.th'):
        super().__init__(start_date, rnday, filename)
        self._source_sink = self._sources

    def __str__(self):
        data = []
        for time in self._source_sink.timevector:
            relative_time = (time - self.start_date).total_seconds()
            if relative_time < 0:
                continue
            line = [f"{relative_time:G}"]
            for element_id in self._source_sink.elements:
                line.append(
                    f'{self._source_sink._data[time][element_id]["temperature"]: .4e}')
            for element_id in self._source_sink.elements:
                line.append(
                    f'{self._source_sink._data[time][element_id]["salinity"]: .4e}')
            data.append(' '.join(line))
        return '\n'.join(data)


class SourceSink:

    _sources = Sources()
    _sinks = Sinks()

    def __init__(
            self,
            filename='source_sink.in'
    ):
        self.filename = filename

    def __str__(self):
        source_id = self._sources.elements
        data = []
        data.append(f'{len(source_id)}')
        for element_id in source_id:
            data.append(f'{element_id}')
        data.append('')
        sink_id = self._sinks.elements
        data.append(f'{len(sink_id)}')
        for element_id in sink_id:
            data.append(f'{element_id}')
        return '\n'.join(data)

    def write(self, path: Union[str, os.PathLike], overwrite: bool = False):
        path = pathlib.Path(path)
        if path.exists() and overwrite is not True:
            raise IOError('File exists and overwrite is not True.')
        open(path / self.filename, 'w').write(str(self))


class Hydrology:

    _sources = Sources()
    _sinks = Sinks()

    def __call__(self, model_driver):
        self.set_start_date(model_driver.param.opt.start_date)
        self.set_rnday(model_driver.param.opt.start_date)

    def set_start_date(self, start_date):
        self._start_date = start_date

    def set_rnday(self, rnday):
        self._rnday = rnday

    def write(
            self,
            path: Union[str, os.PathLike],
            overwrite: bool = False,
            msource: Union[str, bool] = True,
            vsource: Union[str, bool] = True,
            vsink: Union[str, bool] = True,
            source_sink: Union[str, bool] = True,
    ):

        # write source sink
        if source_sink is True:
            fname = 'source_sink.in'
        elif isinstance(source_sink, str):
            fname = source_sink
        if source_sink is not False:
            SourceSink(
                fname
            ).write(path, overwrite)

        if vsource is True:
            fname = 'vsource.th'
        elif isinstance(vsource, str):
            fname = vsource
        if vsource is not False:
            Vsource(
                self._start_date,
                self._rnday,
                fname
            ).write(path, overwrite)

        if msource is True:
            fname = 'msource.th'
        elif isinstance(msource, str):
            fname = msource
        if msource is not False:
            Msource(
                self._start_date,
                self._rnday,
                fname
            ).write(path, overwrite)

        if vsink is True:
            fname = 'vsink.th'
        elif isinstance(vsink, str):
            fname = vsink
        if vsink is not False:
            Vsink(
                self._start_date,
                self._rnday,
                fname
            ).write(path, overwrite)
