from datetime import datetime, timedelta
from functools import partial
import logging
import os
import pathlib
from typing import Union

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
from pyproj import CRS
import pytz
from scipy.spatial import cKDTree
from shapely import ops
from shapely.geometry import Point


from pyschism.dates import localize_datetime


_logger = logging.getLogger(__name__)


class SourceSinkDataset:

    def __init__(self, data):
        self._data = data

    def __len__(self):
        return len(self.data)

    def __str__(self):
        return str(self.data)

    @property
    def elements(self):
        unique_elements = set()
        for elements in self.data.values():
            for element_id in elements.keys():
                unique_elements.add(element_id)
        return list(map(str, sorted(list(map(int, unique_elements)))))

    @property
    def timevector(self):
        return list(sorted(self.data.keys()))

    @property
    def df(self):
        if not hasattr(self, '_df'):
            data = []
            for time, edata in self.data.items():
                for eid, _edata in edata.items():
                    data.append({
                        'time': time,
                        'element_id': eid,
                        **_edata
                    })
            self._df = pd.DataFrame(data)
        return self._df

    @property
    def data(self):
        return self._data


class Sources(SourceSinkDataset):
    pass


class Sinks(SourceSinkDataset):
    pass


class TimeHistoryFile:

    def __init__(self, data, start_date, rnday, filename):
        self.data = data
        self.filename = filename
        self.start_date = start_date
        self.rnday = rnday

    def __str__(self):
        data = []
        for time in self.data.timevector:
            relative_time = (time - self.start_date).total_seconds()
            if relative_time < 0:
                continue
            line = [f"{relative_time:G}"]
            for element_id in self.data.elements:
                line.append(
                    f'{self.data._data[time][element_id]["flow"]:.4e}')
            data.append(' '.join(line))
        return '\n'.join(data)

    def write(self, path: Union[str, os.PathLike], overwrite: bool = False):
        path = pathlib.Path(path)
        if path.exists() and overwrite is not True:
            raise IOError('File exists and overwrite is not True.')
        open(path / self.filename, 'w').write(str(self))


class Vsource(TimeHistoryFile):

    def __init__(self, sources: Sources, start_date, rnday, filename='vsource.th'):
        super().__init__(sources, start_date, rnday, filename)


class Msource(TimeHistoryFile):

    def __init__(self, sources, start_date, rnday, filename='msource.th'):
        super().__init__(sources, start_date, rnday, filename)

    def __str__(self):
        data = []
        for time in self.data.timevector:
            relative_time = (time - self.start_date).total_seconds()
            if relative_time < 0:
                continue
            line = [f"{relative_time:G}"]
            for element_id in self.data.elements:
                line.append(
                    f'{self.data._data[time][element_id]["temperature"]: .4e}')
            for element_id in self.data.elements:
                line.append(
                    f'{self.data._data[time][element_id]["salinity"]: .4e}')
            data.append(' '.join(line))
        return '\n'.join(data)


class Vsink(TimeHistoryFile):

    def __init__(self, sinks: Sinks, start_date, rnday, filename='vsink.th'):
        super().__init__(sinks, start_date, rnday, filename)


class SourceSink:

    def __init__(
            self,
            sources: Sources,
            sinks: Sinks,
            filename='source_sink.in'
    ):
        self.sources = sources
        self.sinks = sinks
        self.filename = filename

    def __str__(self):
        source_id = self.sources.elements
        data = []
        data.append(f'{len(source_id)}')
        for element_id in source_id:
            data.append(f'{element_id}')
        data.append('')
        sink_id = self.sinks.elements
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

    # TODO: This class is missing a time interpolator.
    _data = {}

    def __init__(self, start_date: datetime = None, rnday: timedelta = None):
        self.start_date = start_date
        self.rnday = rnday

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
        if hasattr(self, '_df'):
            del self._df

    def get_element_timeseries(self, element_id):
        data = {}
        element_data = self.df[(self.df['element_id'] == element_id)]
        for row in element_data.sort_values(by=['time']).itertuples():
            data.setdefault(row.time, {}).update({
                'flow': row.flow,
                'temperature': row.temperature,
                'salinity': row.salinity})
        return data

    # def get_interpolated_timeseries(self, element_id):
    #     f = self.get_element_interpolator(element_id)
    #     self.timevector

    def remove_element_timeseries(self, element_id):
        for time in self._data:
            self._data[time].pop(element_id)
        if hasattr(self, '_df'):
            del self._df

    def aggregate_by_radius(self, hgrid, radius, nprocs=-1):

        # --- Generate aggregation mapping

        # gather extreme values
        source_max = {element_id: -float('inf') for element_id
                      in self.sources.elements}
        for element_data in self.sources.data.values():
            for element_id, data in element_data.items():
                source_max[element_id] = np.max(
                    [source_max[element_id], data['flow']])

        sink_max = {element_id: float('inf') for element_id
                    in self.sinks.elements}
        for element_data in self.sinks.data.values():
            for element_id, data in element_data.items():
                sink_max[element_id] = np.min(
                    [sink_max[element_id], data['flow']])

        aggregate_gdf = []
        for element_id, maxflow in {**source_max, **sink_max}.items():
            element_index = hgrid.elements.get_index_by_id(element_id)
            aggregate_gdf.append({
                'element_id': element_id,
                'geometry': hgrid.elements.gdf.loc[element_index].geometry,
                'maxflow': maxflow,

            })
        aggregate_gdf = gpd.GeoDataFrame(
            aggregate_gdf, crs=hgrid.crs).sort_values(
            by='maxflow', key=abs, ascending=False)

        aggregation_mapping = {}
        for row in aggregate_gdf.itertuples():
            if row.element_id in aggregation_mapping:
                continue
            other_sources = aggregate_gdf.loc[
                aggregate_gdf.index.difference([row.Index])]
            point = np.array(row.geometry.centroid)
            circle = get_circle_of_radius(point[0], point[1], radius)
            sources_in_circle = other_sources.loc[other_sources.within(circle)]
            aggregation_mapping[row.element_id] = row.element_id
            for row_in_circle in sources_in_circle.itertuples():
                aggregation_mapping[row_in_circle.element_id] = row.element_id

        # --- move data from one element to the other
        for current, target in aggregation_mapping.items():
            for time, data in self.get_element_timeseries(current).items():
                self.add_data(time, target, **data)
        for current, target in aggregation_mapping.items():
            if current != target:
                self.remove_element_timeseries(current)

        if hasattr(self, '_sources'):
            del self._sources

        if hasattr(self, '_sinks'):
            del self._sinks

        if hasattr(self, '_df'):
            del self._df

    def write(
            self,
            path: Union[str, os.PathLike],
            overwrite: bool = False,
            msource: Union[str, bool] = True,
            vsource: Union[str, bool] = True,
            vsink: Union[str, bool] = True,
            source_sink: Union[str, bool] = True,
    ):

        # unpack sources, sinks
        sources, sinks = self.sources, self.sinks

        # write source sink
        if source_sink is True:
            fname = 'source_sink.in'
        elif isinstance(source_sink, str):
            fname = source_sink
        if source_sink is not False:
            SourceSink(sources, sinks, fname).write(path, overwrite)

        if vsource is True:
            fname = 'vsource.th'
        elif isinstance(vsource, str):
            fname = vsource
        if vsource is not False:
            Vsource(
                sources,
                self.start_date,
                self.rnday,
                fname
            ).write(path, overwrite)

        if msource is True:
            fname = 'msource.th'
        elif isinstance(msource, str):
            fname = msource
        if msource is not False:
            Msource(
                sources,
                self.start_date,
                self.rnday,
                fname
            ).write(path, overwrite)

        if vsink is True:
            fname = 'vsink.th'
        elif isinstance(vsink, str):
            fname = vsink
        if vsink is not False:
            Vsink(
                sinks,
                self.start_date,
                self.rnday,
                fname
            ).write(path, overwrite)

    @property
    def sources(self):
        if not hasattr(self, '_sources'):
            sources = {}
            for element_id in list(map(str, list(
                    sorted(map(int, self.df.element_id.unique()))))):
                element_data = self.df[(self.df['element_id'] == element_id)]
                flow_data = element_data['flow']
                if np.all(flow_data > 0.):
                    # TODO:  Are irregular timeseries allowed?
                    # if not, we need an interpolator here.
                    for row in element_data.sort_values(
                            by=['time']).itertuples():
                        sources.setdefault(row.time, {}).setdefault(
                            element_id, {}).update({
                                'flow': row.flow,
                                'temperature': row.temperature,
                                'salinity': row.salinity})
                # handle elements that are both sources and sinks
                elif not np.all(flow_data < 0) and np.any(flow_data > 0.):
                    for row in element_data.sort_values(
                            by=['time']).itertuples():
                        flow = row.flow if row.flow >= 0. else 0.
                        sources.setdefault(row.time, {}).setdefault(
                            element_id, {}).update({
                                'flow': flow,
                                'temperature': row.temperature,
                                'salinity': row.salinity})
            self._sources = Sources(sources)
        return self._sources

    @property
    def sinks(self):
        if not hasattr(self, '_sinks'):
            sinks = {}
            for element_id in list(map(str, list(
                    sorted(map(int, self.df.element_id.unique()))))):
                element_data = self.df[(self.df['element_id'] == element_id)]
                flow_data = element_data['flow']
                if np.all(flow_data < 0.):
                    # TODO:  Are irregular timeseries allowed?
                    # if not, we need an ingerpolator here.
                    for row in element_data.sort_values(
                            by=['time']).itertuples():
                        sinks.setdefault(row.time, {}).setdefault(
                            element_id, {}).update({
                                'flow': row.flow})
                # handle elements that are both sources and sinks
                elif not np.all(flow_data > 0.) and np.any(flow_data < 0.):
                    for row in element_data.sort_values(
                            by=['time']).itertuples():
                        sinks.setdefault(row.time, {}).setdefault(
                            element_id, {}).update(
                                {'flow': row.flow if row.flow <= 0. else 0.})
            self._sinks = Sinks(sinks)
        return self._sinks

    @property
    def start_date(self):
        if self._start_date is None:
            return self.df.time.min()
        return self._start_date

    @start_date.setter
    def start_date(self, start_date):
        self._start_date = start_date
        if start_date is not None:
            self._start_date = localize_datetime(
                start_date).astimezone(pytz.utc)
        return self._start_date

    @property
    def rnday(self):
        if self._rnday is None:
            return self.df.time.max() - self.df.time.min()
        return self._rnday

    @rnday.setter
    def rnday(self, rnday):
        self._rnday = rnday
        if rnday is not None:
            self._rnday = rnday if isinstance(rnday, timedelta) else \
                timedelta(days=rnday)
        return self._rnday

    @property
    def data(self):
        return self._data

    @property
    def df(self):
        if not hasattr(self, '_df'):
            _data = []
            for time, element_data in self._data.items():
                for element_id, data in element_data.items():
                    _data.append({
                        'time': time,
                        'element_id': element_id,
                        **data
                    })
            self._df = pd.DataFrame(_data)
        return self._df


def get_local_azimutal_projection(lon, lat):
    return CRS(f"+proj=aeqd +R=6371000 +units=m +lat_0={lat} +lon_0={lon}")


def get_circle_of_radius(lon, lat, radius):

    local_azimuthal_projection = "+proj=aeqd +R=6371000 +units=m " \
                                 f"+lat_0={lat} +lon_0={lon}"
    wgs84_to_aeqd = partial(
        pyproj.transform,
        pyproj.Proj("+proj=longlat +datum=WGS84 +no_defs"),
        pyproj.Proj(local_azimuthal_projection),
    )
    aeqd_to_wgs84 = partial(
        pyproj.transform,
        pyproj.Proj(local_azimuthal_projection),
        pyproj.Proj("+proj=longlat +datum=WGS84 +no_defs"),
    )

    center = Point(float(lon), float(lat))
    point_transformed = ops.transform(wgs84_to_aeqd, center)

    return ops.transform(aeqd_to_wgs84, point_transformed.buffer(radius))
