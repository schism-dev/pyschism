import calendar
from datetime import datetime, timedelta, timezone
import json
from typing import Union

import numpy as np
import requests


class CoopsDataCollector:

    url = "https://tidesandcurrents.noaa.gov/api/datagetter?"

    def fetch(self, station_id: str,  variable: str, start_date: datetime,
              rndays: Union[float, timedelta], datum: str = 'NAVD',
              units: str = 'metric'):

        # check start_date
        if not isinstance(start_date, datetime):
            raise TypeError(f'start_date must be a {datetime} instance.')

        # make date GMT
        if start_date.tzinfo is not None and \
                start_date.tzinfo.utcoffset(start_date) is not None:
            tzinfo = start_date.tzinfo
            start_date = start_date.astimezone(timezone(timedelta(0)))
        else:
            tzinfo = False

        # check rndays
        if not isinstance(rndays, timedelta):
            rndays = timedelta(days=rndays)

        end_date = start_date + rndays

        # set datum
        if datum not in ['MHHW', 'MHW', 'MTL', 'MSL', 'MLW', 'MLLW', 'NAVD88',
                         'STND', 'NAVD']:
            raise AttributeError(f'Invalid datum: {datum}.')
        if datum == 'NAVD88':
            datum = 'NAVD'

        responses = list()
        for start_date, end_date in self._iter_datetime_segments(
                start_date, end_date):
            params = self._get_params(station_id, start_date, end_date, datum,
                                      units)
            try:
                r = requests.get(self.url, params=params, timeout=10.)
                r.raise_for_status()
            except requests.exceptions.HTTPError:
                pass  # this means no data avilable
            except requests.exceptions.ConnectionError as errc:
                print("Error Connecting:", errc)
            except requests.exceptions.Timeout as errt:
                print("Timeout Error:", errt)
            except requests.exceptions.RequestException as err:
                print("Unknown error.", err)
            responses.append(r)
        data = dict()
        data['datetime'] = list()
        data['values'] = list()
        for i, response in enumerate(responses):
            json_data = json.loads(response.text)
            if 'error' in json_data.keys():
                _start_date, _end_date = list(
                    self._get_datetime_segments(start_date, end_date))[i]
                data['datetime'].append(_start_date)
                data['values'].append(np.nan)
                data['datetime'].append(_end_date)
                data['values'].append(np.nan)
                continue
            if 'x' not in data.keys():
                data['x'] = float(json_data['metadata']['lon'])
            if 'y' not in data.keys():
                data['y'] = float(json_data['metadata']['lat'])
            if 'name' not in data.keys():
                data['name'] = json_data['metadata']['name']
            for _data in json_data['data']:
                data['datetime'].append(
                    datetime.strptime(_data['t'], '%Y-%m-%d %H:%M'))
                try:
                    data['values'].append(float(_data['v']))
                except ValueError:
                    data['values'].append(np.nan)
        if 'name' not in data.keys():
            data['name'] = ''
        if np.all(data['values'] == np.nan):
            return None
        # change datetimes back to original timezone if tz-aware
        if tzinfo:
            data['datetime'] = [
                dt.replace(tzinfo=timezone.utc).astimezone(tzinfo)
                for dt in data['datetime']]
        return data

    def _get_params(self, station_id, start_date, end_date, datum, units):
        params = {}
        params['station'] = station_id
        params['begin_date'] = start_date.strftime('%Y%m%d %H:%M')
        params['end_date'] = end_date.strftime('%Y%m%d %H:%M')
        params['product'] = 'water_level'
        params['datum'] = datum
        params['units'] = units
        params['time_zone'] = 'gmt'
        params['format'] = 'json'
        params['application'] = 'PySCHISM'
        return params

    def _iter_datetime_segments(self, start_date, end_date):
        """
        https://www.ianwootten.co.uk/2014/07/01/splitting-a-date-range-in-python/
        """

        def get_datespan(interval):
            start_epoch = calendar.timegm(start_date.timetuple())
            end_epoch = calendar.timegm(end_date.timetuple())
            date_diff = end_epoch - start_epoch
            step = date_diff / interval
            delta = timedelta(seconds=step)
            current_date = start_date
            while current_date + delta <= end_date:
                to_date = (current_date + delta)
                yield current_date, to_date
                current_date += delta

        segments = [(start_date, end_date)]
        interval = 2
        while np.any([(_end_date - _start_date).total_seconds()
                      > timedelta(days=31).total_seconds()
                      for _start_date, _end_date in segments]):
            segments = [(from_datetime, to_datetime)
                        for from_datetime, to_datetime
                        in get_datespan(interval)]
            interval += 1
        for _start_date, _end_date in segments:
            yield _start_date, _end_date
