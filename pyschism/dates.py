from datetime import datetime

import numpy as np
import pytz


def pivot_time(input_datetime=None):
    """
    "pivot time" is defined as the nearest floor t00z for any given datetime.
    If this function is called without arguments, it will return the pivot time
    for the current datetime in UTC.
    """
    input_datetime = nearest_cycle_date() if input_datetime is None else \
        localize_datetime(input_datetime).astimezone(pytz.utc)
    return localize_datetime(
        datetime(input_datetime.year, input_datetime.month, input_datetime.day)
        )


def nearest_cycle_date(input_datetime=None, period=6):
    if input_datetime is None:
        input_datetime = localize_datetime(datetime.utcnow())
    current_cycle = int(period * np.floor(input_datetime.hour / period))
    return pytz.timezone('UTC').localize(
        datetime(input_datetime.year, input_datetime.month,
                 input_datetime.day, current_cycle))


def localize_datetime(d):
    # datetime is naïve iff:
    if d.tzinfo is None or d.tzinfo.utcoffset(d) is None:
        return pytz.timezone('UTC').localize(d)
    return d
