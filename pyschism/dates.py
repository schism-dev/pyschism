from datetime import datetime, timedelta
from typing import Union

import numpy as np
import pytz


# def singleton(class_):
#     instances = {}

#     def getinstance(*args, **kwargs):
#         if class_ not in instances:
#             instances[class_] = class_(*args, **kwargs)
#         return instances[class_]
#     return getinstance


# @singleton
class StartDate:

    def __init__(self):
        self.start_date = None

    def __set__(self, obj, val: datetime):
        self.start_date = localize_datetime(val).astimezone(pytz.utc)

    def __get__(self, obj, val) -> datetime:
        return self.start_date

    def __delete__(self, obj):
        self.start_date = None


# @singleton
class EndDate:

    def __init__(self):
        self.end_date = None

    def __set__(self, obj, val: Union[float, timedelta, datetime]):

        if isinstance(val, datetime):
            val = localize_datetime(val).astimezone(pytz.utc)

        elif not isinstance(val, timedelta):
            val = obj.start_date + timedelta(days=float(val))

        elif isinstance(val, timedelta):
            val = obj.start_date + val

        self.end_date = val

    def __get__(self, obj, val) -> datetime:
        return self.end_date

    def __delete__(self, obj):
        self.end_date = None


class SpinupTime:

    def __init__(self):
        self.spinup_time = None

    def __set__(self, obj, val: Union[int, float, timedelta]):
        if not isinstance(val, timedelta):
            val = timedelta(days=float(val))
        self.spinup_time = val

    def __get__(self, obj, val) -> timedelta:
        return self.spinup_time

    def __delete__(self, obj):
        self.spinup_time = None


def nearest_zulu(input_datetime=None):
    """
    "pivot time" is defined as the nearest floor t00z for any given datetime.
    If this function is called without arguments, it will return the pivot time
    for the current datetime in UTC.
    """
    input_datetime = nearest_cycle() if input_datetime is None else \
        localize_datetime(input_datetime).astimezone(pytz.utc)
    return localize_datetime(
        datetime(input_datetime.year, input_datetime.month, input_datetime.day)
    )


def nearest_cycle(input_datetime=None, period=6):
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


def utcnow():
    return localize_datetime(datetime.utcnow())
