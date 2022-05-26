from datetime import datetime
import io
import logging
import os
from os import PathLike
import pathlib
from typing import Union

from matplotlib import pyplot
from matplotlib.axis import Axis
from matplotlib.transforms import Bbox
import numpy as numpy
from pandas import DataFrame
from pyproj import CRS, Transformer
from shapely import ops
from shapely.geometry import Point, Polygon
from stormevents.nhc import VortexTrack
from stormevents.nhc.atcf import ATCF_Mode
import utm

from pyschism.enums import NWSType
from pyschism.forcing.nws.base import NWS
from pyschism.forcing.nws.nws2.sflux import SfluxDataset
from pyschism.mesh import gridgr3


class BestTrackForcing(VortexTrack, NWS):

    def __init__(
        self,
        storm: Union[str, PathLike, DataFrame, io.BytesIO],
        start_date: datetime = None,
        end_date: datetime = None,
        mode: ATCF_Mode = None,
    ):


        VortexTrack.__init__(
            self,
            storm=storm,
            start_date=start_date,
            end_date=end_date,
            file_deck='b',
            advisories=['BEST'],
        )


    def __str__(self):
        """Returns string used in param.nml"""
        return f"{self.dtype.value}"


    @classmethod
    def from_nhc_bdeck(
        cls,
        nhc_bdeck: PathLike,
        start_date: datetime = None,
        end_date: datetime = None,
    ) -> 'NWS':
        return cls.from_file(path=nhc_bdeck, start_date=start_date, end_date=end_date)

    def summary(
        self, output: Union[str, os.PathLike] = None, overwrite: bool = False,
    ):
        min_storm_speed = numpy.min(self.data['speed'])
        max_storm_speed = numpy.max(self.data['speed'])
        track_length = self.distance
        duration = self.duration
        min_central_pressure = numpy.min(self.data['central_pressure'])
        max_wind_speed = numpy.max(self.data['max_sustained_wind_speed'])
        start_loc = (self.data['longitude'][0], self.data['latitude'][0])
        end_loc = (self.data['longitude'].iloc[-1], self.data['latitude'].iloc[-1])
        f = [
            f'Summary of storm: {self.nhc_code}',
            f'min./max. track speed: {min_storm_speed} m/s, {max_storm_speed} m/s',
            f'min. central pressure: {min_central_pressure} hPa',
            f'max. wind speed: {max_wind_speed} kts',
            f'Starting at: {start_loc} and ended at: {end_loc}',
            f'Total track length: {track_length:.2f} km',
            f'Total track duration: {duration:.2f} days',
        ]
        summary = '\n'.join(f)
        if output is not None:
            if not isinstance(output, pathlib.Path):
                output = pathlib.Path(output)
            if overwrite or not output.exists():
                with open(output, 'w+') as fh:
                    fh.write(summary)
            else:
                logging.debug(f'skipping existing file "{output}"')
        return summary

    def write(self, path: PathLike, overwrite: bool = False):
        VortexTrack.to_file(
            self, path=path/'hurricane-track.dat', overwrite=overwrite)

    @property
    def dtype(self) -> NWSType:
        """Returns the datatype of the object"""
        return NWSType(-1)

    def clip_to_bbox(self, bbox, bbox_crs):
        msg = f'bbox must be a {Bbox} instance.'
        assert isinstance(bbox, Bbox), msg
        bbox_pol = Polygon(
            [
                [bbox.xmin, bbox.ymin],
                [bbox.xmax, bbox.ymin],
                [bbox.xmax, bbox.ymax],
                [bbox.xmin, bbox.ymax],
                [bbox.xmin, bbox.ymin],
            ]
        )
        _switch = True
        unique_dates = numpy.unique(self.data['datetime'])
        _found_start_date = False
        for _datetime in unique_dates:
            records = self.data[self.data['datetime'] == _datetime]
            radii = records['radius_of_last_closed_isobar'].iloc[0]
            radii = 1852.0 * radii  # convert to meters
            lon = records['longitude'].iloc[0]
            lat = records['latitude'].iloc[0]
            _, _, number, letter = utm.from_latlon(lat, lon)
            df_crs = CRS.from_epsg(4326)
            utm_crs = CRS.from_epsg(f'326{number}')
            transformer = Transformer.from_crs(df_crs, utm_crs, always_xy=True)
            p = Point(*transformer.transform(lon, lat))
            pol = p.buffer(radii)
            transformer = Transformer.from_crs(utm_crs, bbox_crs, always_xy=True)
            pol = ops.transform(transformer.transform, pol)
            if _switch is True:
                if not pol.intersects(bbox_pol):
                    continue
                else:
                    self.start_date = records['datetime'].iloc[0]
                    _found_start_date = True
                    _switch = False
                    continue

            else:
                if pol.intersects(bbox_pol):
                    continue
                else:
                    self.end_date = records['datetime'].iloc[0]
                    break

        if _found_start_date is False:
            raise Exception(f'No data within mesh bounding box for storm {self.storm_id}.')
