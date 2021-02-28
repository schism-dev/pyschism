from datetime import datetime, timedelta
import pathlib
from typing import Union

import pytz

from pyschism.forcing.tides.tides import Tides
from pyschism.mesh import Hgrid


def datetime_is_naive(d) -> bool:
    return d.tzinfo is None or d.tzinfo.utcoffset(d) is None


def datetime_is_aware(d) -> bool:
    return d.tzinfo is not None and d.tzinfo.utcoffset(d) is not None


class HgridDescriptor:

    def __set__(self, obj, val: Hgrid):
        if not isinstance(val, Hgrid):
            raise TypeError('Argument hgrid must be of type {Hgrid}, not type '
                            f'{type(val)}.')
        obj.__dict__['hgrid'] = val

    def __get__(self, obj, val):
        return obj.__dict__['hgrid']


class StartDateDescriptor:

    def __set__(self, obj, val: datetime):
        if not isinstance(val, datetime):
            raise TypeError(
                    f'Argument start_date must be of type {datetime}, '
                    f'not type {type(val)}.')
        if datetime_is_naive(val):
            val = pytz.timezone('UTC').localize(val)
        obj.__dict__['start_date'] = val

    def __get__(self, obj, val):
        return obj.__dict__['start_date']


class RndayDescriptor:

    def __set__(self, obj, val: Union[int, float, timedelta]):
        if not isinstance(val, (int, float, timedelta)):
            raise TypeError(
                f'Argument rnday must be of type {int}, {float} or '
                f'{timedelta}, not type {type(val)}.')
        if not isinstance(val, timedelta):
            val = timedelta(days=val)
        obj.__dict__['rnday'] = val

    def __get__(self, obj, val) -> timedelta:
        return obj.__dict__['rnday']


class Bctides:

    _hgrid = HgridDescriptor()
    _start_date = StartDateDescriptor()
    _rnday = RndayDescriptor()

    def __init__(
            self,
            hgrid: Hgrid,
            start_date: datetime,
            rnday: timedelta,
            cutoff_depth=50.,
            elevation=True,
            velocity=False,
            temperature=False,
            salinity=False,
            tracers=False,
            tidal_database='hamtide',
            # ts_database='hycom',
    ):
        self._hgrid = hgrid
        self._start_date = start_date
        self._rnday = rnday
        self._cutoff_depth = cutoff_depth
        self._elevation = elevation
        self._velocity = velocity
        self._temperature = temperature
        self._salinity = salinity
        self._tracers = [] if tracers is None or tracers is False else tracers
        self._tides = Tides(database=tidal_database)
        self.tides.use_all()

    def __str__(self):
        f = [
            f'{str(self.start_date)}',
            f'{self.ntip} {self.cutoff_depth}',
        ]
        if self.ntip > 0:
            for constituent in self.tides.get_active_potential_constituents():
                forcing = self.tides(
                    self.start_date, self.rnday, constituent)
                f.append(' '.join([
                    f'{constituent}\n',
                    f'{forcing[0]:G}',
                    f'{forcing[1]:G}',
                    f'{forcing[2]:G}',
                    f'{forcing[3]:G}',
                    f'{forcing[4]:G}']))
        f.append(f'{self.tides.nbfr:d}')
        for constituent in self.tides.get_active_forcing_constituents():
            forcing = self.tides(
                self.start_date, self.rnday, constituent)
            f.append(' '.join([
                f'{constituent}\n',
                f'{forcing[2]:G}',
                f'{forcing[3]:G}',
                f'{forcing[4]:G}']))
        f.append(f'{len(self.hgrid.boundaries.ocean())}')
        for boundary in self.hgrid.boundaries.ocean().itertuples():
            f.append(self.get_forcing(boundary))
        return '\n'.join(f)

    def write(self, path, overwrite: bool = False):
        path = pathlib.Path(path)
        if path.exists() and not overwrite:
            raise IOError('path exists and overwrite is False')
        open(path.resolve(), 'w').write(str(self))

    def get_forcing(self, boundary):

        def forcing_digit(pos):
            if pos == 0:
                return f'{len(boundary.indexes)}'

            elif pos == 1:
                if self.elevation is True:
                    return '3'
                elif self.elevation is False or self.elevation is None:
                    return '0'
                else:
                    raise TypeError('Unhandled elevation source type '
                                    f'{type(self.elevation)}.')

            elif pos == 2:
                if self.velocity is True:
                    return '3'
                elif self.velocity is False or self.velocity is None:
                    return '0'
                else:
                    raise TypeError('Unhandled velocity source type '
                                    f'{type(self.velocity)}.')

            elif pos == 3:
                if self.temperature is False or self.temperature is None:
                    return '0'
                else:
                    raise TypeError('Unhandled temperature source type '
                                    f'{type(self.temperature)}.')

            elif pos == 4:
                if self.salinity is False or self.salinity is None:
                    return '0'
                else:
                    raise TypeError('Unhandled salinity source type '
                                    f'{type(self.salinity)}.')

            else:
                raise ValueError(f'Unhandled argument pos={pos}.')

        def get_forcing(pos, digit):

            f = []

            if pos == 1:  # elevation
                if digit == 3:  # tides
                    if self.elevation is True:  # default tidal elevation
                        for constituent in \
                                self.tides.get_active_forcing_constituents():
                            f.append(f'{constituent}')
                            vertices = self.hgrid.get_xy(
                                crs='EPSG:4326')[boundary.indexes, :]
                            amp, phase = self.tides.get_elevation(
                                constituent, vertices)
                            for i in range(len(amp)):
                                f.append(f'{amp[i]:.8e} {phase[i]:.8e}')
                        return '\n'.join(f)
                    else:
                        raise TypeError(
                            f'Unhandled elevation type {type(self.elevation)}')
                else:
                    raise ValueError(
                        f'Unhandled elevation forcing digit={digit}')
            elif pos == 2:  # velocity
                if self.velocity is True:  # default tidal velocity
                    for constituent in \
                            self.tides.get_active_forcing_constituents():
                        f.append(f'{constituent}')
                        vertices = self.hgrid.get_xy(
                            crs='EPSG:4326')[boundary.indexes, :]
                        uamp, uphase, vamp, vphase = self.tides.get_velocity(
                            constituent, vertices)
                        for i in range(len(vertices)):
                            f.append(f'{uamp[i]:.8e} {uphase[i]:.8e} '
                                     f'{vamp[i]:.8e} {vphase[i]:.8e}')
                    return '\n'.join(f)
                else:
                    raise TypeError(
                        f'Unhandled velocity type {type(self.velocity)}')
            else:
                raise TypeError('Unhandled forcing at pos={pos}.')

        f = [' '.join(list(map(forcing_digit, range(5 + len(self._tracers)))))]

        for i in range(1, 5 + len(self._tracers)):
            if int(forcing_digit(i)) > 0:
                f.append(get_forcing(i, int(forcing_digit(i))))
        return '\n'.join(f)

    @property
    def hgrid(self):
        return self._hgrid

    @property
    def tides(self):
        return self._tides

    @property
    def ntip(self):
        return len(self.tides.get_active_potential_constituents())

    @property
    def cutoff_depth(self):
        return float(self._cutoff_depth)

    @property
    def elevation(self):
        return self._elevation

    @property
    def velocity(self):
        return self._velocity

    @property
    def temperature(self):
        return self._temperature

    @property
    def salinity(self):
        return self._salinity
