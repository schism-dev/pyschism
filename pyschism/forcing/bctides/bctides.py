from datetime import datetime, timedelta
import pathlib
from typing import Union, List, TYPE_CHECKING

# import pytz
# from pyschism import dates
# from pyschism.driver import raise_type_error
from pyschism.forcing.bctides import (
    iettype, ifltype, itetype, isatype, itrtype)
from pyschism.forcing.tides.tides import Tides, TidalDatabase
from pyschism.mesh import Hgrid, Vgrid

if TYPE_CHECKING:
    from pyschism.driver import ModelDriver


def raise_type_error(argname, obj, cls):
    raise TypeError(
        f'Argument {argname} must be of type {cls}, not '
        f'type {type(obj)}.')
# class HgridDescriptor:

#     def __set__(self, obj, val: Hgrid):
#         if not isinstance(val, Hgrid):
#             raise TypeError(
#                 f'Argument hgrid must be of type {Hgrid}, not type '
#                 f'{type(val)}.')
#         obj.__dict__['hgrid'] = val

#     def __get__(self, obj, val):
#         return obj.__dict__['hgrid']


# class StartDateDescriptor:

#     def __set__(self, obj, val: datetime):
#         if not isinstance(val, datetime):
#             raise TypeError(
#                     f'Argument start_date must be of type {datetime}, '
#                     f'not type {type(val)}.')
#         if datetime_is_naive(val):
#             val = pytz.timezone('UTC').localize(val)
#         obj.__dict__['start_date'] = val

#     def __get__(self, obj, val):
#         return obj.__dict__['start_date']


# class RndayDescriptor:

#     def __set__(self, obj, val: Union[int, float, timedelta]):
#         if not isinstance(val, (int, float, timedelta)):
#             raise TypeError(
#                 f'Argument rnday must be of type {int}, {float} or '
#                 f'{timedelta}, not type {type(val)}.')
#         if not isinstance(val, timedelta):
#             val = timedelta(days=val)
#         obj.__dict__['rnday'] = val

#     def __get__(self, obj, val) -> timedelta:
#         return obj.__dict__['rnday']


class Bctides:

    # _hgrid = HgridDescriptor()
    # _start_date = StartDateDescriptor()
    # _rnday = RndayDescriptor()
    # start_date = dates.StartDate()
    # end_date = dates.EndDate()
    # spinunp_time = dates.SpinupTime()
    # rnday = dates.StartDate()

    def __init__(
            self,
            hgrid: Hgrid,
            start_date: datetime,
            rnday: timedelta,
            vgrid: Vgrid = None,
            elevation: iettype.Iettype = None,
            velocity: ifltype.Ifltype = None,
            temperature: itetype.Itetype = None,
            salinity: isatype.Isatype = None,
            tracers: Union[itrtype.Itrtype, List[itrtype.Itrtype]] = None,
            cutoff_depth: float = 50.
    ):
        self.hgrid = hgrid
        self.start_date = start_date
        self.rnday = rnday
        self.vgrid = vgrid
        self.elevation = elevation
        self.velocity = velocity
        self.temperature = temperature
        self.salinity = salinity
        self.tracers = list(tracers) if tracers is not None else []
        self.cutoff_depth = cutoff_depth

    @classmethod
    def from_driver(cls, driver: 'ModelDriver'):
        tides = driver.config.forcings.tides
        baroclinic = driver.config.forcings.baroclinic
        # elevation =
        # velocity =
        # temperature = 
        # salinity = 
        # tracers = 
        return cls(
            driver.config.hgrid,
            driver.param.opt.start_date,
            driver.param.core.rnday,
            driver.config.vgrid,
            elevation,
            velocity,
            temperature,
            salinity,
            tracers
        )

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
            f.append(self.get_forcing_string(boundary))
        return '\n'.join(f)

    def write(
            self,
            output_directory,
            bctides: Union[bool, str] = True,
            elev2D: Union[bool, str] = True,
            uv3D: Union[bool, str] = True,
            flux: Union[bool, str] = True,
            overwrite: bool = False
    ):
        # self.tidal_database.write(path, )
        output_directory = pathlib.Path(output_directory)
        bctides = output_directory / 'bctides.in' if bctides is True \
            else bctides
        if bctides.exists() and not overwrite:
            raise IOError('path exists and overwrite is False')
        with open(bctides, 'w') as f:
            f.write(str(self))
        # write elev2D.th.nc
        if self.elevation is not None:
            if self.elevation.iettype in [4, 5]:
                elev2D = output_directory / 'elev2D.th.nc' if elev2D is True \
                    else elev2D
                self.elevation.write(
                    elev2D,
                    self.hgrid,
                    self.start_date,
                    self.rnday,
                    overwrite
                )

        if self.velocity is not None:
            if self.velocity.ifltype in [4, 5, -4, -5]:
                uv3D = output_directory / 'uv3D.th.nc' if uv3D is True \
                    else uv3D
                self.velocity.write(
                    uv3D,
                    self.hgrid,
                    self.vgrid,
                    self.start_date,
                    self.rnday,
                    overwrite
                )
            elif self.velocity.ifltype == 1:
                flux = output_directory / 'flux.th' if flux is True \
                    else flux
                self.velocity.write(flux, overwrite)

        def write_tracer(tracer):
            tracer.write()

        for tracer in [self.temperature, self.salinity, *self.tracers]:
            if tracer is not None:
                write_tracer(tracer)

    def get_forcing_string(self, boundary):

        def forcing_digit(pos):
            if pos == 0:
                return f'{len(boundary.indexes)}'

            elif pos == 1:
                if self.elevation is None:
                    return '0'
                return str(self.elevation.iettype)

            elif pos == 2:
                if self.velocity is None:
                    return '0'
                return str(self.velocity.ifltype)

            elif pos == 3:
                if self.temperature is None:
                    return '0'
                return str(self.temperature.itetype)

            elif pos == 4:
                if self.salinity is None:
                    return '0'
                return str(self.salinity.isatype)

            elif pos >= 5:
                if len(self.tracers) > 0:
                    if self.tracers[5-pos] is not None:
                        return str(self.tracers[5-pos].itrtype)
                    return '0'
                return '0'

            else:
                raise ValueError(f'Unhandled argument pos={pos}.')

        f = [' '.join(list(map(forcing_digit, range(5 + len(self.tracers)))))]

        for pos in range(1, 5 + len(self.tracers)):
            if pos == 1:
                if self.elevation is not None:
                    f.append(
                        self.elevation.get_boundary_string(
                            self.hgrid, boundary))
            elif pos == 2:
                if self.velocity is not None:
                    f.append(
                        self.velocity.get_boundary_string(self.hgrid, boundary)
                        )
            elif pos == 3:
                if self.temperature is not None:
                    f.append(
                        self.temperature.get_boundary_string(
                            self.hgrid, boundary))
            elif pos == 4:
                if self.salinity is not None:
                    f.append(
                        self.salinity.get_boundary_string(self.hgrid, boundary)
                        )
            elif pos >= 5:
                f.append(self.tracers[5+pos].get_boundary_string(
                    self.hgrid, boundary))
            else:
                raise TypeError(f'Unhandled forcing at pos={pos}.')
        return '\n'.join(f)

    @property
    def hgrid(self):
        return self._hgrid

    @hgrid.setter
    def hgrid(self, hgrid: Hgrid):
        if not isinstance(hgrid, Hgrid):
            raise_type_error('hgrid', hgrid, Hgrid)
        self._hgrid = hgrid

    @property
    def tides(self):
        if not hasattr(self, '_tides'):
            if self.elevation is None and self.velocity is None:
                raise ValueError('Both elevation and velocity are disabled; at least one of the is required.')
            elif self.elevation is not None and self.velocity is None:
                self._tides = Tides(
                    elevation=True,
                    velocity=False,
                    tidal_database=self.elevation.tides.tidal_database
                    )
            elif self.elevation is not None and self.velocity is not None:
                self._tides = Tides(
                    elevation=True,
                    velocity=True,
                    tidal_database=self.elevation.tides.tidal_database
                    )
            else:
                raise ValueError(
                    f'Unhandled combination: self.elevation={self.elvation} and '
                    f'self.velocity={self.velocity}')
        return self._tides

    # @tides.setter
    # def tides(self, tides):
    #     if not isinstance(tides, (Tides, type(None))):
    #         raise_type_error('tides', tides, Tides)
    #     self._tides = tides

    @property
    def ntip(self):
        return len(self.tides.get_active_potential_constituents())

    @property
    def Z0(self):
        if hasattr(self.tides, '_Z0'):
            return self.tides._Z0

    @Z0.setter
    def Z0(self, Z0):
        self.tides.add_Z0(Z0)

    @property
    def cutoff_depth(self):
        return self._cutoff_depth

    @cutoff_depth.setter
    def cutoff_depth(self, cutoff_depth: float):
        self._cutoff_depth = float(cutoff_depth)

    @property
    def tidal_database(self):
        return self._tidal_database

    @tidal_database.setter
    def tidal_database(self, tidal_database: Union[TidalDatabase, str]):
        if tidal_database is not None:
            if not isinstance(tidal_database, TidalDatabase):
                tidal_database = TidalDatabase(tidal_database)
        self._tidal_database = tidal_database

    # @property
    # def subtidal_database(self):
    #     return self._subtidal_database

    # @subtidal_database.setter
    # def subtidal_database(self, subtidal_database: SubTidalDatabase):
    #     if subtidal_database is not None:
    #         # self._subtidal_database = Tides(subtidal_database=subtidal_database)
    #     else:
    #         self._subtidal_database = None

    @property
    def elevation(self):
        return self._elevation

    @elevation.setter
    def elevation(self, elevation):
        if elevation is not None:
            assert isinstance(elevation, iettype.Iettype)
        self._elevation = elevation

    @property
    def velocity(self):
        return self._velocity

    @velocity.setter
    def velocity(self, velocity):
        if velocity is not None:
            assert isinstance(velocity, ifltype.Ifltype)
        self._velocity = velocity

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, temperature: Union[itetype.Itetype, None]):
        if temperature is not None:
            assert isinstance(temperature, itetype.Itetype)
        self._temperature = temperature

    @property
    def salinity(self):
        return self._salinity

    @salinity.setter
    def salinity(self, salinity: Union[isatype.Isatype, None]):
        if salinity is not None:
            assert isinstance(salinity, isatype.Isatype)
        self._salinity = salinity
