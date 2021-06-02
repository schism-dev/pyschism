from datetime import datetime, timedelta
import pathlib
from typing import Any, Dict, Union, List, TYPE_CHECKING

# import pytz
# from pyschism import dates
# from pyschism.driver import raise_type_error
from pyschism.forcing.bctides import (
    iettype, ifltype, itetype, isatype,
    itrtype
)
from pyschism.forcing.tides.tides import Tides, TidalDatabase
# from pyschism.mesh import (
#     Hgrid,
#     Vgrid
# )

if TYPE_CHECKING:
    from pyschism.driver import ModelDriver


def raise_type_error(argname, obj, cls):
    raise TypeError(
        f'Argument {argname} must be of type {cls}, not '
        f'type {type(obj)}.')


class Bctides:

    def __init__(
            self,
            hgrid,  #: Hgrid,
            start_date: datetime,
            rnday: timedelta,
            vgrid=None,
            cutoff_depth: float = 50.
    ):
        self.hgrid = hgrid
        self.start_date = start_date
        self.rnday = rnday
        self.vgrid = vgrid
        self.cutoff_depth = cutoff_depth

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
        f.append(f'{self.nbfr:d}')
        if self.nbfr > 0:
            for constituent in self.tides.get_active_forcing_constituents():
                forcing = self.tides(
                    self.start_date, self.rnday, constituent)
                f.append(' '.join([
                    f'{constituent}\n',
                    f'{forcing[2]:G}',
                    f'{forcing[3]:G}',
                    f'{forcing[4]:G}']))
        f.append(f'{len(self.hgrid.boundaries.open)}')
        for boundary in self.hgrid.boundaries.open.itertuples():
            f.append(self.get_forcing_string(boundary))
        return '\n'.join(f)

    def write(
            self,
            output_directory,
            bctides: Union[bool, str] = True,
            elev2D: Union[bool, str] = True,
            uv3D: Union[bool, str] = True,
            tem3D: Union[bool, str] = True,
            sal3D: Union[bool, str] = True,
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
        elev2D = output_directory / 'elev2D.th.nc' if elev2D \
            is True else elev2D
        self.hgrid.boundaries.elev2d().write(
                    elev2D,
                    self.start_date,
                    self.rnday,
                    timedelta(days=1),
                    overwrite
                )
        # write uv3D.th.nc
        uv3D = output_directory / 'uv3D.th.nc' if uv3D \
            is True else uv3D
        self.hgrid.boundaries.uv3d(self.vgrid).write(
                    uv3D,
                    self.start_date,
                    self.rnday,
                    timedelta(days=1),
                    overwrite
                )
        # write TEM_3D.th.nc
        tem3D = output_directory / 'TEM_3D.th.nc' if tem3D \
            is True else tem3D
        self.hgrid.boundaries.tem3d(self.vgrid).write(
                    tem3D,
                    self.start_date,
                    self.rnday,
                    timedelta(days=1),
                    overwrite
                )
        sal3D = output_directory / 'SAL_3D.th.nc' if sal3D \
            is True else sal3D
        self.hgrid.boundaries.sal3d(self.vgrid).write(
                    sal3D,
                    self.start_date,
                    self.rnday,
                    timedelta(days=1),
                    overwrite
                )

        # for boundary in self.hgrid.boundaries.open.itertuples():
        #     if boundary.iettype is not None:
        #         if hasattr(boundary.iettype, 'write'):
        #             elev2D = output_directory / 'elev2D.th.nc' if elev2D \
        #                 is True else elev2D
        #             boundary.iettype.write(
        #                 elev2D,
        #                 self.hgrid,
        #                 self.start_date,
        #                 self.rnday,
        #                 overwrite
        #             )

        #             boundary.iettype.write()

        # if self.elevation is not None:
        #     if self.elevation.iettype in [4, 5]:
        #         elev2D = output_directory / 'elev2D.th.nc' if elev2D is True \
        #             else elev2D
        #         self.elevation.write(
        #             elev2D,
        #             self.hgrid,
        #             self.start_date,
        #             self.rnday,
        #             overwrite
        #         )

        # if self.velocity is not None:
        #     if self.velocity.ifltype in [4, 5, -4, -5]:
        #         uv3D = output_directory / 'uv3D.th.nc' if uv3D is True \
        #             else uv3D
        #         self.velocity.write(
        #             uv3D,
        #             self.hgrid,
        #             self.vgrid,
        #             self.start_date,
        #             self.rnday,
        #             overwrite
        #         )
        #     elif self.velocity.ifltype == 1:
        #         flux = output_directory / 'flux.th' if flux is True \
        #             else flux
        #         self.velocity.write(flux, overwrite)

        # def write_tracer(tracer):
        #     tracer.write()

        # for tracer in [self.temperature, self.salinity, *self.tracers]:
        #     if tracer is not None:
        #         write_tracer(tracer)

    def get_forcing_string(self, boundary):

        bctypes = [
            boundary.iettype,
            boundary.ifltype,
            boundary.itetype,
            boundary.isatype,
        ]

        def get_focing_digit(bctype):
            if bctype is not None:
                # sensitive to MRO.
                return str(getattr(
                    bctype,
                    f"{bctype.__class__.__bases__[0].__name__.lower()}"))
            return '0'
        line = [
            f'{len(boundary.indexes)}',
            *[digit for digit in map(get_focing_digit, bctypes)]
        ]

        f = [' '.join(line)]
        for bctype in bctypes:
            if bctype is not None:
                f.append(bctype.get_boundary_string(self.hgrid, boundary))
        return '\n'.join(f)

    @property
    def tides(self):
        if not hasattr(self, '_tides'):
            # get the first one you can find, since the Tides object is a
            # singleton.
            tides = None
            for boundary in self.hgrid.boundaries.open.itertuples():
                if boundary.iettype is not None:
                    if hasattr(boundary.iettype, "tides"):
                        tides = boundary.iettype.tides
                        break
                    elif boundary.ifltype is not None:
                        if hasattr(boundary.ifltype, "tides"):
                            tides = boundary.ifltype.tides
                            break
            self._tides = tides
        return self._tides

    @property
    def tracers(self) -> List[Dict[Any, Union[itrtype.Itrtype, None]]]:
        # if not hasattr(self, '_tracers'):
        #     # tracers: List[Dict[Any, Union[itrtype.Itrtype, None]]] = []
        #     boundary_data = {}
        #     for boundary in self.hgrid.boundaries.open.itertuples():
        #         itrtypes = boundary.itrtype
        #         if itrtypes is None:
        #             tracers.append({})
        #         for tracer in boundary.itr
        #             tracers.append()
        #         tracer.setdefault(

        #             )

        # _itrtype = boundary.itrtype
        # return self._tracers
        # TODO: Cheating for now...
        return []

    @property
    def ntip(self):
        if self.tides is None:
            return 0
        return len(self.tides.get_active_potential_constituents())

    @property
    def nbfr(self):
        if self.tides is None:
            return 0
        return self.tides.nbfr

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

    # @property
    # def subtidal_database(self):
    #     return self._subtidal_database

    # @subtidal_database.setter
    # def subtidal_database(self, subtidal_database: SubTidalDatabase):
    #     if subtidal_database is not None:
    #         # self._subtidal_database = Tides(subtidal_database=subtidal_database)
    #     else:
    #         self._subtidal_database = None

    # @property
    # def elevation(self):
    #     return self._elevation

    # @elevation.setter
    # def elevation(self, elevation):
    #     if elevation is not None:
    #         assert isinstance(elevation, iettype.Iettype)
    #     self._elevation = elevation

    # @property
    # def velocity(self):
    #     return self._velocity

    # @velocity.setter
    # def velocity(self, velocity):
    #     if velocity is not None:
    #         assert isinstance(velocity, ifltype.Ifltype)
    #     self._velocity = velocity

    # @property
    # def temperature(self):
    #     return self._temperature

    # @temperature.setter
    # def temperature(self, temperature: Union[itetype.Itetype, None]):
    #     if temperature is not None:
    #         assert isinstance(temperature, itetype.Itetype)
    #     self._temperature = temperature

    # @property
    # def salinity(self):
    #     return self._salinity

    # @salinity.setter
    # def salinity(self, salinity: Union[isatype.Isatype, None]):
    #     if salinity is not None:
    #         assert isinstance(salinity, isatype.Isatype)
    #     self._salinity = salinity


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