from abc import ABC, abstractmethod
from datetime import datetime, timedelta
import logging
import os
import pathlib
from typing import Union
import warnings

import cf
import cftime
from netCDF4 import Dataset
import numpy as np
import pytz

from pyschism.dates import localize_datetime


_logger = logging.getLogger(__name__)


class FieldsDescriptor:

    def __set__(self, obj, fields: cf.FieldList):
        if len(fields) == 0:
            warnings.warn(f"No data with name {obj.name} exists on the input "
                          'dataset.')
        ordered_indexes = np.argsort(
            [field.construct('time').reference_datetime for field
             in fields])
        fields = cf.FieldList([fields[idx] for idx in ordered_indexes])
        obj.__dict__[obj.type] = fields

    def __get__(self, obj, val):
        return obj.__dict__[obj.type]


class Nxgrids:

    def __set__(self, obj, fields: cf.FieldList):
        if len(fields) == 0:
            raise ValueError("No data with name 'lon' exists on the input "
                             'dataset.')
        ordered_indexes = np.argsort(
            [field.construct('time').reference_datetime for field
             in fields])
        fields = cf.FieldList([fields[idx] for idx in ordered_indexes])
        obj.__dict__['nx_grids'] = fields

    def __get__(self, obj, val):
        return obj.__dict__['nx_grids']


class Nygrids:

    def __set__(self, obj, fields: cf.FieldList):
        if len(fields) == 0:
            raise ValueError("No data with name 'lat' exists on the input "
                             'dataset.')
        ordered_indexes = np.argsort(
            [field.construct('time').reference_datetime for field
             in fields])
        fields = cf.FieldList([fields[idx] for idx in ordered_indexes])
        obj.__dict__['ny_grids'] = fields

    def __get__(self, obj, val):
        return obj.__dict__['ny_grids']


class DatetimeArrays:

    def __get__(self, obj, val):
        for field in obj.fields:
            dt_array = []
            for dt in field.construct('time').datetime_array:
                if not isinstance(dt, cftime.DatetimeGregorian):
                    continue
                dt_array.append(pytz.timezone('UTC').localize(
                    datetime(dt.year, dt.month, dt.day, dt.hour, dt.minute)))
            yield dt_array


class DatetimeArray:

    def __get__(self, obj, val):
        dt_array = set()
        for dt_array_x in obj.datetime_arrays:
            for dt in dt_array_x:
                dt_array.add(dt)
        return [dt for dt in sorted(dt_array)]


class ReferenceDatetimes:

    def __get__(self, obj, val):
        for field in obj.fields:
            dt = field.construct('time').reference_datetime
            yield pytz.timezone('UTC').localize(
                datetime(dt.year, dt.month, dt.day, dt.hour, dt.minute)
                )


class Variable:

    fields = FieldsDescriptor()
    reference_datetimes = ReferenceDatetimes()
    datetime_arrays = DatetimeArrays()
    datetime_array = DatetimeArray()

    def __init__(self, fields, var_type, name, long_name,
                 standard_name, units):
        self.type = var_type
        self.name = name
        self.nx_grids = fields.select_by_ncvar('lon')
        self.ny_grids = fields.select_by_ncvar('lat')
        self.fields = fields.select_by_ncvar(self.name)
        self.long_name = long_name
        self.standard_name = standard_name
        self.units = units

    def get_fields(self, start_date: datetime = None,
                   rnday: Union[float, int, timedelta] = None,
                   ) -> cf.FieldList:
        if start_date is None:
            start_date = self.datetime_array[0]

        if start_date < np.min(self.datetime_array) or \
                start_date > np.max(self.datetime_array):
            raise ValueError(f'Requested start date {start_date} is out of '
                             f'range with {np.min(self.datetime_array)} and '
                             f'{np.max(self.datetime_array)}')
        if rnday is None:
            rnday = np.max(self.datetime_array[:-2]) - start_date

        elif isinstance(rnday, (int, float)):
            rnday = timedelta(days=rnday)

        end_date = start_date + rnday  # type: ignore[operator]

        #if end_date <= np.min(self.datetime_array) or \
        #        end_date >= np.max(self.datetime_array):
        #    rnday = end_date - start_date
        #    raise ValueError(f'Requested rnday {rnday} has an end date of '
        #                     f'{end_date} which is out of range with '
        #                     f'start_date={np.min(self.datetime_array)} and '
        #                     f'end_date={np.max(self.datetime_array)}')
        fields = []
        for i, datetime_array in reversed(
                list(enumerate(self.datetime_arrays))):
            field_start_date = np.min(datetime_array)
            if field_start_date >= start_date:
                fields.append(self.fields[i])
            else:
                fields.append(self.fields[i])
                break
        return cf.FieldList(fields)

    def animation(
            self,
            save=False,
            fps=3,
            start_frame=0,
            end_frame=-1,
            figsize=None,
            wireframe=False,
            cmap='jet',
            levels=256,
            show=False,
            xmin=None,
            xmax=None,
            ymin=None,
            ymax=None,
            vmin=None,
            vmax=None,

    ):
        from matplotlib.animation import FuncAnimation
        import matplotlib.pyplot as plt
        # self.current_step = start_frame
        fig = plt.figure(figsize)
        ax = fig.add_subplot(111)

        plt.tight_layout(pad=2)

        # triangulation = self.hgrid.triangulation
        # xmin = np.min(self.hgrid.x) if xmin is None else xmin
        # xmax = np.max(self.hgrid.x) if xmax is None else xmax
        # ymin = np.min(self.hgrid.y) if ymin is None else ymin
        # ymax = np.max(self.hgrid.y) if ymax is None else ymax
        # vmin = np.min(self.values) if vmin is None else vmin
        # vmax = np.max(self.values) if vmax is None else vmax

        # triangulation.set_mask(self.hgrid.elements.get_triangulation_mask(self.parent.wetdry_elem.values))
        
        arrays = self.get_fields()[0].array
        def animate(index):
            # self.current_step = index
            _ax = fig.get_axes()
            ax.clear()
            if len(_ax) > 1:
                cax = _ax[1]
                cax.cla()
            else:
                cax = None

            # triangulation.set_mask(self.hgrid.elements.get_triangulation_mask(self.parent.wetdry_elem.values))

            # if wireframe:
            #     ax.triplot(triangulation, color='k', linewidth=0.7)

            # ax.tricontourf(
            #     triangulation,
            #     self.values,
            #     cmap=cmap,
            #     levels=levels,
            #     vmin=vmin,
            #     vmax=vmax
            #     )
            # print(fields[index].array)
            #print(self.datetime_array[index], arrays[index])
            ax.contourf(arrays[index])

            # ax.set_ylim(ymin, ymax, auto=True)
            # ax.set_xlim(xmin, xmax, auto=True)

            ax.set_xlabel('Longitude (°E)')
            ax.set_ylabel('Latitude (°N)')

            ax.set_title(self.datetime_array[index].strftime('%b %d, %Y %H:%M'))
            # m = plt.cm.ScalarMappable(cmap=cmap)
            # m.set_array(self.values[index])
            # m.set_clim(vmin, vmax)
            # cbar = fig.colorbar(
            #     m, cax=cax, format='%.1f',
            #     boundaries=np.linspace(vmin, vmax, levels)
            # )

            # cbar = fig.colorbar(_ax)
            # cbar.ax.set_ylabel(f'{variable} [{unit}]', rotation=90)

        end_frame = end_frame % len(arrays) if end_frame < 0 else end_frame
        start_frame = start_frame % len(arrays) if start_frame < 0 else start_frame
        frames = range(start_frame, end_frame)
        # print(fields.array[start_frame])
        ax.contourf(arrays[start_frame])
        anim = FuncAnimation(
            fig,
            animate,
            frames,
            blit=False
            )

        if save:
            anim.save(
                pathlib.Path(save),
                writer='imagemagick',
                fps=fps
            )

        if show:
            plt.show()

        return anim


class BaseComponent(ABC):

    name: str = ''

    def write(self, outdir: Union[str, os.PathLike], level: int,
              overwrite: bool = False, start_date: datetime = None,
              rnday: Union[float, int, timedelta] = None):

        assert level in [1, 2]
        outdir = pathlib.Path(outdir)

        if start_date is None:
            for vartype in self.var_types:
                variable = getattr(self, vartype)
                if start_date is None:
                    start_date = np.min(variable.datetime_array)

        if start_date is not None:
            # naive condition
            if start_date.tzinfo is None \
                    or start_date.tzinfo.utcoffset(start_date) is None:
                start_date = pytz.timezone('UTC').localize(start_date)
            timezone = start_date.tzinfo

        stacks = []
        for i, field in enumerate(
            getattr(
                self,
                self.var_types[0]
                ).get_fields(start_date, rnday)
                ):
            stacks.append(f"sflux_{self.name}_{level}.{i+1:04d}.nc")
        for i, filename in enumerate(stacks):
            with Dataset(outdir / filename, 'w',
                         format='NETCDF3_CLASSIC') as dst:
                dst.setncatts({"Conventions": "CF-1.0"})
                # dimensions
                variable = getattr(self, self.var_types[0])
                dst.createDimension('nx_grid',
                                    variable.nx_grids[0].shape[1])
                dst.createDimension('ny_grid',
                                    variable.ny_grids[0].shape[0])
                dst.createDimension('time', None)
                # variables
                # lon
                dst.createVariable('lon', 'f4', ('ny_grid', 'nx_grid'))
                dst['lon'].long_name = "Longitude"
                dst['lon'].standard_name = "longitude"
                dst['lon'].units = "degrees_east"
                dst['lon'][:] = variable.nx_grids[0]
                # lat
                dst.createVariable('lat', 'f4', ('ny_grid', 'nx_grid'))
                dst['lat'].long_name = "Latitude"
                dst['lat'].standard_name = "latitude"
                dst['lat'].units = "degrees_north"
                dst['lat'][:] = variable.ny_grids[0]
                nc_start_date = list(variable.reference_datetimes)[0]
                nc_start_date = nc_start_date.astimezone(timezone)
                dst.createVariable('time', 'f4', ('time',))
                dst['time'].long_name = 'Time'
                dst['time'].standard_name = 'time'
                dst['time'].units = f'days since {nc_start_date.year}-' \
                                    f'{nc_start_date.month}-'\
                                    f'{nc_start_date.day} '\
                                    '00:00:00+' \
                                    f'{nc_start_date.tzinfo}'
                dst['time'].base_date = (
                    nc_start_date.year,
                    nc_start_date.month,
                    nc_start_date.day,
                    0)
                dst['time'][:] = [
                    (localize_datetime(x) - nc_start_date) / timedelta(days=1)
                    for x in variable.datetime_array]
                for vartype in self.var_types:
                    variable = getattr(self, vartype)
                    dst.createVariable(
                        variable.name, 'f4', ('time', 'ny_grid', 'nx_grid'))
                    for field in variable.get_fields(start_date, rnday):
                        dst[variable.name][:] = field
                        setattr(
                            dst[variable.name], "long_name",
                            variable.long_name)
                        setattr(
                            dst[variable.name], "standard_name",
                            variable.standard_name)
                        setattr(
                            dst[variable.name], "units",
                            variable.units)

    @property
    @abstractmethod
    def var_types(self):
        """
        air: ['prmsl', 'spfh', 'stmp', 'uwind', 'vwind']
        prc: ['prate']
        rad: ['dlwrf', 'dswrf']
        """

    @property
    def timevector(self):
        return getattr(self, self.var_types[0]).datetime_array


class AirComponent(BaseComponent):

    name = 'air'
    var_types = ['prmsl', 'spfh', 'stmp', 'uwind', 'vwind']

    def __init__(self, fields: cf.FieldList, prmsl_name='prmsl',
                 spfh_name='spfh', stmp_name='stmp', uwind_name='uwind',
                 vwind_name='vwind'):
        self.prmsl = Variable(fields, 'prmsl', prmsl_name,
                              "Pressure reduced to MSL",
                              "air_pressure_at_sea_level", "Pa")
        self.spfh = Variable(fields, 'spfh', spfh_name,
                             "Surface Specific Humidity (2m AGL)",
                             "specific_humidity", '1')
        self.stmp = Variable(fields, 'stmp', stmp_name,
                             "Surface Air Temperature (2m AGL)",
                             "air_temperature", "K")
        self.uwind = Variable(fields, 'uwind', uwind_name,
                              "Surface Eastward Air Velocity (10m AGL)",
                              "eastward_wind", "m/s")
        self.vwind = Variable(fields, 'vwind', vwind_name,
                              "Surface Northward Air Velocity (10m AGL)",
                              "northward_wind", "m/s")


class PrcComponent(BaseComponent):

    name = 'prc'
    var_types = ['prate']

    def __init__(self, fields: cf.FieldList, prate_name='prate'):
        self.prate = Variable(fields, 'prate', prate_name,
                              "Surface Precipitation Rate",
                              "air_pressure_at_sea_level", "kg/m^2/s")


class RadComponent(BaseComponent):

    name = 'rad'
    var_types = ['dlwrf', 'dswrf']

    def __init__(self, fields: cf.FieldList, dlwrf_name='dlwrf',
                 dswrf_name='dswrf'):
        self.dlwrf = Variable(fields, 'dlwrf', dlwrf_name,
                              "Downward Long Wave Radiation Flux",
                              "surface_downwelling_longwave_flux_in_air",
                              "W/m^2")
        self.dswrf = Variable(fields, 'dswrf', dswrf_name,
                              "Downward Short Wave Radiation Flux",
                              "surface_downwelling_shortwave_flux_in_air",
                              "W/m^2")


class SfluxDataset:

    def __init__(self, resource: Union[str, os.PathLike], prmsl_name='prmsl',
                 spfh_name='spfh', stmp_name='stmp', uwind_name='uwind',
                 vwind_name='vwind', prate_name='prate', dlwrf_name='dlwrf',
                 dswrf_name='dswrf'):
        self.resource = resource
        self.air = AirComponent(self.fields, prmsl_name=prmsl_name,
                                spfh_name=spfh_name, stmp_name=stmp_name,
                                uwind_name=uwind_name, vwind_name=vwind_name)
        self.prc = PrcComponent(self.fields, prate_name=prate_name)
        self.rad = RadComponent(self.fields, dlwrf_name=dlwrf_name,
                                dswrf_name=dswrf_name)

    def write(
            self,
            outdir: Union[str, os.PathLike],
            level: int,
            overwrite: bool = False,
            start_date: datetime = None,
            rnday: Union[float, int, timedelta] = None,
            air=True,
            rad=True,
            prc=True,
            bbox=None,
    ):

        outdir = pathlib.Path(outdir)
        if outdir.name != 'sflux':
            outdir /= 'sflux'
        outdir.mkdir(exist_ok=True)

        if hasattr(self, 'air'):
            if self.air is not None:
                if air is True:
                    self.air.write(outdir, level, overwrite, start_date, rnday)

        if hasattr(self, 'prc'):
            if self.prc is not None:
                if rad is True:
                    self.prc.write(outdir, level, overwrite, start_date, rnday)

        if hasattr(self, 'rad'):
            if self.rad is not None:
                if rad is True:
                    self.rad.write(outdir, level, overwrite, start_date, rnday)

    @property
    def timevector(self):
        for attr in ['air', 'prc', 'rad']:
            if hasattr(self, attr):
                return getattr(self, attr).timevector

    @property
    def resource(self):
        return self._resource

    @resource.setter
    def resource(self, resource):
        self._resource = resource
        if hasattr(self, '_fields'):
            del self._fields

    @property
    def fields(self):
        if not hasattr(self, '_fields'):
            fields = cf.read(self.resource, ignore_read_error=True)
            # check lon
            try:
                lon = fields.select_by_ncvar('lon')[0]
            except IndexError:
                raise ValueError(f"Resource {self.resource} does not contain a "
                                 "'lon' variable.")
            if len(lon.get_data_axes()) != 2:
                raise ValueError("'lon' variable must be a 2-dimensional "
                                 "array")
            fnames = lon.get_filenames()
            lons = fields.select_by_ncvar('lon')
            _logger.info(f'fields.select_by_var() returned {lons}')
            for i in range(len(lons) - 1):
                if not (lon.array == lons[i+1].array).all():
                    raise ValueError(
                        "Invalid sflux dataset. Found two different 'lon' "
                        f"fields on files {fnames} and "
                        f'{lons[i+1].get_filenames()}')
            # check lat
            try:
                lat = fields.select_by_ncvar('lat')[0]
            except IndexError:
                raise ValueError(f"Resource {self.resource} does not contain a "
                                 "'lat' variable.")
            if len(lat.get_data_axes()) != 2:
                raise ValueError("'lat' variable must be a 2-dimensional "
                                 "array")
            fnames = lat.get_filenames()
            lats = fields.select_by_ncvar('lat')
            for i in range(len(lats) - 1):
                if not (lat.array == lats[i+1].array).all():
                    raise ValueError(
                        "Invalid sflux dataset. Found two different 'lat' "
                        f"fields on files {fnames} and "
                        f'{lats[i+1].get_filenames()}')
            self._fields = fields
        return self._fields
