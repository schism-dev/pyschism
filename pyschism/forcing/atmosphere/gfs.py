from datetime import datetime, timedelta
from enum import Enum
import pathlib
import tempfile
from typing import Union
# import warnings

import pytz
from matplotlib.transforms import Bbox  # type: ignore[import]
from netCDF4 import Dataset  # type: ignore[import]
import numpy as np  # type: ignore[import]

from pyschism.enums import GFSProduct
from pyschism.forcing.atmosphere.nws.nws2.sflux import (
    SfluxDataset,
    AirComponent,
    PrcComponent,
    RadComponent,
)

BASE_URL = 'https://nomads.ncep.noaa.gov/dods'


class BaseURL(Enum):
    GFS_0P25 = f'{BASE_URL}/gfs_0p25'
    GFS_0P25_1HR = f'{BASE_URL}/gfs_0p25_1hr'
    GFS_0P50 = f'{BASE_URL}/gfs_0p50'
    GFS_1P00 = f'{BASE_URL}/gfs_1p00'

    @classmethod
    def _missing_(self, name):
        raise ValueError(f'{name} is not a valid GFS product.')


class Tmpdir:

    def __get__(self, obj, val):
        tmpdir = obj.__dict__.get('tmpdir')
        if tmpdir is None:
            tmpdir = tempfile.TemporaryDirectory()
            obj.__dict__['tmpdir'] = tmpdir
            obj.resource = tmpdir.name
        return tmpdir


class StartDate:

    def __set__(self, obj, start_date):
        if start_date is not None:
            obj.__dict__['start_date'] = start_date

    def __get__(self, obj, val):
        start_date = obj.__dict__.get(
            'start_date', obj.utc_timezone.localize(datetime.utcnow()))
        return datetime(start_date.year, start_date.month, start_date.day,
                        int(6 * np.floor(start_date.hour/6)),
                        tzinfo=start_date.tzinfo).astimezone(obj.utc_timezone)


class Inventory:

    def __get__(self, obj, val):
        utc_now = obj.utc_timezone.localize(datetime.utcnow())
        current_cycle = int(6 * np.floor(utc_now.hour/6))
        nomads_now = obj.utc_timezone.localize(datetime(
            utc_now.year, utc_now.month, utc_now.day, current_cycle))
        nomads_start_dates = [nomads_now - i*obj.forecast_interval
                              for i in range(40)]
        inventory = {}
        for reference_date in nomads_start_dates:
            nearest_cycle = int(6 * np.floor(reference_date.hour/6))
            if nearest_cycle != 0:
                # TODO: SCHISM is limited to cycle 00z only.
                continue
            if reference_date + timedelta(hours=0.5) > obj.start_date:
                self.add_dataset_to_inventory(obj, inventory, reference_date)
            else:
                self.add_dataset_to_inventory(obj, inventory, reference_date)
                break
        return inventory

    def add_dataset_to_inventory(self, obj, inventory, reference_date):
        _reference_date = reference_date.strftime('%Y%m%d')
        start_data_url = obj.base_url + f'/gfs{_reference_date}'
        nearest_cycle = int(6 * np.floor(reference_date.hour/6))
        fname = f'{obj.product.value}_{nearest_cycle:02d}z'
        file_url = start_data_url + f'/{fname}'
        try:
            inventory[reference_date] = Dataset(file_url)
        except OSError as e:
            if e.errno == -70:
                cdate = reference_date - obj.forecast_interval
                closest_cycle = int(6 * np.floor(cdate.hour/6))
                nomads_closest = obj.utc_timezone.localize(datetime(
                    cdate.year, cdate.month, cdate.day, closest_cycle))
                if reference_date not in inventory.keys():
                    _reference_date = nomads_closest.strftime('%Y%m%d')
                    start_data_url = obj.base_url + f'/gfs{_reference_date}'
                    fname = f'{obj.product.value}_' \
                            f'{nomads_closest.hour:02d}z'
                    file_url = start_data_url + f'/{fname}'
                    inventory[cdate] = Dataset(file_url)
            else:
                raise e


class GlobalForecastSystem(SfluxDataset):

    _tmpdir = Tmpdir()
    _inventory = Inventory()
    start_date = StartDate()

    def __init__(
        self,
        product: Union[str, GFSProduct] = GFSProduct.GFS_0P25_1HR,
    ):

        if not isinstance(product, GFSProduct):
            product = GFSProduct(product)
        self.product = product
        self.base_url = BaseURL[product.name].value
        self.prmsl_name = 'prmslmsl'
        self.spfh_name = 'spfh2m'
        self.stmp_name = 'tmpsfc'
        self.uwind_name = 'ugrd10m'
        self.vwind_name = 'vgrd10m'
        self.prate_name = 'pratesfc'
        self.dlwrf_name = 'dlwrfsfc'
        self.dswrf_name = 'dswrfsfc'

    def fetch_data(
            self,
            start_date: datetime = None,
            rnday: Union[float, timedelta] = None,
            air: bool = True,
            prc: bool = True,
            rad: bool = True,
            bbox: Bbox = None,
            timezone=None,
    ):
        """Fetches GFS data from NOMADS server. """

        self.start_date = start_date
        self._rnday = rnday
        # inventory = self.get_inventory(start_date, rnday)
        if timezone is None:
            timezone = self.utc_timezone
        if len(self._inventory) == 0:
            raise ValueError("No data. Perhaps the server is down?")
        for date in sorted(self._inventory.keys()):
            date = date.astimezone(timezone)
            nc = self._inventory[date]
            lons = nc.variables['lon'][:]
            lats = nc.variables['lat'][:]
            if bbox is None:
                bbox = Bbox.from_extents(0, -90, 360, 90)
            else:
                xmin = bbox.xmin + 360 if bbox.xmin < 0 else bbox.xmin
                xmax = bbox.xmax + 360 if bbox.xmax < 0 else bbox.xmax
                bbox = Bbox.from_extents(xmin, bbox.ymin, xmax, bbox.ymax)

            lat_idxs = np.where((lats >= bbox.ymin)
                                & (lats <= bbox.ymax))[0]
            lats = lats[lat_idxs]
            lon_idxs = np.where((lons >= bbox.xmin)
                                & (lons <= bbox.xmax))[0]
            _lons = []
            for x in lons[lon_idxs]:
                if x > 180:
                    _lons.append(x-360)
                else:
                    _lons.append(x)
            lons = np.array(_lons)
            x, y = np.meshgrid(lons, lats)
            if air is True:
                with Dataset(pathlib.Path(self._tmpdir.name) /
                             f"air_{self.product.value}_{date}.nc", 'w',
                             format='NETCDF3_CLASSIC') as dst:
                    # global attributes
                    dst.setncatts({"Conventions": "CF-1.0"})
                    # dimensions
                    dst.createDimension('nx_grid', lons.shape[0])
                    dst.createDimension('ny_grid', lats.shape[0])
                    dst.createDimension('time', None)
                    # variables
                    # lon
                    dst.createVariable('lon', 'f4', ('ny_grid', 'nx_grid'))
                    dst['lon'].long_name = "Longitude"
                    dst['lon'].standard_name = "longitude"
                    dst['lon'].units = "degrees_east"
                    dst['lon'][:] = x
                    # lat
                    dst.createVariable('lat', 'f4', ('ny_grid', 'nx_grid'))
                    dst['lat'].long_name = "Latitude"
                    dst['lat'].standard_name = "latitude"
                    dst['lat'].units = "degrees_north"
                    dst['lat'][:] = y
                    # time
                    dst.createVariable('time', 'f4', ('time',))
                    dst['time'].long_name = 'Time'
                    dst['time'].standard_name = 'time'
                    dst['time'].units = f'days since {date.year}-{date.month}'\
                                        f'-{date.day}'
                    dst['time'].base_date = (date.year, date.month, date.day,
                                             date.hour)
                    dst['time'][:] = np.arange(
                        nc['time'].resolution,
                        nc['time'].resolution*(nc['time'].shape[0]+1),
                        step=nc['time'].resolution)

                    # prmsl
                    dst.createVariable(
                        'prmsl', nc.variables[self.prmsl_name].datatype,
                        ('time', 'ny_grid', 'nx_grid'))
                    dst['prmsl'].long_name = "Pressure reduced to MSL"
                    dst['prmsl'].standard_name = "air_pressure_at_sea_level"
                    dst['prmsl'].units = "Pa"
                    dst['prmsl'][:] = nc.variables[
                        self.prmsl_name][:, lat_idxs, lon_idxs]

                    # spfh
                    dst.createVariable(
                        'spfh', nc.variables[self.spfh_name].datatype,
                        ('time', 'ny_grid', 'nx_grid'))
                    dst['spfh'].long_name = "Surface Specific Humidity "\
                                            "(2m AGL)"
                    dst['spfh'].standard_name = "specific_humidity"
                    dst['spfh'].units = "1"
                    dst['spfh'][:] = nc.variables[
                        self.spfh_name][:, lat_idxs, lon_idxs]

                    # stmp
                    dst.createVariable(
                        'stmp', nc.variables[self.stmp_name].datatype,
                        ('time', 'ny_grid', 'nx_grid'))
                    dst['stmp'].long_name = "Surface Air Temperature (2m AGL)"
                    dst['stmp'].standard_name = "air_temperature"
                    dst['stmp'].units = "K"
                    dst['stmp'][:] = nc.variables[
                        self.stmp_name][:, lat_idxs, lon_idxs]

                    # uwind
                    dst.createVariable(
                        'uwind', nc.variables[self.uwind_name].datatype,
                        ('time', 'ny_grid', 'nx_grid'))
                    dst['uwind'].long_name = "Surface Eastward Air Velocity "\
                                             "(10m AGL)"
                    dst['uwind'].standard_name = "eastward_wind"
                    dst['uwind'].units = "m/s"
                    dst['uwind'][:] = nc.variables[
                        self.uwind_name][:, lat_idxs, lon_idxs]

                    # vwind
                    dst.createVariable(
                        'vwind', nc.variables[self.vwind_name].datatype,
                        ('time', 'ny_grid', 'nx_grid'))
                    dst['vwind'].long_name = "Surface Northward Air Velocity "\
                                             "(10m AGL)"
                    dst['vwind'].standard_name = "northward_wind"
                    dst['vwind'].units = "m/s"
                    dst['vwind'][:] = nc.variables[
                        self.vwind_name][:, lat_idxs, lon_idxs]

            if prc is True:
                with Dataset(pathlib.Path(self._tmpdir.name) /
                             f"prc_{self.product.value}_{date}.nc", 'w',
                             format='NETCDF3_CLASSIC') as dst:
                    # global attributes
                    dst.setncatts({"Conventions": "CF-1.0"})
                    # dimensions
                    dst.createDimension('nx_grid', lons.shape[0])
                    dst.createDimension('ny_grid', lats.shape[0])
                    dst.createDimension('time', None)
                    # variables
                    # lon
                    dst.createVariable('lon', 'f4', ('ny_grid', 'nx_grid'))
                    dst['lon'].long_name = "Longitude"
                    dst['lon'].standard_name = "longitude"
                    dst['lon'].units = "degrees_east"
                    dst['lon'][:] = x
                    # lat
                    dst.createVariable('lat', 'f4', ('ny_grid', 'nx_grid'))
                    dst['lat'].long_name = "Latitude"
                    dst['lat'].standard_name = "latitude"
                    dst['lat'].units = "degrees_north"
                    dst['lat'][:] = y
                    # time
                    dst.createVariable('time', 'f4', ('time',))
                    dst['time'].long_name = 'Time'
                    dst['time'].standard_name = 'time'
                    dst['time'].units = f'days since {date.year}-{date.month}'\
                                        f'-{date.day} {date.hour:02d}:'\
                                        f"{date.minute:02d} " \
                                        f'{date.tzinfo}'
                    dst['time'].base_date = (date.year, date.month, date.day,
                                             date.hour)
                    dst['time'][:] = np.arange(
                        nc['time'].resolution,
                        nc['time'].resolution*(len(nc['time'][:])),
                        step=nc['time'].resolution)

                    # prate
                    dst.createVariable(
                        'prate', nc.variables[self.prate_name].datatype,
                        ('time', 'ny_grid', 'nx_grid'))
                    dst['prate'].long_name = "Surface Precipitation Rate"
                    dst['prate'].standard_name = "air_pressure_at_sea_level"
                    dst['prate'].units = "kg/m^2/s"
                    dst['prate'][:] = nc.variables[
                        self.prate_name][:, lat_idxs, lon_idxs]

            if rad is True:
                with Dataset(pathlib.Path(self._tmpdir.name) /
                             f"rad_{self.product.value}_{date}.nc", 'w',
                             format='NETCDF3_CLASSIC') as dst:
                    # global attributes
                    dst.setncatts({"Conventions": "CF-1.0"})
                    # dimensions
                    dst.createDimension('nx_grid', lons.shape[0])
                    dst.createDimension('ny_grid', lats.shape[0])
                    dst.createDimension('time', None)
                    # variables
                    # lon
                    dst.createVariable('lon', 'f4', ('ny_grid', 'nx_grid'))
                    dst['lon'].long_name = "Longitude"
                    dst['lon'].standard_name = "longitude"
                    dst['lon'].units = "degrees_east"
                    dst['lon'][:] = x
                    # lat
                    dst.createVariable('lat', 'f4', ('ny_grid', 'nx_grid'))
                    dst['lat'].long_name = "Latitude"
                    dst['lat'].standard_name = "latitude"
                    dst['lat'].units = "degrees_north"
                    dst['lat'][:] = y
                    # time
                    dst.createVariable('time', 'f4', ('time',))
                    dst['time'].long_name = 'Time'
                    dst['time'].standard_name = 'time'
                    dst['time'].units = f'days since {date.year}-{date.month}'\
                                        f'-{date.day}T{date.hour:02d}:'\
                                        f'{date.minute:02d}'\
                                        f'{date.tzinfo}'
                    dst['time'].base_date = (date.year, date.month, date.day,
                                             date.hour)
                    dst['time'][:] = np.arange(
                        nc['time'].resolution,
                        nc['time'].resolution*(len(nc['time'][:])),
                        step=nc['time'].resolution)

                    # dlwrf
                    dst.createVariable(
                        'dlwrf', nc.variables[self.dlwrf_name].datatype,
                        ('time', 'ny_grid', 'nx_grid'))
                    dst['dlwrf'].long_name = "Downward Long Wave Radiation "\
                                             "Flux"
                    dst['dlwrf'].standard_name = "surface_downwelling_"\
                                                 "longwave_flux_in_air"
                    dst['dlwrf'].units = "W/m^2"
                    dst['dlwrf'][:] = nc.variables[
                        self.dlwrf_name][:, lat_idxs, lon_idxs]

                    # dswrf
                    dst.createVariable(
                        'dswrf', nc.variables[self.dswrf_name].datatype,
                        ('time', 'ny_grid', 'nx_grid'))
                    dst['dswrf'].long_name = "Downward Short Wave Radiation "\
                                             "Flux"
                    dst['dswrf'].standard_name = "surface_downwelling_"\
                                                 "shortwave_flux_in_air"
                    dst['dswrf'].units = "W/m^2"
                    dst['dswrf'][:] = nc.variables[
                        self.dswrf_name][:, lat_idxs, lon_idxs]

        self.resource = self._tmpdir.name
        if air:
            self.air = AirComponent(self.fields)
        if prc:
            self.prc = PrcComponent(self.fields)
        if rad:
            self.rad = RadComponent(self.fields)

    @property
    def forecast_interval(self):
        return timedelta(hours=6)

    @property
    def utc_timezone(self):
        return pytz.timezone('UTC')

    @property
    def nomads_server_ttl(self):
        return timedelta(days=10)
