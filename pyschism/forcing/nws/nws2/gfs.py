import signal
import os
import errno
from functools import wraps
from datetime import datetime, timedelta
from enum import Enum
import pathlib
import tempfile
from typing import Union
import logging

import pytz
from matplotlib.transforms import Bbox
from netCDF4 import Dataset
import numpy as np

from pyschism.enums import GFSProduct
from pyschism.forcing.nws.nws2.sflux import (
    SfluxDataset,
    AirComponent,
    PrcComponent,
    RadComponent,
)
from pyschism.dates import nearest_zulu, localize_datetime, nearest_cycle

BASE_URL = "https://nomads.ncep.noaa.gov/dods"
logger = logging.getLogger(__name__)


class BaseURL(Enum):
    GFS_0P25 = f"{BASE_URL}/gfs_0p25"
    GFS_0P25_1HR = f"{BASE_URL}/gfs_0p25_1hr"
    GFS_0P50 = f"{BASE_URL}/gfs_0p50"
    GFS_1P00 = f"{BASE_URL}/gfs_1p00"

    @classmethod
    def _missing_(self, name):
        raise ValueError(f"{name} is not a valid GFS product.")


class TimeoutError(Exception):
    pass


def timeout(seconds=30, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)

    return decorator


class GFSInventory:
    def __init__(self, product="gfs_0p25_1hr", start_date=None, rnday=4, bbox=None):
        self.product = (
            GFSProduct(product.upper())
            if not isinstance(product, GFSProduct)
            else product
        )
        self.start_date = (
            nearest_cycle()
            if start_date is None
            else localize_datetime(start_date).astimezone(pytz.utc)
        )
        self.rnday = rnday if isinstance(rnday, timedelta) else timedelta(days=rnday)
        # if self.start_date != nearest_cycle(self.start_date):
        #     raise NotImplementedError(
        #         "Argument start_date is does not align with any GFS cycle " "times."
        #     )
        self.check_dates(self.start_date, self.rnday)
        end_date = (self.start_date + self.rnday).astimezone(pytz.utc)
        end_date.replace(tzinfo=None)
        self._files = {
            _: None
            for _ in np.arange(
                self.start_date,
                end_date + self.output_interval,
                self.output_interval,
            ).astype(datetime)
        }
        
        nearest_end_date = nearest_cycle(datetime.utcnow(), period=6)
        nearest_end_date.replace(tzinfo=None)
        for dt in self.nearest_zulus:
            if None not in list(self._files.values()):
                break
            base_url = (
                BASE_URL
                + f"/{self.product.value}"
                + f'/gfs{nearest_zulu(dt).strftime("%Y%m%d")}'
            )
            for cycle in reversed(range(0, 24, 6)):
                if np.datetime64(dt + timedelta(hours=cycle)) > np.datetime64(
                    nearest_end_date
                ):
                    continue
                test_url = f"{base_url}/" + f"{self.product.name.lower()}_{cycle:02d}z"
                nc = self.fetch_nc_by_url(test_url)
                if nc is None:
                    continue
                file_dates = self.get_nc_datevector(nc)
                for _datetime in reversed(list(self._files.keys())):
                    if _datetime in file_dates:
                        if self._files[_datetime] is None:
                            self._files[_datetime] = nc
                if not any(nc is None for nc in self._files.values()):
                    break

        missing_records = [dt for dt, nc in self._files.items() if nc is None]
        if len(missing_records) > 0:
            raise ValueError(f"No GFS data for dates: {missing_records}.")

        self._bbox = self._modified_bbox(bbox)

    def check_dates(self, start_date, rnday):

        end_date = (start_date + rnday).astimezone(pytz.utc)
        end_date.replace(tzinfo=None)
        nc = None
        i = 0
        while nc is None:

            dt = nearest_cycle(
                nearest_cycle(datetime.utcnow()) - (i * self.output_interval)
            )
            cycle = (6 * np.floor(dt.hour / 6)) % 24
            base_url = (
                BASE_URL + f"/{self.product.value}" + f'/gfs{dt.strftime("%Y%m%d")}'
            )
            test_url = f"{base_url}/" + f"{self.product.name.lower()}_{int(cycle):02d}z"
            nc = self.fetch_nc_by_url(test_url)
            i += 1
        datevec = self.get_nc_datevector(nc)
        # print(np.datetime64(np.max(datevec)))
        # print(np.datetime64(end_date - self.output_interval))
        if np.datetime64(np.max(datevec)) < np.datetime64(
            end_date - 2 * self.output_interval
        ):
            raise IOError(
                f"Requested end date at {end_date - self.output_interval} is larger than the max allowed GFS end_date "
                f" {np.max(datevec) - self.output_interval}"
            )

    @staticmethod
    def fetch_nc_by_url(url):
        try:
            logger.info(f"Checking url: {url}")

            @timeout()
            def get_netcdf_timeout():
                return Dataset(url)

            nc = get_netcdf_timeout()
            logger.info("Success!")
            return nc
        except OSError as e:
            if e.errno == -70:
                print()
                logger.info("URL not yet in server.")
                return None
            elif e.errno == -73:
                nc = False

                @timeout()
                def retry():
                    try:
                        return Dataset(url)
                    except Exception:
                        return False

                while not isinstance(nc, Dataset):
                    logger.info("retrying...")
                    nc = retry()

            elif e.errno == -72:
                return None

            else:
                raise e

    def put_sflux_field(self, gfs_varname: str, dst: Dataset, sflux_varname: str):

        lon_idxs, lat_idxs = self._modified_bbox_indexes(self._bbox)
        for i, (dt, nc) in enumerate(self._files.items()):
            time_index = self.get_nc_time_index(nc, dt)
            if time_index == 0:
                time_index += 1
            logger.info(
                f"Putting GFS field {gfs_varname} for time {dt} as "
                f"{sflux_varname} from file "
                f'{nc.filepath().replace(f"{BASE_URL}/", "")}[{time_index}].'
            )
            # TODO: exponential time backoff retry
            def put_nc_field():
                try:
                    dst[sflux_varname][i, :, :] = nc.variables[gfs_varname][
                        time_index, lat_idxs, lon_idxs
                    ]
                    dst.sync()
                    if np.any(dst[sflux_varname][i, :, :].mask) is True:
                        raise Exception(
                            f"Junk value in nc_field {gfs_varname}, index {time_index}"
                        )
                    return True
                except RuntimeError:
                    return False

            success = False
            while success is False:
                success = put_nc_field()

    def get_nc_time_index(self, nc, dt):
        return np.where(np.in1d(self.get_nc_datevector(nc), [dt]))[0][0]

    def get_nc_datevector(self, nc):
        try:
            base_date = localize_datetime(
                datetime.strptime(nc["time"].minimum.split("z")[-1], "%d%b%Y")
            ) + timedelta(hours=float(nc["time"].minimum.split("z")[0]))
            return np.arange(
                base_date,
                base_date + timedelta(len(nc["time"][:])),
                self.output_interval,
            ).astype(datetime)
        except RuntimeError:
            return self.get_nc_datevector(nc)

    def get_sflux_timevector(self):
        timevec = list(self._files.keys())
        _nearest_zulu = nearest_zulu(np.min(timevec))
        return [
            (localize_datetime(x) - _nearest_zulu) / timedelta(days=1) for x in timevec
        ]

    def xy_grid(self):
        lon_idxs, lat_idxs = self._modified_bbox_indexes(self._bbox)
        lon = []
        for x in self.lon[lon_idxs]:
            if x > 180:
                lon.append(x - 360)
            else:
                lon.append(x)
        return np.meshgrid(np.array(lon), self.lat[lat_idxs])

    @property
    def nearest_zulu(self):
        if not hasattr(self, "_nearest_zulu"):
            self._nearest_zulu = nearest_zulu()
        return self._nearest_zulu

    @property
    def nearest_zulus(self):
        return np.arange(
            self.nearest_zulu,
            self.nearest_zulu - timedelta(days=10),
            -timedelta(days=1),
        ).astype(datetime)

    @property
    def output_interval(self):
        if self.product == GFSProduct.GFS_0P25_1HR:
            return timedelta(hours=1)
        return timedelta(hours=6)

    @property
    def lon(self):
        if not hasattr(self, "_lon"):
            nc = self._files[list(self._files.keys())[0]]
            self._lon = nc.variables["lon"][:]
            if not hasattr(self, "_lat"):
                self._lat = nc.variables["lat"][:]
        return self._lon

    @property
    def lat(self):
        if not hasattr(self, "_lat"):
            nc = self._files[list(self._files.keys())[0]]
            self._lat = nc.variables["lat"][:]
            if not hasattr(self, "_lon"):
                self._lon = nc.variables["lon"][:]
        return self._lat

    def _modified_bbox(self, bbox=None):
        if bbox is None:
            return Bbox.from_extents(0, -90, 360, 90)
        else:
            xmin = bbox.xmin + 360 if bbox.xmin < 0 else bbox.xmin
            xmax = bbox.xmax + 360 if bbox.xmax < 0 else bbox.xmax
            return Bbox.from_extents(xmin, bbox.ymin, xmax, bbox.ymax)

    def _modified_bbox_indexes(self, bbox):
        lat_idxs = np.where((self.lat >= bbox.ymin-1.0) & (self.lat <= bbox.ymax+1.0))[0]
        lon_idxs = np.where((self.lon >= bbox.xmin-1.0) & (self.lon <= bbox.xmax+1.0))[0]
        return lon_idxs, lat_idxs


class GlobalForecastSystem(SfluxDataset):
    def __init__(
        self,
        product: Union[str, GFSProduct] = GFSProduct.GFS_0P25_1HR,
    ):
        self.prmsl_name = "prmslmsl"
        self.spfh_name = "spfh2m"
        self.stmp_name = "tmpsfc"
        self.uwind_name = "ugrd10m"
        self.vwind_name = "vgrd10m"
        self.prate_name = "pratesfc"
        self.dlwrf_name = "dlwrfsfc"
        self.dswrf_name = "dswrfsfc"
        self.product = (
            GFSProduct(product) if not isinstance(product, GFSProduct) else product
        )
        self.air = None
        self.prc = None
        self.rad = None

    def write(
        self,
        outdir: Union[str, os.PathLike],
        level,
        start_date: datetime = None,
        rnday: Union[float, timedelta] = 4,
        air: bool = True,
        prc: bool = True,
        rad: bool = True,
        bbox: Bbox = None,
        overwrite: bool = False,
    ):
        """Fetches GFS data from NOMADS server."""
        logger.info("Fetching GFS data.")
        self.start_date = (
            nearest_cycle()
            if start_date is None
            else localize_datetime(start_date).astimezone(pytz.utc)
        )
        self.rnday = rnday if isinstance(rnday, timedelta) else timedelta(days=rnday)

        inventory = GFSInventory(
            self.product, self.start_date, self.rnday + self.output_interval, bbox
        )

        nx_grid, ny_grid = inventory.xy_grid()
        if air is True:
            with Dataset(
                self.tmpdir / f"air_{inventory.product.value}_"
                f"{str(self.start_date)}.nc",
                "w",
                format="NETCDF3_CLASSIC",
            ) as dst:

                # global attributes
                dst.setncatts({"Conventions": "CF-1.0"})
                # dimensions
                dst.createDimension("nx_grid", nx_grid.shape[1])
                dst.createDimension("ny_grid", ny_grid.shape[0])
                dst.createDimension("time", None)
                # variables
                # lon
                dst.createVariable("lon", "f4", ("ny_grid", "nx_grid"))
                dst["lon"].long_name = "Longitude"
                dst["lon"].standard_name = "longitude"
                dst["lon"].units = "degrees_east"
                dst["lon"][:] = nx_grid
                # lat
                dst.createVariable("lat", "f4", ("ny_grid", "nx_grid"))
                dst["lat"].long_name = "Latitude"
                dst["lat"].standard_name = "latitude"
                dst["lat"].units = "degrees_north"
                dst["lat"][:] = ny_grid
                # time
                dst.createVariable("time", "f4", ("time",))
                dst["time"].long_name = "Time"
                dst["time"].standard_name = "time"
                date = nearest_zulu(self.start_date)
                dst["time"].units = (
                    f"days since {date.year}-{date.month}"
                    f"-{date.day} 00:00"
                    f"{date.tzinfo}"
                )
                dst["time"].base_date = (date.year, date.month, date.day, 0)
                dst["time"][:] = inventory.get_sflux_timevector()

                for var in AirComponent.var_types:
                    dst.createVariable(var, "f4", ("time", "ny_grid", "nx_grid"))
                    logger.info(f"Put field {var}")
                    inventory.put_sflux_field(getattr(self, f"{var}_name"), dst, var)

                # prmsl
                dst["prmsl"].long_name = "Pressure reduced to MSL"
                dst["prmsl"].standard_name = "air_pressure_at_sea_level"
                dst["prmsl"].units = "Pa"

                # spfh
                dst["spfh"].long_name = "Surface Specific Humidity " "(2m AGL)"
                dst["spfh"].standard_name = "specific_humidity"
                dst["spfh"].units = "1"

                # stmp
                dst["stmp"].long_name = "Surface Air Temperature (2m AGL)"
                dst["stmp"].standard_name = "air_temperature"
                dst["stmp"].units = "K"

                # uwind
                dst["uwind"].long_name = "Surface Eastward Air Velocity " "(10m AGL)"
                dst["uwind"].standard_name = "eastward_wind"
                dst["uwind"].units = "m/s"

                # vwind
                dst["vwind"].long_name = "Surface Northward Air Velocity " "(10m AGL)"
                dst["vwind"].standard_name = "northward_wind"
                dst["vwind"].units = "m/s"

        if prc is True:
            with Dataset(
                self.tmpdir / f"prc_{inventory.product.value}_"
                f"{str(self.start_date)}.nc",
                "w",
                format="NETCDF3_CLASSIC",
            ) as dst:

                # global attributes
                dst.setncatts({"Conventions": "CF-1.0"})
                # dimensions
                dst.createDimension("nx_grid", nx_grid.shape[1])
                dst.createDimension("ny_grid", ny_grid.shape[0])
                dst.createDimension("time", None)
                # lon
                dst.createVariable("lon", "f4", ("ny_grid", "nx_grid"))
                dst["lon"].long_name = "Longitude"
                dst["lon"].standard_name = "longitude"
                dst["lon"].units = "degrees_east"
                dst["lon"][:] = nx_grid
                # lat
                dst.createVariable("lat", "f4", ("ny_grid", "nx_grid"))
                dst["lat"].long_name = "Latitude"
                dst["lat"].standard_name = "latitude"
                dst["lat"].units = "degrees_north"
                dst["lat"][:] = ny_grid
                # time
                dst.createVariable("time", "f4", ("time",))
                dst["time"].long_name = "Time"
                dst["time"].standard_name = "time"
                date = nearest_zulu(self.start_date)
                dst["time"].units = (
                    f"days since {date.year}-{date.month}"
                    f"-{date.day} 00:00"
                    f"{date.tzinfo}"
                )
                dst["time"].base_date = (date.year, date.month, date.day, 0)
                dst["time"][:] = inventory.get_sflux_timevector()

                for var in PrcComponent.var_types:
                    dst.createVariable(var, float, ("time", "ny_grid", "nx_grid"))
                    logger.info(f"Put field {var}")
                    inventory.put_sflux_field(getattr(self, f"{var}_name"), dst, var)
                # prate
                dst["prate"].long_name = "Surface Precipitation Rate"
                dst["prate"].standard_name = "air_pressure_at_sea_level"
                dst["prate"].units = "kg/m^2/s"

        if rad is True:
            with Dataset(
                self.tmpdir / f"rad_{inventory.product.value}_"
                f"{str(self.start_date)}.nc",
                "w",
                format="NETCDF3_CLASSIC",
            ) as dst:
                # global attributes
                dst.setncatts({"Conventions": "CF-1.0"})
                # dimensions
                dst.createDimension("nx_grid", nx_grid.shape[1])
                dst.createDimension("ny_grid", ny_grid.shape[0])
                dst.createDimension("time", None)
                # lon
                dst.createVariable("lon", "f4", ("ny_grid", "nx_grid"))
                dst["lon"].long_name = "Longitude"
                dst["lon"].standard_name = "longitude"
                dst["lon"].units = "degrees_east"
                dst["lon"][:] = nx_grid
                # lat
                dst.createVariable("lat", "f4", ("ny_grid", "nx_grid"))
                dst["lat"].long_name = "Latitude"
                dst["lat"].standard_name = "latitude"
                dst["lat"].units = "degrees_north"
                dst["lat"][:] = ny_grid
                # time
                dst.createVariable("time", "f4", ("time",))
                dst["time"].long_name = "Time"
                dst["time"].standard_name = "time"
                date = nearest_zulu(self.start_date)
                dst["time"].units = (
                    f"days since {date.year}-{date.month}"
                    f"-{date.day} 00:00"
                    f"{date.tzinfo}"
                )
                dst["time"].base_date = (date.year, date.month, date.day, 0)
                dst["time"][:] = inventory.get_sflux_timevector()

                for var in RadComponent.var_types:
                    dst.createVariable(var, float, ("time", "ny_grid", "nx_grid"))
                    logger.info(f"Put field {var}")
                    inventory.put_sflux_field(getattr(self, f"{var}_name"), dst, var)

                # dlwrf
                dst["dlwrf"].long_name = "Downward Long Wave Radiation " "Flux"
                dst["dlwrf"].standard_name = (
                    "surface_downwelling_" "longwave_flux_in_air"
                )
                dst["dlwrf"].units = "W/m^2"

                # dswrf
                dst["dswrf"].long_name = "Downward Short Wave Radiation " "Flux"
                dst["dswrf"].standard_name = (
                    "surface_downwelling_" "shortwave_flux_in_air"
                )
                dst["dswrf"].units = "W/m^2"

        self.resource = list(self.tmpdir.glob("*.nc"))
        if air is True:
            self.air = AirComponent(self.fields)

        if prc is True:
            self.prc = PrcComponent(self.fields)

        if rad is True:
            self.rad = RadComponent(self.fields)

        super().write(
            outdir,
            level,
            overwrite=overwrite,
            start_date=start_date,
            rnday=rnday,
            air=air,
            rad=rad,
            prc=prc,
        )
        del self._tmpdir

    @property
    def tmpdir(self):
        if not hasattr(self, "_tmpdir"):
            self._tmpdir = tempfile.TemporaryDirectory()
        return pathlib.Path(self._tmpdir.name)

    @property
    def output_interval(self):
        if self.product == GFSProduct.GFS_0P25_1HR:
            return timedelta(hours=1)
        return timedelta(hours=6)
