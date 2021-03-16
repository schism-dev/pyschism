import logging
import os
import pathlib

import appdirs
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata

from pyschism.forcing.tides.base import TidalDataProvider


logger = logging.getLogger(__name__)

TPXO_ELEVATION = 'h_tpxo9.v1.nc'
TPXO_VELOCITY = 'uv_tpxo9.v1.nc'


class TPXO(TidalDataProvider):

    def __init__(self, elevation_file=None, velocity_file=None):
        if elevation_file is None:
            elevation_file = os.getenv('TPXO_ELEVATION')
            if elevation_file is None:
                elevation_file = pathlib.Path(
                    appdirs.user_data_dir('tpxo')) / TPXO_ELEVATION
        if not elevation_file.exists():
            raise FileNotFoundError('\n'.join([
                f'No TPXO file found at "{elevation_file}".',
                'New users will need to register and request a copy of '
                f'the TPXO9 NetCDF file (specifically `{TPXO_ELEVATION}`) '
                'from the authors at https://www.tpxo.net.',
                'Once you obtain `h_tpxo9.v1.nc`, you can follow one of the '
                'following options: ',
                f'1) copy or symlink the file to "{elevation_file}"',
                f'2) set the environment variable `{TPXO_ELEVATION}` to point'
                ' to the file',
            ]))
        self._h = Dataset(elevation_file)

        # if velocity_file is None:
        #     velocity_file = os.getenv('TPXO_VELOCITY')
        #     if velocity_file is None:
        #         velocity_file = pathlib.Path(
        #             appdirs.user_data_dir('tpxo')) / TPXO_VELOCITY
        # if not velocity_file.exists():
        #     raise FileNotFoundError('\n'.join([
        #         f'No TPXO file found at "{velocity_file}".',
        #         'New users will need to register and request a copy of '
        #         f'the TPXO9 NetCDF file (specifically `{TPXO_VELOCITY}`) '
        #         'from the authors at https://www.tpxo.net.',
        #         'Once you obtain `h_tpxo9.v1.nc`, you can follow one of the following options: ',
        #         f'1) copy or symlink the file to "{velocity_file}"',
        #         f'2) set the environment variable `{TPXO_VELOCITY}` to point to the file',
        #     ]))
        # self._uv = pathlib.Path(velocity_file)

    def get_elevation(self, constituent, vertices):
        logger.info('Querying TPXO for elevation constituent '
                    f'{constituent}.')
        amp = 0.01 * self._get_interpolation(
            'elevation', 'ha', constituent, vertices)
        phase = self._get_interpolation(
            'elevation', 'hp', constituent, vertices)
        return amp, phase

    def get_velocity(self, constituent, vertices):
        raise NotImplementedError('Velocity not implemented for TPXO.')
    #     logger.info('Querying TPXO for velocity constituent '
    #                 f'{constituent}.')
    #     uamp = 0.01 * self._get_interpolation(
    #         'velocity', 'UAMP', constituent, vertices)
    #     uphase = self._get_interpolation(
    #         'velocity', 'UPHA', constituent, vertices)
    #     vamp = 0.01 * self._get_interpolation(
    #         'velocity', 'VAMP', constituent, vertices)
    #     vphase = self._get_interpolation(
    #         'velocity', 'VPHA', constituent, vertices)
    #     return uamp, uphase, vamp, vphase

    @property
    def constituents(self):
        if not hasattr(self, '_constituents'):
            self._constituents = [
                c.capitalize() for c in self._h['con'][:].astype(
                    '|S1').tostring().decode('utf-8').split()]
        return self._constituents

    @property
    def x(self) -> np.ndarray:
        return self._h['lon_z'][:, 0].data

    @property
    def y(self) -> np.ndarray:
        return self._h['lat_z'][0, :].data

    def _get_interpolation(self, phys_var, ncvar, constituent, vertices):
        lower_c = [c.lower() for c in self.constituents]
        if phys_var == 'elevation':
            ncarray = self._h
        elif phys_var == 'velocity':
            raise NotImplementedError('What\'s the variable for velocity?')
        array = ncarray[ncvar][
                lower_c.index(constituent.lower()), :, :].flatten()
        _x = np.asarray([x + 360. for x in vertices[:, 0] if x < 0]).flatten()
        _y = vertices[:, 1].flatten()
        x, y = np.meshgrid(self.x, self.y, indexing='ij')
        x = x.flatten()
        y = y.flatten()
        dx = np.mean(np.diff(self.x))
        dy = np.mean(np.diff(self.y))

        # buffer the bbox by 2 difference units
        _idx = np.where(
                np.logical_and(
                        np.logical_and(
                                x >= np.min(_x) - 2 * dx,
                                x <= np.max(_x) + 2 * dx
                        ),
                        np.logical_and(
                                y >= np.min(_y) - 2 * dy,
                                y <= np.max(_y) + 2 * dy
                        )
                )
        )

        # "method" can be 'spline' or any string accepted by griddata()'s method kwarg.
        return griddata(
                (x[_idx], y[_idx]),
                array[_idx],
                (_x, _y),
                method='nearest'
        )
