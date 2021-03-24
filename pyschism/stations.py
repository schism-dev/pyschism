from datetime import timedelta
import pathlib
from typing import Union, Dict, List, Any

from pyproj import Transformer, CRS  # type: ignore[import]
from shapely.geometry import (   # type: ignore[import]
    Polygon,
    MultiPolygon,
    Point
)

from pyschism.enums import (
    StationOutputVariables,
    StationOutputIndex
)


class Stations:

    def __init__(self, nspool_sta: Union[int, timedelta],
                 crs: Union[str, CRS] = None, elev: bool = False,
                 air_pressure: bool = False, windx: bool = False,
                 windy: bool = False, T: bool = False, S: bool = False,
                 u: bool = False, v: bool = False, w: bool = False):
        """Acts as container for requesting point outputs to SCHISM

        Args:
            nspool_sta: If type is <int>, this refers to output skip, else use
                a timedelta object.
            crs: Coordinate reference system for the container. All points
                added using the add_station() method should be in this CRS.
            elev (optional): Request elevation output, default: False
            air_pressure (optional): Request air_pressure output,
                default: False
            windx (optional): Request wind zonal output,
                default: False
            windy (optional): Request wind meridional output, default: False
            T (optional): Request temperature output, default: False
            S (optional): Request salinity output, default: False
            u (optional): Request zonal fluid velocity output, default: False
            v (optional): Request meridional fluid velocity output,
                default: False
            w (optional): Request vertical fluid velocity output,
                default: False
        """

        if not isinstance(nspool_sta, (int, timedelta)):
            raise TypeError('nspool_sta must be an int or timedelta')
        self._nspool_sta = nspool_sta

        if isinstance(crs, str):
            crs = CRS.from_user_input(crs)
        self._crs = crs

        # contains main station data (coordinates and id)
        self._stations: List[Dict[str, Any]] = []

        # init the properties
        self.elev = elev
        self.air_pressure = air_pressure
        self.windx = windx
        self.windy = windy
        self.T = T
        self.S = S
        self.u = u
        self.v = v
        self.w = w

    def __iter__(self):
        for i, s in enumerate(self._stations):
            yield i+1, s['x'], s['y'], s['z'], s['comment']

    def __str__(self):
        f = [f'{self.state}',
             f'{len(self._stations)}']
        for i, x, y, z, comment in self:
            if comment is None:
                comment = ''
            else:
                comment = comment.split('!')
                comment = ' '.join(comment)
            f.append(f'{i} {x} {y} {z} ! {comment}')
        return '\n'.join(f)

    @staticmethod
    def from_file(file: Union[str, pathlib.Path],
                  nspool_sta: Union[int, timedelta],
                  crs: Union[str, CRS] = None, **kwargs):
        """Import station.in file as a Stations instance.

        Args:
            file: Path to the input.in file to import
            nspool_sta: Output frequency request
            crs (optional): Specify CRS of station.in. If no crs is given,
                the program assumes that the file is in the same CRS as the
                Mesh.
            elev (optional): Request elevation output, default: False
            air_pressure (optional): Request air_pressure output,
                default: False
            windx (optional): Request wind zonal output,
                default: False
            windy (optional): Request wind meridional output,
                default: False
            T (optional): Request temperature output, default: False
            S (optional): Request salinity output, default: False
            u (optional): Request zonal fluid velocity output,
                default: False
            v (optional): Request meridional fluid velocity output,
                default: False
            w (optional): Request vertical fluid velocity output,
                default: False

        Returns:
            instance of Stations
        """

        with open(pathlib.Path(file), 'r') as f:
            states = f.readline().split()[:9]
            # override file state with kwargs request
            for i, state in enumerate(states):
                states[i] = kwargs.get(
                    StationOutputIndex(i).name.lower(), False)
            stations = Stations(nspool_sta, crs=crs,
                                **{var.value: bool(states[i])
                                    for i, var in
                                    enumerate(StationOutputVariables)})
            comment: Union[str, None] = None
            for i in range(int(f.readline().split()[0])):
                line = f.readline().split()
                id, x, y, z = line[:4]
                if len(line) > 4:
                    comment = ' '.join(line[5:])
                else:
                    comment = None
                stations.add_station(float(x), float(y), float(z), comment)
        return stations

    def add_station(self, x: float, y: float, z: float = 0,
                    comment: str = None):
        """Adds a station coordinates to the container.

        Args:
            x: x-coordinate
            y: y-coordinate
            z: z-coordinate
            comment (optional): Can be used to store the COOPS id, useful
                during post-processing for validation.
        """
        self._stations.append({'x': x, 'y': y, 'z': z, 'comment': comment})

    def get_active_vars(self) -> List[str]:
        """Returns a list of the names of activated output variables.

        Returns:
            list of strings for each active station output request.
        """
        return [var.value for var in StationOutputVariables
                if getattr(self, var.value) is True]

    def transform_to(self, dst_crs: Union[str, CRS]):
        """Transforms the horizontal coordinates of the stations.

        If :class:`pyschism.Stations` was instantiated with crs=None, then
        this function cannot be used.
        """
        if isinstance(dst_crs, str):
            dst_crs = CRS.from_user_input(dst_crs)
        else:
            if not isinstance(dst_crs, CRS):
                raise TypeError(f'Input must be {str} or {CRS}')
        if dst_crs.equals(self.crs):
            return
        transformer = Transformer.from_crs(self.crs, dst_crs, always_xy=True)
        x = [_['x'] for _ in self._stations]
        y = [_['y'] for _ in self._stations]
        xy = list(zip(*transformer.transform(x, y)))
        for i, (x, y) in enumerate(xy):
            self._stations[i]['x'] = x
            self._stations[i]['y'] = y
        self._crs = dst_crs

    def clip(self, geometry: Union[Polygon, MultiPolygon]):
        """Eliminates any stations not contained within the given geometry.

        It is assumed that the input geometry matches the CRS of the
        Stations instance

        Args:
            geometry: Polygon or MultiPolygon used for clipping the
                stations. Normally this comes from Hgrid.get_multipolygon()
        """
        for i, s in reversed(list(enumerate(self._stations))):
            if not geometry.contains(Point(s['x'], s['y'])):
                self._stations.pop(i)

    def write(self, path: Union[str, pathlib.Path], overwrite: bool = False):
        """Writes the SCHISM station.in file to disk.

        Args:
            path: Path to output on disk (need to include filename).
            overwrite (optional): Allow/disallow file overwrite.
        """
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise IOError("File exists and overwrite is False.")
        with open(path, 'w') as f:
            f.write(str(self))

    @property
    def stations(self) -> List[Dict]:
        """Returns list of currently loaded stations."""
        return self._stations

    @property
    def nspool_sta(self) -> Union[int, timedelta]:
        """Returns output frequency."""
        return self._nspool_sta

    @property
    def crs(self) -> Union[CRS, None]:
        """Returns coordinate reference system of the current instance."""
        return self._crs

    @property
    def state(self) -> str:
        """Returns corresponding string that goes into bctide.in"""
        return ' '.join([str(int(getattr(self, var.value)))
                         for var in StationOutputVariables])

    @property
    def elev(self) -> bool:
        """Returns state (ON/OFF) of variable elev request."""
        return self._elev

    @elev.setter
    def elev(self, elev: bool):
        assert isinstance(elev, bool), 'Argument to "elev" must be boolean'
        self._elev = elev

    @property
    def air_pressure(self) -> bool:
        """Returns state (ON/OFF) of variable air_pressure request."""
        return self._air_pressure

    @air_pressure.setter
    def air_pressure(self, air_pressure: bool):
        assert isinstance(air_pressure, bool), \
            'Argument to "air_pressure" must be boolean'
        self._air_pressure = air_pressure

    @property
    def windx(self) -> bool:
        """Returns state (ON/OFF) of variable windx request."""
        return self._windx

    @windx.setter
    def windx(self, windx: bool):
        assert isinstance(windx, bool), 'Argument to "windx" must be boolean'
        self._windx = windx

    @property
    def windy(self) -> bool:
        """Returns state (ON/OFF) of variable windy request."""
        return self._windy

    @windy.setter
    def windy(self, windy: bool):
        assert isinstance(windy, bool), 'Argument to "windy" must be boolean'
        self._windy = windy

    @property
    def T(self) -> bool:
        """Returns state (ON/OFF) of variable T request."""
        return self._T

    @T.setter
    def T(self, T: bool):
        assert isinstance(T, bool), 'Argument to "T" must be boolean'
        self._T = T

    @property
    def S(self) -> bool:
        """Returns state (ON/OFF) of variable S request."""
        return self._S

    @S.setter
    def S(self, S: bool):
        assert isinstance(S, bool), 'Argument to "S" must be boolean'
        self._S = S

    @property
    def u(self) -> bool:
        """Returns state (ON/OFF) of variable u request."""
        return self._u

    @u.setter
    def u(self, u: bool):
        assert isinstance(u, bool), 'Argument to "u" must be boolean'
        self._u = u

    @property
    def v(self) -> bool:
        """Returns state (ON/OFF) of variable v request."""
        return self._v

    @v.setter
    def v(self, v: bool):
        assert isinstance(v, bool), 'Argument to "v" must be boolean'
        self._v = v

    @property
    def w(self) -> bool:
        """Returns state (ON/OFF) of variable w request."""
        return self._w

    @w.setter
    def w(self, w: bool):
        assert isinstance(w, bool), 'Argument to "w" must be boolean'
        self._w = w
