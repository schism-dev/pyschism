from datetime import timedelta
import pathlib
from typing import Union, Dict, List, Any

from pyproj import Transformer, CRS  # type: ignore[import]
from shapely.geometry import Polygon, MultiPolygon, Point


output_vars = ['elev', 'air_pressure', 'windx', 'windy', 'T', 'S', 'u', 'v',
               'w']


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
            The remaining arguments activate the outputs for the run.
        """

        if not isinstance(nspool_sta, (int, timedelta)):
            raise TypeError('nspool_sta must be an int or timedelta')
        self.__nspool_sta = nspool_sta

        if isinstance(crs, str):
            crs = CRS.from_user_input(crs)
        self.__crs = crs

        # contains main station data (coordinates and id)
        self.__stations: List[Dict[str, Any]] = []

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
        for i, s in enumerate(self.__stations):
            yield i+1, s['x'], s['y'], s['z'], s['comment']

    def __str__(self):
        f = [f'{self.state}',
             f'{len(self.__stations)}']
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
            **kwargs: state of each output variable

        Returns:
            instance of Stations
        """

        with open(pathlib.Path(file), 'r') as f:
            states = f.readline().split()[:9]
            # override file state with kwargs request
            for i, state in enumerate(states):
                states[i] = kwargs.get(output_vars[i], False)
            stations = Stations(nspool_sta, crs=crs,
                                **{var: bool(states[i]) for i, var in
                                    enumerate(output_vars)})
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
        self.__stations.append({'x': x, 'y': y, 'z': z, 'comment': comment})

    def get_active_vars(self):
        return [var for var in output_vars if getattr(self, var) is True]

    def transform_to(self, dst_crs: Union[str, CRS]):
        if isinstance(dst_crs, str):
            dst_crs = CRS.from_user_input(dst_crs)
        else:
            if not isinstance(dst_crs, CRS):
                raise TypeError(f'Input must be {str} or {CRS}')
        if dst_crs.equals(self.crs):
            return
        transformer = Transformer.from_crs(self.crs, dst_crs, always_xy=True)
        x = [_['x'] for _ in self.__stations]
        y = [_['y'] for _ in self.__stations]
        xy = list(zip(*transformer.transform(x, y)))
        for i, (x, y) in enumerate(xy):
            self.__stations[i]['x'] = x
            self.__stations[i]['y'] = y
        self.__crs = dst_crs

    def clip(self, geometry: Union[Polygon, MultiPolygon]):
        """Eliminates any stations not contained within the given geometry.

        It is assumed that the input geometry matches the CRS of the
        Stations instance

        Args:
            geometry: Polygon or MultiPolygon used for clipping the
                stations. Normally this comes from Hgrid.get_multipolygon()
        """
        
        eliminate: List[int] = []
        for id, x, y, _, _ in self:
            if not geometry.contains(Point(x, y)):
                eliminate.insert(id-1, 0)
        for id in eliminate:
            del(self.__stations[id])


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
    def stations(self):
        return self.__stations

    @property
    def nspool_sta(self):
        return self.__nspool_sta

    @property
    def crs(self):
        return self.__crs

    @property
    def state(self):
        return ' '.join([str(int(getattr(self, var))) for var in output_vars])

    @property
    def elev(self):
        return self.__elev

    @elev.setter
    def elev(self, elev: bool):
        assert isinstance(elev, bool), 'Argument to "elev" must be boolean'
        self.__elev = elev

    @property
    def air_pressure(self):
        return self.__air_pressure

    @air_pressure.setter
    def air_pressure(self, air_pressure: bool):
        assert isinstance(air_pressure, bool), \
            'Argument to "air_pressure" must be boolean'
        self.__air_pressure = air_pressure

    @property
    def windx(self):
        return self.__windx

    @windx.setter
    def windx(self, windx: bool):
        assert isinstance(windx, bool), 'Argument to "windx" must be boolean'
        self.__windx = windx

    @property
    def windy(self):
        return self.__windy

    @windy.setter
    def windy(self, windy: bool):
        assert isinstance(windy, bool), 'Argument to "windy" must be boolean'
        self.__windy = windy

    @property
    def T(self):
        return self.__T

    @T.setter
    def T(self, T: bool):
        assert isinstance(T, bool), 'Argument to "T" must be boolean'
        self.__T = T

    @property
    def S(self):
        return self.__S

    @S.setter
    def S(self, S: bool):
        assert isinstance(S, bool), 'Argument to "S" must be boolean'
        self.__S = S

    @property
    def u(self):
        return self.__u

    @u.setter
    def u(self, u: bool):
        assert isinstance(u, bool), 'Argument to "u" must be boolean'
        self.__u = u

    @property
    def v(self):
        return self.__v

    @v.setter
    def v(self, v: bool):
        assert isinstance(v, bool), 'Argument to "v" must be boolean'
        self.__v = v

    @property
    def w(self):
        return self.__w

    @w.setter
    def w(self, w: bool):
        assert isinstance(w, bool), 'Argument to "w" must be boolean'
        self.__w = w
