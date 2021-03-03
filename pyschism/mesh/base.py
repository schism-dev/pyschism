from abc import ABC
from collections import defaultdict
from functools import lru_cache
import logging
from itertools import permutations
import os
import pathlib
from typing import Union, Sequence, Hashable, List, Dict

import geopandas as gpd  # type: ignore[import]
from matplotlib.collections import PolyCollection  # type: ignore[import]
from matplotlib.path import Path  # type: ignore[import]
from matplotlib.tri import Triangulation  # type: ignore[import]
from matplotlib.transforms import Bbox  # type: ignore[import]
import numpy as np  # type: ignore[import]
from pyproj import Transformer, CRS  # type: ignore[import]
from shapely.geometry import (  # type: ignore[import]
    MultiPolygon,
    Polygon,
    LineString,
    LinearRing,
)

from pyschism.mesh.parsers import grd
from pyschism.figures import figure

_logger = logging.getLogger(__name__)

class Description:

    def __set__(self, obj, description):
        if description is None:
            description = ''
        if not isinstance(description, str):
            raise TypeError('Argument description must be a string, not '
                            f'type {type(description)}.')
        obj.__dict__['description'] = description.replace('\n', ' ')

    def __get__(self, obj, val):
        return obj.__dict__['description']


class Nodes:

    def __set__(self, obj: "Gr3", nodes: Dict[Hashable, List[List]]):
        """Setter for the nodes attribute.

        Argument nodes must be of the form:
            {id: [(x0, y0), z0]}
            or
            {id: [(x0, y0), [z0, ..., zn]}

        Gr3 format is assumed to be exclusively a 2D format that can hold
        triangles or quads.

        """
        for coords, _ in nodes.values():
            if len(coords) != 2:
                raise ValueError(
                    'Coordinate vertices for a gr3 type must be 2D, but got '
                    f'coordinates {coords}.')
        obj.__dict__['nodes'] = nodes
        self.gr3 = obj
        # self.ball = NodeBall(self.gr3)

        obj.__dict__['id_to_index'] = {
            self.id()[i]: i for i in range(len(self.id()))}

        obj.__dict__['index_to_id'] = {
            i: self.id()[i] for i in range(len(self.id()))}

    def __call__(self):
        return self.gr3.__dict__['nodes']

    def id(self):
        node_id = self.gr3.__dict__.get('node_id')
        if node_id is None:
            node_id = list(self().keys())
            self.gr3.__dict__['node_id'] = node_id
        return node_id

    def index(self):
        node_index = self.gr3.__dict__.get('node_index')
        if node_index is None:
            node_index = np.arange(len(self()))
            self.gr3.__dict__['node_index'] = node_index
        return node_index

    def coord(self):
        coord = self.gr3.__dict__.get('coord')
        if coord is None:
            coord = np.array([coords for coords, _ in self().values()])
            self.gr3.__dict__['coord'] = coord
        return coord

    def values(self):
        values = self.gr3.__dict__.get('values')
        if values is None:
            values = np.array([val for _, val in self().values()])
            self.gr3.__dict__['values'] = values
        return values

    def get_index_by_id(self, id: Hashable):
        return self.gr3.__dict__['id_to_index'][id]

    def get_id_by_index(self, index: int):
        return self.gr3.__dict__['index_to_id'][index]

    def get_indexes_around_index(self, index):
        indexes_around_index = self.__dict__.get('indexes_around_index')
        if indexes_around_index is None:
            def append(geom):
                for simplex in geom:
                    for i, j in permutations(simplex, 2):
                        indexes_around_index[i].add(j)
            indexes_around_index = defaultdict(set)
            append(self.gr3.elements.triangles())
            append(self.gr3.elements.quads())
            self.__dict__['indexes_around_index'] = indexes_around_index
        return list(indexes_around_index[index])


class Elements:

    def __set__(self, obj: "Gr3", elements: Dict[Hashable, Sequence]):
        if not isinstance(elements, dict):
            raise TypeError('Argument elements must be a dict.')
        vertex_id_set = set(obj.nodes.id())
        for i, element in enumerate(elements):
            if not isinstance(element, Sequence):
                raise TypeError(f'Element with index {i} of the elements '
                                f'argument must be of type {Sequence}, not '
                                f'type {type(element)}.')
            if not set(element).issubset(vertex_id_set):
                ValueError(f'Element with index {i} is not a subset of the '
                           "coordinate id's.")
        obj.__dict__['elements'] = elements
        self.gr3 = obj

    def __call__(self):
        return self.gr3.__dict__["elements"]

    def id(self):
        element_id = self.gr3.__dict__.get('element_id')
        if element_id is None:
            element_id = list(self().keys())
            self.gr3.__dict__['element_id'] = element_id
        return element_id

    def index(self):
        element_index = self.gr3.__dict__.get('element_index')
        if element_index is None:
            element_index = np.arange(len(self()))
            self.gr3.__dict__['element_index'] = element_index
        return element_index

    def array(self):
        element_array = self.gr3.__dict__.get('element_array')
        if element_array is None:
            rank = int(max(map(len, self().values())))
            element_array = np.full((len(self()), rank), -1)
            for i, element in enumerate(self().values()):
                row = np.array(list(map(self.gr3.nodes.get_index_by_id, element)))
                element_array[i, :len(row)] = row
                element_array = np.ma.masked_equal(element_array, -1)
                self.gr3.__dict__['element_array'] = element_array
        return element_array

    def triangles(self):
        triangles = self.gr3.__dict__.get('triangles')
        if triangles is None:
            triangles = np.array(
                [list(map(self.gr3.nodes.get_index_by_id, element))
                for element in self().values()
                if len(element) == 3])
            self.gr3.__dict__['triangles'] = triangles
        return triangles

    def quads(self):
        quads = self.gr3.__dict__.get('quads')
        if quads is None:
            quads = np.array(
                [list(map(self.gr3.nodes.get_index_by_id, element))
                for element in self().values()
                if len(element) == 4])
            self.gr3.__dict__['quads'] = quads
        return quads

    def triangulation(self):
        triangles = self.triangles().tolist()
        for quad in self.quads():
            triangles.append([quad[0], quad[1], quad[3]])
            triangles.append([quad[1], quad[2], quad[3]])
        return Triangulation(
            self.gr3.nodes.coord()[:, 0],
            self.gr3.nodes.coord()[:, 1],
            triangles)

    def geodataframe(self):
        data = []
        for id, element in self().items():
            data.append({
                'geometry': Polygon(
                    self.gr3.nodes.coord()[list(
                        map(self.gr3.nodes.get_index_by_id, element))]),
                'id': id})
        return gpd.GeoDataFrame(data, crs=self.gr3.crs)


class Crs:

    def __set__(self, obj, val):
        if val is not None:
            obj.__dict__['crs'] = CRS.from_user_input(val)

    def __get__(self, obj, val):
        return obj.__dict__.get('crs')


class Edges:

    def __init__(self, grd: "Gr3"):
        self.gr3 = grd

    @lru_cache(maxsize=1)
    def __call__(self) -> gpd.GeoDataFrame:
        data = []
        for ring in self.gr3.hull.rings().itertuples():
            coords = ring.geometry.coords
            for i in range(1, len(coords)):
                data.append({
                    "geometry": LineString([coords[i-1], coords[i]]),
                    "bnd_id": ring.bnd_id,
                    "type": ring.type})
        return gpd.GeoDataFrame(data, crs=self.gr3.crs)

    def exterior(self):
        return self().loc[self()['type'] == 'exterior']

    def interior(self):
        return self().loc[self()['type'] == 'interior']


class Rings:

    def __init__(self, grd: "Gr3"):
        self.gr3 = grd

    @lru_cache(maxsize=1)
    def __call__(self) -> gpd.GeoDataFrame:
        tri = self.gr3.elements.triangulation()
        idxs = np.vstack(list(np.where(tri.neighbors == -1))).T
        boundary_edges = []
        for i, j in idxs:
            boundary_edges.append(
                (tri.triangles[i, j], tri.triangles[i, (j+1) % 3]))
        sorted_rings = sort_rings(edges_to_rings(boundary_edges),
                                  self.gr3.nodes.coord())
        data = []
        for bnd_id, rings in sorted_rings.items():
            coords = self.gr3.nodes.coord()[rings['exterior'][:, 0], :]
            geometry = LinearRing(coords)
            data.append({
                    "geometry": geometry,
                    "bnd_id": bnd_id,
                    "type": 'exterior'
                })
            for interior in rings['interiors']:
                coords = self.gr3.nodes.coord()[interior[:, 0], :]
                geometry = LinearRing(coords)
                data.append({
                    "geometry": geometry,
                    "bnd_id": bnd_id,
                    "type": 'interior'
                })
        return gpd.GeoDataFrame(data, crs=self.gr3.crs)

    def exterior(self):
        return self().loc[self()['type'] == 'exterior']

    def interior(self):
        return self().loc[self()['type'] == 'interior']


class Hull:

    def __init__(self, grd: "Gr3"):
        self.gr3 = grd
        self.edges = Edges(grd)
        self.rings = Rings(grd)

    @lru_cache(maxsize=1)
    def __call__(self) -> gpd.GeoDataFrame:
        data = []
        for bnd_id in np.unique(self.rings()['bnd_id'].tolist()):
            exterior = self.rings().loc[
                (self.rings()['bnd_id'] == bnd_id) &
                (self.rings()['type'] == 'exterior')]
            interiors = self.rings().loc[
                (self.rings()['bnd_id'] == bnd_id) &
                (self.rings()['type'] == 'interior')]
            data.append({
                    "geometry": Polygon(
                        exterior.iloc[0].geometry.coords,
                        [row.geometry.coords for _, row
                            in interiors.iterrows()]),
                    "bnd_id": bnd_id
                })
        return gpd.GeoDataFrame(data, crs=self.gr3.crs)

    @lru_cache(maxsize=1)
    def exterior(self):
        data = []
        for exterior in self.rings().loc[
                self.rings()['type'] == 'exterior'].itertuples():
            data.append({"geometry": Polygon(exterior.geometry.coords)})
        return gpd.GeoDataFrame(data, crs=self.gr3.crs)

    @lru_cache(maxsize=1)
    def interior(self):
        data = []
        for interior in self.rings().loc[
                self.rings()['type'] == 'interior'].itertuples():
            data.append({"geometry": Polygon(interior.geometry.coords)})
        return gpd.GeoDataFrame(data, crs=self.gr3.crs)

    @lru_cache(maxsize=1)
    def implode(self) -> gpd.GeoDataFrame:
        return gpd.GeoDataFrame(
            {"geometry": MultiPolygon([polygon.geometry for polygon
                                       in self().itertuples()])},
            crs=self.gr3.crs)

    def multipolygon(self) -> MultiPolygon:
        return self.implode().iloc[0].geometry


class Gr3(ABC):

    _description = Description()
    _nodes = Nodes()
    _elements = Elements()
    _crs = Crs()

    def __init__(self, nodes, elements=None, description=None, crs=None):
        self._nodes = nodes
        self._elements = elements
        self._description = description
        self._crs = crs
        self._hull = Hull(self)

    def __str__(self):
        return grd.to_string(**self.to_dict())

    def __repr__(self):
        return ", ".join([
            f"{self.nodes()}",
            f"elements={self.elements()}",
            f"description={self.description}",
            f"crs={self.crs}"])

    def to_dict(self):
        return {
            "description": self.description,
            "nodes": self.nodes(),
            "elements": self.elements(),
            "crs": self.crs}

    def write(self, path, overwrite=False):
        grd.write(self.to_dict(), path, overwrite)

    def get_x(self, crs: Union[CRS, str] = None):
        return self.get_xy(crs)[:, 0]

    def get_y(self, crs: Union[CRS, str] = None):
        return self.get_xy(crs)[:, 1]

    def get_xy(self, crs: Union[CRS, str] = None):
        if crs is not None:
            crs = CRS.from_user_input(crs)
            if not crs.equals(self.crs):
                transformer = Transformer.from_crs(
                    self.crs, crs, always_xy=True)
                x, y = transformer.transform(self.x, self.y)
                return np.vstack([x, y]).T
        return np.vstack([self.x, self.y]).T

    def get_bbox(self, crs: Union[CRS, str] = None) -> Bbox:
        vertices = self.get_xy(crs)
        x0 = np.min(vertices[:, 0])
        x1 = np.max(vertices[:, 0])
        y0 = np.min(vertices[:, 1])
        y1 = np.max(vertices[:, 1])
        return Bbox([[x0, y0], [x1, y1]])

    def transform_to(self, dst_crs):
        """Transforms coordinate system of mesh in-place.
        """
        dst_crs = CRS.from_user_input(dst_crs)
        if not self.crs.equals(dst_crs):
            xy = self.get_xy(dst_crs)
            self._nodes = {self.nodes.id()[i]:
                           (coord.tolist(), self.nodes.values()[i])
                           for i, coord in enumerate(xy)}
            self._crs = dst_crs

    def vertices_around_vertex(self, index):
        return self.nodes.vertices_around_vertex(index)

    @classmethod
    def open(cls, file: Union[str, os.PathLike],
             crs: Union[str, CRS] = None):
        return cls(**grd.read(pathlib.Path(file), boundaries=False))

    @figure
    def tricontourf(self, axes=None, show=True, figsize=None, **kwargs):
        if len(self.triangles) > 0:
            axes.tricontourf(self.triangulation, self.values, **kwargs)
        return axes

    @figure
    def tripcolor(self, axes=None, show=True, figsize=None, **kwargs):
        if len(self.triangles) > 0:
            axes.tripcolor(self.triangulation, self.values, **kwargs)
        return axes

    @figure
    def triplot(
        self,
        axes=None,
        show=False,
        figsize=None,
        linewidth=0.07,
        color='black',
        **kwargs
    ):
        if len(self.triangles) > 0:
            kwargs.update({'linewidth': linewidth})
            kwargs.update({'color': color})
            axes.triplot(self.triangulation, **kwargs)
        return axes

    @figure
    def quadplot(
        self,
        axes=None,
        show=False,
        figsize=None,
        facecolor='none',
        edgecolor='k',
        linewidth=0.07,
        **kwargs
    ):
        if len(self.quads) > 0:
            pc = PolyCollection(
                self.coords[self.quads],
                facecolor=facecolor,
                edgecolor=edgecolor,
                linewidth=0.07,
            )
            axes.add_collection(pc)
        return axes

    @figure
    def quadface(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        if len(self.quads) > 0:
            pc = PolyCollection(
                self.coords[self.quads],
                **kwargs
            )
            quad_value = np.mean(self.values[self.quads], axis=1)
            pc.set_array(quad_value)
            axes.add_collection(pc)
        return axes

    @figure
    def plot_wireframe(self, axes=None, show=False, **kwargs):
        axes = self.triplot(axes=axes, **kwargs)
        axes = self.quadplot(axes=axes, **kwargs)
        return axes

    @property
    def nodes(self):
        return self._nodes

    @property
    def coords(self):
        return self._nodes.coord()

    @property
    def vertices(self):
        return self._nodes.coord()

    @property
    def vertex_id(self):
        return self._nodes.id()

    @property
    def elements(self):
        return self._elements

    @property
    def element_id(self):
        return self._elements.id()

    @property
    def values(self):
        return self._nodes.values()

    @property
    def description(self):
        return self._description

    @property
    def crs(self):
        return self._crs

    @property
    def x(self):
        return self._nodes.coord()[:, 0]

    @property
    def y(self):
        return self._nodes.coord()[:, 1]

    @property
    def hull(self):
        return self._hull

    @property
    def triangles(self):
        return self._elements.triangles()

    @property
    def quads(self):
        return self._elements.quads()

    @property
    def triangulation(self):
        return self._elements.triangulation()

    @property
    def bbox(self):
        return self.get_bbox()


def edges_to_rings(edges):
    if len(edges) == 0:
        return edges
    # start ordering the edges into linestrings
    edge_collection = list()
    ordered_edges = [edges.pop(-1)]
    e0, e1 = [list(t) for t in zip(*edges)]
    while len(edges) > 0:
        if ordered_edges[-1][1] in e0:
            idx = e0.index(ordered_edges[-1][1])
            ordered_edges.append(edges.pop(idx))
        elif ordered_edges[0][0] in e1:
            idx = e1.index(ordered_edges[0][0])
            ordered_edges.insert(0, edges.pop(idx))
        elif ordered_edges[-1][1] in e1:
            idx = e1.index(ordered_edges[-1][1])
            ordered_edges.append(
                list(reversed(edges.pop(idx))))
        elif ordered_edges[0][0] in e0:
            idx = e0.index(ordered_edges[0][0])
            ordered_edges.insert(
                0, list(reversed(edges.pop(idx))))
        else:
            edge_collection.append(tuple(ordered_edges))
            idx = -1
            ordered_edges = [edges.pop(idx)]
        e0.pop(idx)
        e1.pop(idx)
    # finalize
    if len(edge_collection) == 0 and len(edges) == 0:
        edge_collection.append(tuple(ordered_edges))
    else:
        edge_collection.append(tuple(ordered_edges))
    return edge_collection


def sort_rings(index_rings, vertices):
    """Sorts a list of index-rings.

    Takes a list of unsorted index rings and sorts them into an "exterior" and
    "interior" components. Any doubly-nested rings are considered exterior
    rings.

    TODO: Refactor and optimize. Calls that use :class:matplotlib.path.Path can
    probably be optimized using shapely.
    """

    # sort index_rings into corresponding "polygons"
    areas = list()
    for index_ring in index_rings:
        e0, e1 = [list(t) for t in zip(*index_ring)]
        areas.append(float(Polygon(vertices[e0, :]).area))

    # maximum area must be main mesh
    idx = areas.index(np.max(areas))
    exterior = index_rings.pop(idx)
    areas.pop(idx)
    _id = 0
    _index_rings = dict()
    _index_rings[_id] = {
        'exterior': np.asarray(exterior),
        'interiors': []
    }
    e0, e1 = [list(t) for t in zip(*exterior)]
    path = Path(vertices[e0 + [e0[0]], :], closed=True)
    while len(index_rings) > 0:
        # find all internal rings
        potential_interiors = list()
        for i, index_ring in enumerate(index_rings):
            e0, e1 = [list(t) for t in zip(*index_ring)]
            if path.contains_point(vertices[e0[0], :]):
                potential_interiors.append(i)
        # filter out nested rings
        real_interiors = list()
        for i, p_interior in reversed(
                list(enumerate(potential_interiors))):
            _p_interior = index_rings[p_interior]
            check = [index_rings[k]
                     for j, k in
                     reversed(list(enumerate(potential_interiors)))
                     if i != j]
            has_parent = False
            for _path in check:
                e0, e1 = [list(t) for t in zip(*_path)]
                _path = Path(vertices[e0 + [e0[0]], :], closed=True)
                if _path.contains_point(vertices[_p_interior[0][0], :]):
                    has_parent = True
            if not has_parent:
                real_interiors.append(p_interior)
        # pop real rings from collection
        for i in reversed(sorted(real_interiors)):
            _index_rings[_id]['interiors'].append(
                np.asarray(index_rings.pop(i)))
            areas.pop(i)
        # if no internal rings found, initialize next polygon
        if len(index_rings) > 0:
            idx = areas.index(np.max(areas))
            exterior = index_rings.pop(idx)
            areas.pop(idx)
            _id += 1
            _index_rings[_id] = {
                'exterior': np.asarray(exterior),
                'interiors': []
            }
            e0, e1 = [list(t) for t in zip(*exterior)]
            path = Path(vertices[e0 + [e0[0]], :], closed=True)
    return _index_rings


def signed_polygon_area(vertices):
    # https://code.activestate.com/recipes/578047-area-of-polygon-using-shoelace-formula/
    n = len(vertices)  # of vertices
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += vertices[i][0] * vertices[j][1]
        area -= vertices[j][0] * vertices[i][1]
        return area / 2.0
