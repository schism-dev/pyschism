from abc import ABC
from collections import defaultdict
import os
import pathlib
from typing import Union, Sequence, Hashable, List

from matplotlib.collections import PolyCollection  # type: ignore[import]
from matplotlib.path import Path  # type: ignore[import]
from matplotlib.tri import Triangulation  # type: ignore[import]
from matplotlib.transforms import Bbox  # type: ignore[import]
import numpy as np  # type: ignore[import]
from pyproj import Transformer, CRS  # type: ignore[import]
from shapely.geometry import MultiPolygon, Polygon  # type: ignore[import]

from pyschism.mesh.parsers import grd
from pyschism import figures as fig


class Vertices:

    def __set__(self, obj, vertices):
        vertices = np.array(vertices)
        if vertices.shape[1] != 2:
            raise ValueError('Argument vertices must be castable to a Nx2 '
                             'matrix, that is, only 2-dimensional meshes are '
                             f'supported. Got shape {vertices.shape}')
        obj.__dict__['vertices'] = vertices

    def __get__(self, obj, val):
        return obj.__dict__['vertices']


class VertexId:

    def __set__(self, obj, vertex_id):
        vertex_id = np.array(vertex_id).flatten()
        if np.max(vertex_id.shape) != np.max(obj.vertices.shape):
            raise ValueError(f'Argument vertex_id has dimension mismatch with '
                             f'vertices attribute. Got shape {vertex_id.shape}'
                             ' but vertices attribute has shape '
                             f'{self.vertices.shape}.')
        obj.__dict__['vertex_id'] = vertex_id.tolist()

    def __get__(self, obj, val):
        return obj.__dict__['vertex_id']


class Elements:

    def __set__(self, obj, elements: Sequence):
        if not isinstance(elements, Sequence):
            raise TypeError('Argument elements must be iterable.')
        vertex_id_set = set(obj.vertex_id)
        for i, element in enumerate(elements):
            if not isinstance(element, Sequence):
                raise TypeError(f'Element with index {i} of the elements '
                                f'argument must be of type {Sequence}, not '
                                f'type {type(element)}.')
            if not set(element).issubset(vertex_id_set):
                ValueError(f'Element with index {i} is not a subset of the '
                           "coordinate id's.")
        obj.__dict__['elements'] = elements

    def __get__(self, obj, val):
        return obj.__dict__['elements']


class ElementId:

    def __set__(self, obj, element_id: Union[Sequence, None]):
        if element_id is None:
            element_id = list(range(len(obj._elements)))
        else:
            if not isinstance(element_id, Sequence):
                raise TypeError('Argument element_id must be iterable.')
            if len(element_id) != len(obj.elements):
                ValueError(
                    "Argument element_id has dimension mismatch with elements "
                    f'table. There are {obj.elements.shape[0]} elements, but '
                    f"{len(element_id)} element id's were given.")
        obj.__dict__['element_id'] = element_id

    def __get__(self, obj, val):
        return obj.__dict__['element_id']


class Values:

    def __set__(self, obj, values):
        if values is not None:
            if not isinstance(values, Sequence):
                raise TypeError('Argument values must be a sequence type, not '
                                f'type {type(values)}.')
            values = np.array(values)
            if obj._vertices.shape[0] != values.shape[0]:
                raise TypeError('Dimension mismatch: First dimension of values '
                                'argument must match the first dimension of the '
                                'vertices.')
        else:
            values = np.full(obj._vertices.shape, np.nan)
        obj.__dict__['values'] = np.ma.masked_equal(values, np.nan)

    def __get__(self, obj, val):
        return obj.__dict__['values']


class Header:

    def __set__(self, obj, description):
        if description is None:
            description = ''
        if not isinstance(description, str):
            raise TypeError('Argument description must be a string, not '
                            f'type {type(description)}.')
        obj.__dict__['description'] = description.replace('\n', ' ')

    def __get__(self, obj, val):
        return obj.__dict__['description']


class Crs:

    def __set__(self, obj, val):
        if val is not None:
            obj.__dict__['crs'] = CRS.from_user_input(val)

    def __get__(self, obj, val):
        return obj.__dict__.get('crs')


class Tria3:

    def __get__(self, obj, val):
        tria3 = obj.__dict__.get("tria3")
        if tria3 is None:
            tria3 = np.array(
                [list(map(obj.get_vertex_index_by_id, element)) for element
                 in obj.elements if len(element) == 3])
            obj.__dict__['tria3'] = tria3
        return tria3


class Quad4:

    def __get__(self, obj, val):
        if obj.__dict__.get("quad4") is None:
            obj.__dict__["quad4"] = np.array(
                [list(map(obj.get_vertex_index_by_id, element)) for element
                 in obj.elements if len(element) == 4])
        return obj.__dict__["quad4"]


class IndexRingCollection:

    def __get__(self, obj, val):

        if obj.__dict__.get("index_ring_collection") is not None:
            return obj.__dict__.get("index_ring_collection")

        boundary_edges = list()
        tri = obj.triangulation
        idxs = np.vstack(
            list(np.where(tri.neighbors == -1))).T
        for i, j in idxs:
            boundary_edges.append(
                (int(tri.triangles[i, j]),
                    int(tri.triangles[i, (j+1) % 3])))
        index_ring_collection = sort_edges(boundary_edges)
        # sort index_rings into corresponding "polygons"
        areas = list()
        for index_ring in index_ring_collection:
            e0, e1 = [list(t) for t in zip(*index_ring)]
            areas.append(float(Polygon(obj.vertices[e0, :]).area))

        # maximum area must be main mesh
        idx = areas.index(np.max(areas))
        exterior = index_ring_collection.pop(idx)
        areas.pop(idx)
        _id = 0
        _index_ring_collection = dict()
        _index_ring_collection[_id] = {
            'exterior': np.asarray(exterior),
            'interiors': []
            }
        e0, e1 = [list(t) for t in zip(*exterior)]
        path = Path(obj.vertices[e0 + [e0[0]], :], closed=True)
        while len(index_ring_collection) > 0:
            # find all internal rings
            potential_interiors = list()
            for i, index_ring in enumerate(index_ring_collection):
                e0, e1 = [list(t) for t in zip(*index_ring)]
                if path.contains_point(obj.vertices[e0[0], :]):
                    potential_interiors.append(i)
            # filter out nested rings
            real_interiors = list()
            for i, p_interior in reversed(
                    list(enumerate(potential_interiors))):
                _p_interior = index_ring_collection[p_interior]
                check = [index_ring_collection[_]
                         for j, _ in reversed(
                            list(enumerate(potential_interiors)))
                         if i != j]
                has_parent = False
                for _path in check:
                    e0, e1 = [list(t) for t in zip(*_path)]
                    _path = Path(obj.vertices[e0 + [e0[0]], :], closed=True)
                    if _path.contains_point(
                            obj.vertices[_p_interior[0][0], :]):
                        has_parent = True
                if not has_parent:
                    real_interiors.append(p_interior)
            # pop real rings from collection
            for i in reversed(sorted(real_interiors)):
                _index_ring_collection[_id]['interiors'].append(
                    np.asarray(index_ring_collection.pop(i)))
                areas.pop(i)
            # if no internal rings found, initialize next polygon
            if len(index_ring_collection) > 0:
                idx = areas.index(np.max(areas))
                exterior = index_ring_collection.pop(idx)
                areas.pop(idx)
                _id += 1
                _index_ring_collection[_id] = {
                    'exterior': np.asarray(exterior),
                    'interiors': []
                    }
                e0, e1 = [list(t) for t in zip(*exterior)]
                path = Path(obj.vertices[e0 + [e0[0]], :], closed=True)
        obj.__dict__['index_ring_collection'] = _index_ring_collection
        return obj.__dict__['index_ring_collection']


class OuterRingCollection:

    def __get__(self, obj, val):
        if obj.__dict__.get('outer_ring_collection') is None:
            obj.__dict__['outer_ring_collection'] = defaultdict()
            for key, ring in obj.index_ring_collection.items():
                obj.__dict__['outer_ring_collection'][key] = ring['exterior']
        return obj.__dict__['outer_ring_collection']


class InnerRingCollection:

    def __get__(self, obj, val):
        if obj.__dict__.get('inner_ring_collection') is None:
            obj.__dict__['inner_ring_collection'] = defaultdict()
            for key, rings in obj.index_ring_collection.items():
                obj.__dict__['inner_ring_collection'][key] = rings['interiors']
        return obj.__dict__['inner_ring_collection']


class TriangulationDescriptor:

    def __get__(self, obj, val):
        if obj.__dict__.get('triangulation') is None:
            obj.__dict__['triangulation'] = Triangulation(
                obj.vertices[:, 0], obj.vertices[:, 1],
                triangles=obj.tria3)
        return obj.__dict__['triangulation']


class Gr3(ABC):

    _vertices = Vertices()
    _vertex_id = VertexId()
    _elements = Elements()
    _element_id = ElementId()
    _values = Values()
    _description = Header()
    _crs = Crs()
    _tria3 = Tria3()
    _quad4 = Quad4()
    _index_ring_collection = IndexRingCollection()
    _outer_ring_collection = OuterRingCollection()
    _inner_ring_collection = InnerRingCollection()
    _triangulation = TriangulationDescriptor()

    def __init__(self, vertices, elements=None, values=None,
                 vertex_id=None, element_id=None, description=None, crs=None):
        self._vertices = vertices
        self._vertex_id = vertex_id
        self._elements = elements
        self._values = values
        self._element_id = element_id
        self._description = description
        self._crs = crs
        self._vertex_id_to_index = {self.vertex_id[i]: i for i in
                                    range(len(self.vertex_id))}
        self._vertex_index_to_id = {i: self.vertex_id[i] for i in
                                    range(len(self.vertex_id))}

    def __str__(self):
        return grd.dict_to_string(self.to_dict())

    def to_dict(self):
        return {
            "description": self.description,
            "vertices": self.vertices,
            "elements": self.elements,
            "vertex_id": self.vertex_id,
            "element_id": self.element_id,
            "values": self.values}

    def write(self, path, overwrite=False):
        grd.write(self.to_dict(), path, overwrite)

    def get_vertex_index_by_id(self, id: Hashable):
        return self._vertex_id_to_index[id]

    def get_vertex_id_by_index(self, index: int):
        return self._vertex_index_to_id[index]

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

    def get_multipolygon(self, crs: Union[CRS, str] = None) -> MultiPolygon:
        pol_col: List[Polygon] = []
        for id, ring in self.index_ring_collection.items():
            outer = self.get_xy(crs)[ring['exterior'][:, 0], :]
            inner: List[Polygon] = []
            for inner_ring in ring['interiors']:
                inner.append(self.get_xy(crs)[inner_ring[:, 0], :])
            pol_col.append(Polygon(outer, inner))
        return MultiPolygon(pol_col)

    def get_bbox(self, crs: Union[CRS, str] = None) -> Bbox:
        vertices = self.get_xy(crs)
        x0 = np.min(vertices[:, 0])
        x1 = np.max(vertices[:, 0])
        y0 = np.min(vertices[:, 1])
        y1 = np.max(vertices[:, 1])
        return Bbox([[x0, y0], [x1, y1]])

    @classmethod
    def open(cls, file: Union[str, os.PathLike],
             crs: Union[str, CRS] = None):
        return cls(**grd.read(pathlib.Path(file), boundaries=False))

    @fig._figure
    def tricontourf(self, axes=None, show=True, figsize=None, **kwargs):
        if len(self.tria3) > 0:
            axes.tricontourf(self.triangulation, self.values, **kwargs)
        return axes

    @fig._figure
    def tripcolor(self, axes=None, show=True, figsize=None, **kwargs):
        if len(self.tria3) > 0:
            axes.tripcolor(self.triangulation, self.values, **kwargs)
        return axes

    @fig._figure
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

    @fig._figure
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

    @fig._figure
    def quadface(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        if len(self.quad4) > 0:
            pc = PolyCollection(
                self.coords[self.quad4],
                **kwargs
                )
            quad_value = np.mean(self.values[self.quad4], axis=1)
            pc.set_array(quad_value)
            axes.add_collection(pc)
        return axes

    @fig._figure
    def plot_wireframe(self, axes=None, show=False, **kwargs):
        axes = self.triplot(axes=axes, **kwargs)
        axes = self.quadplot(axes=axes, **kwargs)
        return axes

    @property
    def coords(self):
        return self._vertices

    @property
    def vertices(self):
        return self._vertices

    @property
    def vertex_id(self):
        return self._vertex_id

    @property
    def elements(self):
        return self._elements

    @property
    def element_id(self):
        return self._element_id

    @property
    def values(self):
        return self._values

    @property
    def description(self):
        return self._description

    @property
    def crs(self):
        return self._crs

    @property
    def x(self):
        return self._vertices[:, 0]

    @property
    def y(self):
        return self._vertices[:, 1]

    @property
    def tria3(self):
        return self._tria3

    @property
    def quad4(self):
        return self._quad4

    @property
    def index_ring_collection(self):
        return self._index_ring_collection

    @property
    def outer_ring_collection(self):
        return self._outer_ring_collection

    @property
    def inner_ring_collection(self):
        return self._inner_ring_collection

    @property
    def triangulation(self):
        return self._triangulation

    @property
    def bbox(self):
        return self.get_bbox()


def sort_edges(edges):
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


def signed_polygon_area(vertices):
    # https://code.activestate.com/recipes/578047-area-of-polygon-using-shoelace-formula/
    n = len(vertices)  # of vertices
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += vertices[i][0] * vertices[j][1]
        area -= vertices[j][0] * vertices[i][1]
        return area / 2.0
