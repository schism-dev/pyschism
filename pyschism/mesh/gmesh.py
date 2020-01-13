import uuid
import numpy as np
from pyschism.mesh.figures import _figure
from collections.abc import Iterable, Mapping
from matplotlib.tri import Triangulation
from matplotlib.collections import PolyCollection
from pyproj import CRS, Transformer
from functools import lru_cache


class Gmesh:
    """
    Base class for all geographic unstructured mesh.
    """
    def __init__(
        self,
        nodes,
        elements,
        crs=None,
        description=None,
    ):
        self._nodes = nodes
        self._elements = elements
        self._crs = crs
        self._description = description

    def transform_to(self, dst_crs):
        dst_crs = CRS.from_user_input(dst_crs)
        if self.srs != dst_crs.srs:
            transformer = Transformer.from_crs(
                self.crs, dst_crs, always_xy=True)
            self._vertices = list(zip(*transformer.transform(self.x, self.y)))
            self._crs = dst_crs

    @_figure
    def tricontourf(self, axes=None, show=True, figsize=None, **kwargs):
        if len(self.triangles) > 0:
            axes.tricontourf(self.triangulation, self.values, **kwargs)
        return axes

    @_figure
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

    @_figure
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
                self.vertices[self.quads],
                facecolor=facecolor,
                edgecolor=edgecolor,
                linewidth=0.07,
                )
            axes.add_collection(pc)
        return axes

    @_figure
    def quadface(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        if len(self.quads) > 0:
            pc = PolyCollection(
                self.vertices[self.quads],
                **kwargs
                )
            quad_values = np.mean(self.values[self.quads], axis=1)
            pc.set_array(quad_values)
            axes.add_collection(pc)
        return axes

    @_figure
    def plot_wireframe(self, axes=None, **kwargs):
        self.triplot(axes=axes, **kwargs)
        self.quadplot(axes=axes, **kwargs)
        return axes

    @property
    def nodes(self):
        return self._nodes

    @property
    def elements(self):
        return self._elements

    @property
    @lru_cache
    def vertices(self):
        return np.array(
            [(x, y) for x, y, _ in self.nodes.values()]
            )

    @property
    @lru_cache
    def values(self):
        return np.array(
            [z for _, _, z in self.nodes.values()]
            )

    @property
    def xy(self):
        return self.vertices

    @property
    def x(self):
        return self.vertices[:, 0]

    @property
    def y(self):
        return self.vertices[:, 1]

    @property
    @lru_cache
    def triangulation(self):
        if len(self.triangles) > 0:
            return Triangulation(self.x, self.y, self.triangles)

    @property
    @lru_cache
    def triangles(self):
        return self._get_geom_by_type(3)

    @property
    @lru_cache
    def quads(self):
        return self._get_geom_by_type(4)

    @property
    def description(self):
        return self._description

    @property
    def crs(self):
        return self._crs

    @description.setter
    def description(self, description):
        self._description = description

    @property
    def _nodes(self):
        return self.__nodes

    @property
    def _elements(self):
        return self.__elements

    @property
    def _crs(self):
        return self.__crs

    @property
    def _description(self):
        return self.__description

    @_nodes.setter
    def _nodes(self, nodes):
        msg = "nodes argument must be a dictionary of the form "
        msg += "\\{node_id:  (x, y, z)\\}"
        assert isinstance(nodes, Mapping), msg
        for node in nodes.values():
            assert len(node) == 3, msg
            for coord in node:
                assert isinstance(coord, (float, int)), msg
        self.__nodes = nodes

    @_elements.setter
    def _elements(self, elements):
        msg = "elemens argument must be a dictionary of the form "
        msg += "\\{element_id:  (e0, ..., en)\\} where n==2 or n==3."
        assert isinstance(elements, Mapping), msg
        for geom in elements.values():
            msg = f"Found an element with {len(geom)} sides. "
            msg += "Only triangles and quadrilaterals are supported."
            assert len(geom) in [3, 4], msg
            # msg = "Boundary indexes must be a subset of the node id's."
            # assert set(geom).issubset(self.nodes.keys()), msg
        self.__elements = elements

    @_crs.setter
    def _crs(self, crs):
        self.__crs = crs

    @_description.setter
    def _description(self, description):
        if description is None:
            description = str(uuid.uuid4())[:8]
        self.__description = description

    def _get_geom_by_type(self, lenght):
        geom = list()
        for _geom in self.elements.values():
            _geom = list(map(int, _geom))
            if len(_geom) == lenght:
                geom.append([idx-1 for idx in _geom])
        return np.array(geom)
