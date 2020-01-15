import uuid
import numpy as np
from pyschism.mesh.figures import _figure
from collections.abc import Mapping
from matplotlib.tri import Triangulation
from matplotlib.collections import PolyCollection
from pyproj import CRS, Transformer, Proj
from functools import lru_cache


class Geomesh:
    """
    Base class for all geographic unstructured mesh.
    """
    def __init__(
        self,
        coords,
        triangles=None,
        quads=None,
        values=None,
        crs=None,
        description=None,
    ):
        self._coords = coords
        self._triangles = triangles
        self._quads = quads
        self._values = values
        self._crs = crs
        self._description = description

    def transform_to(self, dst_crs):
        dst_crs = CRS.from_user_input(dst_crs)
        if self.srs != dst_crs.srs:
            transformer = Transformer.from_crs(
                self.crs, dst_crs, always_xy=True)
            self._vertices = list(zip(*transformer.transform(self.x, self.y)))
            self._crs = dst_crs

    def get_index(self, id):
        return self.node_indexes[id]

    @classmethod
    def from_gr3(cls, nodes, elements):
        return cls(*cls._gr3_to_mesh(nodes, elements))

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
                self.coords[self.quads],
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
                self.coords[self.quads],
                **kwargs
                )
            quad_values = np.mean(self.values[self.quads], axis=1)
            pc.set_array(quad_values)
            axes.add_collection(pc)
        return axes

    @_figure
    def plot_wireframe(self, axes=None, show=False, **kwargs):
        axes = self.triplot(axes=axes, **kwargs)
        axes = self.quadplot(axes=axes, **kwargs)
        return axes

    @property
    @lru_cache
    def coords(self):
        return np.array(
            [coord for coord in self._coords.values()]
            )

    @property
    @lru_cache
    def node_indexes(self):
        return {id: index for index, id in enumerate(self._coords)}

    @property
    @lru_cache
    def triangles(self):
        return np.array(
            [list(map(self.get_index, index))
             for index in self._triangles.values()])

    @property
    @lru_cache
    def quads(self):
        return np.array(
            [list(map(self.get_index, index))
             for index in self._quads.values()])

    @property
    def xy(self):
        return self.coords

    @property
    def x(self):
        return self.coords[:, 0]

    @property
    def y(self):
        return self.coords[:, 1]

    @property
    def values(self):
        return self._values

    @property
    @lru_cache
    def triangulation(self):
        if len(self.triangles) > 0:
            return Triangulation(self.x, self.y, self.triangles)

    @property
    def description(self):
        return self._description

    @property
    def proj(self):
        return Proj(self.crs)

    @property
    def srs(self):
        return self.proj.srs

    @property
    def crs(self):
        return self._crs

    @description.setter
    def description(self, description):
        self._description = description

    @property
    def _coords(self):
        return self.__coords

    @property
    def _values(self):
        return self.__values

    @property
    def _crs(self):
        return self.__crs

    @property
    def _description(self):
        return self.__description

    @staticmethod
    def _gr3_to_mesh(nodes, elements):
        # cast gr3 inputs into a geomesh structure format
        coords = {id: (x, y) for id, ((x, y), value) in nodes.items()}
        triangles = {id: geom for id, geom in elements.items()
                     if len(geom) == 3}
        quads = {id: geom for id, geom in elements.items()
                 if len(geom) == 4}
        values = [value for coord, value in nodes.values()]
        return coords, triangles, quads, values

    @_coords.setter
    def _coords(self, coords):
        msg = "coords argument must be a dictionary of the form "
        msg += "\\{coord_id:  (x, y)\\}"
        assert isinstance(coords, Mapping), msg
        for coord in coords.values():
            assert len(coord) == 2, msg
            for coord in coord:
                assert isinstance(coord, (float, int)), msg
        self.__coords = coords

    @_values.setter
    def _values(self, values):
        if values is None:
            values = np.full(self.coords.shape[0], np.nan)
        values = np.asarray(values)
        assert values.shape[0] == self.coords.shape[0]
        self.__values = values

    @_crs.setter
    def _crs(self, crs):
        if crs is not None:
            crs = CRS.from_user_input(crs)
        self.__crs = crs

    @_description.setter
    def _description(self, description):
        if description is None:
            description = str(uuid.uuid4())[:8]
        self.__description = description


Gmesh = Geomesh  # Namespace alias
