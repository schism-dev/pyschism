import numpy as np
from pyproj import Proj, CRS, Transformer
from functools import lru_cache
import logging
import uuid
import pathlib
from pyschism.mesh import gr3
from pyschism.figures import _figure
from collections.abc import Mapping
from matplotlib.tri import Triangulation
from matplotlib.collections import PolyCollection
from matplotlib.transforms import Bbox


class EuclideanMesh2D:

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

    @classmethod
    def open(cls, path, crs=None, fmt="grd"):
        assert fmt.lower() in ['grd', 'gr3', 'adcirc', 'schism', '2dm', 'sms', 'msh']
        
        if fmt.lower() in ['grd', 'gr3', 'adcirc', 'schism']:
            return cls.open_gr3(path, crs)
        
        elif fmt.lower() in ['msh']:
            return cls.open_msh(path, crs)

        elif fmt.lower() in ['2dm', 'sms']:
            return cls.open_2dm(path, crs)

    @classmethod
    def open_2dm(cls, path, crs=None):
        raise NotImplementedError

    @classmethod
    def open_gr3(cls, path, crs=None):
        grd = gr3.reader(path)
        grd.pop('boundaries')
        return cls.from_grd(grd, crs)

    @classmethod
    def from_grd(cls, grd, crs=None):
        grd = gr3.to_gmesh(grd)
        if 'boundaries' in grd:
            grd.pop('boundaries')
        return cls(**gr3.to_gmesh(grd), crs=crs)

    def transform_to(self, dst_crs):
        dst_crs = CRS.from_user_input(dst_crs)
        if self.srs != dst_crs.srs:
            transformer = Transformer.from_crs(
                self.crs, dst_crs,
                always_xy=True
                )
            xy = list(zip(*transformer.transform(self.x, self.y)))
            ids = list(self._coords.keys())
            self._coords = {
                ids[i]: coord for i, coord in enumerate(xy)
            }
            self._crs = dst_crs
    
    def get_node_index(self, id):
        return self.node_index[id]

    def get_node_id(self, index):
        return self.node_id[index]

    def get_element_index(self, id):
        return self.element_index[id]

    def get_element_id(self, index):
        return self.element_id[index]
    
    def get_x(self, crs=None):
        return self.get_xy(crs)[:, 0]

    def get_y(self, crs=None):
        return self.get_xy(crs)[:, 1]

    def get_xy(self, crs=None):
        if crs is not None:
            crs = CRS.from_user_input(crs)
            transformer = Transformer.from_crs(self.crs, crs, always_xy=True)
            x, y = transformer.transform(self.x, self.y)
            return np.vstack([x, y]).T
        else:
            return np.vstack([self.x, self.y]).T

    def write(self, path, overwrite=False, fmt='gr3'):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            msg = 'File exists, use overwrite=True to allow overwrite.'
            raise Exception(msg)
        with open(path, 'w') as f:
            f.write(self.ascii_string(fmt))


    def ascii_string(self, fmt):
        if fmt.lower() in ['grd', 'gr3', 'adcirc', 'schism']:
            return self.gr3
        if fmt.lower() in ['sms', '2dm']:
            return self.sms

    @_figure
    def tricontourf(self, axes=None, show=True, figsize=None, **kwargs):
        if len(self.tria3) > 0:
            axes.tricontourf(self.triangulation, self.values, **kwargs)
        return axes

    @_figure
    def tripcolor(self, axes=None, show=True, figsize=None, **kwargs):
        if len(self.tria3) > 0:
            axes.tripcolor(self.triangulation, self.values, **kwargs)
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
        if len(self.quad4) > 0:
            pc = PolyCollection(
                self.coords[self.quad4],
                **kwargs
                )
            quad_value = np.mean(self.values[self.quad4], axis=1)
            pc.set_array(quad_value)
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
    def xy(self):
        return self.coords

    @property
    def x(self):
        return self.coords[:, 0]

    @property
    def y(self):
        return self.coords[:, 1]

    @property
    @lru_cache
    def values(self):
        return self._values

    @property
    def triangulation(self):
        if len(self.tria3) > 0:
            return Triangulation(self.x, self.y, self.tria3)

    @property
    def triangles(self):
        return self.tria3

    @property
    def quads(self):
        return self.quad4

    @property
    @lru_cache
    def nodes(self):
        return list(self._nodes.items())

    @property
    @lru_cache
    def elements(self):
        return list(self._elements.values())

    @property
    @lru_cache
    def coords_id(self):
        return self._coords.keys()

    @property
    @lru_cache
    def triangles_id(self):
        return self._triangles.keys()

    @property
    @lru_cache
    def quads_id(self):
        return self._quads.keys()

    @property
    @lru_cache
    def node_id(self):
        return {index: id for index, id in enumerate(self._coords)}

        
    @property
    @lru_cache
    def element_id(self):
        return {index: id for index, id in enumerate(self._elements)}


    @property
    @lru_cache
    def node_index(self):
        return {id: index for index, id in enumerate(self._coords)}


    @property
    @lru_cache
    def element_index(self):
        return {id: index for index, id in enumerate(self._elements)}

    @property
    def bbox(self):
        x0 = np.min(self.x)
        x1 = np.max(self.x)
        y0 = np.min(self.y)
        y1 = np.max(self.y)
        return Bbox([[x0, y0], [x1, y1]])

    @staticmethod
    def parse_2dm(path):
        grd = dict()
        with open(pathlib.Path(path), 'r') as f:
            f.readline()
            while 1:
                line = f.readline().split()
                if len(line) == 0:
                    break
                if line[0] in ['E3T', 'E4Q']:
                    grd[line[0]].update({
                        line[1]: line[2:]
                        })
                if line[0] == 'ND':
                    grd[line[0]].update({
                        line[1]: (
                            list(map(float, line[2:-1])), float(line[-1])
                            )
                        })
        return grd

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

    @property
    def gr3(self):
        f = f"{self.description}\n"
        f += f"{len(self.elements)} "
        f += f"{len(self.nodes)}\n"
        # TODO: Make faster using np.array2string
        for id, ((x, y), value) in self._nodes.items():
            f += f"{id} "
            f += f"{x:<.16E} "
            f += f" {y:<.16E} "
            f += f"{value:<.16E}\n"
        for id, element in self._elements.items():
            f += f"{id} "
            f += f"{len(element):d} "
            for tag in element:
                f += f"{tag} "
            f += "\n"
        return f

    @property
    def sms(self):
        def geom_string(geom_type, geom):
            f = ''
            if geom is not None:
                for i in range(geom.shape[0]):
                    f += f"{geom_type} {i + 1} "
                    for j in range(len(geom[i, :])):
                        f += f"{geom[i, j]+1} "
                    f += "\n"
            return f

        f = "MESH2D\n"
        # TODO: Make faster using np.array2string
        f += geom_string("E3T", self.triangles)
        f += geom_string("E4Q", self.quads)
        for i in range(self.coords.shape[0]):
            f += f"ND {i + 1} "
            f += f"{self.x[i]:<.16E} "
            f += f"{self.y[i]:<.16E} "
            f += f"{-self.values[i]:<.16E}\n"
        return f

    @property
    @lru_cache
    def grd(self):
        description = self.description
        if self.crs is not None:
            description += f"; {self.crs.srs}"
        return {
            "nodes": self.nodes,
            "elements": self.elements,
            "description": description,
        }
    
    @property
    @lru_cache
    def tria3(self):
        return np.array(
            [list(map(self.get_node_index, index))
             for id, index in self._triangles.items()])

    @property
    @lru_cache
    def quad4(self):
        return np.array(
            [list(map(self.get_node_index, index))
             for id, index in self._quads.items()])
   
    @property
    def logger(self):
        try:
            return self.__logger
        except AttributeError:
            self.__logger = logging.getLogger(
                __name__ + '.' + self.__class__.__name__)
            return self.__logger

    @description.setter
    def description(self, description):
        self._description = description

    def _certify_input_geom(self, geom_type, geom):
        geom_types = {
            "triangles": 3,
            "quads": 4,
        }
        _geom = np.asarray(list(geom.values()))
        if len(_geom) > 0:
            msg = f'Invalid shape for geom {_geom.shape}.'
            assert _geom.shape[1] == geom_types[geom_type], msg
            msg = f"{geom_type} members must be a subset of coords keys"
            for IDtags in geom.values():
                for IDtag in IDtags:
                    assert IDtag in self.coords_id, msg

    @property
    def _coords(self):
        return self.__coords

    @property
    def _triangles(self):
        return self.__triangles

    @property
    def _quads(self):
        return self.__quads

    @property
    def _values(self):
        return self.__values

    @property
    def _crs(self):
        return self.__crs

    @property
    def _description(self):
        return self.__description

    @property
    @lru_cache
    def _nodes(self):
        return {id: ((x, y), -self.values[i]) for i, (id, (x, y))
                in enumerate(self._coords.items())}

    @property
    @lru_cache
    def _elements(self):
        elements_id = list()
        elements_id.extend(list(self._triangles.keys()))
        elements_id.extend(list(self._quads.keys()))
        elements_id = range(1, len(elements_id)+1) \
            if len(set(elements_id)) != len(elements_id) else elements_id
        elements = list()
        elements.extend(list(self._triangles.values()))
        elements.extend(list(self._quads.values()))
        elements = {
            elements_id[i]: indexes for i, indexes in enumerate(elements)}
        return elements

    @_coords.setter
    def _coords(self, coords):
        msg = "coord argument must be a dictionary of the form "
        msg += "\\{coord_id:  (x, y)\\}"
        assert isinstance(coords, Mapping), msg
        for coord in coords.values():
            assert len(coord) == 2, msg
            assert isinstance(coord[0], (float, int)), msg
            assert isinstance(coord[1], (float, int)), msg
        self.__coords = coords
        type(self).coords.fget.cache_clear()

    @_values.setter
    def _values(self, values):
        if values is None:
            values = []
        values = np.asarray(values)
        if len(values) > 0:
            if len(values.shape) in [1, 2]:
                assert values.shape[0] == self.coords.shape[0]
            elif len(values.shape) == 3:
                assert values.shape[1] == self.coords.shape[0]
            else:
                msg = f"input values has invalid shape: {values.shape}"
                raise NotImplementedError(msg)

        else:
            values = np.full(self.coords.shape[0], np.nan)
        self.__values = values
        type(self).values.fget.cache_clear()

    @_triangles.setter
    def _triangles(self, triangles):
        if triangles is None:
            triangles = {}
        self._certify_input_geom("triangles", triangles)
        self.__triangles = triangles

    @_quads.setter
    def _quads(self, quads):
        if quads is None:
            quads = {}
        self._certify_input_geom("quads", quads)
        self.__quads = quads


    @_crs.setter
    def _crs(self, crs):
        if crs is not None:
            crs = CRS.from_user_input(crs)
        self.__crs = crs

    @_description.setter
    def _description(self, description):
        if description is None:
            description = uuid.uuid4().hex[:8]
        assert isinstance(description, str)
        self.__description = description
