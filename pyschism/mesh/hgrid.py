import numpy as np
from functools import lru_cache
from collections import defaultdict
from collections.abc import Iterable, Mapping
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.collections import PolyCollection
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pyproj import CRS, Transformer
from pyschism.mesh.friction import ManningsN, DragCoefficient, RoughnessLength
from pyschism.mesh.gr3 import parse_gr3, write_gr3
import pyschism.mesh.figures as fig


class Hgrid:

    def __init__(
        self,
        nodes,
        elements,
        boundaries=None,
        crs=None,
        description=None,
    ):
        self._nodes = nodes
        self._elements = elements
        self._boundaries = boundaries
        self._crs = crs
        self._description = description

    @classmethod
    def open(cls, hgrid, crs=None):
        kwargs = cls.parse(hgrid)
        kwargs.update({"crs": crs})
        return cls(**kwargs)

    def set_friction(self, value, ftype='manning'):

        # certify ftype
        _ftypes = {
            'manning': ManningsN,
            'drag': DragCoefficient,
            'rough': RoughnessLength
        }
        msg = f"ftype argument must be one of {_ftypes.keys()}"
        assert ftype in _ftypes, msg

        # certify value
        msg = "value argument must be an instance of type "
        msg += f"{int}, {float} or an iterable ."
        assert isinstance(value, (Iterable, int, float)), msg

        _ftypes = {}

        if isinstance(value, (int, float)):
            if ftype == 'manning':
                self._fgrid = ManningsN.constant(self, value)

        return self.fgrid

    def add_boundary_type(self, ibtype):
        if ibtype not in self.boundaries:
            self._boundaries[ibtype] = defaultdict()

    def add_boundary_data(self, ibtype, id, indexes, **properties):
        assert set(indexes).issubset(set(self.vertices.keys()))
        self._boundary[ibtype] = {
            'indexes': indexes,
            **properties
        }

    def transform_to(self, dst_crs):
        dst_crs = CRS.from_user_input(dst_crs)
        if self.srs != dst_crs.srs:
            transformer = Transformer.from_crs(
                self.crs, dst_crs, always_xy=True)
            self._vertices = list(zip(*transformer.transform(self.x, self.y)))
            self._crs = dst_crs

    def dump(self, path, overwrite=False):
        grd = {
            'description': self.description,
            'nodes': self.nodes,
            'elements': self.elements,
            'boundaries': self.boundaries,
        }
        write_gr3(grd, path, overwrite)

    def _figure(f):
        def decorator(*argv, **kwargs):
            axes = fig.get_axes(
                kwargs.get('axes', None),
                kwargs.get('figsize', None)
                )
            kwargs.update({'axes': axes})
            axes = f(*argv, **kwargs)
            if kwargs.get('show', False):
                axes.axis('scaled')
                plt.show()
            return axes
        return decorator

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

    @_figure
    def make_plot(
        self,
        axes=None,
        vmin=None,
        vmax=None,
        show=False,
        title=None,
        figsize=rcParams["figure.figsize"],
        extent=None,
        cbar_label=None,
        **kwargs
    ):
        if vmin is None:
            vmin = np.min(self.values)
        if vmax is None:
            vmax = np.max(self.values)
        kwargs.update(**fig.get_topobathy_kwargs(self.values, vmin, vmax))
        col_val = kwargs.pop('col_val')
        self.tricontourf(axes=axes, vmin=vmin, vmax=vmax, **kwargs)
        kwargs.pop('levels')
        self.quadface(axes=axes, **kwargs)
        axes.axis('scaled')
        if extent is not None:
            axes.axis(extent)
        if title is not None:
            axes.set_title(title)
        mappable = ScalarMappable(cmap=kwargs['cmap'])
        mappable.set_array([])
        mappable.set_clim(vmin, vmax)
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("bottom", size="2%", pad=0.5)
        cbar = plt.colorbar(
            mappable,
            cax=cax,
            # extend=cmap_extend,
            orientation='horizontal'
        )
        if col_val != 0:
            cbar.set_ticks([vmin, vmin + col_val * (vmax-vmin), vmax])
            cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
        else:
            cbar.set_ticks([vmin, vmax])
            cbar.set_ticklabels([np.around(vmin, 2), np.around(vmax, 2)])
        if cbar_label is not None:
            cbar.set_label(cbar_label)
        return axes

    @_figure
    def plot_boundary(
        self,
        ibtype,
        id,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        boundary = [int(idx)-1 for idx in self.boundaries[ibtype][id]]
        axes.plot(self.x[boundary], self.y[boundary], **kwargs)
        return axes

    @_figure
    def plot_boundaries(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        kwargs.update({'axes': axes})
        for ibtype, bnds in self.boundaries.items():
            for bnd in bnds:
                axes = self.plot_boundary(ibtype, bnd, **kwargs)
                kwargs.update({'axes': axes})
        return kwargs['axes']

    @staticmethod
    def parse(path):
        return parse_gr3(path)

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
    def triangles(self):
        try:
            return self.__triangles
        except AttributeError:
            self.__triangles = self._get_geom_by_length(3)
            return self.__triangles

    @property
    @lru_cache
    def quads(self):
        try:
            return self.__quads
        except AttributeError:
            self.__quads = self._get_geom_by_length(4)
            return self.__quads

    @property
    def boundaries(self):
        return self._boundaries.copy()

    @property
    def triangulation(self):
        return Triangulation(self.x, self.y, self.triangles)

    @property
    def crs(self):
        return self._crs

    @property
    def description(self):
        return self._description

    @property
    def fgrid(self):
        return self._fgrid

    def _get_geom_by_length(self, lenght):
        geom = list()
        for _geom in self.elements.values():
            _geom = list(map(int, _geom))
            if len(_geom) == lenght:
                geom.append([idx-1 for idx in _geom])
        return np.array(geom)

    @property
    def _nodes(self):
        return self.__nodes

    @property
    def _elements(self):
        return self.__elements

    @property
    def _boundaries(self):
        return self.__boundaries

    @property
    def _crs(self):
        return self.__crs

    @property
    def _description(self):
        return self.__description

    @property
    def _fgrid(self):
        return self.__fgrid

    @description.setter
    def description(self, description):
        self._description = description

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

    @_boundaries.setter
    def _boundaries(self, boundaries):
        if boundaries is not None:
            msg = "elemens argument must be a dictionary of the form "
            msg += "\\{element_id:  (e0, ..., en)\\} where n==2 or n==3."
            assert isinstance(boundaries, Mapping), msg
            # for bnds in boundaries.values():
            #     for geom in bnds.values():
            #         msg = "Boundary indexes must be a subset of the node id's."
            #         assert set(geom).issubset(self.nodes), msg
        else:
            boundaries = {}
        # ocean boundaries
        if None not in boundaries:
            boundaries[None] = {}
        # land boundaries
        if 0 not in boundaries:
            boundaries[0] = {}
        # interior boundaries
        if 1 not in boundaries:
            boundaries[1] = {}

        self.__boundaries = boundaries

    @_crs.setter
    def _crs(self, crs):
        self.__crs = crs

    @_description.setter
    def _description(self, description):
        self.__description = description

    @_fgrid.setter
    def _fgrid(self, fgrid):
        self.__fgrid = fgrid
