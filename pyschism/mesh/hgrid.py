import numpy as np
from functools import lru_cache
from collections.abc import Iterable
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.tri import Triangulation
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pyproj import CRS, Transformer
from pyschism.mesh.friction import ManningsN, DragCoefficient, RoughnessLength
from pyschism.mesh.gr3 import parse_gr3, write_gr3
import pyschism.mesh.figures as fig


class Hgrid:

    def __init__(
        self,
        vertices,
        values,
        triangles=None,
        quads=None,
        ocean_boundaries=None,
        land_boundaries=None,
        interior_boundaries=None,
        crs=None,
        description=None,
    ):
        self._vertices = vertices
        self._values = values
        self._triangles = triangles
        self._quads = quads
        self._ocean_boundaries = ocean_boundaries
        self._land_boundaries = land_boundaries
        self._interior_boundaries = interior_boundaries
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
            'vertices': self.vertices,
            'values': self.values,
            'triangles': self.triangles,
            'quads': self.quads,
            'ocean_boundaries': self.ocean_boundaries,
            'land_boundaries': self.land_boundaries,
            'interior_boundaries': self.interior_boundaries,
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
        if len(self.triangles) > 0:
            self.tricontourf(axes=axes, vmin=vmin, vmax=vmax, **kwargs)
        kwargs.pop('levels')
        if len(self.quads) > 0:
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
    def plot_ocean_boundary(
        self,
        id,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        kwargs.update({'axes': axes})
        return self._make_boundary_plot(self.ocean_boundaries[id], **kwargs)

    @_figure
    def plot_land_boundary(
        self,
        id,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        kwargs.update({'axes': axes})
        return self._make_boundary_plot(self.land_boundaries[id], **kwargs)

    @_figure
    def plot_interior_boundary(
        self,
        id,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        kwargs.update({'axes': axes})
        return self._make_boundary_plot(self.interior_boundaries[id], **kwargs)

    @_figure
    def plot_boundaries(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        kwargs.update({'axes': axes})
        kwargs.update({'axes': self.plot_ocean_boundaries(**kwargs)})
        kwargs.update({'axes': self.plot_land_boundaries(**kwargs)})
        kwargs.update({'axes': self.plot_interior_boundaries(**kwargs)})
        return kwargs['axes']

    @_figure
    def plot_ocean_boundaries(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        kwargs.update({'axes': axes})
        for id in self.ocean_boundaries:
            axes = self.plot_ocean_boundary(id, **kwargs)
            kwargs.update({'axes': axes})
        return axes

    @_figure
    def plot_land_boundaries(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        kwargs.update({'axes': axes})
        for id in self.ocean_boundaries:
            axes = self.plot_land_boundary(id, **kwargs)
            kwargs.update({'axes': axes})
        return axes

    @_figure
    def plot_interior_boundaries(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        kwargs.update({'axes': fig.get_axes(axes, figsize)})
        for id in self.interior_boundaries:
            axes = self.plot_interior_boundary(id, **kwargs)
            kwargs.update({'axes': axes})
        if show:
            plt.show()
        return axes

    @staticmethod
    def parse(path):
        return parse_gr3(path)

    @property
    @lru_cache
    def vertices(self):
        return np.array(self._vertices)

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
    def values(self):
        return np.array(self._values)

    @property
    @lru_cache
    def triangles(self):
        return np.array(self._triangles)

    @property
    @lru_cache
    def quads(self):
        return np.array(self._quads)

    @property
    def ocean_boundaries(self):
        return dict(self._ocean_boundaries)

    @property
    def land_boundaries(self):
        return dict(self._land_boundaries)

    @property
    def interior_boundaries(self):
        return dict(self._interior_boundaries)

    @property
    def crs(self):
        return self._crs

    @property
    def description(self):
        return self._description

    @property
    @lru_cache
    def triangulation(self):
        return Triangulation(self.x, self.y, self.triangles)

    @property
    def fgrid(self):
        return self._fgrid

    def _make_boundary_plot(
        self,
        boundary,
        axes=None,
        show=False,
        title=None,
        **kwargs
    ):
        axes = fig.get_axes(axes)
        axes.plot(self.x[boundary], self.y[boundary], **kwargs)
        if show:
            plt.show()
        return axes

    @property
    def _vertices(self):
        return self.__vertices

    @property
    def _values(self):
        return self.__values

    @property
    def _triangles(self):
        return self.__triangles

    @property
    def _quads(self):
        return self.__quads

    @property
    def _ocean_boundaries(self):
        return self.__ocean_boundaries

    @property
    def _land_boundaries(self):
        return self.__land_boundaries

    @property
    def _interior_boundaries(self):
        return self.__interior_boundaries

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

    @_vertices.setter
    def _vertices(self, vertices):
        self.__vertices = vertices

    @_values.setter
    def _values(self, values):
        self.__values = values

    @_triangles.setter
    def _triangles(self, triangles):
        self.__triangles = triangles

    @_quads.setter
    def _quads(self, quads):
        self.__quads = quads

    @_ocean_boundaries.setter
    def _ocean_boundaries(self, ocean_boundaries):
        self.__ocean_boundaries = ocean_boundaries

    @_land_boundaries.setter
    def _land_boundaries(self, land_boundaries):
        self.__land_boundaries = land_boundaries

    @_interior_boundaries.setter
    def _interior_boundaries(self, interior_boundaries):
        self.__interior_boundaries = interior_boundaries

    @_crs.setter
    def _crs(self, crs):
        self.__crs = crs

    @_description.setter
    def _description(self, description):
        self.__description = description

    @_fgrid.setter
    def _fgrid(self, fgrid):
        self.__fgrid = fgrid
