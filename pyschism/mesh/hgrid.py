import numpy as np
from collections import defaultdict
from collections.abc import Iterable, Mapping
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from functools import lru_cache
from pyschism.mesh import gr3
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pyschism.mesh.gmesh import Gmesh
from pyschism.mesh.friction import ManningsN, DragCoefficient, RoughnessLength
import pyschism.mesh.figures as fig


class Hgrid(Gmesh):
    """
    Class that represents the unstructured planar mesh used by SCHISM.
    """

    def __init__(
        self,
        nodes,
        elements,
        boundaries=None,
        crs=None,
        description=None,
    ):
        coords = {id: (x, y) for id, ((x, y), value) in nodes.items()}
        values = [-value for coord, value in nodes.values()]
        triangles = {id: geom for id, geom in elements.items()
                     if len(geom) == 3}
        quads = {id: geom for id, geom in elements.items()
                 if len(geom) == 4}
        super().__init__(coords, triangles, quads, values, crs, description)
        self._boundaries = boundaries

    @staticmethod
    def open(hgrid, crs=None):
        kwargs = gr3.reader(hgrid)
        kwargs.update({"crs": crs})
        return Hgrid(**kwargs)

    def set_friction(self, value, ftype='manning'):

        # certify ftype
        ftypes = {
            'manning': ManningsN,
            'drag': DragCoefficient,
            'rough': RoughnessLength
        }
        msg = f"ftype argument must be one of {ftypes.keys()}"
        assert ftype.lower() in ftypes, msg

        # certify value
        msg = "value argument must be an instance of type "
        msg += f"{int}, {float} or an iterable ."
        assert isinstance(value, (Iterable, int, float)), msg

        if isinstance(value, (int, float)):
            if ftype == 'manning':
                self._fgrid = ftypes[ftype].constant(self, value)

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

    def write(self, path, overwrite=False):
        grd = {
            'description': self.description,
            'nodes': self.nodes,
            'elements': self.elements,
            'boundaries': self.boundaries,
        }
        gr3.writer(grd, path, overwrite)

    @fig._figure
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

    @fig._figure
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

    @fig._figure
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

    @property
    @lru_cache
    def nodes(self):
        return {id: ((x, y), self.values[i]) for i, (id, (x, y))
                in enumerate(self._coords.items())}

    @property
    @lru_cache
    def elements(self):
        keys = [id for id in self._triangles]
        keys.extend([id for id in self._quads])
        keys.sort(key=int)
        geom = dict(self._triangles.items())
        geom.update(dict(self._quads.items()))
        elements = dict()
        for i, id in enumerate(keys):
            elements[id] = geom[id]
        return elements

    @property
    def boundaries(self):
        return self._boundaries.copy()

    @property
    def fgrid(self):
        return self._fgrid

    @property
    def _boundaries(self):
        return self.__boundaries

    @property
    def _fgrid(self):
        return self.__fgrid

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

    @_fgrid.setter
    def _fgrid(self, fgrid):
        self.__fgrid = fgrid

