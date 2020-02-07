from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from pyschism.mesh import grd, sms2dm
import numpy as np
import pathlib
import logging
from collections import defaultdict
import fiona
from functools import lru_cache
from shapely.geometry import LineString, mapping
from pyschism import figures as fig
from pyschism.mesh.base import EuclideanMesh2D


class Hgrid(EuclideanMesh2D):
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
        super().__init__(**grd.euclidean_mesh({
            'nodes': nodes,
            'elements': elements,
            'description': description,
            'crs': crs
            }))
        self._boundaries = boundaries

    @staticmethod
    def open(path, crs=None):
        msh = grd.reader(path)
        msh.update({"crs": crs})
        return Hgrid(**msh)

    def add_boundary_type(self, ibtype):
        if ibtype not in self.boundaries:
            self.__boundaries[ibtype] = defaultdict()
        else:
            msg = f"Cannot add boundary_type={ibtype}: boundary type already "
            msg += "exists."
            raise Exception(msg)

    def set_boundary_data(self, ibtype, id, indexes, **properties):
        msg = "Indexes must be subset of node id's."
        for idx in np.asarray(indexes).flatten():
            assert idx in self.node_index.keys(), msg
        self.__boundaries[ibtype][id] = {
            'indexes': indexes,
            **properties
        }

    def clear_boundaries(self):
        self.__boundaries = {}
        self.__boundaries[None] = {}

    def delete_boundary_type(self, ibtype):
        del self.__boundaries[ibtype]

    def delete_boundary_data(self, ibtype, id):
        del self.__boundaries[ibtype][id]

    def write_boundaries(self, path, overwrite=False):
        path = pathlib.Path(path)
        if path.exists() and not overwrite:
            msg = "Destination path exists and overwrite=False"
            raise IOError(msg)
        with fiona.open(
                    path.absolute(),
                    'w',
                    driver='ESRI Shapefile',
                    crs=self.crs.srs,
                    schema={
                        'geometry': 'LineString',
                        'properties': {
                            'id': 'int',
                            'ibtype': 'str',
                            'bnd_id': 'str'
                            }
                        }) as dst:
            _cnt = 0
            for ibtype, bnds in self.boundaries.items():
                for id, bnd in bnds.items():
                    idxs = list(map(self.get_node_index, bnd['indexes']))
                    linear_ring = LineString(self.xy[idxs].tolist())
                    dst.write({
                            "geometry": mapping(linear_ring),
                            "properties": {
                                "id": _cnt,
                                "ibtype": ibtype,
                                "bnd_id": f"{ibtype}:{id}"
                                }
                            })
                    _cnt += 1

    def generate_boundaries(
        self,
        threshold=0.,
        land_ibtype=0,
        interior_ibtype=1,
    ):
        if np.any(np.isnan(self.values)):
            msg = "Mesh contains invalid values. Raster values must "
            msg += "be interpolated to the mesh before generating "
            msg += "boundaries."
            raise Exception(msg)

        # generate exterior boundaries
        for ring in self.outer_ring_collection.values():
            # find boundary edges
            edge_tag = np.full(ring.shape, 0)
            edge_tag[np.where(self.values[ring[:, 0]] < threshold)[0], 0] = -1
            edge_tag[np.where(self.values[ring[:, 1]] < threshold)[0], 1] = -1
            edge_tag[np.where(self.values[ring[:, 0]] >= threshold)[0], 0] = 1
            edge_tag[np.where(self.values[ring[:, 1]] >= threshold)[0], 1] = 1
            # sort boundary edges
            ocean_boundary = list()
            land_boundary = list()
            for i, (e0, e1) in enumerate(edge_tag):
                if np.any(np.asarray((e0, e1)) == -1):
                    ocean_boundary.append(tuple(ring[i, :]))
                elif np.any(np.asarray((e0, e1)) == 1):
                    land_boundary.append(tuple(ring[i, :]))
            ocean_boundaries = self.sort_edges(ocean_boundary)
            land_boundaries = self.sort_edges(land_boundary)
            _bnd_id = len(self.boundaries[None])
            for bnd in ocean_boundaries:
                e0, e1 = [list(t) for t in zip(*bnd)]
                e0 = list(map(self.get_node_id, e0))
                data = e0 + [self.get_node_id(e1[-1])]
                self.set_boundary_data(None, _bnd_id, data)
                _bnd_id += 1
            # add land boundaries
            if land_ibtype not in self.boundaries:
                self.add_boundary_type(land_ibtype)
            _bnd_id = len(self._boundaries[land_ibtype])
            for bnd in land_boundaries:
                e0, e1 = [list(t) for t in zip(*bnd)]
                e0 = list(map(self.get_node_id, e0))
                data = e0 + [self.get_node_id(e1[-1])]
                self.set_boundary_data(land_ibtype, _bnd_id, data)
                _bnd_id += 1
        # generate interior boundaries
        _bnd_id = 0
        _interior_boundaries = defaultdict()
        for interiors in self.inner_ring_collection.values():
            for interior in interiors:
                e0, e1 = [list(t) for t in zip(*interior)]
                if self.signed_polygon_area(self.coords[e0, :]) < 0:
                    e0 = list(reversed(e0))
                    e1 = list(reversed(e1))
                e0 = list(map(self.get_node_id, e0))
                e0 += [e0[0]]
                _interior_boundaries[_bnd_id] = e0
                _bnd_id += 1
        self.add_boundary_type(interior_ibtype)
        for bnd_id, data in _interior_boundaries.items():
            self.set_boundary_data(interior_ibtype, bnd_id, data)

    @fig._figure
    def make_plot(
        self,
        axes=None,
        vmin=None,
        vmax=None,
        show=False,
        title=None,
        # figsize=rcParams["figure.figsize"],
        extent=None,
        cbar_label=None,
        **kwargs
    ):
        if vmin is None:
            vmin = np.min(self.values)
        if vmax is None:
            vmax = np.max(self.values)
        kwargs.update(**fig.get_topobathy_kwargs(self.values, vmin, vmax))
        kwargs.pop('col_val')
        levels = kwargs.pop('levels')
        if vmin != vmax:
            self.tricontourf(
                axes=axes,
                levels=levels,
                vmin=vmin,
                vmax=vmax,
                **kwargs
            )
        else:
            self.tripcolor(axes=axes, **kwargs)
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
            orientation='horizontal'
        )
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
        tags=True,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):

        boundary = list(map(
            self.get_node_index, self.boundaries[ibtype][id]['indexes']))
        p = axes.plot(self.x[boundary], self.y[boundary], **kwargs)
        if tags:
            axes.text(
                self.x[boundary[len(boundary)//2]],
                self.y[boundary[len(boundary)//2]],
                f"ibtype={ibtype}\nid={id}",
                color=p[-1].get_color()
                )
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
            for id in bnds:
                axes = self.plot_boundary(ibtype, id, **kwargs)
                kwargs.update({'axes': axes})
        return kwargs['axes']

    @property
    def boundaries(self):
        return self._boundaries

    @property
    @lru_cache
    def logger(self):
        return logging.getLogger(__name__ + '.' + self.__class__.__name__)

    @property
    def _grd(self):
        grd = super()._grd
        grd.update({'boundaries': self.boundaries})
        return grd

    @property
    def _sms2dm(self):
        _sms2dm = super()._sms2dm
        _sms2dm.update({'nodestrings': sms2dm.nodestrings(self.boundaries)})
        return _sms2dm

    @property
    def _boundaries(self):
        return self.__boundaries

    @_boundaries.setter
    def _boundaries(self, boundaries):
        self.clear_boundaries() # clear
        if boundaries is not None:
            for ibtype, bnds in boundaries.items():
                if ibtype is not None:
                    self.add_boundary_type(ibtype)
                for id, bnd in bnds.items():
                    if 'properties' in bnd.keys():
                        properties = bnd['properties']
                    else:
                        properties = {}
                    self.set_boundary_data(
                        ibtype,
                        id,
                        bnd['indexes'],
                        **properties
                    )
