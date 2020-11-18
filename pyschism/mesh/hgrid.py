from collections import defaultdict, namedtuple
import pathlib

import fiona  # type: ignore[import]
from matplotlib.cm import ScalarMappable  # type: ignore[import]
import matplotlib.pyplot as plt  # type: ignore[import]
from mpl_toolkits.axes_grid1 import make_axes_locatable  # type: ignore[import]
import numpy as np  # type: ignore[import]
from shapely.geometry import LineString, mapping  # type: ignore[import]

from pyschism import figures as fig
from pyschism.mesh.parsers import grd
from pyschism.mesh.base import Gr3, sort_edges, signed_polygon_area


class Boundaries:

    def __set__(self, obj, boundaries):
        obj.__dict__['boundaries'] = {}
        obj.__dict__['boundaries'][None] = {}
        self._obj = obj
        if boundaries is not None:
            for ibtype, bnds in boundaries.items():
                if ibtype is not None:
                    self._add_boundary_type(ibtype)
                for id, bnd in bnds.items():
                    if 'properties' in bnd.keys():
                        properties = bnd['properties']
                    else:
                        properties = {}
                    self._set_boundary_data(
                        obj.vertex_id,
                        ibtype,
                        id,
                        bnd['indexes'],
                        **properties
                    )

    def __get__(self, obj, val):
        return obj.__dict__['boundaries']

    def __iter__(self):
        for ibtype, bnd in self._obj.boundaries:
            yield ibtype, bnd

    def _add_boundary_type(self, ibtype):
        if ibtype not in self._obj.boundaries:
            self._obj.boundaries[ibtype] = defaultdict()
        else:
            raise ValueError(f'Boundary with ibtype {ibtype} already exists.')

    def _set_boundary_data(self, vertex_id, ibtype, id, indexes, **properties):
        if not set(indexes).issubset(set(vertex_id)):
            raise ValueError("Indexes must be subset of node id's.")
        self._obj.boundaries[ibtype][id] = namedtuple(
            'boundary', ['indexes', 'properties'])(
            indexes=indexes, properties=properties)

    def to_shapefile(self, path, overwrite=False):
        path = pathlib.Path(path)
        if path.exists() and not overwrite:
            raise IOError("Destination path exists and overwrite=False")
        with fiona.open(path, 'w', driver='ESRI Shapefile',
                        crs=self._obj.crs.srs,
                        schema={
                            'geometry': 'LineString',
                            'properties': {
                                'id': 'int',
                                'ibtype': 'str',
                                'bnd_id': 'str'}}) as dst:
            _cnt = 0
            for ibtype, bnds in self._obj.boundaries:
                for id, bnd in bnds.items():
                    idxs = list(map(self._obj.get_vertex_by_id, bnd.indexes))
                    linear_ring = LineString(
                        self._obj.vertices[idxs].tolist())
                    dst.write({
                            "geometry": mapping(linear_ring),
                            "properties": {
                                "id": _cnt,
                                "ibtype": ibtype,
                                "bnd_id": f"{ibtype}:{id}"}})
                    _cnt += 1

    def auto_generate(
        self,
        threshold=0.,
        land_ibtype=0,
        interior_ibtype=1,
    ):
        values = self._obj.values
        if np.any(np.isnan(values)):
            raise Exception("Mesh contains invalid values. Raster values must"
                            "be interpolated to the mesh before generating "
                            "boundaries.")

        # drop current boundaries
        self._obj.__dict__['boundaries'] = {}
        self._obj.__dict__['boundaries'][None] = {}

        # generate exterior boundaries
        for ring in self._obj.outer_ring_collection.values():
            # find boundary edges
            edge_tag = np.full(ring.shape, 0)
            edge_tag[np.where(values[ring[:, 0]] < threshold)[0], 0] = -1
            edge_tag[np.where(values[ring[:, 1]] < threshold)[0], 1] = -1
            edge_tag[np.where(values[ring[:, 0]] >= threshold)[0], 0] = 1
            edge_tag[np.where(values[ring[:, 1]] >= threshold)[0], 1] = 1
            # sort boundary edges
            ocean_boundary = list()
            land_boundary = list()
            for i, (e0, e1) in enumerate(edge_tag):
                if np.any(np.asarray((e0, e1)) == -1):
                    ocean_boundary.append(tuple(ring[i, :]))
                elif np.any(np.asarray((e0, e1)) == 1):
                    land_boundary.append(tuple(ring[i, :]))
            ocean_boundaries = sort_edges(ocean_boundary)
            land_boundaries = sort_edges(land_boundary)
            _bnd_id = len(self._obj.boundaries[None])
            for bnd in ocean_boundaries:
                e0, e1 = [list(t) for t in zip(*bnd)]
                e0 = list(map(self._obj.get_vertex_id_by_index, e0))
                data = e0 + [self._obj.get_vertex_id_by_index(e1[-1])]
                self._set_boundary_data(self._obj.vertex_id,
                                        None, _bnd_id, data)
                _bnd_id += 1
            # add land boundaries
            if land_ibtype not in self._obj.boundaries:
                self._add_boundary_type(land_ibtype)
            _bnd_id = len(self._obj.boundaries[land_ibtype])
            for bnd in land_boundaries:
                e0, e1 = [list(t) for t in zip(*bnd)]
                e0 = list(map(self._obj.get_vertex_id_by_index, e0))
                data = e0 + [self._obj.get_vertex_id_by_index(e1[-1])]
                self._set_boundary_data(self._obj.vertex_id, land_ibtype,
                                        _bnd_id, data)
                _bnd_id += 1
        # generate interior boundaries
        _bnd_id = 0
        _interior_boundaries = defaultdict()
        for interiors in self._obj.inner_ring_collection.values():
            for interior in interiors:
                e0, e1 = [list(t) for t in zip(*interior)]
                if signed_polygon_area(self._obj._vertices[e0, :]) < 0:
                    e0 = list(reversed(e0))
                    e1 = list(reversed(e1))
                e0 = list(map(self._obj.get_vertex_id_by_index, e0))
                e0.append(e0[0])
                _interior_boundaries[_bnd_id] = e0
                _bnd_id += 1
        self._add_boundary_type(interior_ibtype)
        for bnd_id, data in _interior_boundaries.items():
            self._set_boundary_data(self._obj._vertex_id, interior_ibtype,
                                    bnd_id, data)

    @fig._figure
    def plot(
        self,
        ibtype,
        id,
        tags=True,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        boundary = list(map(self._obj.get_vertex_index_by_id,
                            self._obj.boundaries[ibtype][id].indexes))
        p = axes.plot(self._obj.x[boundary], self._obj.y[boundary], **kwargs)
        if tags:
            axes.text(
                self._obj.x[boundary[len(boundary)//2]],
                self._obj.y[boundary[len(boundary)//2]],
                f"ibtype={ibtype}\nid={id}",
                color=p[-1].get_color()
                )
        return axes

    @fig._figure
    def make_plot(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        kwargs.update({'axes': axes})
        for ibtype, bnds in self._obj.boundaries.items():
            for id in bnds:
                axes = self.plot(ibtype, id, **kwargs)
                kwargs.update({'axes': axes})
        return kwargs['axes']


class Hgrid(Gr3):
    """
    Class that represents the unstructured planar mesh used by SCHISM.
    """
    _boundaries = Boundaries()

    def __init__(self, *args, boundaries=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._boundaries = boundaries

    @staticmethod
    def open(path, crs=None):
        _grd = grd.read(path, crs=crs)
        _grd['values'] = [-val for val in _grd['values']]
        return Hgrid(**_grd)

    def to_dict(self):
        return {
            "description": self.description,
            "vertices": self.vertices,
            "elements": self.elements,
            "vertex_id": self.vertex_id,
            "element_id": self.element_id,
            "values": [-val for val in self.values],
            "boundaries": self.boundaries}

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

    @property
    def boundaries(self):
        return self.__dict__['boundaries']

    @property
    def open_boundaries(self):
        return self.__dict__['boundaries'][None]
