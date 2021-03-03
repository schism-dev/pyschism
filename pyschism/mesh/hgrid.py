from functools import lru_cache
import logging
from typing import Union

import geopandas as gpd
from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from shapely.geometry import LineString

from pyschism.figures import figure, get_topobathy_kwargs
from pyschism.mesh.parsers import grd
from pyschism.mesh.base import Gr3  # , sort_edges, signed_polygon_area

_logger = logging.getLogger(__name__)


class HgridBoundaries:

    def __init__(self, hgrid: "Hgrid", boundaries: Union[dict, None]):
        ocean_boundaries = []
        land_boundaries = []
        interior_boundaries = []
        if boundaries is not None:
            for ibtype, bnds in boundaries.items():
                if ibtype is None:
                    for id, data in bnds.items():
                        indexes = list(map(hgrid.nodes.get_index_by_id,
                                       data['indexes']))
                        ocean_boundaries.append({
                            'id': id,
                            "index_id": data['indexes'],
                            "indexes": indexes,
                            'geometry': LineString(hgrid.vertices[indexes])
                            })

                elif str(ibtype).endswith('1'):
                    for id, data in bnds.items():
                        indexes = list(map(hgrid.nodes.get_index_by_id,
                                       data['indexes']))
                        interior_boundaries.append({
                            'id': id,
                            'ibtype': ibtype,
                            "index_id": data['indexes'],
                            "indexes": indexes,
                            'geometry': LineString(hgrid.vertices[indexes])
                            })
                else:
                    for id, data in bnds.items():
                        _indexes = np.array(data['indexes'])
                        if _indexes.ndim > 1:
                            # ndim > 1 implies we're dealing with an ADCIRC
                            # mesh that includes boundary pairs, such as weir
                            new_indexes = []
                            for i, line in enumerate(_indexes.T):
                                if i % 2 != 0:
                                    new_indexes.extend(np.flip(line))
                                else:
                                    new_indexes.extend(line)
                            _indexes = np.array(new_indexes).flatten()
                        else:
                            _indexes = _indexes.flatten()
                        indexes = list(map(hgrid.nodes.get_index_by_id,
                                       _indexes))

                        land_boundaries.append({
                            'id': id,
                            'ibtype': ibtype,
                            "index_id": data['indexes'],
                            "indexes": indexes,
                            'geometry': LineString(hgrid.vertices[indexes])
                            })

        self._ocean = gpd.GeoDataFrame(ocean_boundaries)
        self._land = gpd.GeoDataFrame(land_boundaries)
        self._interior = gpd.GeoDataFrame(interior_boundaries)
        self._hgrid = hgrid
        self._data = boundaries

    def ocean(self):
        return self._ocean

    def land(self):
        return self._land

    def interior(self):
        return self._interior

    @property
    def data(self):
        return self._data

    @lru_cache(maxsize=1)
    def __call__(self):
        data = []
        for bnd in self.ocean().itertuples():
            data.append({
                'id': bnd.id,
                'ibtype': None,
                "index_id": bnd.index_id,
                "indexes": bnd.indexes,
                'geometry': bnd.geometry})

        for bnd in self.land().itertuples():
            data.append({
                'id': bnd.id,
                'ibtype': bnd.ibtype,
                "index_id": bnd.index_id,
                "indexes": bnd.indexes,
                'geometry': bnd.geometry})

        for bnd in self.interior().itertuples():
            data.append({
                'id': bnd.id,
                'ibtype': bnd.ibtype,
                "index_id": bnd.index_id,
                "indexes": bnd.indexes,
                'geometry': bnd.geometry})

        return gpd.GeoDataFrame(data, crs=self._hgrid.crs)

    def __len__(self):
        return len(self())


class Hgrid(Gr3):
    """
    Class that represents the unstructured planar mesh used by SCHISM.
    """
    # _boundaries = BoundariesDescriptor()

    def __init__(self, *args, boundaries=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._boundaries = HgridBoundaries(self, boundaries)

    @staticmethod
    def open(path, crs=None):
        if str(path).endswith('.ll') and crs is None:
            crs = 'epsg:4326'
        _grd = grd.read(path, crs=crs)
        _grd['nodes'] = {id: (coords, -val) for id, (coords, val)
                         in _grd['nodes'].items()}
        return Hgrid(**_grd)

    def to_dict(self):
        _grd = super().to_dict()
        _grd.update({
            "nodes": {id: (coord, -val) for id, (coord, val)
                      in self.nodes().items()},
            "boundaries": self.boundaries.data})
        return _grd

    @figure
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
        kwargs.update(**get_topobathy_kwargs(self.values, vmin, vmax))
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
        return self._boundaries

    @property
    def ocean_boundaries(self):
        return self.boundaries.ocean

    @property
    def land_boundaries(self):
        return self.boundaries.land

    @property
    def interior_boundaries(self):
        return self.boundaries.interior


#     def auto_generate(
#         self,
#         threshold=0.,
#         land_ibtype=0,
#         interior_ibtype=1,
#     ):
#         values = self._obj.values
#         if np.any(np.isnan(values)):
#             raise Exception("Mesh contains invalid values. Raster values must"
#                             "be interpolated to the mesh before generating "
#                             "boundaries.")

#         # drop current boundaries
#         self._obj.__dict__['boundaries'] = {}
#         self._obj.__dict__['boundaries'][None] = {}

#         # generate exterior boundaries
#         for ring in self._obj.outer_ring_collection.values():
#             # find boundary edges
#             edge_tag = np.full(ring.shape, 0)
#             edge_tag[np.where(values[ring[:, 0]] < threshold)[0], 0] = -1
#             edge_tag[np.where(values[ring[:, 1]] < threshold)[0], 1] = -1
#             edge_tag[np.where(values[ring[:, 0]] >= threshold)[0], 0] = 1
#             edge_tag[np.where(values[ring[:, 1]] >= threshold)[0], 1] = 1
#             # sort boundary edges
#             ocean_boundary = list()
#             land_boundary = list()
#             for i, (e0, e1) in enumerate(edge_tag):
#                 if np.any(np.asarray((e0, e1)) == -1):
#                     ocean_boundary.append(tuple(ring[i, :]))
#                 elif np.any(np.asarray((e0, e1)) == 1):
#                     land_boundary.append(tuple(ring[i, :]))
#             ocean_boundaries = sort_edges(ocean_boundary)
#             land_boundaries = sort_edges(land_boundary)
#             _bnd_id = len(self._obj.boundaries[None])
#             for bnd in ocean_boundaries:
#                 e0, e1 = [list(t) for t in zip(*bnd)]
#                 e0 = list(map(self._obj.get_vertex_id_by_index, e0))
#                 data = e0 + [self._obj.get_vertex_id_by_index(e1[-1])]
#                 self._set_boundary_data(self._obj.vertex_id,
#                                         None, _bnd_id, data)
#                 _bnd_id += 1
#             # add land boundaries
#             if land_ibtype not in self._obj.boundaries:
#                 self._add_boundary_type(land_ibtype)
#             _bnd_id = len(self._obj.boundaries[land_ibtype])
#             for bnd in land_boundaries:
#                 e0, e1 = [list(t) for t in zip(*bnd)]
#                 e0 = list(map(self._obj.get_vertex_id_by_index, e0))
#                 data = e0 + [self._obj.get_vertex_id_by_index(e1[-1])]
#                 self._set_boundary_data(self._obj.vertex_id, land_ibtype,
#                                         _bnd_id, data)
#                 _bnd_id += 1
#         # generate interior boundaries
#         _bnd_id = 0
#         _interior_boundaries = defaultdict()
#         for interiors in self._obj.inner_ring_collection.values():
#             for interior in interiors:
#                 e0, e1 = [list(t) for t in zip(*interior)]
#                 if signed_polygon_area(self._obj._vertices[e0, :]) < 0:
#                     e0 = list(reversed(e0))
#                     e1 = list(reversed(e1))
#                 e0 = list(map(self._obj.get_vertex_id_by_index, e0))
#                 e0.append(e0[0])
#                 _interior_boundaries[_bnd_id] = e0
#                 _bnd_id += 1
#         self._add_boundary_type(interior_ibtype)
#         for bnd_id, data in _interior_boundaries.items():
#             self._set_boundary_data(self._obj._vertex_id, interior_ibtype,
#                                     bnd_id, data)

#     @figure
#     def plot(
#         self,
#         ibtype,
#         id,
#         tags=True,
#         axes=None,
#         show=False,
#         figsize=None,
#         **kwargs
#     ):
#         boundary = list(map(self._obj.get_vertex_index_by_id,
#                             self._obj.boundaries[ibtype][id].indexes))
#         p = axes.plot(self._obj.x[boundary], self._obj.y[boundary], **kwargs)
#         if tags:
#             axes.text(
#                 self._obj.x[boundary[len(boundary)//2]],
#                 self._obj.y[boundary[len(boundary)//2]],
#                 f"ibtype={ibtype}\nid={id}",
#                 color=p[-1].get_color()
#                 )
#         return axes

#     @figure
#     def make_plot(
#         self,
#         axes=None,
#         show=False,
#         figsize=None,
#         **kwargs
#     ):
#         kwargs.update({'axes': axes})
#         for ibtype, bnds in self._obj.boundaries.items():
#             for id in bnds:
#                 axes = self.plot(ibtype, id, **kwargs)
#                 kwargs.update({'axes': axes})
#         return kwargs['axes']

