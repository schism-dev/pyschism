from typing import List, Union

import numpy as np
import geopandas as gpd
from shapely.geometry import LineString

# from pyschism.forcing.bctides.mod3d import TEM_3D, SAL_3D
# from pyschism.forcing.bctides.nudge import TEM_Nudge, SAL_Nudge
# from pyschism.forcing.bctides.elev2d import Elev2D
# from pyschism.forcing.bctides.uv3d import UV3D
# from pyschism.forcing.bctides.iettype import Iettype
# from pyschism.forcing.bctides.ifltype import Ifltype
# from pyschism.forcing.bctides.itetype import Itetype
# from pyschism.forcing.bctides.isatype import Isatype
# from pyschism.forcing.bctides.itrtype import Itrtype


class Boundaries:
    def __init__(self, hgrid, boundaries: Union[dict, None]):

        ocean_boundaries = []
        land_boundaries = []
        interior_boundaries = []
        if boundaries is not None:
            for ibtype, bnds in boundaries.items():
                if ibtype is None:
                    for id, data in bnds.items():
                        indexes = list(
                            map(hgrid.nodes.get_index_by_id, data["indexes"])
                        )
                        ocean_boundaries.append(
                            {
                                "id": str(id + 1),  # hacking it
                                "index_id": data["indexes"],
                                "indexes": indexes,
                                "geometry": LineString(hgrid.vertices[indexes]),
                                # "iettype": None,
                                # "ifltype": None,
                                # "itetype": None,
                                # "isatype": None,
                                # "itrtype": {},
                                # 'nudge_temperature': None,
                                # 'nudge_salinity': None,
                            }
                        )

                elif str(ibtype).endswith("1"):
                    for id, data in bnds.items():
                        indexes = list(
                            map(hgrid.nodes.get_index_by_id, data["indexes"])
                        )
                        interior_boundaries.append(
                            {
                                "id": str(id + 1),
                                "ibtype": ibtype,
                                "index_id": data["indexes"],
                                "indexes": indexes,
                                "geometry": LineString(hgrid.vertices[indexes]),
                            }
                        )
                else:
                    for id, data in bnds.items():
                        _indexes = np.array(data["indexes"])
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
                        indexes = list(map(hgrid.nodes.get_index_by_id, _indexes))

                        land_boundaries.append(
                            {
                                "id": str(id + 1),
                                "ibtype": ibtype,
                                "index_id": data["indexes"],
                                "indexes": indexes,
                                "geometry": LineString(hgrid.vertices[indexes]),
                            }
                        )

        self.open = gpd.GeoDataFrame(
            ocean_boundaries, crs=hgrid.crs if len(ocean_boundaries) > 0 else None
        )
        self.land = gpd.GeoDataFrame(
            land_boundaries, crs=hgrid.crs if len(land_boundaries) > 0 else None
        )
        self.interior = gpd.GeoDataFrame(
            interior_boundaries, crs=hgrid.crs if len(interior_boundaries) > 0 else None
        )
        self.hgrid = hgrid
        self.data = boundaries

    # def elev2d(self):
    #     return Elev2D(self.hgrid)

    # def uv3d(self, vgrid):
    #     return UV3D(self.hgrid, vgrid)

    # def tem3d(self, vgrid):
    #     return TEM_3D(self.hgrid, vgrid)

    # def sal3d(self, vgrid):
    #     return SAL_3D(self.hgrid, vgrid)

    # def TEM_nudge(self, vgrid, data_source, rlmax=1.5, rnu_day=0.25):
    #     return TEM_Nudge(self.hgrid, vgrid, data_source, rlmax, rnu_day)

    # def SAL_nudge(self, vgrid, data_source, rlmax=1.5, rnu_day=0.25):
    #     return SAL_Nudge(self.hgrid, vgrid, data_source, rlmax, rnu_day)

    # def set_forcing(
    #         self,
    #         boundary_id: int,
    #         iettype: Iettype = None,
    #         ifltype: Ifltype = None,
    #         itetype: Itetype = None,
    #         isatype: Isatype = None,
    #         itrtype: Union[Itrtype, List[Itrtype]] = None,
    # ):

    #     def check_input_type(argname, obj, cls):
    #         if not isinstance(obj, cls):
    #             raise TypeError(
    #                 f'Argument {argname} must be of type {cls}, not '
    #                 f'type {type(obj)}.')

    #     def modify_dataframe(argname, obj, cls):
    #         if obj is not None:
    #             check_input_type(argname, obj, cls)
    #             idxs = self.open[self.open['id'] == boundary_id].index.values
    #             for idx in idxs:
    #                 self.open.at[idx, argname] = obj

    #     modify_dataframe('iettype', iettype, Iettype)
    #     modify_dataframe('ifltype', ifltype, Ifltype)
    #     modify_dataframe('itetype', itetype, Itetype)
    #     modify_dataframe('isatype', isatype, Isatype)
    #     modify_dataframe('itrtype', itrtype, Itrtype)


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
