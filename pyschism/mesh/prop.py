import pathlib
from typing import Union

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, MultiPolygon

from pyschism.mesh.base import Gr3


class PropField:

    def __init__(self, gr3: Gr3, element_values: np.array):

        if not isinstance(gr3, Gr3):
            raise TypeError(
                f'Argument gr3 must be an instance of type {Gr3}, not type '
                f'{type(gr3)}.')

        values = np.array(element_values).flatten().shape[0]
        if gr3.elements.shape[0] != values:
            raise ValueError(
                'Shape mismatch between element_values and hgrid.')

        self.elements = gr3.elements
        self.values = values

    def __str__(self):
        f = []
        for i, (iele, element) in enumerate(self.elements.items()):
            f.append(f'{iele:d} {self.values[i]:G}')
        return '\n'.join(f)

    @classmethod
    def constant(cls, gr3: Gr3, value: np.array):
        return cls(gr3, np.full((gr3.elements.shape[0],), value))

    def write(self, path, overwrite: bool = False):
        path = pathlib.Path(path)
        if path.exists() and not overwrite:
            raise IOError('path exists and overwrite is False')
        with open(path, 'w') as f:
            f.write(str(self))

    @classmethod
    def by_region(
            cls,
            gr3: Gr3,
            region: Union[Polygon, MultiPolygon],
            inner_value: float,
            outer_value: float,
            op='touches',
    ):

        if not isinstance(gr3, Gr3):
            raise TypeError(
                f'Argument gr3 must be an instance of type {Gr3}, not type '
                f'{type(gr3)}.')

        if not isinstance(region, (Polygon, MultiPolygon)):
            raise TypeError(
                f'Argument region must be an instance of types {Polygon} or '
                f'{MultiPolygon}, not type {type(region)}.')

        if isinstance(region, Polygon):
            region = [region]

        obj = cls.constant(gr3, np.nan)
        # gdf_in = gpd.sjoin(
        #     gr3.elements.gdf,
        #     gpd.GeoDataFrame(
        #         {'geometry': region},
        #         crs=gr3.crs
        #     ),
        #     op=op
        # )
        gr3_gdf = gr3.elements.gdf
        gdf_in = gr3_gdf.geometry.touches(region)
        inner_idxs = [i.index for i in gdf_in.itertuples()]
        obj.values[inner_idxs] = inner_value
        outer_idxs = gr3_gdf.loc[gr3_gdf.index.difference(gdf_in)]
        obj.values[outer_idxs] = outer_value
        return obj


class Fluxflag(PropField):
    """ Class for writing fluxflag.prop file, which is parameter for
        checking volume and salt conservation.
    """
    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, -1)

    @classmethod
    def by_region(cls, region: Union[Polygon, MultiPolygon],
                  value: int):
        if value not in [1, -1]:
            raise ValueError('Argument value must be 1 or -1.')
        return super(cls).by_region(
            region, inner_value=value, outer_value=-value)


class Tvdflag(PropField):
    """Class for writing tvd.prop file, which specify horizontal regions 
       where upwind or TVD/TVD^2 is used based on the element property values
       (0: upwind; 1: TVD/TVD^2).
    """
    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 1)

    @classmethod
    def by_region(cls, region: Union[Polygon, MultiPolygon],
                  value: int):
        if value not in [1, -1]:
            raise ValueError('Argument value must be 1 or -1.')
        return super(cls).by_region(
            region, inner_value=value, outer_value=-value)







        # hgrid = hgrid.to_dict()
        # Get lon/lat of nodes
        # nodes = hgrid['nodes']
        # lon = []
        # lat = []
        # for id, (coords, values) in nodes.items():
        #     lon.append(coords[0])
        #     lat.append(coords[1])

        # Get centroid of elements and check if it is in the region
        # elements = hgrid['elements']
        # out = []
        # for id, element in elements.items():
        #     i34 = len(element)
        #     if i34 == 3:
        #         v1 = int(element[0])
        #         v2 = int(element[1])
        #         v3 = int(element[2])
        #         xtmp = (lon[v1-1] + lon[v2-1] + lon[v3-1])/3
        #         ytmp = (lat[v1-1] + lat[v2-1] + lat[v3-1])/3
        #     else:
        #         v1 = int(element[0])
        #         v2 = int(element[1])
        #         v3 = int(element[2])
        #         v4 = int(element[3])
        #         xtmp = (lon[v1-1] + lon[v2-1] + lon[v3-1] + lon[v4-1])/4
        #         ytmp = (lat[v1-1] + lat[v2-1] + lat[v3-1] + lat[v4-1])/4
        #     p = Point((xtmp, ytmp))
        #     if p.within(region):
        #         value = 0
        #     line = [f'{id}']
        #     line.extend([f'{value}'])
        #     line.extend(['\n'])
        #     out.append(' '.join(line))
        # self.out = out
        # return self.out



        # hgrid = hgrid.to_dict()
        # Get lon/lat of nodes
        # nodes = hgrid['nodes']
        # lon = []
        # lat = []
        # for id, (coords, values) in nodes.items():
        #     lon.append(coords[0])
        #     lat.append(coords[1])

        # # Get centroid of elements and check if it is in the region
        # elements = hgrid['elements']
        # out = []
        # for id, element in elements.items():
        #     i34 = len(element)
        #     if i34 == 3:
        #         v1 = int(element[0])
        #         v2 = int(element[1])
        #         v3 = int(element[2])
        #         xtmp = (lon[v1-1] + lon[v2-1] + lon[v3-1])/3
        #         ytmp = (lat[v1-1] + lat[v2-1] + lat[v3-1])/3
        #     else:
        #         v1 = int(element[0])
        #         v2 = int(element[1])
        #         v3 = int(element[2])
        #         v4 = int(element[3])
        #         xtmp = (lon[v1-1] + lon[v2-1] + lon[v3-1] + lon[v4-1])/4
        #         ytmp = (lat[v1-1] + lat[v2-1] + lat[v3-1] + lat[v4-1])/4
        #     p = Point((xtmp, ytmp))
        #     if p.within(region):
        #         value = 0
        #     line = [f'{id}']
        #     line.extend([f'{value}'])
        #     line.extend(['\n'])
        #     out.append(' '.join(line))
        # self.out = out
        # return self.out