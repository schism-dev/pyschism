import pathlib
from typing import List, Tuple, Union

from matplotlib.collections import PolyCollection
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, LinearRing

from pyschism.mesh.base import Gr3
from pyschism.figures import figure


class Prop:

    def __init__(self, gr3: Gr3, element_values: np.array):

        if not isinstance(gr3, Gr3):
            raise TypeError(
                f'Argument gr3 must be an instance of type {Gr3}, not type '
                f'{type(gr3)}.')

        values = np.array(element_values).flatten()
        if len(gr3.elements) != values.shape[0]:
            raise ValueError(
                'Shape mismatch between element_values and hgrid.')

        self.gr3 = gr3
        self.values = values

    def __str__(self):
        f = []
        for i, (iele, element) in enumerate(self.gr3.elements.to_dict().items()):
            f.append(f'{iele} {self.values[i]:G}')
        return '\n'.join(f)

    @classmethod
    def constant(cls, gr3: Gr3, value: np.array):
        return cls(gr3, np.full((len(gr3.elements),), value))

    def write(self, path, overwrite: bool = False):
        path = pathlib.Path(path)
        if path.exists() and not overwrite:
            raise IOError('path exists and overwrite is False')
        with open(path, 'w') as f:
            f.write(str(self))

    @figure
    def make_plot(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):

        tria = PolyCollection(self.gr3.coords[self.gr3.elements.triangles])
        quad = PolyCollection(self.gr3.coords[self.gr3.elements.quads])
        tria.set_array(self.values[self.gr3.elements.tri_idxs])
        quad.set_array(self.values[self.gr3.elements.qua_idxs])
        axes.add_collection(tria)
        axes.add_collection(quad)
        return axes


class Fluxflag(Prop):
    """ Class for writing fluxflag.prop file, which is parameter for
        checking volume and salt conservation.
    """
    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, -1)

    def add_region(self, upstream: LinearRing, downstream: LinearRing):

        def get_next_value_pair():
            next_value = np.max(self.values) + 1
            return next_value + 1, next_value
        upstream_val, downstream_val = get_next_value_pair()

        upstream_node_mask = self.gr3.nodes.gdf.within(Polygon(upstream)).to_numpy()
        upstream_elem_idxs = np.where(np.any(upstream_node_mask[self.gr3.elements.array], axis=1))
        self.values[upstream_elem_idxs] = upstream_val

        downstream_node_mask = self.gr3.nodes.gdf.within(Polygon(downstream)).to_numpy()
        downstream_elem_idxs = np.where(np.any(downstream_node_mask[self.gr3.elements.array], axis=1))
        self.values[downstream_elem_idxs] = downstream_val

    @classmethod
    def from_prop_table(cls, hgrid, prop_table):
        obj = cls.default(hgrid)
        for regions in prop_table_to_list_of_tuples(prop_table):
            obj.add_region(*regions)
        return obj


class Tvdflag(Prop):
    """Class for writing tvd.prop file, which specify horizontal regions 
       where upwind or TVD/TVD^2 is used based on the element property values
       (0: upwind; 1: TVD/TVD^2).
    """

    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 1)

    @classmethod
    def from_geometry(
            cls,
            gr3: Gr3,
            region: Union[Polygon, MultiPolygon],
            inner_value: float = 0,
            outer_value: float = 1,
    ):
        '''inner_value == 0 implies the use of updwind, TVD is the default
        everywhere'''

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
        elements_gdf = gr3.elements.gdf
        gdf_in = elements_gdf.geometry.intersects(region)
        obj.values[gdf_in] = inner_value
        outer_indexes = np.setdiff1d(elements_gdf.index, np.where(gdf_in))
        obj.values[outer_indexes] = outer_value
        return obj


def reg2multipoly(file):
    with open(file) as f:
        f.readline()
        npoly = int(f.readline())
        polygons = []
        for _ in range(npoly):
            line = f.readline().split()
            if len(line) == 0:
                break
            nvrt = int(line[0])
            tflag = int(line[1])
            exterior = []
            for _ in range(nvrt):
                exterior.append(tuple(list(map(float, f.readline().split()))))
            interiors = []
            line = f.readline().split()
            if len(line) != 0:
                while tflag == 1:
                    nvrt = int(line[0])
                    tflag = int(line[1])
                    interior = []
                    for _ in range(nvrt):
                        interior.append(tuple(list(map(float, f.readline().split()))))
                    interiors.append(interior)
            polygons.append(Polygon(exterior, interiors))
    return MultiPolygon(polygons)


def prop_table_to_list_of_tuples(file) -> List[Tuple[Polygon, Polygon]]:
    output = []
    file = pathlib.Path(file)
    with open(file) as f:
        lines = f.read()
    for line in lines.split('\n'):
        if ';' not in line:
            break
        reg_upstream, reg_downstream = line.split(';')
        reg_downstream_file = reg_downstream.split(':')[0]
        reg_upstream_file = reg_upstream.split(':')[0]
        reg_upstream_file = reg_upstream_file.strip()
        reg_downstream_file = reg_downstream_file.strip()
        downstream = reg2multipoly(file.parent / reg_downstream_file).geoms[0].exterior
        upstream = reg2multipoly(file.parent / reg_upstream_file).geoms[0].exterior
        output.append((upstream, downstream))
    return output
