from typing import Union

from shapely.geometry import Polygon, MultiPolygon, Point

from pyschism.mesh.hgrid import Hgrid

#from pyschism.mesh.base import Gr3

class PropField:

#    def __init__(self, hgrid):
#        self.hgrid = Hgrid 
    def __init__():
        pass

    @classmethod
    def constant(self, hgrid, value):

        hgrid = hgrid.to_dict()
        elements = hgrid['elements']
        out = []
        for iele, element in elements.items():
            line = [f'{iele}']
            line.extend([f'{value}'])
            line.extend([f'\n'])
            out.append(' '.join(line))

        return out

    def define_by_region( 
        hgrid, 
        region: Union[Polygon, MultiPolygon],
        value):

        hgrid = hgrid.to_dict()
        #Get lon/lat of nodes
        nodes = hgrid['nodes']
        lon = []
        lat = []
        for id, (coords, values) in nodes.items():
            lon.append(coords[0])
            lat.append(coords[1])
        
        #Get centroid of elements and check if it is in the region
        elements = hgrid['elements']
        out = []
        for id, element in elements.items():
            i34 = len(element)
            if i34 == 3:
                v1 = int(element[0])
                v2 = int(element[1])
                v3 = int(element[2])
                xtmp = (lon[v1-1] + lon[v2-1] + lon[v3-1])/3
                ytmp = (lat[v1-1] + lat[v2-1] + lat[v3-1])/3
            else:
                v1 = int(element[0])
                v2 = int(element[1])
                v3 = int(element[2])
                v4 = int(element[3])
                xtmp = (lon[v1-1] + lon[v2-1] + lon[v3-1] + \
                    lon[v4-1])/4
                ytmp = (lat[v1-1] + lat[v2-1] + lat[v3-1] + \
                    lat[v4-1])/4
            p = Point((xtmp, ytmp))
            if p.within(region):
               value = 0 
            line = [f'{id}'] 
            line.extend([f'{value}'])
            line.extend([f'\n'])
            out.append(' '.join(line))
        return out
 
class Fluxflag(PropField):
    """ Class for writing fluxflag.prop file, which is parameter for
        checking volume and salt conservation."""

    pass

class Tvdflag(PropField):
    """Class for writing tvd.prop file, which specify horizontal regions 
       where upwind or TVD/TVD^2 is used based on the element property values
       (0: upwind; 1: TVD/TVD^2). """

    pass
