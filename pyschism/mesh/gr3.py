import pathlib
import numpy as np
from pyschism.mesh.gmesh import Gmesh
from functools import lru_cache
from collections import defaultdict


def writer(grd, path, overwrite=False):
    path = pathlib.Path(path)
    if path.is_file() and not overwrite:
        msg = 'File exists, pass overwrite=True to allow overwrite.'
        raise Exception(msg)
    with open(path, 'w') as f:
        write_gr3(f, grd)
    return 0  # for unittests


def reader(path):
    grd = dict()
    grd['nodes'] = defaultdict(list)
    grd['elements'] = defaultdict(list)
    with open(pathlib.Path(path), 'r') as f:
        grd['description'] = f.readline().strip('\n')
        NE, NP = map(int, f.readline().split())
        for i in range(NP):
            id, x, y, z = f.readline().split()
            grd['nodes'][id] = ((float(x), float(y)), float(z))
        for i in range(NE):
            geom = f.readline().split()
            grd['elements'][geom[0]] = [x for x in geom[2:]]
        # Assume EOF if NOPE is empty.
        try:
            NOPE = int(f.readline().split()[0])
        except IndexError:
            return grd
        # let NOPE=-1 mean an ellipsoidal-mesh
        # reassigning NOPE to 0 until further implementation is applied.
        grd['boundaries'] = defaultdict(dict)
        _bnd_id = 0
        f.readline()  # Not used.
        while _bnd_id < NOPE:
            NETA = int(f.readline().split()[0])
            _cnt = 0
            grd['boundaries'][None][_bnd_id] = list()
            while _cnt < NETA:
                grd['boundaries'][None][_bnd_id].append(
                    f.readline().split()[0].strip())
                _cnt += 1
            _bnd_id += 1
        NBOU = int(f.readline().split()[0])
        _nbnd_cnt = 0
        f.readline()  # not used
        while _nbnd_cnt < NBOU:
            npts, ibtype = map(int, f.readline().split()[:2])
            _pnt_cnt = 0
            if ibtype not in grd['boundaries']:
                _bnd_id = 0
            else:
                _bnd_id = len(grd['boundaries'][ibtype])
            grd['boundaries'][ibtype][_bnd_id] = list()
            while _pnt_cnt < npts:
                line = f.readline().split()
                grd['boundaries'][ibtype][_bnd_id].append(line[0])
                _pnt_cnt += 1
            _nbnd_cnt += 1
    # typecast defaultdict to regular dict
    for key, value in grd.items():
        if key != 'description':
            grd[key] = dict(value)
    return grd


def write_gr3(f, grd):
    write_graph(f, grd)
    if 'boundaries' in grd.keys():
        write_boundaries(f, grd)


def write_graph(f, grd):
    f.write(f"{grd['description']}\n")
    f.write(f"{len(grd['elements'])} ")
    f.write(f"{len(grd['nodes'])}\n")
    # TODO: Make faster using np.array2string
    for id, ((x, y), z) in grd['nodes'].items():
        f.write(f"{id} ")
        f.write(f"{x:<.16E} ")
        f.write(f"{y:<.16E} ")
        f.write(f"{z:<.16E}\n")
    for id, geom in grd['elements'].items():
        f.write(f"{id} ")
        f.write(f"{len(geom):d} ")
        for idx in geom:
            f.write(f"{idx} ")
        f.write("\n")


def write_boundaries(f, grd):
    # ocean boundaries
    f.write(f"{len(grd['boundaries'][None]):d} ")
    f.write("! total number of ocean boundaries\n")
    # count total number of ocean boundaries
    _sum = np.sum(
        [len(boundary) for boundary in grd['boundaries'][None].values()]
        )
    f.write(f"{int(_sum):d} ! total number of ocean boundary nodes\n")
    # write ocean boundary indexes
    for i, boundary in grd['boundaries'][None].items():
        f.write(f"{len(boundary):d}")
        f.write(f" ! number of nodes for ocean_boundary_{i}\n")
        for idx in boundary:
            f.write(f"{idx}\n")
    # count remaining boundaries
    _cnt = 0
    for key in grd['boundaries']:
        if key is not None:
            _cnt += len(grd['boundaries'][key])
    f.write(f"{_cnt:d}  ! total number of non-ocean boundaries\n")
    # count remaining boundary nodes
    _cnt = 0
    for ibtype in grd['boundaries']:
        if ibtype is not None:
            for bnd in grd['boundaries'][ibtype].values():
                _cnt += np.asarray(bnd).size
    f.write(f"{_cnt:d} ! Total number of non-ocean boundary nodes\n")
    # all additional boundaries
    for ibtype, boundaries in grd['boundaries'].items():
        if ibtype is None:
            continue
        for id, boundary in boundaries.items():
            f.write(f"{len(boundary):d} ")
            f.write(f"{ibtype} ")
            f.write(f"! boundary {ibtype}:{id}\n")
            for idx in boundary:
                f.write(f"{idx}\n")


def to_mesh(nodes, elements):
    # cast gr3 inputs into a geomesh structure format
    coords = {id: (x, y) for id, ((x, y), value) in nodes.items()}
    triangles = {id: geom for id, geom in elements.items()
                 if len(geom) == 3}
    quads = {id: geom for id, geom in elements.items()
             if len(geom) == 4}
    values = [value for coord, value in nodes.values()]
    return coords, triangles, quads, values


class Gr3(Gmesh):

    def __init__(
        self,
        nodes,
        elements,
        crs=None,
        description=None,
    ):
        super().__init__(*to_mesh(nodes, elements), crs, description)

    @classmethod
    def open(cls, gr3, crs=None):
        kwargs = reader(gr3)
        kwargs.update({"crs": crs})
        return cls(**kwargs)

    def write(self, path, overwrite=False):
        writer(self.grd, path, overwrite)

    @property
    @lru_cache
    def grd(self):
        description = self.description
        if self.crs is not None:
            description += f" CRS: {self.crs.srs}"
        return {
            "nodes": self.nodes,
            "elements": self.elements,
            "description": description,
        }

    @property
    @lru_cache
    def nodes(self):
        return {id: ((x, y), -self.values[i]) for i, (id, (x, y))
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
