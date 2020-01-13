import pathlib
import numpy as np
from collections import defaultdict


def reader(path):
    grd = dict()
    grd['nodes'] = defaultdict(list)
    grd['elements'] = defaultdict(list)
    grd['boundaries'] = defaultdict(dict)
    with open(pathlib.Path(path), 'r') as f:
        grd['description'] = f.readline().strip('\n')
        NE, NP = map(int, f.readline().split())
        for i in range(NP):
            id, x, y, z = f.readline().split()
            grd['nodes'][id] = (float(x), float(y), -float(z))
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
                _bnd_id = len(grd['boundaries'][ibtype]) + 1
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


def get_gr3_graph(grd):
    f = f"{grd['description']}\n"
    f += f"{len(grd['elements'])} "
    f += f"{len(grd['nodes'])}\n"
    # TODO: Make faster using np.array2string
    for id, (x, y, z) in grd['nodes'].items():
        f += f"{id} "
        f += f"{x:<.16E} "
        f += f"{y:<.16E} "
        f += f"{z:<.16E}\n"
    for id, geom in grd['elements'].items():
        f += f"{id} "
        f += f"{len(geom):d} "
        for idx in geom:
            f += f"{idx} "
        f += "\n"
    return f


def get_gr3_boundaries(grd):
    # ocean boundaries
    f = ""
    f += f"{len(grd['boundaries'][None]):d} "
    f += "! total number of ocean boundaries\n"
    # count total number of ocean boundaries
    _sum = np.sum(
        [len(boundary) for boundary in grd['boundaries'][None].values()]
        )
    f += f"{int(_sum):d} ! total number of ocean boundary nodes\n"
    # write ocean boundary indexes
    for i, boundary in grd['boundaries'][None].items():
        f += f"{len(boundary):d}"
        f += f" ! number of nodes for ocean_boundary_{i}\n"
        for idx in boundary:
            f += f"{idx}\n"
    # count remaining boundaries
    _cnt = 0
    for key in grd['boundaries']:
        if key is not None:
            _cnt += len(grd['boundaries'][key])
    f += f"{_cnt:d}  ! total number of non-ocean boundaries\n"
    # count remaining boundary nodes
    _cnt = 0
    for ibtype in grd['boundaries']:
        if ibtype is not None:
            for bnd in grd['boundaries'][ibtype].values():
                _cnt += np.asarray(bnd).size
    f += f"{_cnt:d} ! Total number of non-ocean boundary nodes\n"
    # all additional boundaries
    for ibtype, boundaries in grd['boundaries'].items():
        if ibtype is None:
            continue
        for id, boundary in boundaries.items():
            f += f"{len(boundary):d} "
            f += f"{ibtype} "
            f += f"! boundary {ibtype}:{id}\n"
            for idx in boundary:
                f += f"{idx}\n"
    return f


def gr3(grd):
    f = get_gr3_graph(grd)
    if 'boundaries' in grd.keys():
        f += get_gr3_boundaries(grd)
    return f


def writer(grd, path, overwrite=False):
    path = pathlib.Path(path)
    if path.is_file() and not overwrite:
        msg = 'File exists, pass overwrite=True to allow overwrite.'
        raise Exception(msg)
    with open(path, 'w') as f:
        f.write(gr3(grd))
