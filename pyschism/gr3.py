import pathlib
import numpy as np
from collections import defaultdict


def parse_gr3(path):
    grd = defaultdict()
    grd['vertices'] = list()
    grd['values'] = list()
    grd['triangles'] = list()
    grd['quads'] = list()
    grd['ocean_boundaries'] = defaultdict(list)
    grd['land_boundaries'] = defaultdict(list)
    grd['interior_boundaries'] = defaultdict(list)
    with open(pathlib.Path(path), 'r') as f:
        grd['description'] = f.readline().strip('\n')
        NE, NP = map(int, f.readline().split())
        for i in range(NP):
            _, x, y, z = f.readline().split()
            grd['vertices'].append((float(x), float(y)))
            grd['values'].append(-float(z))
        for i in range(NE):
            line = f.readline().split()
            geom_t = int(line[1])
            geom = [int(x)-1 for x in line[2:]]
            if geom_t == 3:
                grd['triangles'].append(geom)
            elif geom_t == 4:
                grd['quads'].append(geom)
            else:
                msg = 'ERROR: only supports triangular and quad meshes.'
                raise NotImplementedError(msg)
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
            while _cnt < NETA:
                grd['ocean_boundaries'][_bnd_id].append(
                    int(f.readline().split()[0])-1)
                _cnt += 1
            _bnd_id += 1
        NBOU = int(f.readline().split()[0])
        _nbnd_cnt = 0
        _land_bnd_id = 0
        _inte_bnd_id = 0
        f.readline()  # not used
        while _nbnd_cnt < NBOU:
            npts, ibtype = map(int, f.readline().split()[:2])
            msg = f"ERROR: Found ibtype {ibtype} which is an unknown "
            msg += "boundary type for SCHISM."
            assert ibtype in [0, 1]
            _pnt_cnt = 0
            while _pnt_cnt < npts:
                line = f.readline().split()
                if ibtype == 0:
                    grd['land_boundaries'][_land_bnd_id].append(
                        int(line[0])-1)
                elif ibtype == 1:
                    grd['interior_boundaries'][_inte_bnd_id].append(
                        int(line[0])-1)
                _pnt_cnt += 1
            if ibtype == 0:
                _land_bnd_id += 1
            elif ibtype == 1:
                _inte_bnd_id += 1
            _nbnd_cnt += 1
    return grd


def get_gr3(grd):
    f = f"{grd['description']}\n"
    f += f"{len(grd['triangles']) + len(grd['quads'])} "
    f += f"{grd['values'].shape[0]}\n"
    # TODO: Make faster using np.array2string
    for i in range(grd['values'].shape[0]):
        f += f"{i + 1} "
        f += f"{grd['vertices'][i][0]:<.16E} "
        f += f" {grd['vertices'][i][1]:<.16E} "
        f += f"{-grd['values'][i]:<.16E}\n"
    _cnt = 0
    for geom_t in ['triangles', 'quads']:
        for i in range(len(grd[geom_t])):
            f += f"{_cnt + 1} "
            f += f"{len(grd[geom_t][i]):d} "
            for idx in grd[geom_t][i]:
                f += f"{idx+1:d} "
            f += "\n"
            _cnt += 1
    # ocean boundaries
    f += f"{len(grd['ocean_boundaries'].keys()):d} "
    f += "! total number of ocean boundaries\n"
    # count total number of ocean boundaries
    _sum = np.sum(
        [len(boundary) for boundary in grd['ocean_boundaries'].values()]
        )
    f += f"{int(_sum):d} ! total number of ocean boundary nodes\n"
    # write ocean boundary indexes
    for i, boundary in grd['ocean_boundaries'].items():
        f += f"{len(boundary):d}"
        f += f" ! number of nodes for ocean_boundary_{i}\n"
        for idx in boundary:
            f += f"{idx+1:d}\n"
    # count remaining boundaries
    _cnt = len(grd['land_boundaries'].keys())
    _cnt += len(grd['interior_boundaries'].keys())
    f += f"{_cnt:d}  ! total number of non-ocean boundaries\n"
    _cnt = int(np.sum([len(x) for x in grd['land_boundaries'].values()]))
    _cnt += int(
        np.sum([len(x) for x in grd['interior_boundaries'].values()]))
    f += f"{_cnt:d} ! Total number of non-ocean boundary nodes\n"
    # write land boundaries
    for i, boundary in grd['land_boundaries'].items():
        f += f"{len(boundary):d} "
        f += f"0 "
        f += "! number of nodes for land_boundary_"
        f += f"{i}\n"
        for idx in boundary:
            f += f"{idx+1:d}\n"
    # write inner boundaries
    for i, boundary in grd['interior_boundaries'].items():
        f += f"{len(boundary):d} "
        f += f"1 "
        f += "! number of nodes for inner_boundary_"
        f += f"{i + len(grd['land_boundaries'].keys())}\n"
        for idx in boundary:
            f += f"{idx+1:d}\n"
    return f


def write_gr3(grd, path, overwrite=False):
    path = pathlib.Path(path)
    if path.is_file() and not overwrite:
        msg = 'File exists, pass overwrite=True to allow overwrite.'
        raise Exception(msg)
    with open(path, 'w') as f:
        f.write(get_gr3(grd))
