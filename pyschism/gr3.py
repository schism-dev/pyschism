import pathlib
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
        grd['description'] = f"{f.readline()}"
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
