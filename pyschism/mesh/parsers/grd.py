from collections import defaultdict
import os
import pathlib
from typing import Union, Dict, TextIO, Sequence
import warnings

import numpy as np  # type: ignore[import]
from pyproj import CRS  # type: ignore[import]
from pyproj.exceptions import CRSError  # type: ignore[import]


def buffer_to_dict(buf: TextIO):
    description = buf.readline()
    NE, NP = map(int, buf.readline().split())
    vertex_id = []
    vertices = []
    values = []
    for _ in range(NP):
        line = buf.readline().split()
        vertex_id.append(line[0])
        vertices.append((float(line[1]), float(line[2])))
        values.append(float(line[3]))
    element_id = []
    elements = []
    for _ in range(NE):
        line = buf.readline().split()
        element_id.append(line[0])
        elements.append(line[2:])
    # Assume EOF if NOPE is empty.
    try:
        NOPE = int(buf.readline().split()[0])
    except IndexError:
        return {'description': description,
                'vertices': vertices,
                'elements': elements,
                'values': values,
                'vertex_id': vertex_id,
                'element_id': element_id}
    # let NOPE=-1 mean an ellipsoidal-mesh
    # reassigning NOPE to 0 until further implementation is applied.
    boundaries: Dict = defaultdict(dict)
    _bnd_id = 0
    buf.readline()  # Not used.
    while _bnd_id < NOPE:
        NETA = int(buf.readline().split()[0])
        _cnt = 0
        boundaries[None][_bnd_id] = dict()
        boundaries[None][_bnd_id]['indexes'] = list()
        while _cnt < NETA:
            boundaries[None][_bnd_id]['indexes'].append(
                buf.readline().split()[0].strip())
            _cnt += 1
        _bnd_id += 1
    NBOU = int(buf.readline().split()[0])
    _nbnd_cnt = 0
    buf.readline()  # not used
    while _nbnd_cnt < NBOU:
        npts, ibtype = map(int, buf.readline().split()[:2])
        _pnt_cnt = 0
        if ibtype not in boundaries:
            _bnd_id = 0
        else:
            _bnd_id = len(boundaries[ibtype])
        boundaries[ibtype][_bnd_id] = dict()
        boundaries[ibtype][_bnd_id]['indexes'] = list()
        while _pnt_cnt < npts:
            line = buf.readline().split()
            boundaries[ibtype][_bnd_id]['indexes'].append(line[0])
            _pnt_cnt += 1
        _nbnd_cnt += 1
    return {'description': description,
            'vertices': vertices,
            'elements': elements,
            'values': values,
            'vertex_id': vertex_id,
            'element_id': element_id,
            'boundaries': boundaries}


def dict_to_string(grd):
    """
    must contain keys:
        description
        vertices
        elements
        vertex_id
        element_id
        values
        boundaries (optional)
            indexes
    """
    NE, NP = len(grd['elements']), len(grd['vertices'])
    out = [
        f"{grd['description']}",
        f"{NE} {NP}"
    ]
    for i in range(NP):
        vertices = ' '.join([f'{axis:<.16E}' for axis in grd['vertices'][i]])
        if isinstance(grd['values'][i], Sequence):
            values = ' '.join([f'{value:<.16E}' for value in grd['values'][i]])
        else:
            values = grd['values'][i]
        out.append(f"{grd['vertex_id'][i]} {vertices} {values}")
    for i in range(NE):
        elements = ' '.join([f'{element}' for element in grd['elements'][i]])
        out.append(f"{grd['element_id'][i]} {len(grd['elements'][i])} "
                   f"{elements}")
    if 'boundaries' not in grd:
        return "\n".join(out)
    ocean_boundaries = grd['boundaries'].get(None, {})
    out.append(f"{len(grd['boundaries'][None]):d} ! total number of ocean "
               "boundaries")
    # count total number of ocean boundaries
    if len(ocean_boundaries) > 0:
        cnt = 0
        for id, bnd in ocean_boundaries.items():
            cnt += len(bnd.indexes)
    out.append(f"{cnt:d} ! total number of ocean boundary nodes")
    # write ocean boundary indexes
    for i, boundary in ocean_boundaries.items():
        out.append(
                f"{len(boundary.indexes):d} ! number of nodes for "
                f"ocean_boundary_{i}")
        out.extend([idx for idx in boundary.indexes])
    # remaining boundaries
    cnt = 0
    for key in grd['boundaries']:
        if key is not None:
            for bnd in grd['boundaries'][key]:
                cnt += 1
    out.append(f"{cnt:d}  ! total number of non-ocean boundaries")
    # count remaining boundary nodes
    cnt = 0
    for ibtype in grd['boundaries']:
        if ibtype is not None:
            for bnd in grd['boundaries'][ibtype].values():
                cnt += np.asarray(bnd.indexes).size
    out.append(f"{cnt:d} ! Total number of non-ocean boundary nodes")
    # all additional boundaries
    for ibtype, boundaries in grd['boundaries'].items():
        if ibtype is None:
            continue
        for id, boundary in boundaries.items():
            out.append(f"{len(boundary.indexes):d} {ibtype} ! boundary "
                       f"{ibtype}:{id}")
            # indexes = []
            # for idx in boundary.indexes:
            #     indexes.append(f"{idx}")
            # indexes = ' '.join(indexes)
            # out.append(f"{indexes}")
            out.extend([idx for idx in boundary.indexes])
    return "\n".join(out)


def read(resource: Union[str, os.PathLike], boundaries: bool = True, crs=True):
    """Converts a file-like object representing a grd-formatted unstructured
    mesh into a python dictionary:

    Args:
        resource: Path to file on disk or file-like object such as
            :class:`io.StringIO`
    """
    resource = pathlib.Path(resource)
    with open(resource, 'r') as stream:
        grd = buffer_to_dict(stream)
    if boundaries is False:
        grd.pop('boundaries', None)
    if crs is True:
        crs = None
    if crs is None:
        for try_crs in grd['description'].split():
            try:
                crs = CRS.from_user_input(try_crs)
                break
            except CRSError:
                pass
    if crs is None:
        warnings.warn(f'File {str(resource)} does not contain CRS '
                      'information.')
    if crs is not False:
        grd.update({'crs': crs})
    return grd


def write(grd, path, overwrite=False):
    path = pathlib.Path(path)
    if path.is_file() and not overwrite:
        raise Exception('File exists, pass overwrite=True to allow overwrite.')
    with open(path, 'w') as f:
        f.write(dict_to_string(grd))
