from collections import defaultdict
import os
import numbers
import pathlib
from typing import Union, Dict, TextIO
import warnings

import numpy as np  # type: ignore[import]
from pyproj import CRS  # type: ignore[import]
from pyproj.exceptions import CRSError  # type: ignore[import]


def buffer_to_dict(buf: TextIO):
    description = buf.readline().strip()
    NE, NP = map(int, buf.readline().split())
    nodes = {}
    for _ in range(NP):
        line = buf.readline().strip('\n').split()
        # Gr3/fort.14 format cannot distinguish between a 2D mesh with one
        # vector value (e.g. velocity, which uses 2 columns) or a 3D mesh with
        # one scalar value. This is a design problem of the mesh format, which
        # renders it ambiguous, and the main reason why the use of fort.14/grd
        # formats is discouraged, in favor of UGRID.
        # Here, we assume the input mesh is strictly a 2D mesh, and the data
        # that follows is an array of values.
        if len(line[3:]) == 1:
            nodes[line[0]] = [
                (float(line[1]), float(line[2])), float(line[3])]
        else:
            nodes[line[0]] = [
                (float(line[1]), float(line[2])),
                [float(line[i]) for i in range(3, len(line[3:]))]
            ]
    elements = {}
    for _ in range(NE):
        line = buf.readline().split()
        elements[line[0]] = line[2:]
    # Assume EOF if NOPE is empty.
    try:
        NOPE = int(buf.readline().split()[0])
    except IndexError:
        return {'description': description,
                'nodes': nodes,
                'elements': elements}
    # let NOPE=-1 mean an ellipsoidal-mesh
    # reassigning NOPE to 0 until further implementation is applied.
    boundaries: Dict = defaultdict(dict)
    _bnd_id = 0
    buf.readline()
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
    buf.readline()
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
            if len(line) == 1:
                boundaries[ibtype][_bnd_id]['indexes'].append(line[0])
            else:
                index_construct = []
                for val in line:
                    if '.' in val:
                        continue
                    index_construct.append(val)
                boundaries[ibtype][_bnd_id]['indexes'].append(index_construct)
            _pnt_cnt += 1
        _nbnd_cnt += 1
    return {'description': description,
            'nodes': nodes,
            'elements': elements,
            'boundaries': boundaries}


def to_string(description, nodes, elements, boundaries=None, crs=None):
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
    NE, NP = len(elements), len(nodes)
    out = [f"{description}", f"{NE} {NP}"]
    # TODO: Probably faster if using np.array2string
    for id, (coords, values) in nodes.items():
        if isinstance(values, numbers.Number):
            values = [values]
        line = [f"{id}"]
        line.extend([f"{x:<.8f}" for x in coords])
        line.extend([f"{x:<.8f}" for x in values])
        out.append(" ".join(line))

    for id, element in elements.items():
        line = [f"{id}"]
        line.append(f"{len(element)}")
        line.extend([f"{e}" for e in element])
        out.append(" ".join(line))
    if boundaries is None:
        out.append('')
        return "\n".join(out)
    if boundaries is not None:
        out.append(f"{len(boundaries[None]):d} "
                   "! total number of ocean boundaries")
        # count total number of ocean boundaries
        _sum = 0
        for bnd in boundaries[None].values():
            _sum += len(bnd['indexes'])
        out.append(f"{int(_sum):d} ! total number of ocean boundary nodes")
        # write ocean boundary indexes
        for i, boundary in boundaries[None].items():
            out.append(f"{len(boundary['indexes']):d}"
                       f" ! number of nodes for ocean_boundary_{i}")
            for idx in boundary['indexes']:
                out.append(f"{idx}")
    else:
        out.append("0 ! total number of ocean boundaries")
        out.append("0 ! total number of ocean boundary nodes")
    # remaining boundaries
    boundaries = {} if boundaries is None else boundaries
    _cnt = 0
    for key in boundaries:
        if key is not None:
            for bnd in boundaries[key]:
                _cnt += 1
    out.append(f"{_cnt:d}  ! total number of non-ocean boundaries")
    # count remaining boundary nodes
    _cnt = 0
    for ibtype in boundaries:
        if ibtype is not None:
            for bnd in boundaries[ibtype].values():
                _cnt += np.asarray(bnd['indexes']).size
    out.append(f"{_cnt:d} ! Total number of non-ocean boundary nodes")
    # all additional boundaries
    for ibtype, boundaries in boundaries.items():
        if ibtype is None:
            continue
        for id, boundary in boundaries.items():
            line = [
                f"{len(boundary['indexes']):d}",
                f"{ibtype}",
                f"! boundary {ibtype}:{id}"]
            out.append(' '.join(line))
            for idx in boundary['indexes']:
                out.append(f"{idx}")
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
                      'information and no CRS was given.')
    if crs is not False:
        grd.update({'crs': crs})
    return grd


def write(grd, path, overwrite=False):
    path = pathlib.Path(path)
    if path.is_file() and not overwrite:
        raise Exception('File exists, pass overwrite=True to allow overwrite.')
    with open(path, 'w') as f:
        f.write(to_string(**grd))
