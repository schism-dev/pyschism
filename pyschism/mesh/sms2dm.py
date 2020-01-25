import pathlib

def reader(path):
    sms2dm = dict()
    with open(pathlib.Path(path), 'r') as f:
        f.readline()
        while 1:
            line = f.readline().split()
            if len(line) == 0:
                break
            if line[0] in ['E3T', 'E4Q']:
                sms2dm[line[0]].update({
                    line[1]: line[2:]
                    })
            if line[0] == 'ND':
                sms2dm[line[0]].update({
                    line[1]: (
                        list(map(float, line[2:-1])), float(line[-1])
                        )
                    })
    return sms2dm

def writer(sms2dm, path, overwrite=False):
    path = pathlib.Path(path)
    if path.is_file() and not overwrite:
        msg = 'File exists, pass overwrite=True to allow overwrite.'
        raise Exception(msg)
    with open(path, 'w') as f:
        f.write(string(sms2dm))
    return 0  # for unittests


def string(sms2dm):
    f = graph(sms2dm)
    f += boundaries(sms2dm)
    return f

def boundaries(sms2dm):
    f = ''
    if 'boundaries' in sms2dm.keys():
        for ibtype, bnds in sms2dm['boundaries'].items():
            for id, bnd in bnds.items():
                f += nodestring(bnd['indexes'])
    return f

def geom_string(geom_type, geom):
    assert geom_type in ['E3T', 'E4Q', 'E6T', 'E8Q', 'E9Q']
    f = ''
    for i in range(len(geom)):
        f += f"{geom_type} {i + 1} "
        for j in range(len(geom[i, :])):
            f += f"{geom[i, j]+1} "
        f += "\n"
    return f


def nodestring(geom):
    raise NotImplementedError


def triangular_elements(geom):
    f = ''
    if geom is not None:
        f += geom_string("E3T", geom)
    return f

def quadrilateral_elements(geom):
    f = ''
    if geom is not None:
        for i in range(len(geom)):
            f += f"E4Q {i + 1} "
            for j in range(len(geom[i, :])):
                f += f"{geom[i, j]+1} "
            f += "\n"
    return f

def graph(sms2dm):
    f = "MESH2D\n"
    # TODO: Make faster using np.array2string
    if 'triangles' in sms2dm:
        f += geom_string("E3T", sms2dm['triangles'])
    if 'quads' in sms2dm:
        f += geom_string("E4Q", sms2dm['quads'])
    for i in range(len(sms2dm['coords'])):
        f += f"ND {i + 1} "
        f += f"{sms2dm['coords'][i][0]:<.16E} "
        f += f"{sms2dm['coords'][i][1]:<.16E} "
        f += f"{-sms2dm['values'][i]:<.16E}\n"
    return f