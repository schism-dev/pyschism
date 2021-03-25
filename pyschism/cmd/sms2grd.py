from argparse import Namespace
import pathlib

from pyschism.mesh.parsers import sms2dm
from pyschism.mesh import Gr3


class Sms2grdCli:

    def __init__(self, args: Namespace):

        if args.output_path.exists() and args.overwrite is not True:
            raise IOError('File exists and overwrite is not True.')

        Gr3(**sms2dm_to_gr3(sms2dm.read(args.sms2dm))).write(
            args.output_path, args.overwrite)


def sms2dm_to_gr3(sms2dm):
    nodes = {id: (coords, -value) for id, (coords, value)
             in sms2dm['ND'].items()}
    elements = {}
    if 'E3T' in sms2dm:
        elements.update({i+1: geom for i, geom
                         in enumerate(sms2dm['E3T'].values())})

    if 'E4Q' in sms2dm:
        offset = len(elements)
        elements.update({i+1+offset: geom for i, geom
                         in enumerate(sms2dm['Q4T'].values())})
    return {
        'nodes': nodes,
        'elements': elements,
        # 'boundaries': boundaries,
        'description': '',
        'crs': 'epsg:4326',
    }


def add_sms2grd(subparser):
    # creating manning
    sms2grd = subparser.add_parser('sms2grd')
    sms2grd.add_argument('sms2dm', type=pathlib.Path)
    sms2grd.add_argument('output_path', type=pathlib.Path)
    sms2grd.add_argument('--overwrite', action='store_true')
