from argparse import Namespace

from pyschism.mesh.parsers import sms2dm
from pyschism.mesh import Gr3


def sms2dm_to_gr3(sms2dm):

    elements = {}
    if 'E3T' in sms2dm:
        elements.update({i+1: geom for i, geom
                         in enumerate(sms2dm['E3T'].values())})

    if 'E4Q' in sms2dm:
        offset = len(elements)
        elements.update({i+1+offset: geom for i, geom
                         in enumerate(sms2dm['Q4T'].values())})
    return {
        'nodes': sms2dm['ND'],
        'elements': elements,
        # 'boundaries': boundaries,
        'description': '',
        'crs': 'epsg:4326',
    }


def add_sms2grd(subparser):
    # creating manning
    sms2grd = subparser.add_parser('sms2grd')
    sms2grd.add_argument('sms2dm')
    sms2grd.add_argument('output_path')
    sms2grd.add_argument('--overwrite', action='store_true')


def sms2grd(args: Namespace):
    Gr3(**sms2dm_to_gr3(sms2dm.read(args.sms2dm))).write(
        args.output_path, args.overwrite)
