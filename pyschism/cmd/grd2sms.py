from argparse import Namespace
import pathlib

from pyschism.mesh import Hgrid


class Grd2SmsCli:

    def __init__(self, args: Namespace):

        if args.output_path.exists() and args.overwrite is not True:
            raise IOError('File exists and overwrite is not True.')

        Hgrid.open(args.hgrid, crs=args.hgrid_crs).write(
            args.output_path, args.overwrite, format='2dm')


def add_grd2sms(subparser):
    grd2sms = subparser.add_parser('grd2sms')
    grd2sms.add_argument('hgrid', type=pathlib.Path)
    grd2sms.add_argument('output_path', type=pathlib.Path)
    grd2sms.add_argument('--overwrite', action='store_true')
    grd2sms.add_argument('--hgrid-crs')
