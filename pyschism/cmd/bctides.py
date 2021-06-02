import argparse
from datetime import datetime, timedelta
import json
import logging
import pathlib
import warnings

from pyproj import CRS

from pyschism.forcing.baroclinic import GOFS  # , RTOFS
# from pyschism.forcing.baroclinic.gofs import GOFSElevation
# from pyschism.forcing.baroclinic.rtofs import RTOFSElevation
from pyschism.forcing.tides import Tides
from pyschism.forcing.bctides import Bctides, iettype, ifltype, itetype, isatype
from pyschism.mesh import Hgrid, Vgrid


logger = logging.getLogger(__name__)


baroclinic_databases = {
    'gofs': GOFS,
    # 'rtofs': RTOFSElevation,
}


def get_tides(args: argparse.Namespace):
    def get_elevation():
        return True if args.tidal_database is not None else False

    def get_velocity():
        vgrid = Vgrid.default() if args.vgrid is None else args.vgrid
        return (
            True
            if args.include_velocity is True or vgrid.is3D() is True
            else False
        )

    return Tides(
        elevation=get_elevation(),
        velocity=get_velocity(),
        tidal_database=args.tidal_database,
    )


class BctidesCli:
    def __init__(self, args: argparse.Namespace):

        # for boundary in args.hgrid.boundaries.open.itertuples():
        #     args.hgrid.boundaries.set_forcing(
        #         boundary.id,
        #         iettype=args.iettype,
        #         ifltype=args.ifltype,
        #         itetype=args.itetype,
        #         isatype=args.isatype,
        #         # itrtype=args.itrtype,
        #     )

        bctides = Bctides(
                args.hgrid,
                args.start_date,
                args.run_days,
                vgrid=args.vgrid,
            )

        if args.Z0 is not None:
            bctides.Z0 = args.Z0

        bctides.write(args.output_directory, overwrite=args.overwrite)


class HgridAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hgrid = Hgrid.open(values)
        if len(hgrid.boundaries.open) == 0:
            raise TypeError(f"Hgrid provided {values} contains no open boundaries.")
        setattr(namespace, self.dest, hgrid)


class HgridCrsAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values is not None:
            namespace.hgrid.nodes._crs = CRS.from_user_input(values)
        setattr(namespace, self.dest, values)


class VgridAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, Vgrid.open(values))


class CustomBoundaryAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if isinstance(values, str):
            if "iettype" in option_string:
                values = json.loads(values)
                for bnd_id, tval in values.items():
                    if not isinstance(tval, int):
                        raise TypeError(
                            'Boundary type must be an integer, not type '
                            f'{type(tval)}.')
                    if tval == 3:
                        namespace.hgrid.boundaries.set_forcing(
                            bnd_id,
                            iettype=iettype.Iettype3(get_tides(namespace))
                            )
                    else:
                        raise ValueError
        setattr(namespace, self.dest, values)


class Iettype3Action(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        for boundary in namespace.hgrid.boundaries.open.itertuples():
            namespace.hgrid.boundaries.set_forcing(
                boundary.id,
                iettype=self.const(get_tides(namespace)))
        setattr(namespace, self.dest, self.const)


class Iettype4Action(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        for boundary in namespace.hgrid.boundaries.open.itertuples():
            namespace.hgrid.boundaries.set_forcing(
                boundary.id,
                iettype=self.const(baroclinic_databases[
                    namespace.baroclinic_database]()))
        setattr(namespace, self.dest, self.const)


class Iettype5Action(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        for boundary in namespace.hgrid.boundaries.open.itertuples():
            namespace.hgrid.boundaries.set_forcing(
                boundary.id,
                iettype=self.const(
                    get_tides(namespace),
                    baroclinic_databases[namespace.baroclinic_database]()
                    )
            )
        setattr(namespace, self.dest, self.const)


class Ifltype3Action(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        for boundary in namespace.hgrid.boundaries.open.itertuples():
            namespace.hgrid.boundaries.set_forcing(
                boundary.id,
                ifltype=values(get_tides(namespace)))
        setattr(namespace, self.dest, values)


class Ifltype4Action(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        for boundary in namespace.hgrid.boundaries.open.itertuples():
            namespace.hgrid.boundaries.set_forcing(
                boundary.id,
                ifltype=self.const(baroclinic_databases[
                    namespace.baroclinic_database]()))
        setattr(namespace, self.dest, values)


class Ifltype5Action(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        for boundary in namespace.hgrid.boundaries.open.itertuples():
            namespace.hgrid.boundaries.set_forcing(
                boundary.id,
                ifltype=self.const(
                    get_tides(namespace),
                    baroclinic_databases[namespace.baroclinic_database]()
                    )
            )
        setattr(namespace, self.dest, values)


def add_bctypes(bctides):
    _iettype = bctides.add_mutually_exclusive_group()
    _iettype.add_argument(
        "--iettype",
        dest="iettype",
        action=CustomBoundaryAction,
        const=iettype.Iettype,
        help='Per-boundary specification option for elevation variable. '
             'The format required is a json-style string. For example: '
             '--iettype=\'{"1": 3}\' would apply iettype-3 to boundary with '
             'id equal to "1".'
    )
    _iettype.add_argument(
        "--elev-th",
        "--iettype-1",
        dest="iettype",
        type=iettype.Iettype1,
        help='Global elevation options for time history.',
    )
    _iettype.add_argument(
        '--elev-val',
        '--iettype-2',
        dest='iettype',
        type=iettype.Iettype2,
        help='Global constant elevation option.',
    )

    _iettype.add_argument(
        "--elev-tides",
        "--iettype-3",
        dest="iettype",
        nargs=0,
        const=iettype.Iettype3,  # the option string is present but not followed by a command-line argument. In this case the value from const will be produced.
        action=Iettype3Action,
    )

    _iettype.add_argument(
        '--elev-subtides',
        '--elev-2d',
        '--iettype-4',
        dest='iettype',
        nargs=0,
        const=iettype.Iettype4,
        action=Iettype4Action,
    )

    _iettype.add_argument(
        '--elev-tides-subtides',
        '--tides-elev-2d',
        '--elev-2d-tides',
        '--iettype-5',
        dest='iettype',
        nargs=0,
        const=iettype.Iettype5,
        action=Iettype5Action,
    )

    _iettype.add_argument(
        '--elev-zero',
        '--iettype-_1',
        dest='iettype',
        type=iettype.Iettype_1,
    )
    _ifltype = bctides.add_mutually_exclusive_group()
    _ifltype.add_argument(
        "--ifltype",
        dest="ifltype",
        nargs=1,  # 0 or more
        const={},
        action=CustomBoundaryAction,
    )
    _ifltype.add_argument(
        '--flux-th',
        '--ifltype-1',
        dest='ifltype',
        type=ifltype.Ifltype1
    )
    _ifltype.add_argument(
        '--flux-val',
        '--ifltype-2',
        dest='ifltype',
        type=ifltype.Ifltype2
    )

    _ifltype.add_argument(
        '--uv-tides',
        '--ifltype-3',
        nargs="?",  # 0 or more
        action=Ifltype3Action,
        dest='ifltype',
        const=ifltype.Ifltype3,
    )

    _ifltype.add_argument(
        '--uv-subtides',
        '--uv-3d',
        '--uv3D',
        '--ifltype-4',
        action=Ifltype4Action,
        dest='ifltype',
        nargs=0,
        const=ifltype.Ifltype4,
    )

    _ifltype.add_argument(
        '--uv-tides-subtides',
        '--uv-3d-tides',
        '--uv3d-tides',
        '--uv3D-tides',
        '--uv3D-tides',
        '--ifltype-5',
        action=Ifltype5Action,
        nargs=0,
        dest='ifltype',
        const=ifltype.Ifltype5,
    )

    _ifltype.add_argument(
        '--uv-zero',
        '--flather',
        action='store_const',
        dest='ifltype',
        const=ifltype.Ifltype_1,
    )
    _itetype = bctides.add_mutually_exclusive_group()
    _itetype.add_argument(
        "--itetype",
        dest="itetype",
        nargs=1,  # 0 or more
        const={},
        action=CustomBoundaryAction,
    )
    _itetype.add_argument(
        '--temp-th',
        '--itetype-1',
        dest='itetype',
        type=itetype.Itetype1
    )
    _itetype.add_argument(
        '--temp-val',
        '--itetype-2',
        dest='itetype',
        type=itetype.Itetype2
    )

    _itetype.add_argument(
        '--temp-ic',
        '--itetype-3',
        dest='itetype',
        type=itetype.Itetype3,
    )

    _itetype.add_argument(
        '--temp-3d',
        '--itetype-4',
        action='store_const',
        dest='itetype',
        const=itetype.Itetype4,
    )
    _isatype = bctides.add_mutually_exclusive_group()
    _isatype.add_argument(
        "--isatype",
        dest="isatype",
        nargs=1,  # 0 or more
        const={},
        action=CustomBoundaryAction,
    )
    _isatype.add_argument(
        '--salt-th',
        '--isatype-1',
        dest='isatype',
        type=isatype.Isatype1
    )
    _isatype.add_argument(
        '--salt-val',
        '--isatype-2',
        dest='isatype',
        type=isatype.Isatype2
    )

    _isatype.add_argument(
        '--salt-ic',
        '--isatype-3',
        dest='isatype',
        type=isatype.Isatype3,
    )

    _isatype.add_argument(
        '--salt-3d',
        '--isatype-4',
        action='store_const',
        dest='isatype',
        const=isatype.Isatype4,
    )


def add_bctides(subparsers):
    bctides = subparsers.add_parser("bctides")
    bctides.add_argument(
        "hgrid",
        action=HgridAction,
    )
    bctides.add_argument(
        '--hgrid-crs',
        action=HgridCrsAction,
    )
    bctides.add_argument(
        "--overwrite", action="store_true", help="Allow overwrite of output file."
    )
    bctides.add_argument(
        "--log-level",
        choices=[name.lower() for name in logging._nameToLevel],
        default="warning",
    )
    bctides.add_argument(
        "--run-days", type=lambda x: timedelta(days=float(x)), required=True
    )
    bctides.add_argument(
        "--start-date",
        type=lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S"),
        # help=''
        # action=DateAction,
    )
    # bctides.add_argument(
    #     "--date-formatting",
    #     type=str,
    #     type=lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S"),
    #     # help=''
    # )
    bctides.add_argument("--output-directory", "-o", type=pathlib.Path, required=True)
    bctides.add_argument(
        "--vgrid",
        action=VgridAction,
        default=Vgrid.default(),
    )
    bctides.add_argument(
        '--tides',
        '--tidal-database',
        choices=['hamtide', 'tpxo'],
        default='hamtide',
        dest='tidal_database'
    )
    bctides.add_argument(
        '--hycom',
        '--baroclinic-database',
        choices=['rtofs', 'gofs'],
        dest='baroclinic_database',
        default='gofs',
    )
    bctides.add_argument(
        '--include-velocity',
        action='store_true'
    )
    bctides.add_argument('--Z0', type=float)
    bctides.add_argument("--cutoff-depth", type=float, default=50.0)
    add_bctypes(bctides)


# --- drafts

# @property
# def logger(self):
#     if not hasattr(self, "_logger"):
#         root_logger = logging.getLogger()
#         root_logger.setLevel(logging._nameToLevel[self.args.log_level.upper()])
#         self._logger = logging.getLogger(__name__)
#         self._logger.setLevel(logging._nameToLevel[self.args.log_level.upper()])
#     return self._logger

# if args.output_file is not None:
    #     bctides.write(args.output_file, overwrite=args.overwrite)
    # else:
    #     p rint(str(bctides))
    # @property
    # def args(self):
    #     return self._args

    # @args.setter
    # def args(self, args: argparse.Namespace):

    #     vgrid = Vgrid.default() if args.vgrid is None else Vgrid.open(args.vgrid)

    #     def get_tides():
    #         def get_elevation():
    #             return True if args.tidal_database is not None else False

    #         def get_velocity():
    #             return (
    #                 True
    #                 if args.include_velocity is True or vgrid.is3D() is True
    #                 else False
    #             )

    #         return Tides(
    #             elevation=get_elevation(),
    #             velocity=get_velocity(),
    #             tidal_database=args.tidal_database,
    #         )

    #     tides = get_tides()

    #     baroclinic = {"gofs": GOFS, "rtofs": RTOFS}[args.baroclinic_database]()

    #     if args.iettype is None:
    #         args.iettype = (
    #             iettype.Iettype3(tides)
    #             if vgrid.is2D()
    #             else iettype.Iettype5(tides, baroclinic)
    #         )

    #     elif issubclass(args.iettype, iettype.Iettype):
    #         if issubclass(args.iettype, iettype.Iettype5):
    #             args.iettype = args.iettype(tides, baroclinic)
    #         elif issubclass(args.iettype, iettype.Iettype3):
    #             args.iettype = args.iettype(tides)
    #         elif issubclass(args.iettype, iettype.Iettype4):
    #             args.iettype = args.iettype(baroclinic)

    #     if args.ifltype is None:
    #         if vgrid.is3D() is True:
    #             args.ifltype = ifltype.Ifltype5(tides, baroclinic)
    #         elif args.include_velocity is True:
    #             args.ifltype = ifltype.Ifltype3(tides)
    #     elif issubclass(args.ifltype, ifltype.Ifltype):
    #         if issubclass(args.ifltype, ifltype.Ifltype4):
    #             args.ifltype = args.ifltype(baroclinic)
    #         elif issubclass(args.ifltype, ifltype.Ifltype5):
    #             args.ifltype = args.ifltype(tides, baroclinic)

    #     def assert_type(name, obj, cls):
    #         if obj is not None:
    #             assert isinstance(obj, cls), (
    #                 f"args.{name} must be an instance of type {cls}, but "
    #                 f"got an object {obj} of type {type(obj)}."
    #             )

    #     assert_type("iettype", args.iettype, iettype.Iettype)
    #     assert_type("ifltype", args.ifltype, ifltype.Ifltype)
    #     assert_type("itetype", args.itetype, itetype.Itetype)
    #     assert_type("isatype", args.isatype, isatype.Isatype)

    #     self._args = args

    # @property
    # def bctides(self):
    #     if not hasattr(self, "_bctides"):
    #         self._bctides = Bctides(
    #             self.hgrid,
    #             self.start_date,
    #             self.args.run_days,
    #             vgrid=self.vgrid,
    #             # elevation=self.args.iettype,
    #             # velocity=self.args.ifltype,
    #             # temperature=self.args.itetype,
    #             # salinity=self.args.isatype,
    #         )
    #         if self.args.Z0 is not None:
    #             self._bctides.Z0 = self.args.Z0

    #     return self._bctides

    # @property
    # def hgrid(self):
    #     if not hasattr(self, "_hgrid"):
    #         self.logger.info(f"Open hgrid file: {self.args.hgrid}")
    #         start = time()
    #         self._hgrid = Hgrid.open(self.args.hgrid, crs=self.args.hgrid_crs)
    #         self.logger.info(f"Reading hgrid file took {time()-start} seconds.")
    #         # init boundary
    #         # for boundary in
    #     return self._hgrid

    # @property
    # def start_date(self):
    #     if self.args.start_date is None:
    #         self.args.start_date = (
    #             dates.nearest_cycle()
    #             if self.args.start_date is None
    #             else self.args.start_date
    #         )
    #     return self.args.start_date

    # inp = args.inputstring
    # add_bctypes(bctides, bnd_id=1)
    # for i in range()

        # add_bctypes(bctides, bnd_id=1)
    # boundaries = bctides.add_subparsers(dest='boundaries_subparsers')
    # iettype_bctypes = boundaries.add_parsers('iettype_bc')
    # boundaries = bctypes.add_parser('--boundaries')
    # add_bctypes(boundaries, bnd_id=1)
    # bctides.add_argument(
        # '--iettype-3-1',
        # '-bnd',
        # nargs=1,
        # action=BoundaryAction
        # action='append',
        # default=[],
        # const=[],
    # )

    # bctides.add_argument(
    #     "--boundary",
    #     "-bnd",
    #     dest="custom_boundaries",
    #     action=CustomBoundariesAction,
    # )
    # import sys
    # known_args = sys.argv[sys.argv.index('bctides')+1:]
    # h_pop_ = 
    # known_args.pop()
    # namespace, _ = bctides.parse_known_args()
    # hgrid = Hgrid.open(namespace.hgrid, crs=namespace.hgrid_crs)
        # boundaries = parser.add_subparsers(dest="boundaries")
        # parsers = []
        # iettypes = 
        # add_bctypes(parser)
    # for i in range(len(hgrid.boundaries.open)):
        # parsers.append(boundaries.add_parser(f"boundary_{i+1}"))
        # add_bctypes(bctides, bnd_id=f'-{i+1}')
            # parsers[i].add_mutually_exclusive_group()
            # parsers[i]s.add_argument(
            #     f"--boundary-{i}",
            #     "-bnd",
            #     required=True,
            # )
        # add a required argument for each boundary
    # _iettype = bctides.add_mutually_exclusive_group(required=True)
    # _iettype.add_argument(
    #     "--elev-th",
    #     "--iettype-1",
    #     dest="iettype",
    #     nargs="?",  # 0 or more
    #     const=iettype.Iettype1,
    #     help='Global elevation options for '
    #     # action=IettypeAction,
    #     # type=iettype.Iettype1,
    # )
    # # hgrid = Hgrid.open()
    # breakpoint()
    # for i in range(len(hgrid.boundaries.open)):
    #     pass
    # _iettype.add_argument(
    #     '--elev-val',
    #     '--iettype-2',
    #     dest='iettype',
    #     type=iettype.Iettype2
    # )

    # _iettype.add_argument(
    #     "--elev-tides",
    #     "--iettype-3",
    #     dest="iettype",
    #     nargs="?",  # 0 or more
    #     const=iettype.Iettype3,  # the option string is present but not followed by a command-line argument. In this case the value from const will be produced.
    #     action=Iettype3Action,
    # )

    # _iettype.add_argument(
    #     '--elev-subtides',
    #     '--elev-2d',
    #     '--iettype-4',
    #     action='store_const',
    #     dest='iettype',
    #     const=iettype.Iettype4,
    # )

    # _iettype.add_argument(
    #     '--elev-tides-subtides',
    #     '--tides-elev-2d',
    #     '--elev-2d-tides',
    #     '--iettype-5',
    #     action='store_const',
    #     dest='iettype',
    #     const=iettype.Iettype5,
    # )

    # _iettype.add_argument(
    #     '--elev-zero',
    #     '--iettype-_1',
    #     dest='iettype',
    #     type=iettype.Iettype_1,
    # )
    # _ifltype = bctides.add_mutually_exclusive_group()
    # _ifltype.add_argument(
    #     '--flux-th',
    #     '--ifltype-1',
    #     dest='ifltype',
    #     type=ifltype.Ifltype1
    # )
    # _ifltype.add_argument(
    #     '--flux-val',
    #     '--ifltype-2',
    #     dest='ifltype',
    #     type=ifltype.Ifltype2
    # )

    # _ifltype.add_argument(
    #     '--uv-tides',
    #     '--ifltype-3',
    #     action='store_const',
    #     dest='ifltype',
    #     const=ifltype.Ifltype3,
    # )

    # _ifltype.add_argument(
    #     '--uv-subtides',
    #     '--uv-3d',
    #     '--uv3D',
    #     '--ifltype-4',
    #     action='store_const',
    #     dest='ifltype',
    #     const=ifltype.Ifltype4,
    # )

    # _ifltype.add_argument(
    #     '--uv-tides-subtides',
    #     '--uv-3d-tides',
    #     '--uv3d-tides',
    #     '--uv3D-tides',
    #     '--uv3D-tides',
    #     '--ifltype-5',
    #     action='store_const',
    #     dest='ifltype',
    #     const=ifltype.Ifltype5,
    # )

    # _ifltype.add_argument(
    #     '--uv-zero',
    #     '--flather',
    #     action='store_const',
    #     dest='ifltype',
    #     const=ifltype.Ifltype_1,
    # )
    # _itetype = bctides.add_mutually_exclusive_group()
    # _itetype.add_argument(
    #     '--temp-th',
    #     '--itetype-1',
    #     dest='itetype',
    #     type=itetype.Itetype1
    # )
    # _itetype.add_argument(
    #     '--temp-val',
    #     '--itetype-2',
    #     dest='itetype',
    #     type=itetype.Itetype2
    # )

    # _itetype.add_argument(
    #     '--temp-ic',
    #     '--itetype-3',
    #     dest='itetype',
    #     type=itetype.Itetype3,
    # )

    # _itetype.add_argument(
    #     '--temp-3d',
    #     '--itetype-4',
    #     action='store_const',
    #     dest='itetype',
    #     const=itetype.Itetype4,
    # )
    # _isatype = bctides.add_mutually_exclusive_group()
    # _isatype.add_argument(
    #     '--salt-th',
    #     '--isatype-1',
    #     dest='isatype',
    #     type=isatype.Isatype1
    # )
    # _isatype.add_argument(
    #     '--salt-val',
    #     '--isatype-2',
    #     dest='isatype',
    #     type=isatype.Isatype2
    # )

    # _isatype.add_argument(
    #     '--salt-ic',
    #     '--isatype-3',
    #     dest='isatype',
    #     type=isatype.Isatype3,
    # )

    # _isatype.add_argument(
    #     '--salt-3d',
    #     '--isatype-4',
    #     action='store_const',
    #     dest='isatype',
    #     const=isatype.Isatype4,
    # )
    # bctides.add_argument(
    #     '-b',
    #     '--boundary',
    #     type=str,
    #     nargs='+',
    #     action='append',
    #     help='file list'
    # )
    # boundaries = bctides.add_subparsers()
    # custom_bounds = boundaries.add_parser('custom_bounds')
    # custom_bounds.add_argument(
    #     '-b',
    #     '--boundary',
    #     type=str,
    #     nargs='+',
    #     action='append',
    #     help='file list'
    # )
