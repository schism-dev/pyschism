from argparse import Namespace
from datetime import datetime, timedelta
import logging
import pathlib
from time import time

from pyschism import dates
from pyschism.forcing.baroclinic import GOFS, RTOFS
from pyschism.forcing.tides import Tides
from pyschism.forcing.bctides import (
    Bctides, iettype, ifltype, itetype, isatype
    )
from pyschism.mesh import Hgrid, Vgrid


class BctidesCli:

    def __init__(self, args: Namespace):
        self.args = args
        # if args.output_file is not None:
        #     bctides.write(args.output_file, overwrite=args.overwrite)
        # else:
        #     print(str(bctides))
        self.bctides.write(self.args.output_directory,
                           overwrite=self.args.overwrite)

    @property
    def args(self):
        return self._args

    @args.setter
    def args(self, args: Namespace):

        vgrid = Vgrid.default() if args.vgrid is None else Vgrid.open(args.vgrid)

        def get_tides():

            def get_elevation():
                return True if args.tidal_database is not None else False

            def get_velocity():
                return True if args.include_velocity is True \
                    or vgrid.is3D() is True else False
            return Tides(
                elevation=get_elevation(),
                velocity=get_velocity(),
                tidal_database=args.tidal_database
            )

        tides = get_tides()

        baroclinic = {
                    'gofs': GOFS,
                    'rtofs': RTOFS
                }[args.baroclinic_database]()

        if args.iettype is None:
            args.iettype = iettype.Iettype3(tides) \
                if vgrid.is2D() else iettype.Iettype5(tides, baroclinic)

        elif issubclass(args.iettype, iettype.Iettype):
            if issubclass(args.iettype, iettype.Iettype5):
                args.iettype = args.iettype(tides, baroclinic)
            elif issubclass(args.iettype, iettype.Iettype3):
                args.iettype = args.iettype(tides)
            elif issubclass(args.iettype, iettype.Iettype4):
                args.iettype = args.iettype(baroclinic)

        if args.ifltype is None:
            if vgrid.is3D() is True:
                args.ifltype = ifltype.Ifltype5(tides, baroclinic)
            elif args.include_velocity is True:
                args.ifltype = ifltype.Ifltype3(tides)
        elif issubclass(args.ifltype, ifltype.Ifltype):
            if issubclass(args.ifltype, ifltype.Ifltype4):
                args.ifltype = args.ifltype(baroclinic)
            elif issubclass(args.ifltype, ifltype.Ifltype5):
                args.ifltype = args.ifltype(tides, baroclinic)

        def assert_type(name, obj, cls):
            if obj is not None:
                assert isinstance(obj, cls), \
                    f'args.{name} must be an instance of type {cls}, but ' \
                    f'got an object {obj} of type {type(obj)}.'

        assert_type('iettype', args.iettype, iettype.Iettype)
        assert_type('ifltype', args.ifltype, ifltype.Ifltype)
        assert_type('itetype', args.itetype, itetype.Itetype)
        assert_type('isatype', args.isatype, isatype.Isatype)

        self._args = args

    @property
    def logger(self):
        if not hasattr(self, '_logger'):
            root_logger = logging.getLogger()
            root_logger.setLevel(logging._nameToLevel[self.args.log_level.upper()])
            self._logger = logging.getLogger(__name__)
            self._logger.setLevel(logging._nameToLevel[self.args.log_level.upper()])
        return self._logger

    @property
    def bctides(self):
        if not hasattr(self, '_bctides'):
            self._bctides = Bctides(
                self.hgrid,
                self.start_date,
                self.args.run_days,
                vgrid=self.vgrid,
                elevation=self.args.iettype,
                velocity=self.args.ifltype,
                temperature=self.args.itetype,
                salinity=self.args.isatype,
            )
            if self.args.Z0 is not None:
                self._bctides.Z0 = self.args.Z0

        return self._bctides

    @property
    def hgrid(self):
        if not hasattr(self, '_hgrid'):
            self.logger.info(f'Open hgrid file: {self.args.hgrid}')
            start = time()
            self._hgrid = Hgrid.open(self.args.hgrid, crs=self.args.hgrid_crs)
            self.logger.info(f'Reading hgrid file took {time()-start} seconds.')
        return self._hgrid

    @property
    def start_date(self):
        if self.args.start_date is None:
            self.args.start_date = dates.nearest_cycle() if \
                self.args.start_date is None else self.args.start_date
        return self.args.start_date

    @property
    def vgrid(self):
        if not hasattr(self, '_vgrid'):
            self._vgrid = Vgrid.default() if self.args.vgrid is None else \
                Vgrid.open(self.args.vgrid)
        return self._vgrid

    @vgrid.setter
    def vgrid(self, vgrid):
        if vgrid is None:
            vgrid = Vgrid.default()
        assert isinstance(vgrid, Vgrid)
        self._vgrid = vgrid

    # @property
    # def iettype(self):
    #     return self.args.iettype
    #     if not hasattr(self, '_iettype'):
    #         if self.args.iettype is None:
    #             if self.baroclinic is None:
    #                 self._iettype = iettype.Iettype3(self.tides)
    #             else:
    #                 self._iettype = iettype.Iettype5(
    #                     self.tides,
    #                     self.baroclinic
    #                 )
    #         else:
    #             self._iettype = self.args.iettype
    #     # print(self._iettype)
    #     # exit()
    #     return self._iettype

    # @property
    # def ifltype(self):
    #     if not hasattr(self, '_ifltype'):
    #         if self.args.ifltype is None:
    #             # if self.vgrid.is3D():
    #             #     self.args.include_velocity = True
    #         #     if self.baroclinic is None:
    #         #         self._ifltype = ifltype.Ifltype3(self.tides)
    #         #     else:
    #         #         self._ifltype = ifltype.Ifltype5(
    #         #             self.tides,
    #         #             self.baroclinic
    #         #         )
    #         # else:
    #         #     self._ifltype = self.args.ifltype
    #     return self._ifltype

    # @property
    # def itetype(self):
    #     if not hasattr(self, '_itetype'):
    #         if self.args.itetype is None:
    #             if self.baroclinic is not None:
    #                 self._itetype = itetype.Itetype4(self.baroclinic)
    #             else:
    #                 self._itetype = None
    #         else:
    #             self._itetype = self.args.itetype
    #     return self._itetype

    # @property
    # def isatype(self):
    #     if not hasattr(self, '_isatype'):
    #         if self.args.isatype is None:
    #             if self.baroclinic is not None:
    #                 self._isatype = isatype.Isatype4(self.baroclinic)
    #             else:
    #                 self._isatype = None
    #         else:
    #             self._isatype = self.args.isatype
    #     return self._isatype

    # @property
    # def elevation(self):
    #     if not hasattr(self, '_elevation'):

    #         if self.tides is not None and self.baroclinic is None:
    #             self._elevation = iettype.Iettype3(self.tides)

    #         elif self.tides is None and self.baroclinic is not None:
    #             self._elevation = iettype.Iettype4(self.baroclinic)

    #         elif self.tides is not None and self.baroclinic is not None:
    #             self._elevation = iettype.Iettype5(
    #                     self.tides,
    #                     self.baroclinic
    #                 )

    #         elif self.tides is None and self.baroclinic is None:
    #             raise NotImplementedError('No tides or baroclinic forcing.')

    #     return self._elevation

    # @property
    # def velocity(self):
    #     if not hasattr(self, '_velocity'):

    #         if self.tides is not None and self.baroclinic is None:
    #             self._velocity = ifltype.Ifltype3(self.tides)

    #         elif self.tides is None and self.baroclinic is not None:
    #             self._velocity = ifltype.Ifltype4(self.baroclinic)

    #         elif self.tides is not None and self.baroclinic is not None:
    #             self._velocity = ifltype.Ifltype5(
    #                     self.tides,
    #                     self.baroclinic
    #                 )

    #         elif self.tides is None and self.baroclinic is None:
    #             raise NotImplementedError('No tides or baroclinic forcing.')

    #     return self._velocity

    # @property
    # def tides(self):
    #     if not hasattr(self, '_tides'):
    #         if self.args.iettype is None or isinstance(
    #                 self.args.iettype, (iettype.Iettype3, iettype.Iettype5)):
    #             if self.vgrid.is2D() is True:
    #                 if self.baroclinic is not None:
    #                     self._tides = Tides(
    #                         elevation=True,
    #                         velocity=True,
    #                         tidal_database=self.args.tidal_database
    #                     )
    #                 else:
    #                     self._tides = Tides(
    #                         elevation=True,
    #                         velocity=False,
    #                         tidal_database=self.args.tidal_database
    #                     )
    #             else:
    #                 self._tides = Tides(
    #                     elevation=True,
    #                     velocity=True,
    #                     tidal_database=self.args.tidal_database
    #                 )
    #         else:
    #             self._tides = None
    #     return self._tides

    # @property
    # def baroclinic(self):
    #     if not hasattr(self, '_baroclinic'):
    #         if self.baroclinic_database is not None:
    #             self._baroclinic = {
    #                 'gofs': GOFS,
    #                 'rtofs': RTOFS
    #             }[self.baroclinic_database]()
    #         else:
    #             self._baroclinic = None
    #     return self._baroclinic

    # @property
    # def temperature(self):
    #     if not hasattr(self, '_temperature'):
    #         self._temperature = None
    #     return self._temperature

    # @property
    # def salinity(self):
    #     if not hasattr(self, '_salinity'):
    #         self._salinity = None
    #     return self._salinity

    # @property
    # def baroclinic_database(self):
    #     if self.args.baroclinic_database is not None:
    #         return self.args.baroclinic_database
    #     if self.vgrid.is3D() is True:
    #         return 'gofs'

    # def baroclinic_database(self):
    #     if self.args.baroclinic_database is not None:
    #         return self.args.baroclinic_database
    #     if self.vgrid.is3D():
    #         return 'gofs'
    # @property
    # def elevation(self):
    #     return bctypes.InitialElevation


def add_bctides(subparsers):
    bctides = subparsers.add_parser('bctides')
    bctides.add_argument('hgrid')
    bctides.add_argument(
        '--run-days',
        type=lambda x: timedelta(days=float(x)),
        required=True
    )
    bctides.add_argument(
        '--start-date',
        type=lambda x: datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'),
        # help=''
    )
    bctides.add_argument('--output-directory', '-o', type=pathlib.Path,
                         required=True)
    bctides.add_argument('--vgrid')
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
    bctides.add_argument(
        '--cutoff-depth', type=float, default=50.)
    bctides.add_argument('--hgrid-crs')
    bctides.add_argument(
        "--overwrite", action="store_true",
        help="Allow overwrite of output file.")
    bctides.add_argument(
        "--log-level",
        choices=[name.lower() for name in logging._nameToLevel],
        default='warning'
    )
    _iettype = bctides.add_mutually_exclusive_group()
    _iettype.add_argument(
        '--elev-th',
        '--iettype-1',
        dest='iettype',
        type=iettype.Iettype1,
    )
    _iettype.add_argument(
        '--elev-val',
        '--iettype-2',
        dest='iettype',
        type=iettype.Iettype2
    )

    _iettype.add_argument(
        '--elev-tides',
        '--iettype-3',
        action='store_const',
        dest='iettype',
        const=iettype.Iettype3,
    )

    _iettype.add_argument(
        '--elev-subtides',
        '--elev-2d',
        '--iettype-4',
        action='store_const',
        dest='iettype',
        const=iettype.Iettype4,
    )

    _iettype.add_argument(
        '--elev-tides-subtides',
        '--tides-elev-2d',
        '--elev-2d-tides',
        '--iettype-5',
        action='store_const',
        dest='iettype',
        const=iettype.Iettype5,
    )

    _iettype.add_argument(
        '--elev-zero',
        '--iettype-_1',
        dest='iettype',
        type=iettype.Iettype_1,
    )
    _ifltype = bctides.add_mutually_exclusive_group()
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
        action='store_const',
        dest='ifltype',
        const=ifltype.Ifltype3,
    )

    _ifltype.add_argument(
        '--uv-subtides',
        '--uv-3d',
        '--uv3D',
        '--ifltype-4',
        action='store_const',
        dest='ifltype',
        const=ifltype.Ifltype4,
    )

    _ifltype.add_argument(
        '--uv-tides-subtides',
        '--uv-3d-tides',
        '--uv3d-tides',
        '--uv3D-tides',
        '--uv3D-tides',
        '--ifltype-5',
        action='store_const',
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
