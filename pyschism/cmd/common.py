import argparse
from datetime import datetime, timedelta
from enum import Enum
import json
import logging
import sys
from time import time

import numpy as np

from pyschism import dates
from pyschism.enums import Stratification, Sflux1Types, Sflux2Types
from pyschism.mesh import Hgrid, Vgrid, gridgr3
from pyschism.forcing import nws, bctides, hycom, source_sink
from pyschism.utils.timedeltatype import TimeDeltaType


logger = logging.getLogger(__name__)


def add_log_level_to_parser(parser):
    parser.add_argument(
        "--log-level",
        choices=["warning", "info", "debug"],
        default="warning",
    )


class HgridAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if self.const is None and not bool(
            set(sys.argv).intersection(["-h", "--help"])
        ):
            tmp_parser = argparse.ArgumentParser(add_help=False)
            tmp_parser.add_argument("--hgrid-crs")
            tmp_args, _ = tmp_parser.parse_known_args()
            logger.info(f"Opening hgrid from {values}...")
            start = time()
            hgrid = Hgrid.open(values, crs=tmp_args.hgrid_crs)
            logger.info(f"Reading hgrid took {time()-start}...")
            setattr(namespace, self.dest, hgrid)
        else:
            setattr(namespace, self.dest, self.const)


def add_hgrid_to_parser(parser, const=None):
    hgrid = parser.add_argument_group("Horizontal grid options")
    hgrid.add_argument(
        "hgrid",
        action=HgridAction,  # ,
        const=const,
        help="Path to the SCHISM hgrid file.",
    )
    hgrid.add_argument(
        "--hgrid-crs", help="Coordinate reference system string of hgrid."
    )


class VgridAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if self.const is None and not bool(
            set(sys.argv).intersection(["-h", "--help"])
        ):  # This is the init
            logger.info(f"Opening vgrid from {values}")
            start = time()
            vgrid = Vgrid.open(values)
            logger.info(f"Opening vgrid took {time()-start}")
            setattr(namespace, self.dest, vgrid)
        else:
            setattr(namespace, self.dest, self.const)


class VgridBinAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, Vgrid.from_binary(namespace.hgrid, values))


def add_vgrid_to_parser(parser, const=None):
    vgrid_group = parser.add_argument_group("Vertical grid options")
    vgrid = vgrid_group.add_mutually_exclusive_group()
    vgrid.add_argument(
        "--vgrid",
        action=VgridAction,  # if not bool(set(sys.argv).intersection(['-h', '--help'])) else None,
        default=Vgrid.default(),
        const=const,
        metavar="PATH_TO_FILE",
        help="Path to vgrid.in file to use in model.",
    )
    vgrid.add_argument(
        "--vgrid-bin",
        action=VgridBinAction,
        nargs="?",
        metavar="PATH_TO_BINARY",
        help="Path to the vgrid generation binary (must exist in PATH). "
        "Defaults to gen_vqs.",
    )


def add_fgrid_to_parser(parser):
    fgrid_group = parser.add_argument_group("Friction grid options.")
    fgrid_group.add_argument("--fgrid", help="Friction grid file.")
    fgrid_group.add_argument("--fgrid-crs")
    fgrid_group.add_argument(
        "--fgrid-type",
        choices=["auto", "manning", "drag", "rough"],
        default="auto",
        help="Can be used to specify which friction type the file is.",
    )


def add_stratification_to_parser(parser):

    stratification_group = parser.add_argument_group(
        "Stratification options."
    ).add_mutually_exclusive_group()
    stratification_group.add_argument(
        "--baroclinic",
        dest="stratification",
        action="store_const",
        const=Stratification.BAROCLINIC,
    )
    stratification_group.add_argument(
        "--barotropic",
        dest="stratification",
        action="store_const",
        const=Stratification.BAROTROPIC,
    )

    # args = parser.parse_known_args()[0]
    # stratification = (
    #     Stratification.BAROTROPIC
    #     if args.vgrid.is2D() is True
    #     else Stratification.BAROCLINIC
    # )
    # stratification_group.set_default(stratification=None)


def add_albedo_to_parser(parser):
    parser.add_argument("--albedo")


def add_diffmin_to_parser(parser):
    parser.add_argument("--diffmin")


def add_diffmax_to_parser(parser):
    parser.add_argument("--diffmax")


def add_watertype_to_parser(parser):
    parser.add_argument("--watertype")


def add_gridgr3_to_parser(parser):
    gridgr3_group = parser.add_argument_group("Auxiliary fields.")
    add_albedo_to_parser(gridgr3_group)
    add_diffmin_to_parser(gridgr3_group)
    add_diffmax_to_parser(gridgr3_group)
    add_watertype_to_parser(gridgr3_group)


def add_tvdflag_to_parser(parser):
    parser.add_argument("--tvdflag")


def add_fluxflag_to_parser(parser):
    parser.add_argument("--fluxflag")


def add_prop_to_parser(parser):
    prop_group = parser.add_argument_group("Prop fields.")
    add_tvdflag_to_parser(prop_group)
    add_fluxflag_to_parser(prop_group)


def add_tidal_constituents_to_parser(parser):

    # TODO: Get list directly from objects
    choices = ["K1", "O1", "P1", "Q1", "MM", "Mf", "M4", "MN4", "MS4", "2N2", "S1"]

    # class AllConstituentsAction(argparse.Action):
    #     def __call__(self, parser, namespace, values, option_string=None):
    #         tmp_parser = argparse.ArgumentParser(add_help=False)
    #         add_tidal_database_to_parser(tmp_parser)
    #         tmp_args = tmp_parser.parse_known_args()[0]
    #         tides = bctides.Tides(
    #             tidal_database=tmp_args.tidal_database,
    #             elevation=False,
    #             velocity=False,
    #         )
    #         tides.use_all()
    #         setattr(namespace, self.dest, tides)

    tide_group = parser.add_argument_group(
        "Tidal constituent options"
    ).add_mutually_exclusive_group()
    tide_group.add_argument(
        "--all-constituents",
        action="store_const",
        const="all",
        dest="constituents",
        help="Enables the use of all available tidal constituents. (default)",
    )
    tide_group.add_argument(
        "--major-constituents",
        action="store_const",
        const="major",
        dest="constituents",
        help="Enables the use of major group  of constituents only.",
    )
    tide_group.add_argument(
        "-c",
        "--constituents",
        action="append",
        choices=choices,
        dest="constituents",
        help="Enables specific tidal constituent to be forced in the model. "
        "(case-insensitive)",
    )
    tide_group.set_defaults(constituents="all")
    # tmp_parser = argparse.ArgumentParser(add_help=False)
    # add_tidal_database_to_parser(tmp_parser)
    # tmp_args = tmp_parser.parse_known_args()[0]
    # tides = bctides.Tides(
    #     tidal_database=tmp_args.tidal_database,
    #     elevation=False,
    #     velocity=False,
    # )
    # tides.use_all()
    # tide_group.set_defaults(tides=tides)


def add_bctides_options_to_parser(parser):
    bctides_options_group = parser.add_argument_group("Bctides options")
    bctides_options_group.add_argument("--cutoff-depth", type=float, default=50.0)
    bctides_options_group.add_argument(
        "--Z0",
        type=float,
    )


def add_tidal_database_to_parser(parser):

    # class TidalDatabaseAction(argparse.Action):
    #     def __call__(self, parser, namespace, values, option_string=None):
    #         tmp_parser = argparse.ArgumentParser(add_help=False)
    #         add_tidal_constituents_to_parser(tmp_parser)
    #         tmp_args = tmp_parser.parse_known_args()[0]
    #         print(tmp_args)
    #         exit()

    tidal_db_group = parser.add_argument_group("Tidal harmonics database options")
    tidal_db_group.add_argument(
        "--tidal-database",
        choices=["tpxo", "hamtide"],
        default="tpxo",
        # action=TidalDatabaseAction,
        help="Selects which tidal harmonics database to use. (default=tpxo)",
    )


def add_baroclinic_database_to_parser(parser):
    class BaroclinicDatabases(Enum):
        RTOFS = hycom.RTOFS
        GOFS = hycom.GOFS

    class BaroclinicDatabaseAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):

            setattr(namespace, self.dest, values)

    hycom_group = parser.add_argument_group("HYCOM options")
    # *************************************************************************************************
    hycom_group.add_argument(
        "--hycom",
        # "--baroclinic-database",
        choices=[x.name.lower() for x in BaroclinicDatabases],
        dest="baroclinic_database",
        default="gofs",
        action=BaroclinicDatabaseAction,
        help="Options for obtaining HYCOM data. Defaults to GOFS.",
    )


# class CustomBoundaryAction(argparse.Action):
#     def __call__(self, parser, namespace, values, option_string=None):
#         if isinstance(values, str):
#             if "iettype" in option_string:
#                 values = json.loads(values)
#                 for bnd_id, tval in values.items():
#                     if not isinstance(tval, int):
#                         raise TypeError(
#                             'Boundary type must be an integer, not type '
#                             f'{type(tval)}.')
#                     if tval == 3:
#                         namespace.hgrid.boundaries.set_forcing(
#                             bnd_id,
#                             iettype=iettype.Iettype3(get_tides(namespace))
#                             )
#                     else:
#                         raise ValueError
#         setattr(namespace, self.dest, values)


def get_per_boundary_help_message(varname, ibctype):
    return (
        f"Per-boundary specification option for {varname} variable. "
        "The format required is a json-style string. For example: "
        f"--{ibctype}=" + "'{\"1\": 3}' "
        f" would apply --{ibctype}-3 to boundary "
        'with id equal to "1".'
    )


def add_iettype_to_parser(parser):
    class IettypeAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            values = json.loads(values)
            for key, val in values.items():
                values.update({key: self.iettype(int(key), val)})
            setattr(namespace, self.dest, values)

        def iettype(self, key, value):
            values = [1, 2, 3, 4, 5, -1]
            assert value in values, (
                "per-boundary iettype values must be " f"one of {values}, not {value}."
            )
            raise NotImplementedError(
                "Per-boundary specification not yet " "supported."
            )

    class Iettype3Action(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            tmp_parser = argparse.ArgumentParser(add_help=False)
            add_tidal_constituents_to_parser(tmp_parser)
            add_tidal_database_to_parser(tmp_parser)
            tmp_args = tmp_parser.parse_known_args()[0]
            setattr(
                namespace,
                self.dest,
                self.const(
                    constituents=tmp_args.constituents, database=tmp_args.tidal_database
                ),
            )

    class Iettype4Action(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if not bool(set(sys.argv).intersection(["-h", "--help"])):
                if namespace.vgrid.is2D() is True:
                    raise NotImplementedError("--iettype-4 not available for 2D model.")
                tmp_parser = argparse.ArgumentParser(add_help=False)
                add_baroclinic_database_to_parser(tmp_parser)
                tmp_args = tmp_parser.parse_known_args()[0]
                setattr(namespace, self.dest, self.const(tmp_args.baroclinic_database))

    class Iettype5Action(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if not bool(set(sys.argv).intersection(["-h", "--help"])):
                if namespace.vgrid.is2D() is True:
                    raise NotImplementedError("--iettype-5 not available for 2D model.")
                tmp_parser = argparse.ArgumentParser(add_help=False)
                add_tidal_constituents_to_parser(tmp_parser)
                add_tidal_database_to_parser(tmp_parser)
                add_baroclinic_database_to_parser(tmp_parser)
                tmp_args = tmp_parser.parse_known_args()[0]
                iettype3 = bctides.iettype.Iettype3(
                    constituents=tmp_args.constituents, database=tmp_args.tidal_database
                )
                iettype4 = bctides.iettype.Iettype4(
                    data_source=tmp_args.baroclinic_database
                )
                setattr(namespace, self.dest, self.const(iettype3, iettype4))

    iettype_group = parser.add_argument_group(
        "Elevation boundary condition options"
    ).add_mutually_exclusive_group()
    iettype_group.add_argument(
        "--iettype",
        "--per-boundary-elevation",
        dest="iettype",
        metavar="JSON_STYLE_STRING",
        action=IettypeAction,
        # const=bctides.iettype.Iettype,
        help=get_per_boundary_help_message("elevation", "iettype"),
    )
    iettype_group.add_argument(
        "--iettype-1",
        "--elevation-time-history",
        dest="iettype",
        metavar="ELEVATION_TIME_HISTORY",
        type=bctides.iettype.Iettype1,
        help="Global elevation options for time history.  (Not implemented)",
    )
    iettype_group.add_argument(
        "--iettype-2",
        "--constant-elevation-value",
        dest="iettype",
        metavar="VALUE",
        type=bctides.iettype.Iettype2,
        help="Global constant elevation option.",
    )

    iettype_group.add_argument(
        "--iettype-3",
        "--harmonic-tidal-elevation",
        dest="iettype",
        nargs=0,
        const=bctides.iettype.Iettype3,
        action=Iettype3Action,
        help="Enables tidal elevation harmonic forcing.",
    )

    iettype_group.add_argument(
        "--iettype-4",
        "--subtidal-elevation",
        "--elev-2d",
        dest="iettype",
        nargs=0,
        const=bctides.iettype.Iettype4,
        action=Iettype4Action,
        help="Enables subtidal elevation forcing.",
    )

    iettype_group.add_argument(
        "--iettype-5",
        "--subtidal-and-tidal-elevation",
        dest="iettype",
        nargs=0,
        const=bctides.iettype.Iettype5,
        action=Iettype5Action,
        help="Enables the combined elevation tidal harmonic forcing and "
        "subtidal elevation forcing.",
    )

    iettype_group.add_argument(
        "--iettype--1",
        "--zero-elevation",
        dest="iettype",
        const=bctides.iettype.Iettype_1,
        action="store_const",
        help="Zero elevation at boundaries. (Not implemented)",
    )


def add_ifltype_to_parser(parser):
    class IfltypeAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            values = json.loads(values)
            for key, val in values.items():
                values.update({key: self.iettype(int(key), val)})
            setattr(namespace, self.dest, values)

        def ifltype(self, key, value):
            values = [1, 2, 3, 4, 5, -1]
            assert value in values, (
                "per-boundary ifltype values must be " f"one of {values}, not {value}."
            )
            raise NotImplementedError(
                "Per-boundary specification not yet " "supported."
            )

    class Ifltype3Action(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            tmp_parser = argparse.ArgumentParser(add_help=False)
            add_tidal_constituents_to_parser(tmp_parser)
            add_tidal_database_to_parser(tmp_parser)
            tmp_args = tmp_parser.parse_known_args()[0]
            setattr(
                namespace,
                self.dest,
                self.const(
                    constituents=tmp_args.constituents,
                    database=tmp_args.tidal_database,
                ),
            )

    class Ifltype4Action(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if not bool(set(sys.argv).intersection(["-h", "--help"])):
                if namespace.vgrid.is2D() is True:
                    raise NotImplementedError("--ifltype-4 not available for 2D model.")
                tmp_parser = argparse.ArgumentParser(add_help=False)
                add_baroclinic_database_to_parser(tmp_parser)
                tmp_args = tmp_parser.parse_known_args()[0]
                setattr(namespace, self.dest, self.const(tmp_args.baroclinic_database))

    class Ifltype5Action(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if not bool(set(sys.argv).intersection(["-h", "--help"])):
                if namespace.vgrid.is2D() is True:
                    raise NotImplementedError("--ifltype-5 not available for 2D model.")
                tmp_parser = argparse.ArgumentParser(add_help=False)
                add_tidal_constituents_to_parser(tmp_parser)
                add_tidal_database_to_parser(tmp_parser)
                add_baroclinic_database_to_parser(tmp_parser)
                tmp_args = tmp_parser.parse_known_args()[0]
                ifltype3 = bctides.ifltype.Ifltype3(
                    constituents=tmp_args.constituents, database=tmp_args.tidal_database
                )
                ifltype4 = bctides.ifltype.Ifltype4(
                    data_source=tmp_args.baroclinic_database
                )
                setattr(namespace, self.dest, self.const(ifltype3, ifltype4))

    ifltype_group = parser.add_argument_group(
        "Velocity boundary condition options"
    ).add_mutually_exclusive_group()
    # ifltype_group.add_argument(
    #     "--ifltype",
    #     '--per-boundary-velocity',
    #     dest="ifltype",
    #     metavar='JSON_STYLE_STRING',
    #     action=IfltypeAction,
    #     help=get_per_boundary_help_message('velocity', 'ifltype')
    # )
    ifltype_group.add_argument(
        "--ifltype-1",
        "--flux-time-history",
        metavar="VELOCITY_TIME_HISTORY",
        dest="ifltype",
        # type=ifltype.Ifltype1
        help="Global boundary velocity time history file. (Not implemented)",
    )
    ifltype_group.add_argument(
        "--ifltype-2",
        "--constant-flux-value",
        dest="ifltype",
        metavar="VALUE",
        # type=ifltype.Ifltype2,
        help="Global constant flow option.",
    )
    ifltype_group.add_argument(
        "--ifltype-3",
        "--harmonic-tidal-velocity",
        dest="ifltype",
        nargs=0,
        const=bctides.ifltype.Ifltype3,
        action=Ifltype3Action,
        help="Enables tidal velocity harmonic forcing.",
    )

    ifltype_group.add_argument(
        "--ifltype-4",
        "--subtidal-velocity",
        dest="ifltype",
        nargs=0,
        const=bctides.ifltype.Ifltype4,
        action=Ifltype4Action,
        help="Enables subtidal velocity forcing.",
    )

    ifltype_group.add_argument(
        "--ifltype-5",
        "--subtidal-and-tidal-velocity",
        action=Ifltype5Action,
        nargs=0,
        dest="ifltype",
        const=bctides.ifltype.Ifltype5,
        help="Enables the combined velocity tidal harmonic forcing and "
        "subtidal velocity forcing.",
    )

    ifltype_group.add_argument(
        "--ifltype--1",
        "--uv-zero",
        "--flather",
        # action='store_const',
        dest="ifltype",
        # const=ifltype.Ifltype_1,
    )


def add_itetype_to_parser(parser):
    def get_nudge_bnds(namespace):
        tmp_parser = argparse.ArgumentParser(add_help=False)
        add_baroclinic_database_to_parser(tmp_parser)
        add_nudge_to_parser(tmp_parser)
        tmp_args = tmp_parser.parse_known_args()[0]
        bnd_ids = namespace.hgrid.boundaries.open["id"].tolist()
        if tmp_args.nudge_temp is None:
            if namespace.vgrid.is2D():
                tmp_args.nudge_temp = {bnd_id: False for bnd_id in bnd_ids}
            else:
                tmp_args.nudge_temp = {bnd_id: True for bnd_id in bnd_ids}
        elif tmp_args.nudge_temp is True:
            tmp_args.nudge_temp = {bnd_id: True for bnd_id in bnd_ids}
        elif tmp_args.nudge_temp is False:
            tmp_args.nudge_temp = {bnd_id: False for bnd_id in bnd_ids}
        else:
            remaining_bounds = list(
                set(list(tmp_args.nudge_temp.keys())) - set(bnd_ids)
            )
            if np.all(list(tmp_args.nudge_temp.values())) is True:
                for bnd_id in remaining_bounds:
                    tmp_args.nudge_temp.setdefault(bnd_id, False)
            else:
                for bnd_id in remaining_bounds:
                    tmp_args.nudge_temp.setdefault(bnd_id, True)
        remaining_bounds = list(set(tmp_args.nudge_temp.keys()) - set(bnd_ids))
        if len(remaining_bounds) > 0:
            f = ["No nudging boundaries found with ids:"]
            for bnd in remaining_bounds:
                f.append(bnd)
            raise ValueError(" ".join(f))
        return tmp_args.nudge_temp

    class Itetype4Action(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            tmp_parser = argparse.ArgumentParser(add_help=False)
            add_nudge_to_parser(tmp_parser)
            tmp_args = tmp_parser.parse_known_args()[0]
            if namespace.hgrid is not None:
                nudge_bounds = get_nudge_bnds(namespace)
                itetype4 = bctides.itetype.Itetype4(
                    nudge=nudge_bounds["1"],
                    rlmax=tmp_args.rlmax_temp,
                    rnu_day=tmp_args.rnu_day_temp,
                )
                setattr(namespace, self.dest, itetype4)

    itetype_group = parser.add_argument_group(
        "Temperature boundary condition options"
    ).add_mutually_exclusive_group()
    # itetype_group.add_argument(
    #     "--itetype",
    #     dest="itetype",
    #     # nargs=1,  # 0 or more
    #     # const={},
    #     # action=CustomBoundaryAction,
    # )
    itetype_group.add_argument(
        "--itetype-1",
        "--temperature-time-history",
        dest="itetype",
        # type=itetype.Itetype1
    )
    itetype_group.add_argument(
        "--itetype-2",
        "--constant-temperature-value",
        dest="itetype",
        # type=itetype.Itetype2
    )

    itetype_group.add_argument(
        "--itetype-3",
        "--temperature-initial-condition",
        dest="itetype",
        # type=itetype.Itetype3,
    )

    itetype_group.add_argument(
        "--itetype-4",
        "--temp-3d",
        action=Itetype4Action,
        dest="itetype",
        nargs=0,
        const=bctides.itetype.Itetype4,
    )


def add_isatype_to_parser(parser):
    def get_nudge_bnds(namespace):
        tmp_parser = argparse.ArgumentParser(add_help=False)
        add_baroclinic_database_to_parser(tmp_parser)
        add_nudge_to_parser(tmp_parser)
        tmp_args = tmp_parser.parse_known_args()[0]
        bnd_ids = namespace.hgrid.boundaries.open["id"].tolist()
        if tmp_args.nudge_salt is None:
            if namespace.vgrid.is2D():
                tmp_args.nudge_salt = {bnd_id: False for bnd_id in bnd_ids}
            else:
                tmp_args.nudge_salt = {bnd_id: True for bnd_id in bnd_ids}
        elif tmp_args.nudge_salt is True:
            tmp_args.nudge_salt = {bnd_id: True for bnd_id in bnd_ids}
        elif tmp_args.nudge_salt is False:
            tmp_args.nudge_salt = {bnd_id: False for bnd_id in bnd_ids}
        else:
            remaining_bounds = list(
                set(list(tmp_args.nudge_salt.keys())) - set(bnd_ids)
            )
            if np.all(list(tmp_args.nudge_salt.values())) is True:
                for bnd_id in remaining_bounds:
                    tmp_args.nudge_salt.setdefault(bnd_id, False)
            else:
                for bnd_id in remaining_bounds:
                    tmp_args.nudge_salt.setdefault(bnd_id, True)
        remaining_bounds = list(set(tmp_args.nudge_salt.keys()) - set(bnd_ids))
        if len(remaining_bounds) > 0:
            f = ["No nudging boundaries found with ids:"]
            for bnd in remaining_bounds:
                f.append(bnd)
            raise ValueError(" ".join(f))
        return tmp_args.nudge_salt

    class Isatype4Action(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            tmp_parser = argparse.ArgumentParser(add_help=False)
            add_nudge_to_parser(tmp_parser)
            tmp_args = tmp_parser.parse_known_args()[0]
            if namespace.hgrid is not None:
                nudge_bounds = get_nudge_bnds(namespace)
                isatype4 = bctides.isatype.Isatype4(
                    nudge=nudge_bounds["1"],
                    rlmax=tmp_args.rlmax_salt,
                    rnu_day=tmp_args.rnu_day_salt,
                )
                setattr(namespace, self.dest, isatype4)

    isatype_group = parser.add_argument_group(
        "Salinity boundary condition options"
    ).add_mutually_exclusive_group()
    # isatype_group.add_argument(
    #     "--isatype",
    #     dest="isatype",
    #     # nargs=1,  # 0 or more
    #     # const={},
    #     # action=CustomBoundaryAction,
    # )
    isatype_group.add_argument(
        "--isatype-1",
        "--salt-th",
        dest="isatype",
        # type=isatype.Isatype1
    )
    isatype_group.add_argument(
        "--isatype-2",
        "--salt-val",
        dest="isatype",
        # type=isatype.Isatype2
    )

    isatype_group.add_argument(
        "--isatype-3",
        # "--salt-ic",
        dest="isatype",
        # type=isatype.Isatype3,
    )

    isatype_group.add_argument(
        "--isatype-4",
        "--salt-3d",
        action=Isatype4Action,
        dest="isatype",
        nargs=0,
        const=bctides.isatype.Isatype4,
    )


class NudgeAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        nudge = {}
        if "disable" in option_string:
            for disable_bnd in values:
                nudge.setdefault(disable_bnd, False)
            if len(nudge) == 0:
                nudge = False
        else:
            for enable_bnd in values:
                nudge.setdefault(enable_bnd, True)
            if len(nudge) == 0:
                nudge = True
        setattr(namespace, self.dest, nudge)


def add_temperature_nudge_to_parser(parser):
    temp_group = parser.add_argument_group("Nudge options for temperature")
    nudge_parser = temp_group.add_mutually_exclusive_group()
    nudge_parser.add_argument(
        "--nudge-temp",
        "--nudge-temperature",
        "--nudge-t",
        dest="nudge_temp",
        nargs="*",
        action=NudgeAction,
        help="Enables temperature nudging. If no arguments are given, the nudging "
        "is enabled for all the boundaries. If a per-boundary nudging is desired, "
        "you can pass the boundary id's of the boundaries where nudging should "
        "be enabled, for example passing --nudge-temp 1 2 4 will enable nudging "
        "for boundaries with id's 1, 2 and 4, and will disable the nudging "
        "for the remaining boundaries. Note: This option is mutually exclusive "
        "to --disable-nudge-temperature.",
    )
    nudge_parser.add_argument(
        "--disable-nudge-temp",
        "--disable-nudge-temperature",
        "--disable-nudge-t",
        dest="nudge_temp",
        nargs="*",
        action=NudgeAction,
        help="Disables temperature nudging. If no arguments are given, the nudging "
        "is disabled for all the boundaries. If per-boundary nudging disabling is desired, "
        "you can pass the boundary id's of the boundaries where nudging should "
        "be disabled, for example passing --disable-nudge-temp 1 2 4 will disable nudging "
        "for boundaries with id's 1, 2 and 4, and will enable the nudging "
        "for the remaining boundaries. Note: This option is mutually exclusive "
        "to --nudge-temperature.",
    )
    temp_group.set_defaults(nudge_temp=None)
    temp_group.add_argument("--rlmax-temp", default=1.5, type=float)
    temp_group.add_argument("--rnu_day-temp", default=0.25, type=float)


def add_salinity_nudge_to_parser(parser):
    salt_group = parser.add_argument_group("Nudge options for salinity")
    salt_group.add_argument(
        "--nudge-salt",
        "--nudge-salinity",
        "--nudge-s",
        dest="nudge_salt",
        nargs="?",
        action=NudgeAction,
    )
    salt_group.add_argument(
        "--disable-nudge-salt",
        "--disable-nudge-salinity",
        "--disable-nudge-s",
        dest="nudge_salt",
        nargs="?",
        action=NudgeAction,
    )
    salt_group.set_defaults(nudge_salt=None)
    salt_group.add_argument("--rlmax-salt", default=1.5, type=float)
    salt_group.add_argument("--rnu_day-salt", default=0.25, type=float)


def add_nudge_to_parser(parser):
    add_temperature_nudge_to_parser(parser)
    add_salinity_nudge_to_parser(parser)


def add_ibctype_to_parser(parser):
    add_iettype_to_parser(parser)
    add_ifltype_to_parser(parser)
    add_itetype_to_parser(parser)
    add_isatype_to_parser(parser)
    # add_itrtype_to_parser(parser)
    add_nudge_to_parser(parser)


def add_windrot_to_parser(parser):
    parser.add_argument(
        "--windrot",
        # action=WindrotAction,
    )


def add_sflux_options_to_parser(parser):
    parser.add_argument(
        "--file-interval-sflux-1", type=lambda x: timedelta(hours=int(x))
    )
    parser.add_argument(
        "--file-interval-sflux-2", type=lambda x: timedelta(hours=int(x))
    )


class NWS2Action(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) > 2:
            raise ValueError(
                "pyschism forecast init: error: argument --nws-2/--sflux: "
                "expected at most two arguments."
            )
        # if len(values) == 1:
        #     values.append(None)
        sflux_1 = Sflux1Types[values[0].upper()].value(
            product="gfs_0p25_1hr" if values[0].lower() == "gfs" else values[0].lower()
        )
        if len(values) == 2:
            sflux_2 = Sflux2Types[values[1].upper()].value()
        else:
            sflux_2 = None

        tmp_parser = argparse.ArgumentParser(add_help=False)
        if not bool(set(sys.argv).intersection(["-h", "--help"])):
            add_windrot_to_parser(tmp_parser)
            tmp_args = tmp_parser.parse_known_args()[0]
            windrot = (
                gridgr3.Windrot.default(namespace.hgrid)
                if tmp_args.windrot is None
                else tmp_args.windrot
            )
        else:
            windrot = None
        setattr(namespace, self.dest, nws.NWS2(sflux_1, sflux_2, windrot=windrot))


def add_sflux_to_parser(parser):
    parser.add_argument(
        "sflux",
        nargs="+",
        action=NWS2Action,
        metavar="sflux_level",
    )
    add_sflux_options_to_parser(parser)
    add_windrot_to_parser(parser)


def add_nws_to_parser(parser):

    nws_group = parser.add_argument_group("Atmospheric forcing options")
    nws_parser = nws_group.add_mutually_exclusive_group()
    nws_parser.add_argument(
        "--nws-1",
        dest="nws",
    )
    nws_parser.add_argument(
        "--nws-2",
        "--sflux",
        dest="nws",
        nargs="+",
        action=NWS2Action,
        metavar="sflux_level",
    )
    add_sflux_options_to_parser(parser)
    add_windrot_to_parser(parser)
    nws_parser.add_argument(
        "--nws-3",
        dest="nws",
    )
    nws_parser.add_argument(
        "--nws-4",
        dest="nws",
    )
    # nws_group.add_argument('--air')
    # nws_group.add_argument('--prc')
    # nws_group.add_argument('--rad')


def add_source_sink_to_parser(parser):
    class SourceSinkAction(argparse.Action):
        class SourceSinkTypes(Enum):
            NWM = source_sink.NWM

            @classmethod
            def _missing_(cls, name):
                f = [
                    f"{name} is not a valid source_sink type. Valid values are: ",
                ]
                for enum_type in cls:
                    f.append(enum_type.name.lower())
                f.append(".")
                raise ValueError("".join(f))

        # def __init__(self, option_strings, dest, nargs=None, **kwargs):
        #     if nargs is not None:
        #         raise ValueError("nargs not allowed")
        #     super().__init__(option_strings, dest, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            ss = source_sink.SourceSink()

            for source_sink_request in values:
                ss.add_forcing(
                    self.SourceSinkTypes[source_sink_request.upper()].value()
                )

            setattr(namespace, self.dest, ss)

    ss_parser = parser.add_argument_group(
        "Source/Sink forcing options"
    ).add_mutually_exclusive_group()
    ss_parser.add_argument(
        "--source-sink",
        dest="source_sink",
        nargs="+",
        action=SourceSinkAction,
        help="Add source/sink terms to model.",
    )


def add_waves_to_parser(parser):
    raise NotImplementedError("add_waves_to_parser")


def add_generic_ic_to_parser(parser, name):

    # class GenericIcAction(argparse.Action):
    #     def __call__(self, parser, namespace, values, option_string=None):
    #         raise NotImplementedError(f'{option_string} is not implemented')
    #         setattr(namespace, self.dest, values)

    ic = parser.add_argument_group(f"{name.capitalize()} initial condition options.")

    # tmp_parser = argparse.ArgumentParser(add_help=False)
    # tmp_parser.add_argument(f'--{name}_ic')
    # add_source_sink_to_parser(tmp_parser)
    # tmp_args = tmp_parser.parse_known_args()[0]
    # if tmp_args.args.elev_ic is None and tmp_parser.source_sink is not None:
    #     default = gridgr3.ElevIc.default(self.hgrid)

    # if self.args.temp_ic is None and self.bctides.itetype is not None:
    #     self.args.temp_ic = gridgr3.TempIc.from_hycom(
    #         self.hgrid,
    #         # self.bctides.itetype.data_source,
    #         self.start_date,
    #     )

    # if self.args.salt_ic is None and self.bctides.isatype is not None:
    #     self.args.salt_ic = gridgr3.SaltIc.from_hycom(
    #         self.hgrid,
    #         # self.bctides.isatype.data_source,
    #         self.start_date,
    #     )

    ic.add_argument(
        f"--{name}-ic",
        # type=gr3type.open,
        # action=GenericIcAction,
        help=f"{name.capitalize()} initial condition file.",
        # default=default,
    )
    ic.add_argument(f"--{name}-ic-crs")


def add_dates_to_parser(parser):
    def add_datetime_format_to_parser(parser):
        parser.add_argument(
            "--datetime-format", "--strptime", default=r"%Y-%m-%dT%H:%M:%S"
        )

    def add_start_date_to_parser(parser):
        class StartDateAction(argparse.Action):
            def __call__(self, parser, namespace, values, option_string=None):
                tmp_parser = argparse.ArgumentParser(add_help=False)
                add_datetime_format_to_parser(tmp_parser)
                tmp_args = tmp_parser.parse_known_args()[0]
                start_date = datetime.strptime(values, tmp_args.datetime_format)
                setattr(namespace, self.dest, start_date)

        parser.add_argument(
            "--start-date",
            action=StartDateAction,
            default=dates.nearest_cycle(),
        )

    def add_end_date_to_parser(parser):
        class EndDateAction(argparse.Action):
            def __call__(self, parser, namespace, values, option_string=None):
                if not isinstance(values, timedelta):
                    tmp_parser = argparse.ArgumentParser(add_help=False)
                    add_datetime_format_to_parser(tmp_parser)
                    tmp_args = tmp_parser.parse_known_args()[0]
                    values = datetime.strptime(values, tmp_args.datetime_format)
                setattr(namespace, self.dest, values)

        end_time_group = parser.add_mutually_exclusive_group(required=True)

        end_time_group.add_argument("--end-date", action=EndDateAction, dest="end_date")

        end_time_group.add_argument(
            "--run-days",
            action=EndDateAction,
            dest="end_date",
            type=lambda x: timedelta(days=float(x)),
        )

    add_start_date_to_parser(parser)
    add_end_date_to_parser(parser)


def add_elev_ic_to_parser(parser):
    class ElevIcAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):

            if "no" in option_string:
                setattr(namespace, self.dest, False)
                return
            # if user passes argument, it will be opened.
            elif values is not None:
                tmp_parser = argparse.ArgumentParser(add_help=False)
                tmp_parser.add_argument("--elev-ic-crs")
                tmp_parser.add_argument("--hgrid-crs")
                tmp_args = tmp_parser.parse_known_args()[0]
                tmp_args.elev_ic_crs = (
                    tmp_args.hgrid_crs
                    if tmp_args.elev_ic_crs is None and tmp_args.hgrid_crs is not None
                    else tmp_args.elev_ic_crs
                )
                elev_ic = gridgr3.ElevIc.open(values, crs=tmp_args.elev_ic_crs)
            else:
                elev_ic = gridgr3.ElevIc.default(namespace.hgrid)
            setattr(namespace, self.dest, elev_ic)

    elev_group = parser.add_argument_group("Elevation initial condition options.")
    elev_ic = elev_group.add_mutually_exclusive_group()
    elev_ic.add_argument(
        "--elev-ic",
        # type=Elev,
        action=ElevIcAction,
        help="Elevation initial condition file. If omitted, one will be created.",
        nargs="?",
        dest="elev_ic",
    )
    elev_ic.add_argument(
        "--no-elev-ic",
        # type=Elev,
        action=ElevIcAction,
        help="Elevation initial condition file. If omitted, one will be created.",
        nargs=0,
        dest="elev_ic",
    )
    elev_group.add_argument(
        "--elev-ic-crs",
        help="Elevation initial condition CRS for grid file. If omitted, it will be assumed that this file "
        "is in the same CRS as the hgrid.",
    )


def add_temp_ic_to_parser(parser):
    add_generic_ic_to_parser(parser, "temp")


def add_salt_ic_to_parser(parser):
    add_generic_ic_to_parser(parser, "salt")


def add_ic_to_parser(parser):
    add_elev_ic_to_parser(parser)
    add_temp_ic_to_parser(parser)
    add_salt_ic_to_parser(parser)


def add_surface_outputs_to_parser(parser):
    def hydro():
        """
        hydro output options
        """
        return {
            1: ("elev", "0: off; 1: on - elev. [m]"),
            2: ("air_pressure", "air pressure [Pa]"),
            3: ("air_temperature", "air temperature [C]"),
            4: ("specific_humidity", "Specific humidity [-]"),
            5: ("solar_radiation", "solar (shortwave) radiation [W/m/m]"),
            6: ("sensible_flux", "sensible flux (positive upward) [W/m/m]"),
            7: ("latent_heat", "latent heat flux (positive upward) [W/m/m]"),
            8: (
                "upward_longwave",
                "upward longwave radiation (positive upward) [W/m/m]",
            ),
            9: (
                "downward_longwave",
                "downward longwave radiation (positive downward) [W/m/m]",
            ),
            10: ("total_heat_flux", "total flux=-flsu-fllu-(radu-radd) [W/m/m]"),
            11: ("evaporation", "evaporation rate [kg/m/m/s]"),
            12: ("precipitation", "precipitation rate [kg/m/m/s]"),
            13: ("bottom_stress", "Bottom stress vector [kg/m/s^2(Pa)]"),
            14: ("wind_speed", "wind velocity vector [m/s]"),
            15: ("wind_stress", "wind stress vector [m^2/s/s]"),
            16: ("dahv", "depth-averaged vel vector [m/s]"),
            17: ("vertical_velocity", "vertical velocity [m/s]"),
            18: ("temp", "water temperature [C]"),
            19: ("salt", "water salinity [PSU]"),
            20: ("water_density", "water density [kg/m^3]"),
            21: ("diffusivity", "eddy diffusivity [m^2/s]"),
            22: ("viscosity", "eddy viscosity [m^2/s]"),
            23: ("TKE", "turbulent kinetic energy"),
            24: ("mixing-lenght", "turbulent mixing length [m]"),
            25: ("hvel", "horizontal vel vector [m/s]"),
            26: ("hvel_side", "horizontal vel vector defined @side [m/s]"),
            27: ("wvel_elem", "vertical vel. @elem [m/s]"),
            28: ("temp_elem", "T @prism centers [C]"),
            29: ("salt_elem", "S @prism centers [PSU]"),
            30: (
                "pressure_gradient",
                "Barotropic pressure gradient force vector (m.s-2) @side " "centers",
            ),
        }

    def wwm():
        """
        WWM output options
        """
        return {
            1: ("WWM_1", "sig. height (m)"),
            2: ("WWM_2", "Mean average period (sec) - TM01"),
            3: (
                "WWM_3",
                "Zero down crossing period for comparison with buoy (s) - TM02",
            ),
            4: ("WWM_4", "Average period of wave runup/overtopping - TM10"),
            5: ("WWM_5", "Mean wave number (1/m)"),
            6: ("WWM_6", "Mean wave length (m)"),
            7: (
                "WWM_9",
                "Mean average energy transport direction (degr) - MWD in NDBC?",
            ),
            8: ("WWM_10", "Mean directional spreading (degr)"),
            9: ("WWM_11", "Discrete peak period (sec) - Tp"),
            10: (
                "WWM_12",
                "Continuous peak period based on higher order moments (sec)",
            ),
            11: ("WWM_13", "Peak phase vel. (m/s)"),
            12: ("WWM_14", "Peak n-factor."),
            13: ("WWM_15", "Peak group vel. (m/s)"),
            14: ("WWM_16", "Peak wave number"),
            15: ("WWM_17", "Peak wave length"),
            16: ("WWM_18", "Peak (dominant) direction (degr)"),
            17: ("WWM_19", "Peak directional spreading"),
            18: ("WWM_20", "Discrete peak direction (radian?) "),
            19: ("WWM_21", "Orbital vel. (m/s) "),
            20: ("WWM_22", "RMS Orbital vel. (m/s) "),
            21: ("WWM_23", "Bottom excursion period (sec?) "),
            22: ("WWM_24", "Bottom wave period (sec) "),
            23: ("WWM_25", "Uresell number based on peak period "),
            24: ("WWM_26", "Friction velocity (m/s?) "),
            25: ("WWM_27", "Charnock coefficient "),
            26: ("WWM_28", "Rougness length "),
            27: ("WWM_energy_dir", "WWM_energy vector"),
            28: (
                "wave-force",
                "Wave force vector (m.s-2) computed by wwm @side centers and "
                "whole levels",
            ),
        }

    def gen():
        """
        gen output options
        """
        return {
            1: ("GEN_1", "1st tracer"),
            2: ("GEN_2", "2nd tracer"),
        }

    def age():
        """
        age output options
        """
        return {
            1: ("AGE_1", 'Indices from "1" to "ntracer_age/2"; [days]'),
            2: ("AGE_2", 'Indices from "1" to "ntracer_age/2"; [days]'),
        }

    def sed():
        """
        sed output options
        """
        return {
            1: ("SED_depth_change", "bottom depth _change_ from init. condition (m)"),
            2: ("SED_D50", " Bed median grain size in the active layer (mm)"),
            3: ("SED_bed_stress", " Bottom shear stress (Pa)"),
            4: ("SED_bed_roughness", " Bottom roughness lenghth (mm)"),
            5: ("SED_TSC", "total suspended concentration (g/L)"),
            6: ("bed_thickness", " total bed thickness @elem (m)"),
            7: ("bed_age", " total bed age over all layers @elem (sec)"),
            8: ("z0st", " Sediment transport roughness length @elem (m) (z0st_elem)"),
            9: ("z0cr", "current-ripples roughness length @elem (m) (z0cr_elem)"),
            10: ("z0sw", "sand-waves roughness length (m) @elem (z0sw_elem)"),
            11: ("z0wr", "wave-ripples roughness length @elem (m) (z0wr_elem)"),
            12: ("SED3D_1", "conc. of 1st class (one output need by each class) [g/L]"),
            13: (
                "SED_bdld_1",
                "Bedload transport rate vector (kg.m-1.s-1) for 1st tracer (one "
                "output need by tracer)",
            ),
            14: (
                "SED_bedfrac_1",
                "Bed fraction 1st tracer (one output need by each class) [-]",
            ),
            15: ("SED3D_2", "conc. of 2nd class"),
            16: ("SED_bdld_2", "Bedload transport of 2nd class"),
            17: ("SED_bedfrac_3", "Bed fraction of 2nd class"),
        }

    def eco():
        """
        EcoSim output options
        """
        return {1: ("ECO_1", "EcoSim outputs")}

    def icm():
        """
        ICM output options
        """
        return {
            1: ("ICM_Chl", "Chlorophyll"),
            2: ("ICM_pH", "PH values (ICM_PH on)"),
            3: ("ICM_PrmPrdt", "ICM primary production @elem [gC/m^3/day]"),
            4: ("ICM_DIN", "ICM totoal inorganic nitrogen (DIN) @elem [gN/m^3]"),
            5: ("ICM_PON", "ICM paticulate organic nitrogen (PON) @elem [gN/m^3]"),
            6: (
                "ICM_SED_BENDOC",
                "ICM bed sediment flux arrays: SED_BENDOC (output "
                "name:ICM_SED_BENDOC) @elem [gC/(m^2 day)]",
            ),
            7: (
                "ICM_SED_BENNH4",
                "ICM bed sediment flux arrays: SED_BENNH4 (output "
                "name:ICM_SED_BENNH4) @elem [gC/(m^2 day)]",
            ),
            8: (
                "ICM_SED_BENNO3",
                "ICM bed sediment flux arrays: SED_BENNO3 (output "
                "name:ICM_SED_BENNO3)@elem [gC/(m^2 day)]",
            ),
            9: (
                "ICM_SED_BENPO4",
                "ICM bed sediment flux arrays: SED_BENPO4 (output "
                "name:ICM_SED_BENPO4) @elem [gC/(m^2 day)]",
            ),
            10: (
                "ICM_SED_BENCOD",
                "ICM bed sediment flux arrays: SED_BENCOD (output "
                "name:ICM_SED_BENCOD) @elem [gC/(m^2 day)]",
            ),
            11: (
                "ICM_SED_BENDO",
                "ICM bed sediment flux arrays: SED_BENDO (output "
                "name:ICM_SED_BENDO) @elem [gC/(m^2 day)]",
            ),
            12: (
                "ICM_SED_BENSA",
                "ICM bed sediment flux arrays: SED_BENSA (output "
                "name:ICM_SED_BENSA) @elem [gC/(m^2 day)]",
            ),
            13: ("ICM_lfsav", "ICM SAV leaf biomass @elem [gC/m^3] (k=1 is surface)"),
            14: ("ICM_stsav", "ICM SAV stem biomass @elem [gC/m^3]"),
            15: ("ICM_rtsav", "ICM SAV root biomass @elem [gC/m^3]"),
            16: ("ICM_tlfsav", "ICM SAV total leaf biomass @elem [gC/m^2]"),
            17: ("ICM_tstsav", "ICM SAV total stem biomass @elem [gC/m^2]"),
            18: ("ICM_trtsav", "ICM SAV total root biomass @elem [gC/m^2]"),
            19: ("ICM_hcansav", "ICM SAV canopy height @elem [m]"),
            20: ("ICM_CNH4", "bottom NH4 conc"),
            21: ("ICM_CNH3", "bottom NO3 conc"),
            22: ("ICM_CPIP", "bottom P conc"),
            23: ("ICM_CPOS", "bottom Si conc"),
            24: ("ICM_CCH4", "bottom CH4 conc"),
            25: ("ICM_CSO4", "bottom SO4 conc"),
            26: ("ICM_CH2S", "bottom H2S conc"),
            27: ("ICM_SEDPON1", "bottom PON g1 conc"),
            28: ("ICM_SEDPON2", "bottom PON g2 conc"),
            29: ("ICM_SEDPON3", "bottom PON g3 conc"),
            30: ("ICM_SEDPOP1", "bottom POP g1 conc"),
            31: ("ICM_SEDPOP2", "bottom POP g2 conc"),
            32: ("ICM_SEDPOP3", "bottom POP g3 conc"),
            33: ("ICM_SEDPOC1", "bottom POC g1 conc"),
            34: ("ICM_SEDPOC2", "bottom POC g2 conc"),
            35: ("ICM_SEDPOC3", "bottom POC g3 conc"),
            36: ("ICM_EROH2S", "erosion flux H2S"),
            37: ("ICM_EROLPOC", "erosion flux LPOC"),
            38: ("ICM_ERORPOC", "erosion flux RPOC"),
            39: ("ICM_DO_consumption", "DO consumption"),
            40: ("ICM_GP1", "PB growth #1"),
            41: ("ICM_GP2", "PB growth #2"),
            42: ("ICM_GP3", "PB growth #3"),
            43: ("ICM_1", "Zoo. #1"),
            44: ("ICM_2", "Zoo. #2"),
            45: ("ICM_3", "phyto #1"),
            46: ("ICM_4", "phyto #2"),
            47: ("ICM_5", "phyto #3"),
            48: ("ICM_6", "RPOC"),
            49: ("ICM_7", "LPOC"),
            50: ("ICM_8", "DOC"),
            51: ("ICM_9", "RPON"),
            52: ("ICM_10", "LPON"),
            53: ("ICM_11", "DON"),
            54: ("ICM_12", "NH4"),
            55: ("ICM_13", "NO3"),
            56: ("ICM_14", "RPOP"),
            57: ("ICM_15", "LPOP"),
            58: ("ICM_16", "DOP"),
            59: ("ICM_17", "PO4t"),
            60: ("ICM_18", "Si- biogenic"),
            61: ("ICM_19", "available Si"),
            62: ("ICM_20", "COD: Chemical oxygen demand"),
            63: ("ICM_21", "DO"),
            64: ("ICM_22", "TIC"),
            65: ("ICM_23", "ALK"),
            66: ("ICM_24", "CA"),
            67: ("ICM_25", "CACO3"),
        }

    def cos():
        """
        CoSINE output options
        """
        return {
            1: ("COS_1", "COS_1"),
            2: ("COS_2", "COS_2"),
            3: ("COS_3", "COS_3"),
            4: ("COS_4", "COS_4"),
            5: ("COS_5", "COS_5"),
            6: ("COS_6", "COS_6"),
            7: ("COS_7", "COS_7"),
            8: ("COS_8", "COS_8"),
            9: ("COS_9", "COS_9"),
            10: ("COS_10", "COS_10"),
            11: ("COS_11", "COS_11"),
            12: ("COS_12", "COS_12"),
            13: ("COS_13", "COS_13"),
        }

    def fib():
        """
        Fecal indicating bacteria output options
        """
        return {1: ("FIB_1", "FIB_1")}

    def sed2d():
        """
        SED2D output options
        """
        return {
            1: ("SED2D_depth_change", "bottom depth _change_ from init. condition (m)"),
            2: ("SED2D_Cd", "drag coefficient used in transport formulae"),
            3: ("SED2D_cflsed", "Courant number (b.qtot.dt / h.dx)"),
            4: ("SED2D_d50", "Top layer d50 (m)"),
            5: ("SED2D_total_transport", "total transport rate vector (kg/m/s)"),
            6: ("SED2D_susp_load", "suspended tranport rate vector (kg/m/s)"),
            7: ("SED2D_bed_load", "bedload transport rate vector (kg/m/s)"),
            8: (
                "SED2D_average_transport",
                "time averaged total transport rate vector (kg/m/s)",
            ),
            9: ("SED2D_bottom_slope", "bottom slope vector (m/m); negative uphill"),
            10: ("z0eq2d", "Total roughness length @elem (m) (z0eq)"),
            11: ("z0cr2d", "current-ripples roughness length @elem (m) (z0cr)"),
            12: ("z0sw2d", "sand-waves roughness length @elem (m) (z0sw)"),
            13: ("z0wr2d", "wave-ripples roughness length @elem (m) (z0wr)"),
        }

    def mar():
        """
        marsh output options
        """
        return {
            1: ("marsh_flag", "marsh_flag"),
        }

    def ice():
        """
        ice output options
        """
        return {
            1: ("ICE_velocity", "ice advective velcoity vector [m/s]"),
            2: ("ICE_strain_rate", "strain rate @ elem [1/sec]"),
            3: ("ICE_net_heat_flux", "net heat flux to ocean (>0 warm up SST) [W/m/m]"),
            4: (
                "ICE_fresh_water_flux",
                "net fresh water flux to ocean (>0 freshens up SSS) [kg/s/m/m]",
            ),
            5: ("ICE_top_T", "ice temperature [C] at air-ice interface"),
            6: ("ICE_tracer_1", "ice volume [m]"),
            7: ("ICE_tracer_2", "ice concentration [-]"),
            8: ("ICE_tracer_3", "snow volume [m]"),
        }

    def ana():
        return {
            1: ("ANA_air_pres_grad_x", "x-component of 𝛁air_pres/ρ0 [m/s/s]"),
            2: ("ANA_air_pres_grad_y", "y-component of 𝛁air_pres/ρ0 [m/s/s]"),
            3: ("ANA_tide_pot_grad_x", "α*g*𝛁Ψ [m/s/s] (gradient of tidal potential)"),
            4: ("ANA_tide_pot_grad_y", "α*g*𝛁Ψ [m/s/s]"),
            5: ("ANA_hor_viscosity_x", "𝛁·(μ𝛁u) [m/s/s] (horizontal viscosity)"),
            6: ("ANA_hor_viscosity_y", "𝛁·(μ𝛁u) [m/s/s]"),
            7: (
                "ANA_bclinic_force_x",
                "-g/rho0* ∫_z^η dr_dx dz  [m/s/s] (b-clinic gradient)",
            ),
            8: ("ANA_bclinic_force_y", "-g/rho0* ∫_z^η dr_dy dz  [m/s/s]"),
            9: (
                "ANA_vert_viscosity_x",
                "d (ν du/dz)/dz [m/s/s] - no vegetation effects (vertical "
                "viscosity)",
            ),
            10: (
                "ANA_vert_viscosity_y",
                "d (ν dv/dz)/dz [m/s/s] - no vegetation effects",
            ),
            11: ("ANA_mom_advection_x", "(u·𝛁) u [m/s/s] (momentum advection)"),
            12: ("ANA_mom_advection_y", "(u·𝛁) u [m/s/s]"),
            13: ("ANA_Richardson", "gradient Richardson number [-]"),
            14: (
                "ANA_transport_min_dt_elem",
                "min time step at each element over all subcycles in horizontal "
                "transport solver [s]  ",
            ),
        }

    surface_outputs_group = parser.add_argument_group("Surface output options")
    surface_outputs_group.add_argument(
        "--nspool",
        help="If passing an integer, it is interpreted as timestep, "
        "if passing a float, it is interpreted as hours.",
    )

    outputs = {
        "hyd": hydro(),
        "wwm": wwm(),
        "gen": gen(),
        "age": age(),
        "sed": sed(),
        "eco": eco(),
        "icm": icm(),
        "cos": cos(),
        "fib": fib(),
        "sed2d": sed2d(),
        "mar": mar(),
        "ice": ice(),
        "ana": ana(),
    }
    for short_name, output in outputs.items():
        for id, (long_name, help_msg) in output.items():
            surface_outputs_group.add_argument(
                f"--{long_name.lower().replace('_', '-')}",
                f"-{short_name}{id}",
                help=help_msg,
                action="store_true",
            )


def add_stations_outputs_to_parser(parser):
    stations_output_group = parser.add_argument_group("Stations output options")
    stations_output_group.add_argument(
        "--stations-file",
        metavar="PATH",
    )
    stations_output_group.add_argument("--stations-file-crs")
    stations_output_group.add_argument(
        "--nspool-sta",
        help="Output interval for stations. If a station file is provided, "
        "and this option is omitted, the stations outputs will be at a six-minute "
        "iterval. If this parameter is provided as a float, the float will be interpreted in minutes, "
        "if an integer is provided, it will be interpreted as iteration step interval.",
    )
    stations_output_group.add_argument(
        "--stations-elev",
        action="store_true",
        help="Include water level output for stations.",
    )
    stations_output_group.add_argument(
        "--stations-prmsl",
        action="store_true",
        help="Include pressure at mean sea level output for stations.",
    )
    stations_output_group.add_argument(
        "--stations-uwind",
        action="store_true",
        help="Include wind u-component output for stations.",
    )
    stations_output_group.add_argument(
        "--stations-vwind",
        action="store_true",
        help="Include wind v-component output for stations.",
    )
    stations_output_group.add_argument(
        "--stations-temp",
        action="store_true",
        help="Include water temperature output for stations.",
    )
    stations_output_group.add_argument(
        "--stations-sal",
        action="store_true",
        help="Include water salinity output for stations.",
    )
    stations_output_group.add_argument(
        "--stations-uvel",
        action="store_true",
        help="Include water u-component output for stations.",
    )
    stations_output_group.add_argument(
        "--stations-vvel",
        action="store_true",
        help="Include water v-component output for stations.",
    )
    stations_output_group.add_argument(
        "--stations-wvel",
        action="store_true",
        help="Include water level output for stations.",
    )


def add_vmin_to_parser(parser):
    parser.add_argument("--vmin", type=float)


def add_vmax_to_parser(parser):
    parser.add_argument("--vmax", type=float)


def add_schism_binary_to_parser(parser):
    parser.add_argument("--schism-binary", default="pschism_TVD-VL")


def add_modules_to_parser(parser):
    parser.add_argument("--modules-init")
    parser.add_argument("--modulepath")
    parser.add_argument("--module", action="append", dest="modules")


def add_slurm_to_parser(parser):
    parser.add_argument("--account")
    parser.add_argument("--partition")
    parser.add_argument(
        "--walltime",
        type=lambda x: TimeDeltaType()(x),  #   help="In hours, float."
    )
    parser.add_argument("--ntasks", required=True, type=int)
    parser.add_argument("--filename")
    parser.add_argument("--rundir")
    parser.add_argument("--run-name")
    parser.add_argument("--mail-type")
    parser.add_argument("--mail-user")
    parser.add_argument("--log-filename")
    parser.add_argument("--path-prepend")
    parser.add_argument("--ld_library_path-prepend")
    parser.add_argument("--nodes")
    parser.add_argument("--launcher", default="srun")
    parser.add_argument("--extra-commands", action="append")


def add_workload_manager_options_to_parser(parser):
    wlm = parser.add_subparsers(dest="workload_manager")
    add_slurm_to_parser(wlm.add_parser("slurm"))

    # --- Parser design #1
    # class WorkloadManagerAction(argparse.Action):

    #     def __call__(self, parser, namespace, values, option_string=None):
    #         setattr(namespace, self.dest, values[2:])
    # options = parser.add_argument_group('Workload manager options')
    # manager = options.add_mutually_exclusive_group()
    # manager.add_argument(
    #     '--torque',
    #     '--pbs',
    #     action=WorkloadManagerAction,
    #     dest='workload_manager',
    #     nargs=0,
    # )
    # manager.add_argument(
    #     '--slurm',
    #     action=WorkloadManagerAction,
    #     dest='workload_manager',
    #     nargs=0,
    # )

    # options.add_argument('--account')
    # options.add_argument('--partition')
    # options.add_argument('--walltime', type=float, help="In hours, float.")
    # options.add_argument('--filename')

    # --- Parser design # 2: Using subparsers
    # options = parser.add_argument_group('Workload manager options')
    # add_torque_to_parser(wlm.add_parser("torque"))
    # add_pbs_to_parser(wlm.add_parser("pbs"))

    # torque = wlm.add_parser("torque")

    # server_config = init.add_subparsers(dest="server_config")
    # slurm = server_config.add_parser(
    #     'slurm', help="Add options for slurm run configuration.")
    # # slurm.add_argument('--account')
    # # slurm.add_argument('--partition')
    # # slurm.add_argument('--walltime', type=float, help="In hours, float.")
    # # slurm.add_argument('--slurm-filename')
    # # slurm.add_argument('--slurm-rundir')
    # # slurm.add_argument('--run-name')
    # # slurm.add_argument('--mail-type')
    # # slurm.add_argument('--mail-user')
    # # slurm.add_argument('--log-filename')
    # slurm.add_argument('--path-prefix')
    # # slurm.add_argument('--slurm-nodes')
    # # slurm.add_argument('--slurm-launcher', default='srun')
    # # slurm.add_argument('--extra-commands', action='append')

    # # add server options
    # parser.add_argument('--hostname')
    # parser.add_argument('--port', type=int)
    # parser.add_argument('--wdir', required=True if '--hostname' in sys.argv else False)
    # parser.add_argument('--keep-wdir', action='store_true')
    # parser.add_argument('--binaries-path', '--binaries-prefix', dest='binaries_prefix')
    # parser.add_argument('--source-script')
    # parser.add_argument('--additional-mpi-options')

    # # make nproc required when using ssh
    # # args = parser.parse_known_args()[0]
    # if '--hostname' in sys.argv:
    #     parser.add_argument('--nproc', '--ncpu', type=int, required=True)
    # else:
    #     parser.add_argument('--nproc', '--ncpu', type=int, default=-1)

    # # add resource manager option
    # manager = parser.add_mutually_exclusive_group()
    # manager.add_argument('--use-torque', action='store_true')
    # manager.add_argument('--use-pbs', action='store_true')
    # manager.add_argument('--use-slurm', action='store_true')

    # # flag some options as required when a resource manager is enabled
    # _required = '--use-torque' in sys.argv
    # _required = _required | ('--use-pbs' in sys.argv)
    # _required = _required | ('--use-slurm' in sys.argv)

    # # resource manager specific options
    # # parser.add_argument('--account', required=_required)
    # parser.add_argument('--slurm-ntasks', required=_required, type=int)
    # # parser.add_argument('--walltime', required=_required, type=float)
    # # parser.add_argument('--partition')
    # # parser.add_argument('--slurm-filename')
    # # parser.add_argument('--slurm-rundir')
    # # parser.add_argument('--run-name')
    # # parser.add_argument('--mail-type')
    # # parser.add_argument('--mail-user')
    # # parser.add_argument('--log-filename')
    # # parser.add_argument('--slurm-nodes')
    # # parser.add_argument('--slurm-launcher', default='srun')
    # # parser.add_argument('--extra-commands', action='append')
    # # parser.add_argument('--module', default=list(), action='append', dest='modules')
