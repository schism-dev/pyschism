import argparse
from datetime import datetime, timedelta
import logging
import pathlib

from pyschism.cmd import common
from pyschism.forcing.bctides import Bctides


logger = logging.getLogger(__name__)


class BctidesCli:
    def __init__(self, args: argparse.Namespace):
        Bctides(
                args.hgrid,
                vgrid=args.vgrid,
                iettype=args.iettype,
                ifltype=args.ifltype,
                isatype=args.isatype,
                itetype=args.itetype,
                # itrtype=args.itrtype,
                cutoff_depth=args.cutoff_depth,
            ).write(
                args.output_directory,
                args.start_date,
                args.end_date,
                overwrite=args.overwrite,
                parallel_download=args.parallel_download
        )

    @staticmethod
    def add_subparser_action(subparsers):
        add_bctides_options_to_parser(subparsers.add_parser('bctides'))


def add_bctides_options_to_parser(parser):
    common.add_hgrid_to_parser(parser)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwrite of output files."
    )
    common.add_log_level_to_parser(parser)
    common.add_dates_to_parser(parser)
    # parser.add_argument(
    #     "--run-days",
    #     type=lambda x: timedelta(days=float(x)),
    #     required=True,
    #     help='Number of simulation days.'

    # )
    # parser.add_argument(
    #     "--start-date",
    #     type=lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S"),
    #     # help=''
    #     # action=DateAction,
    # )
    # bctides.add_argument(
    #     "--date-formatting",
    #     type=str,
    #     type=lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S"),
    #     # help=''
    # )
    parser.add_argument(
        "--output-directory",
        "-o",
        type=pathlib.Path,
        required=True
    )
    common.add_vgrid_to_parser(parser)
    common.add_tidal_constituents_to_parser(parser)
    common.add_tidal_database_to_parser(parser)
    common.add_baroclinic_database_to_parser(parser)
    common.add_bctides_options_to_parser(parser)
    common.add_ibctype_to_parser(parser)
    parser.add_argument("--parallel-download", action='store_true')


#     add_nudge_to_parser(parser)


# class NudgeAction(argparse.Action):
#     def __call__(self, parser, namespace, values, option_string=None):
        
#         if len(values) == 0:
#             values = None

#         if values is None:
#             if namespace.vgrid.is2D():
#                 values = False
#             else:
#                 values = True
#             if 'disable' in option_string:
#                 values = False

#         setattr(namespace, self.dest, values)


# def add_nudge_to_parser(parser):
#     parser.add_argument(
#         '--nudge-temp',
#         '--nudge-temperature',
#         '--nudge-t',
#         dest='nudge_temp',
#         nargs='*',
#         action=NudgeAction,
#     )
#     parser.add_argument(
#         '--disable-nudge-temp',
#         '--disable-nudge-temperature',
#         '--disable-nudge-t',
#         dest='nudge_temp',
#         nargs=0,
#         action=NudgeAction,
#     )
#     parser.set_defaults(nudge_temp=None)
#     parser.add_argument('--rlmax-temp', default=1.5, type=float)
#     parser.add_argument('--rnu_day-temp', default=0.25, type=float)
#     parser.add_argument(
#         '--nudge-salt',
#         '--nudge-salinity',
#         '--nudge-s',
#         dest='nudge_salt',
#         nargs='?',
#         action=NudgeAction,
#     )
#     parser.add_argument(
#         '--disable-nudge-salt',
#         '--disable-nudge-salinity',
#         '--disable-nudge-s',
#         dest='nudge_salt',
#         nargs='?',
#         action=NudgeAction,
#     )
#     parser.set_defaults(nudge_salt=None)
#     parser.add_argument('--rlmax-salt', default=1.5, type=float)
#     parser.add_argument('--rnu_day-salt', default=0.25, type=float)


# def get_tides(args: argparse.Namespace):
#     def get_elevation():
#         return True if args.tidal_database is not None else False

#     def get_velocity():
#         return (
#             True
#             if args.include_velocity is True or args.vgrid.is3D() is True
#             else False
#         )

#     return Tides(
#         elevation=get_elevation(),
#         velocity=get_velocity(),
#         tidal_database=args.tidal_database,
#     )


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


# class Iettype3Action(argparse.Action):

#     def __call__(self, parser, namespace, values, option_string=None):
#         for boundary in namespace.hgrid.boundaries.open.itertuples():
#             namespace.hgrid.boundaries.set_forcing(
#                 boundary.id,
#                 iettype=self.const(get_tides(namespace)))
#         setattr(namespace, self.dest, self.const)


# class Iettype4Action(argparse.Action):

#     def __call__(self, parser, namespace, values, option_string=None):
#         for boundary in namespace.hgrid.boundaries.open.itertuples():
#             namespace.hgrid.boundaries.set_forcing(
#                 boundary.id,
#                 iettype=self.const(baroclinic_databases[
#                     namespace.baroclinic_database]()))
#         setattr(namespace, self.dest, self.const)


# class Iettype5Action(argparse.Action):

#     def __call__(self, parser, namespace, values, option_string=None):
#         for boundary in namespace.hgrid.boundaries.open.itertuples():
#             namespace.hgrid.boundaries.set_forcing(
#                 boundary.id,
#                 iettype=self.const(
#                     get_tides(namespace),
#                     baroclinic_databases[namespace.baroclinic_database]()
#                     )
#             )
#         setattr(namespace, self.dest, self.const)


# class Ifltype3Action(argparse.Action):

#     def __call__(self, parser, namespace, values, option_string=None):
#         for boundary in namespace.hgrid.boundaries.open.itertuples():
#             namespace.hgrid.boundaries.set_forcing(
#                 boundary.id,
#                 ifltype=values(get_tides(namespace)))
#         setattr(namespace, self.dest, values)


# class Ifltype4Action(argparse.Action):

#     def __call__(self, parser, namespace, values, option_string=None):
#         for boundary in namespace.hgrid.boundaries.open.itertuples():
#             namespace.hgrid.boundaries.set_forcing(
#                 boundary.id,
#                 ifltype=self.const(baroclinic_databases[
#                     namespace.baroclinic_database]()))
#         setattr(namespace, self.dest, values)


# class Ifltype5Action(argparse.Action):
#     def __call__(self, parser, namespace, values, option_string=None):
#         for boundary in namespace.hgrid.boundaries.open.itertuples():
#             namespace.hgrid.boundaries.set_forcing(
#                 boundary.id,
#                 ifltype=self.const(
#                     get_tides(namespace),
#                     baroclinic_databases[namespace.baroclinic_database]()
#                     )
#             )
#         setattr(namespace, self.dest, values)


# class Itetype4Action(argparse.Action):
#     def __call__(self, parser, namespace, values, option_string=None):
#         tmp_parser = argparse.ArgumentParser(add_help=False)
#         add_nudge_to_parser(tmp_parser)
#         tmp_args, _ = tmp_parser.parse_known_args()
#         bnd_ids = namespace.hgrid.boundaries.open['id'].tolist()
#         if tmp_args.nudge_temp is None:
#             if namespace.vgrid.is2D():
#                 tmp_args.nudge_temp = []
#             else:
#                 tmp_args.nudge_temp = bnd_ids
#         elif tmp_args.nudge_temp is True:
#             tmp_args.nudge_temp = bnd_ids
#         elif tmp_args.nudge_temp is False:
#             tmp_args.nudge_temp = []
#         remaining_bounds = list(set(tmp_args.nudge_temp) - set(bnd_ids))
#         if len(remaining_bounds) > 0:
#             raise ValueError(f"No boundary with id's {remaining_bounds}.")
#         for boundary in namespace.hgrid.boundaries.open.itertuples():
#             namespace.hgrid.boundaries.set_forcing(
#                 boundary.id,
#                 itetype=self.const(
#                     baroclinic_databases[namespace.baroclinic_database](),
#                     ))
#             if boundary.id in tmp_args.nudge_temp:
#                 namespace.hgrid.boundaries.open.loc[boundary.Index].itetype.nudge = True
#                 namespace.hgrid.boundaries.open.loc[boundary.Index].itetype.rlmax = tmp_args.rlmax_temp               
#                 namespace.hgrid.boundaries.open.loc[boundary.Index].itetype.rnu_day = tmp_args.rnu_day_temp
#             else:
#                 namespace.hgrid.boundaries.open.loc[boundary.Index].itetype.nudge = False

#         setattr(namespace, self.dest, self.const)


# class Isatype4Action(argparse.Action):
#     def __call__(self, parser, namespace, values, option_string=None):
#         # logger.debug('Init isatype-4')
#         tmp_parser = argparse.ArgumentParser(add_help=False)
#         add_nudge_to_parser(tmp_parser)
#         tmp_args, _ = tmp_parser.parse_known_args()
#         bnd_ids = namespace.hgrid.boundaries.open['id'].tolist()
#         if tmp_args.nudge_salt is None:
#             if namespace.vgrid.is2D():
#                 tmp_args.nudge_salt = []
#             else:
#                 tmp_args.nudge_salt = bnd_ids
#         elif tmp_args.nudge_salt is True:
#             tmp_args.nudge_salt = bnd_ids
#         elif tmp_args.nudge_salt is False:
#             tmp_args.nudge_salt = []
#         # logger.debug('Init isatype-4')
#         remaining_bounds = list(set(tmp_args.nudge_salt) - set(bnd_ids))
#         if len(remaining_bounds) > 0:
#             raise ValueError(f"No boundary with id's {remaining_bounds}.")
#         for boundary in namespace.hgrid.boundaries.open.itertuples():
#             namespace.hgrid.boundaries.set_forcing(
#                 boundary.id,
#                 isatype=self.const(baroclinic_databases[
#                     namespace.baroclinic_database]()))
#             if boundary.id in tmp_args.nudge_salt:
#                 namespace.hgrid.boundaries.open.loc[boundary.Index].isatype.nudge = True
#                 namespace.hgrid.boundaries.open.loc[boundary.Index].isatype.rlmax = tmp_args.rlmax_salt               
#                 namespace.hgrid.boundaries.open.loc[boundary.Index].isatype.rnu_day = tmp_args.rnu_day_salt
#             else:
#                 namespace.hgrid.boundaries.open.loc[boundary.Index].isatype.nudge = False
#         setattr(namespace, self.dest, self.const)
