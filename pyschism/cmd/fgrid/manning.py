from pathlib import Path
from argparse import Namespace

from shapely.geometry import Polygon, MultiPolygon
import geopandas as gpd

from pyschism.mesh.fgrid import ManningsN
from pyschism.mesh import Hgrid


class ManningsNCli:

    def __init__(self, args: Namespace):

        if args.sub_action == 'generate':

            outdir = args.out_dir
            if not outdir:
                outdir = Path(args.hgrid).parent
            outdir = Path(outdir)

            hgrid = Hgrid.open(args.hgrid, crs=args.hgrid_crs)

            if args.constant is not None:
                mann_obj = ManningsN.constant(hgrid, args.constant)

            elif args.linear is not None:
                # NOTE: args.linear can be empty list in which case
                # we need to use defaults
                keys = ["min_value", "max_value", "min_depth", "max_depth"]
                linear_args = {
                    keys[i]: val for i, val in enumerate(args.linear)}

                mann_obj = ManningsN.linear_with_depth(
                    hgrid, **linear_args)

            for i, (shpfile, value) in enumerate(args.region_specific):
                try:
                    shpfile = Path(shpfile)
                except TypeError:
                    raise TypeError(
                        f"Invalid input type ({type(value)}) for"
                        f" region shapefile path")

                try:
                    value = float(value)
                except ValueError:
                    raise TypeError(
                        f"Invalid input type ({type(value)}) for"
                        f" constant value across region")

                gdf = gpd.read_file(shpfile)

                if gdf.crs and mann_obj.crs and not gdf.crs.equals(mann_obj.crs):
                    gdf = gdf.to_crs(mann_obj.crs)
                multipolygon = gdf.unary_union
                if not isinstance(multipolygon, (Polygon, MultiPolygon)):
                    raise ValueError(
                        f"Invalid input shape for region {i}")

                mann_obj.add_region(multipolygon, value)

            mann_obj.write(outdir/'manning.gr3', overwrite=args.overwrite)

            return

        raise NotImplementedError(f'Unhandled CLI action: {args.action}.')


def manning_subparser(subparsers):

    # creating manning
    manning = subparsers.add_parser('manning')

    sub_actions = manning.add_subparsers(dest='sub_action')

    gen_mann = sub_actions.add_parser('generate')

    gen_mann.add_argument('--hgrid', required=True)
    gen_mann.add_argument('--hgrid-crs')
    gen_mann.add_argument('--out-dir', required=False)
    gen_mann.add_argument('--overwrite', action='store_true')

    gen_mode_grp = gen_mann.add_mutually_exclusive_group(required=True)
    gen_mode_grp.add_argument('--constant', type=float)
    gen_mode_grp.add_argument('--linear', nargs='*', type=float)

    gen_mann.add_argument(
        '--region-specific', nargs=2, action='append',
        default=list())


add_manning = manning_subparser
