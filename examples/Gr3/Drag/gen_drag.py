from pyschism.mesh.hgrid import Hgrid
from pyschism.mesh.fgrid import DragCoefficient
import argparse
import logging

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.getLogger("pyschism").setLevel(logging.DEBUG)


def generate_drag_coefficient(
    hgrid_file: str,
    depth1: float,
    depth2: float,
    bfric_river: float,
    bfric_land: float,
    regions: list[str] = None,
    values: list[float] = None,
    flags: list[int] = None,
    output_file: str = "drag.gr3",
) -> None:
    """
    Generate drag coefficient file for SCHISM model.
    
    Args:
        hgrid_file: Path to the hgrid.gr3 file
        depth1: First depth threshold
        depth2: Second depth threshold
        bfric_river: Bottom friction coefficient for rivers
        bfric_land: Bottom friction coefficient for land
        regions: List of region file paths
        values: List of values to apply to regions
        flags: List of flags for region modifications (0: reset, 1: add)
        output_file: Output file path for drag coefficients
    """
    hgrid = Hgrid.open(hgrid_file, crs="epsg:4326")
    fgrid = DragCoefficient.linear_with_depth(hgrid, depth1, depth2,
                                              bfric_river, bfric_land)

    if regions and values and flags:
        for reg, value, flag in zip(regions, values, flags):
            fgrid.modify_by_region(hgrid, reg, value, depth1, flag)

    fgrid.write(output_file, overwrite=True)


def list_of_strings(arg: str) -> list[str]:
    """Convert comma-separated string to list of strings."""
    return arg.split(",")


def list_of_floats(arg: str) -> list[float]:
    """Convert comma-separated string to list of floats."""
    return [float(x) for x in arg.split(",")]


def list_of_ints(arg: str) -> list[int]:
    """Convert comma-separated string to list of integers."""
    return [int(x) for x in arg.split(",")]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate drag coefficient file for SCHISM model")

    parser.add_argument("hgrid_file", type=str, help="Path to hgrid.gr3 file")
    parser.add_argument("--depth1",
                        type=float,
                        default=-1.0,
                        help="First depth threshold")
    parser.add_argument("--depth2",
                        type=float,
                        default=-3.0,
                        help="Second depth threshold")
    parser.add_argument("--bfric-river",
                        type=float,
                        default=0.0025,
                        help="Bottom friction coefficient for rivers")
    parser.add_argument("--bfric-land",
                        type=float,
                        default=0.025,
                        help="Bottom friction coefficient for land")
    parser.add_argument("--regions",
                        type=list_of_strings,
                        help="Comma-separated list of region file paths")
    parser.add_argument(
        "--values",
        type=list_of_floats,
        help="Comma-separated list of values to apply to regions")
    parser.add_argument(
        "--flags",
        type=list_of_ints,
        help=
        "Comma-separated list of flags for region modifications (0: reset, 1: add)"
    )
    parser.add_argument("--output",
                        type=str,
                        default="drag.gr3",
                        help="Output file path for drag coefficients")

    args = parser.parse_args()

    generate_drag_coefficient(
        hgrid_file=args.hgrid_file,
        depth1=args.depth1,
        depth2=args.depth2,
        bfric_river=args.bfric_river,
        bfric_land=args.bfric_land,
        regions=args.regions,
        values=args.values,
        flags=args.flags,
        output_file=args.output,
    )
