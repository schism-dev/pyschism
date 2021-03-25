import matplotlib.pyplot as plt

from pyschism.mesh import Hgrid


def plot_hgrid(args):
    hgrid = Hgrid.open(args.path, crs=args.crs)
    ax = None
    if not args.no_topobathy:
        ax = hgrid.make_plot(
                vmin=args.vmin,
                vmax=args.vmax,
            )
    if args.show_elements:
        ax = hgrid.triplot(axes=ax)

    # if args.plot_boundaries:
    #     hgrid.plot_boundaries(axes=ax)

    plt.show()


def add_plot_hgrid(subparser):
    plot = subparser.add_parser('plot')
    plot.add_argument('path')
    plot.add_argument('--crs')
    plot.add_argument('--vmin')
    plot.add_argument('--vmax')
    plot.add_argument('--show-elements', action='store_true')
    plot.add_argument('--no-topobathy')
