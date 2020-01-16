#!/usr/bin/env python
import sys
import argparse
import matplotlib.pyplot as plt
from pyschism.mesh import Mesh


class PlotMeshCommand:

    def __init__(self, args):
        self._args = args

    def run(self):
        self._make_plot()
        self._make_triplot()
        self._make_boundary_plot()
        self._save_fig()
        self._show_fig()
        return 0

    def _make_plot(self):
        if not self.args.no_topobathy:
            self.mesh.make_plot(
                axes=self.ax,
                vmin=self.args.vmin,
                vmax=self.args.vmax,
                # levels=self.args.levels
                )

    def _make_triplot(self):
        if self.args.plot_elements:
            self.mesh.hgrid.plot_wireframe(axes=self.ax)

    def _make_boundary_plot(self):
        if self.args.plot_boundaries:
            self.mesh.hgrid.plot_boundaries(axes=self.ax)

    def _save_fig(self):
        if self.args.save_path:
            self.fig.savefig(self.args.save_path, bbox_inches='tight')

    def _show_fig(self):
        if not self.args.no_show:
            plt.show()

    @property
    def args(self):
        return self._args

    @property
    def mesh(self):
        try:
            return self.__mesh
        except AttributeError:
            self.__mesh = Mesh.open(
                self.args.hgrid,
                vgrid=self.args.vgrid,
                crs=self.args.crs,
                )
            return self.__mesh

    @property
    def fig(self):
        try:
            return self.__fig
        except AttributeError:
            self.__fig = plt.figure()
            return self.__fig

    @property
    def ax(self):
        try:
            return self.__ax
        except AttributeError:
            self.__ax = self.fig.add_subplot(111)
            return self.__ax

    @property
    def _args(self):
        return self.__args

    @_args.setter
    def _args(self, args):
        self.__args = args


def parse_args():
    parser = argparse.ArgumentParser(
            description="Program to see a quick plot of an SCHISM mesh.")
    parser.add_argument('hgrid')
    parser.add_argument('--vgrid')
    parser.add_argument('--crs')
    parser.add_argument("--vmin", type=float)
    parser.add_argument("--vmax", type=float)
    parser.add_argument("--no-topobathy", action="store_true",)
    parser.add_argument("--plot-elements", action="store_true")
    parser.add_argument("--plot-boundaries", action="store_true")
    parser.add_argument("--save-path", "--save")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-show", action='store_true')
    return parser.parse_args()


def main():
    return PlotMeshCommand(parse_args()).run()


# https://medium.com/opsops/how-to-test-if-name-main-1928367290cb
def init():
    if __name__ == "__main__":
        sys.exit(main())


init()
