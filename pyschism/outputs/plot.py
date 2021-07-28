import pathlib

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from netCDF4 import Dataset
import numpy as np
from pyugrid import UGrid

from pyschism.enums import OutputVariableUnit, OutputVariableShortName


class PlotOutputCombined:

    def __init__(self, path):
        self.nc = Dataset(pathlib.Path(path))

    def plot(self, variable, show=False, index=None):
        if index is None:
            self.animation(
                variable,
                # show=True,
                save='/home/jreniel/pyschism/examples/example_1/test.gif',
                vmin=0,
                vmax=3,
                # start_frame=200,
                # end_frame=300,
                )
        else:
            var = self.nc[variable]
            ugrid = UGrid.from_nc_dataset(self.nc)
            x = ugrid.nodes[:, 0]
            y = ugrid.nodes[:, 1]
            triangulation = Triangulation(x, y, ugrid.faces[:, :3])
            triangulation.set_mask(self.nc['wetdry_elem'][index])
            plt.tricontourf(
                triangulation, var[index, :], levels=256, cmap='jet')
            plt.gca().axis('scaled')
            if show:
                plt.show()

    def animation(
            self,
            variable,
            save=False,
            fps=3,
            start_frame=0,
            end_frame=-1,
            figsize=None,
            wireframe=False,
            cmap='jet',
            levels=256,
            show=False,
            xmin=None,
            xmax=None,
            ymin=None,
            ymax=None,
            vmin=None,
            vmax=None,

    ):

        fig = plt.figure(figsize)
        ax = fig.add_subplot(111)

        plt.tight_layout(pad=2)

        ugrid = UGrid.from_nc_dataset(self.nc)
        x = ugrid.nodes[:, 0]
        y = ugrid.nodes[:, 1]
        triangulation = Triangulation(x, y, ugrid.faces[:, :3])
        xmin = np.min(x) if xmin is None else xmin
        xmax = np.max(x) if xmax is None else xmax
        ymin = np.min(y) if ymin is None else ymin
        ymax = np.max(y) if ymax is None else ymax
        vmin = np.min(self.nc[variable]) if vmin is None else vmin
        vmax = np.max(self.nc[variable]) if vmax is None else vmax
        unit = OutputVariableUnit[OutputVariableShortName(variable).name].value

        def animate(index):
            _ax = fig.get_axes()
            ax.clear()
            if len(_ax) > 1:
                cax = _ax[1]
                cax.cla()
            else:
                cax = None

            triangulation.set_mask(self.nc['wetdry_elem'][index])

            if wireframe:
                ax.triplot(triangulation, color='k', linewidth=0.7)

            ax.tricontourf(
                triangulation,
                self.nc[variable][index, :],
                cmap=cmap,
                levels=levels,
                vmin=vmin,
                vmax=vmax
                )

            ax.set_ylim(ymin, ymax, auto=True)
            ax.set_xlim(xmin, xmax, auto=True)

            ax.set_xlabel('Longitude (°E)')
            ax.set_ylabel('Latitude (°N)')

            # ax.set_title(dates[i].strftime('%b %d, %Y %H:%M'))
            m = plt.cm.ScalarMappable(cmap=cmap)
            m.set_array(self.nc[variable][index, :])
            m.set_clim(vmin, vmax)
            cbar = fig.colorbar(m, cax=cax, format='%.1f', boundaries=np.linspace(vmin, vmax, levels))

            # cbar = fig.colorbar(_ax)
            cbar.ax.set_ylabel(f'{variable} [{unit}]', rotation=90)

        end_frame = end_frame % self.nc[variable].shape[0] \
            if end_frame < 0 else end_frame
        start_frame = start_frame % self.nc[variable].shape[0] \
            if start_frame < 0 else start_frame
        frames = range(start_frame, end_frame)
        anim = FuncAnimation(
            fig,
            animate,
            frames,
            blit=False
            )

        if save:
            anim.save(
                pathlib.Path(save),
                writer='ffmpeg',
                fps=fps
            )

        if show:
            plt.show()

        return anim
