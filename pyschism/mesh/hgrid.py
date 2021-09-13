import pathlib
import tempfile

from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import requests

from pyschism.figures import figure, get_topobathy_kwargs
from pyschism.mesh.parsers import grd
from pyschism.mesh.base import Gr3  # , sort_edges, signed_polygon_area
from pyschism.mesh.boundaries import Boundaries


class Hgrid(Gr3):
    """
    Class that represents the unstructured planar mesh used by SCHISM.
    """

    def __init__(self, *args, boundaries=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._boundaries = Boundaries(self, boundaries)

    @staticmethod
    def open(path, crs=None):
        if str(path).endswith(".ll") and crs is None:
            crs = "epsg:4326"

        try:
            response = requests.get(path)
            response.raise_for_status()
            tmpfile = tempfile.NamedTemporaryFile()
            with open(tmpfile.name, "w") as fh:
                fh.write(response.text)
            _grd = grd.read(pathlib.Path(tmpfile.name), crs=crs)
        except Exception:
            _grd = grd.read(path, crs=crs)

        _grd["nodes"] = {
            id: (coords, -val) for id, (coords, val) in _grd["nodes"].items()
        }

        return Hgrid(**_grd)

    def to_dict(self, boundaries=True):
        _grd = super().to_dict()
        if boundaries is True:
            _grd.update(
                {
                    "nodes": {
                        id: (coord, -val)
                        for id, (coord, val) in self.nodes.to_dict().items()
                    },
                    "boundaries": self.boundaries.data,
                }
            )
            return _grd
        return _grd

    def copy(self):
        return self.__class__(**super().to_dict())

    @figure
    def make_plot(
        self,
        axes=None,
        vmin=None,
        vmax=None,
        show=False,
        title=None,
        # figsize=rcParams["figure.figsize"],
        extent=None,
        cbar_label=None,
        **kwargs
    ):
        if vmin is None:
            vmin = np.min(self.values)
        if vmax is None:
            vmax = np.max(self.values)
        kwargs.update(**get_topobathy_kwargs(self.values, vmin, vmax))
        kwargs.pop("col_val")
        levels = kwargs.pop("levels")
        if vmin != vmax:
            self.tricontourf(axes=axes, levels=levels, vmin=vmin, vmax=vmax, **kwargs)
        else:
            self.tripcolor(axes=axes, **kwargs)
        self.quadface(axes=axes, **kwargs)
        axes.axis("scaled")
        if extent is not None:
            axes.axis(extent)
        if title is not None:
            axes.set_title(title)
        mappable = ScalarMappable(cmap=kwargs["cmap"])
        mappable.set_array([])
        mappable.set_clim(vmin, vmax)
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("bottom", size="2%", pad=0.5)
        cbar = plt.colorbar(mappable, cax=cax, orientation="horizontal")
        cbar.set_ticks([vmin, vmax])
        cbar.set_ticklabels([np.around(vmin, 2), np.around(vmax, 2)])
        if cbar_label is not None:
            cbar.set_label(cbar_label)
        return axes

    @property
    def boundaries(self):
        return self._boundaries

    @property
    def ocean_boundaries(self):
        return self.boundaries.ocean

    @property
    def land_boundaries(self):
        return self.boundaries.land

    @property
    def interior_boundaries(self):
        return self.boundaries.interior
