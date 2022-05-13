import os
import pathlib
from typing import Union

import f90nml  # type: ignore[import]

from pyschism.enums import NWSType
from pyschism.forcing.nws.base import NWS
from pyschism.forcing.nws.nws2.sflux import SfluxDataset
from pyschism.mesh import gridgr3


SFLUX_DEFAULTS = f90nml.read(pathlib.Path(__file__).parent / "sflux_inputs.txt")


class NWS2(NWS):
    def __init__(
        self,
        sflux_1: SfluxDataset = None,
        sflux_2: SfluxDataset = None,
        windrot: gridgr3.Windrot = None,
    ):
        """Loads SfluxDataset to use as NWS2 input."""

        self.sflux_1 = sflux_1
        self.sflux_2 = sflux_2
        self.windrot = windrot

    def __str__(self):
        data = []
        data = "\n".join(data)
        return f"&sflux_inputs\n{data}/\n"

    @classmethod
    def read(cls, path, sflux_1_glob="*_1.*", sflux_2_glob="*_2.*"):
        path = pathlib.Path(path)
        sflux_2 = list(path.glob(sflux_2_glob))
        return cls(
            sflux_1=SfluxDataset(list(path.glob(sflux_1_glob))),
            sflux_2=SfluxDataset(sflux_2) if len(sflux_2) > 0 else None,
        )

    def write(
        self,
        path: Union[str, os.PathLike],
        start_date=None,
        end_date=None,
        bbox=None,
        overwrite: bool = False,
        windrot: bool = True,
        air=True,
        rad=True,
        prc=True,
    ):

        # write sflux namelist
        path = pathlib.Path(path)
        if path.name != "sflux":
            path /= "sflux"
        path.mkdir(parents=True, exist_ok=overwrite)
        with open(path / "sflux_inputs.txt", "w") as f:
            f.write(str(self))
        # write sflux data
        if self.sflux_1 is not None:
            self.sflux_1.write(
                path,
                1,
                start_date=start_date,
                rnday=end_date,
                air=air,
                rad=rad,
                prc=prc,
                bbox=bbox,
                overwrite=overwrite,
            )
        if self.sflux_2 is not None:
            self.sflux_2.write(
                path,
                2,
                start_date=start_date,
                rnday=end_date,
                air=air,
                rad=rad,
                prc=prc,
                bbox=bbox,
                overwrite=overwrite,
            )
        # # write windrot data
        if windrot is not False and self.windrot is not None:
            windrot = "windrot_geo2proj.gr3" if windrot is True else windrot
            self.windrot.write(path.parent / "windrot_geo2proj.gr3", overwrite)

    @property
    def dtype(self) -> NWSType:
        """Returns the datatype of the object"""
        return NWSType(2)

    @property
    def sflux_1(self) -> SfluxDataset:
        return self._sflux_1

    @sflux_1.setter
    def sflux_1(self, sflux_1):
        if not isinstance(sflux_1, SfluxDataset):
            raise TypeError(
                f"Argument sflux_1 must be of type {SfluxDataset},"
                f" not type {type(sflux_1)}"
            )
        self._sflux_1 = sflux_1

    @property
    def sflux_2(self) -> SfluxDataset:
        return self._sflux_2

    @sflux_2.setter
    def sflux_2(self, sflux_2):
        if sflux_2 is not None:
            if not isinstance(sflux_2, SfluxDataset):
                raise TypeError(
                    f"Argument sflux_2 must be of type {SfluxDataset}, not "
                    f"type {type(sflux_2)}."
                )
        self._sflux_2 = sflux_2

    @property
    def windrot(self):
        return self._windrot

    @windrot.setter
    def windrot(self, windrot: Union[gridgr3.Windrot, None]):
        if not isinstance(windrot, (gridgr3.Windrot, type(None))):
            raise TypeError(
                f"Argument windrot must be of type {gridgr3.Windrot} or None, "
                f"not type {type(windrot)}."
            )
        self._windrot = windrot
