from functools import lru_cache
import os
from typing import Union, List

import numpy as np  # type: ignore[import]
from pyproj import CRS  # type: ignore[import]

from pyschism.forcing.tides.bctypes import BoundaryCondition
from pyschism.forcing.tides.tides import Tides
from pyschism.forcing.atmosphere.nws import NWS
from pyschism.forcing.hydrology.base import Hydrology
# from pyschism.forcing.atmosphere.nws.nws2 import NWS2
from pyschism.mesh import Hgrid, Vgrid, Fgrid
from pyschism.enums import Coriolis


# class HgridDescriptor:

#     def __set__(self, obj, hgrid: Hgrid):
#         if not isinstance(hgrid, Hgrid):
#             raise TypeError(f'Argument hgrid must be of type {Hgrid}, '
#                             f'not {type(hgrid)}')
#         obj.__dict__['hgrid'] = hgrid

#     def __get__(self, obj, val):
#         return obj.__dict__['hgrid']


# class VgridDescriptor:

#     def __set__(self, obj, vgrid: Union[Vgrid, None]):
#         if not isinstance(vgrid, (Vgrid, type(None))):
#             raise TypeError(f'Argument vgrid must be of type {Vgrid} or None, '
#                             f'not {type(vgrid)}')
#         if vgrid is None:
#             vgrid = Vgrid()
#         obj.__dict__['vgrid'] = vgrid

#     def __get__(self, obj, val):
#         return obj.__dict__['vgrid']


# class FgridDescriptor:

#     def __set__(self, obj, fgrid: Fgrid):
#         if not isinstance(fgrid, Fgrid):
#             raise TypeError(f'Argument fgrid must be of type {Fgrid} or None, '
#                             f'not {type(fgrid)}')
#         obj.__dict__['fgrid'] = fgrid

#     def __get__(self, obj, val):
#         return obj.__dict__['fgrid']


# class OpenBoundariesDescriptor:

#     def __get__(self, obj, val):
#         open_boudaries = obj.__dict__.get('open_boudaries')
#         if open_boudaries is None:
#             open_boudaries = {}
#             for bnd in obj.hgrid.boundaries.ocean.itertuples():
#                 open_boudaries[bnd.id] = {
#                     'indexes': bnd.indexes, 'forcing': None}
#             obj.__dict__['open_boudaries'] = open_boudaries
#         return open_boudaries


# class NwsDescriptor:

#     def __set__(self, obj, nws: NWS):
#         if not isinstance(nws, NWS):
#             raise TypeError(f"Argument nws must be of type {NWS}, not "
#                             f"type {type(nws)}.")
#         obj.__dict__['nws'] = nws

#     def __get__(self, obj, val):
#         return obj.__dict__.get('nws')




# class NcorDescriptor:

#     def __set__(self, obj, ncor: Coriolis):
#         if not isinstance(ncor, Coriolis):
#             raise TypeError(f"ncor must be of type {Coriolis}, not type "
#                             f"{type(ncor)}.")
#         if ncor == Coriolis.AUTO:
#             if obj.hgrid.ics == 1:
#                 obj.sfea0 = np.median(obj.hgrid.get_y("EPSG:4326"))
#             elif obj.hgrid.ics == 2:
#                 pass  # nothing to do for ics=2
#             else:
#                 raise ValueError(
#                     f'Unknown hgrid.ics parameter {obj.hgrid.ics}')

#         elif ncor == Coriolis.CORICOEFF:
#             obj.coricoef = 0.

#         elif ncor == Coriolis.RLATITUDE:
#             obj.rlatitude = 46.

#         else:
#             raise NotImplementedError(
#                 f"Unknown value for Coriolis enum type {ncor}.")
#         obj.__dict__['ncor'] = ncor

#     def __get__(self, obj, val):
#         return obj.__dict__.get('ncor', Coriolis.AUTO)


# class Sfea0Descriptor:

#     def __set__(self, obj, sfea0: float):
#         obj.__dict__['sfea0'] = sfea0

#     def __get__(self, obj, val):
#         return obj.__dict__.get('sfea0')


# class CoricoeffDescriptor:

#     def __set__(self, obj, coricoeff: float):
#         obj.__dict__['coricoeff'] = coricoeff

#     def __get__(self, obj, val):
#         return obj.__dict__.get('coricoeff', 0.)


# class RlatitudeDescriptor:

#     def __set__(self, obj, rlatitude: float):
#         obj.__dict__['rlatitude'] = rlatitude

#     def __get__(self, obj, val):
#         return obj.__dict__.get('rlatitude', 46.)


class OpenBoundaries:

    def __init__(self, hgrid: Hgrid):
        open_boundaries = {}
        for bnd in hgrid.boundaries.ocean().itertuples():
            open_boundaries[bnd.id] = {
                'indexes': bnd.indexes, 'forcing': None}
        self._hgrid = hgrid
        self._open_boundaries = open_boundaries

    def __call__(self):
        return self._open_boundaries

    def __len__(self):
        return len(self._hgrid.boundaries.ocean())

    def __getitem__(self, id):
        return self._open_boundaries[id]

    def __iter__(self):
        for id, data in self._open_boundaries.items():
            yield id, data


class ModelDomain:

    # _hgrid = HgridDescriptor()
    # _vgrid = VgridDescriptor()
    # _fgrid = FgridDescriptor()
    # _open_boundaries = OpenBoundariesDescriptor()
    # _nws = NwsDescriptor()
    # _ncor = NcorDescriptor()
    # _ics = IcsDescriptor()
    # sfea0 = Sfea0Descriptor()
    # coricoef = CoricoeffDescriptor()
    # rlatitude = RlatitudeDescriptor()

    def __init__(self, hgrid: Hgrid, vgrid: Vgrid, fgrid: Fgrid):
        """Class representing a SCHISM computational domain.

        This class combines the horizontal grid (hgrid), vertical grid (vgrid)
        and friction/drag grids (fgrid). Additionally, this class holds
        information about forcings.
        Args:
            hgrid: :class:`pyschism.mesh.Hgrid` instance.
            vgrid: :class:`pyschism.mesh.Vgrid` instance.
            fgrid: :class:`pyschism.mesh.Fgrid` derived instance.
        """
        self._hgrid = hgrid
        self._vgrid = vgrid
        self._fgrid = fgrid
        self._open_boundaries = OpenBoundaries(hgrid)
        self._ncor = Coriolis.AUTO
        self._nws: Union[NWS, None] = None
        self._hydrology: List[Hydrology] = []

    @staticmethod
    def open(hgrid: Union[str, os.PathLike], fgrid: Union[str, os.PathLike],
             vgrid: os.PathLike = None, hgrid_crs: Union[str, CRS] = None,
             fgrid_crs: Union[str, CRS] = None):
        """Open files from disk"""
        return ModelDomain(
            Hgrid.open(hgrid, hgrid_crs),
            Vgrid.open(vgrid) if vgrid is not None else Vgrid(),
            Fgrid.open(fgrid, fgrid_crs))

    def add_boundary_condition(self, forcing: BoundaryCondition, id=None):
        if id is None:
            for i in range(len(self.open_boundaries)):
                self.add_boundary_condition(forcing, i)
        else:
            if not isinstance(forcing, BoundaryCondition):
                raise TypeError("Argument must be of type "
                                f"{BoundaryCondition} but got type "
                                f"{type(forcing)}")
            self.open_boundaries[id]['forcing'] = forcing

    def set_atmospheric_forcing(self, atmospheric_forcing: NWS):
        self._nws = atmospheric_forcing

    def set_coriolis(self, ncor: Coriolis):
        self._ncor = ncor

    def add_hydrology(self, hydrology: Hydrology):
        assert isinstance(hydrology, Hydrology), \
            f"Argument hydrology must be of type {Hydrology}, " \
            f"not type {type(hydrology)}."
        self._hydrology.append(hydrology)

    @lru_cache(maxsize=1)
    def get_active_potential_constituents(self):
        # PySCHISM allows the user to input the tidal potentials individually
        # for each boundary, however, SCHISM supports only a global
        # specification. Here, we collect all the activated tidal potentials
        # on each boundary and activate them all globally
        const = dict()
        for id in self.open_boundaries():
            forcing = self.open_boundaries[id]['forcing']
            if isinstance(forcing, Tides):
                for active in forcing.get_active_potential_constituents():
                    const[active] = True
        return tuple(const.keys())

    @lru_cache(maxsize=1)
    def get_active_forcing_constituents(self):
        # PySCHISM allows the user to input the tidal forcings individually
        # for each boundary, however, SCHISM supports only a global
        # specification. Here, we collect all the activated tidal forcings
        # on each boundary and activate them all globally
        const = dict()
        for id in self.open_boundaries():
            forcing = self.open_boundaries[id]['forcing']
            if isinstance(forcing, Tides):
                for active in forcing.get_active_forcing_constituents():
                    const[active] = True
        return tuple(const.keys())

    def make_plot(self, **kwargs):
        if self.vgrid.is3D():
            raise NotImplementedError(
                "Plotting not yet supported for 3D meshes.")
        elif self.vgrid.is2D():
            self.hgrid.make_plot(**kwargs)

    @property
    def hgrid(self):
        return self._hgrid

    @property
    def vgrid(self):
        return self._vgrid

    @property
    def fgrid(self):
        return self._fgrid

    @property
    def nws(self):
        return self._nws

    @property
    def hydrology(self):
        return self._hydrology

    @property
    def ics(self):
        if self.hgrid.crs is None:
            return None
        elif self.hgrid.crs.is_geographic:
            return 2
        else:
            return 1

    @property
    def ncor(self):
        return self._ncor

    @property
    def open_boundaries(self):
        return self._open_boundaries

    @property
    def bctides(self):
        return self._bctides

    @property
    def _hgrid(self):
        return self.__hgrid

    @_hgrid.setter
    def _hgrid(self, hgrid: Hgrid):
        assert isinstance(hgrid, Hgrid), \
            f"Argument hgrid must be of type {Hgrid}, not type {type(hgrid)}."

        self.__hgrid = hgrid

    @property
    def _ncor(self):
        return self.__ncor

    @_ncor.setter
    def _ncor(self, ncor: Coriolis):
        if not isinstance(ncor, Coriolis):
            raise TypeError(f"ncor must be of type {Coriolis}, not type "
                            f"{type(ncor)}.")
        if ncor == Coriolis.AUTO:
            if self.ics == 1:
                self.sfea0 = np.median(self.hgrid.get_y("EPSG:4326"))
            elif self.ics == 2:
                pass  # nothing to do for ics=2
            else:
                raise ValueError(
                    f'Unknown hgrid.ics parameter {self.ics}')

        elif ncor == Coriolis.CORICOEFF:
            self.coricoef = 0.

        elif ncor == Coriolis.RLATITUDE:
            self.rlatitude = 46.

        else:
            raise NotImplementedError(
                f"Unknown value for Coriolis enum type {ncor}.")
        self.__ncor = ncor
