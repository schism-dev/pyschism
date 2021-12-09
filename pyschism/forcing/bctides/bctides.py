from datetime import datetime, timedelta
import pathlib
from typing import Dict, Union
import logging

from pyschism import dates

from pyschism.mesh.vgrid import Vgrid
from pyschism.forcing.bctides import iettype, ifltype, isatype, itetype, itrtype, Tides
from pyschism.forcing.bctides.elev2d import Elev2D
from pyschism.forcing.bctides.uv3d import UV3D
from pyschism.forcing.bctides.mod3d import TEM_3D, SAL_3D

logger = logging.getLogger(__name__)


class IbctypeDescriptor:
    def __init__(self, name, bctype):
        self.name = name
        self.bctype = bctype

    def __get__(self, obj, val):
        return obj.gdf[self.name]

    def __set__(self, obj, val):
        if val is not None:
            if isinstance(val, dict):
                for bnd_id, ibctype in val.items():
                    if not isinstance(val, (self.bctype, type(None))):
                        raise TypeError(
                            f"Argument {val} must be of type {self.bctype} "
                            f" or None, not type {type(ibctype)}."
                        )
                # TODO
                raise NotImplementedError("Need to find the column name")
                idxs = obj.gdf[
                    (obj.gdf["id"] == bnd_id & obj.gdf["id"] == bnd_id)
                ].index.values
                for idx in idxs:
                    obj.gdf.at[idx, val] = obj
            else:
                if not isinstance(val, self.bctype):
                    raise TypeError(
                        f"Argument {self.name} must be of type "
                        f"{self.bctype}, not type {type(val)}."
                    )
                obj.gdf[self.name] = val


class BctidesMeta(type):
    def __new__(meta, name, bases, attrs):
        bctypes = {
            "iettype": iettype.Iettype,
            "ifltype": ifltype.Ifltype,
            "isatype": isatype.Isatype,
            "itetype": itetype.Itetype,
            "itrtype": itrtype.Itrtype,
        }
        for name, ibctype in bctypes.items():
            attrs[name] = IbctypeDescriptor(name, ibctype)
        return type(name, bases, attrs)


class Bctides(metaclass=BctidesMeta):

    start_date = dates.StartDate()
    end_date = dates.EndDate()

    def __init__(
        self,
        hgrid,
        vgrid=None,
        iettype: Union[Dict, iettype.Iettype] = None,
        ifltype: Union[Dict, ifltype.Ifltype] = None,
        isatype: Union[Dict, isatype.Isatype] = None,
        itetype: Union[Dict, itetype.Itetype] = None,
        itrtype: Union[Dict, itrtype.Itrtype] = None,
        cutoff_depth: float = 50.0,
    ):
        self.hgrid = hgrid
        self.vgrid = Vgrid.default() if vgrid is None else vgrid
        self.cutoff_depth = cutoff_depth
        self.iettype = iettype
        self.ifltype = ifltype
        self.isatype = isatype
        self.itetype = itetype
        self.itrtype = itrtype

    def __str__(self):
        f = [
            f"{str(self.start_date)}",
            f"{self.ntip} {self.cutoff_depth}",
        ]
        if self.ntip > 0:
            for constituent in self.tides.get_active_potential_constituents():
                forcing = self.tides(self.start_date, self.rnday, constituent)
                f.append(
                    " ".join(
                        [
                            f"{constituent}\n",
                            f"{forcing[0]:G}",
                            f"{forcing[1]:G}",
                            f"{forcing[2]:G}",
                            f"{forcing[3]:G}",
                            f"{forcing[4]:G}",
                        ]
                    )
                )
        f.append(f"{self.nbfr:d}")
        if self.nbfr > 0:
            for constituent in self.tides.get_active_forcing_constituents():
                forcing = self.tides(self.start_date, self.rnday, constituent)
                f.append(
                    " ".join(
                        [
                            f"{constituent}\n",
                            f"{forcing[2]:G}",
                            f"{forcing[3]:G}",
                            f"{forcing[4]:G}",
                        ]
                    )
                )
        global_constituents = self.tides.get_active_constituents()
        f.append(f"{len(self.gdf)}")
        for boundary in self.gdf.itertuples():
            f.append(self.get_forcing_string(boundary, global_constituents))
        return "\n".join(f)

    def write(
        self,
        output_directory,
        start_date: datetime = None,
        end_date: Union[datetime, timedelta] = None,
        bctides: Union[bool, str] = True,
        elev2D: Union[bool, str] = True,
        uv3D: Union[bool, str] = True,
        tem3D: Union[bool, str] = True,
        sal3D: Union[bool, str] = True,
        overwrite: bool = False,
        parallel_download=False,
        progress_bar=True,
    ):
        if start_date is not None:
            self.start_date = start_date
        if end_date is not None:
            self.end_date = end_date
        # self.tidal_database.write(path, )
        output_directory = pathlib.Path(output_directory)
        output_directory.mkdir(exist_ok=overwrite, parents=True)
        bctides = output_directory / "bctides.in" if bctides is True else bctides
        if bctides.exists() and not overwrite:
            raise IOError("path exists and overwrite is False")
        with open(bctides, "w") as f:
            f.write(str(self))
        ## write nudge
        #for bctype, tracer in {"itetype": "TEM", "isatype": "SAL"}.items():
        #    for boundary in self.gdf.itertuples():
        #        data_source = getattr(boundary, bctype)
        #        if data_source is not None:

        #            # I admit this exec is hacky.
        #            # pros: works well, it's simple, we don't need a return value
        #            # cons: might be confusing to read.
        #            # This generates all the nudges and writes the nudge files.
        #            exec(
        #                f"from pyschism.forcing.bctides.nudge import {tracer}_Nudge;"
        #                f"_tracer = output_directory / f'{tracer}_nudge.gr3' if {tracer.lower()}3D is True else {tracer};"
        #                f"_tr={tracer}_Nudge(self, data_source, rlmax=data_source.rlmax, rnu_day=data_source.rnu_day);"
        #                f'logger.info(f"Writing {tracer} nudge to file '
        #                + r'{_tracer}");'
        #                "_tr.write(_tracer, overwrite=overwrite)"
        #            )
        #            break

        def write_elev2D():
            _elev2D = output_directory / "elev2D.th.nc" if elev2D is True else elev2D
            Elev2D(self).write(
                _elev2D,
                self.start_date,
                self.rnday,
                timedelta(days=1),
                overwrite,
                progress_bar=progress_bar,
            )

        def write_uv3D():
            # write uv3D.th.nc
            _uv3D = output_directory / "uv3D.th.nc" if uv3D is True else uv3D
            UV3D(self).write(
                _uv3D,
                self.start_date,
                self.rnday,
                timedelta(days=1),
                overwrite,
                progress_bar=progress_bar,
            )

        def write_tem3D():
            # write TEM_3D.th.nc
            _tem3D = output_directory / "TEM_3D.th.nc" if tem3D is True else tem3D
            TEM_3D(self).write(
                _tem3D,
                self.start_date,
                self.rnday,
                timedelta(days=1),
                overwrite,
                progress_bar=progress_bar,
            )

        def write_sal3D():
            _sal3D = output_directory / "SAL_3D.th.nc" if sal3D is True else sal3D
            SAL_3D(self).write(
                _sal3D,
                self.start_date,
                self.rnday,
                timedelta(days=1),
                overwrite,
                progress_bar=progress_bar,
            )

        if parallel_download is True:
            from multiprocessing import Process

            jobs = [
                Process(target=f)
                for f in (write_elev2D, write_uv3D, write_tem3D, write_sal3D)
            ]
            for job in jobs:
                job.start()
            for job in jobs:
                job.join()
        else:
            if elev2D:
                write_elev2D()
            if uv3D:
                write_uv3D()
            if tem3D:
                write_tem3D()
            if sal3D:
                write_sal3D()

        # def write_tracer(tracer):
        #     tracer.write()

        # for tracer in [self.temperature, self.salinity, *self.tracers]:
        #     if tracer is not None:
        #         write_tracer(tracer)

    def get_forcing_string(self, boundary, global_constituents):

        bctypes = [
            boundary.iettype,
            boundary.ifltype,
            boundary.itetype,
            boundary.isatype,
        ]

        def get_focing_digit(bctype):
            if bctype is not None:
                # sensitive to MRO.
                return str(
                    getattr(
                        bctype, f"{bctype.__class__.__bases__[0].__name__.lower()}")
                )
            return "0"

        line = [
            f"{len(boundary.indexes)}",
            *[digit for digit in map(get_focing_digit, bctypes)],
        ]

        f = [" ".join(line)]
        for bctype in bctypes:
            if bctype is not None:
                f.append(
                    bctype.get_boundary_string(
                        self.hgrid, boundary, global_constituents=global_constituents
                    )
                )
        return "\n".join(f)

    @property
    def gdf(self):
        if not hasattr(self, "_gdf"):
            self._gdf = self.hgrid.boundaries.open.copy()
            self._gdf["iettype"] = None
            self._gdf["ifltype"] = None
            self._gdf["isatype"] = None
            self._gdf["itetype"] = None
            self._gdf["itrtype"] = None
        return self._gdf

    @property
    def ntip(self):
        return len(self.tides.get_active_potential_constituents())

    @property
    def nbfr(self):
        return len(self.tides.get_active_forcing_constituents())

    @property
    def rnday(self):
        return self.end_date - self.start_date

    @property
    def tides(self):

        if not hasattr(self, "_tides"):

            class TidalConstituentCombiner(Tides):
                def __init__(self, gdf):
                    self.gdf = gdf
                    afc = self.get_active_forcing_constituents()
                    apc = self.get_active_potential_constituents()
                    for constituent in set([*afc, *apc]):
                        self.use_constituent(
                            constituent,
                            forcing=True if constituent in afc else False,
                            potential=True if constituent in apc else False,
                        )

                def get_active_forcing_constituents(self):
                    active_constituents = set()
                    for row in self.gdf.itertuples():
                        if row.iettype is not None:
                            if row.iettype.iettype in [3, 5]:
                                [
                                    active_constituents.add(x)
                                    for x in row.iettype.tides.get_active_constituents()
                                ]
                        if row.ifltype is not None:
                            if row.ifltype.ifltype in [3, 5]:
                                [
                                    active_constituents.add(x)
                                    for x in row.ifltype.tides.get_active_constituents()
                                ]

                    return list(active_constituents)

                def get_active_potential_constituents(self):
                    active_constituents = set()
                    for row in self.gdf.itertuples():
                        if row.iettype is not None:
                            if row.iettype.iettype in [3, 5]:
                                [
                                    active_constituents.add(x)
                                    for x in row.iettype.tides.get_active_potential_constituents()
                                ]
                        if row.ifltype is not None:
                            if row.ifltype.ifltype in [3, 5]:
                                [
                                    active_constituents.add(x)
                                    for x in row.ifltype.tides.get_active_potential_constituents()
                                ]

                    return list(active_constituents)

                @property
                def constituents(self):
                    if not hasattr(self, "_constituents"):
                        self._constituents = sorted(
                            list(
                                set(
                                    [
                                        *self.get_active_potential_constituents(),
                                        *self.get_active_forcing_constituents(),
                                    ]
                                )
                            )
                        )
                    return self._constituents

            self._tides = TidalConstituentCombiner(self.gdf)

        return self._tides


class TidesCombiner(Tides):
    def __init__(self, bctides):
        self.bctides = bctides

    def get_active_potential_constituents(self):
        const = dict()
        # for row in self.bctides.gdf.itertuples():

    #         forcing = data['forcing']
    #         if isinstance(forcing, Tides):
    #             for active in forcing.get_active_potential_constituents():
    #                 const[active] = True
    #     return tuple(const.keys())

    # def get_active_forcing_constituents(self):
    #     # set active tidal forcing constituents
    #     const = dict()
    #     for id, data in self._model_domain.open_boundaries:
    #         forcing = data['forcing']
    #         if isinstance(forcing, Tides):
    #             for active in forcing.get_active_forcing_constituents():
    #                 const[active] = True
    #     return tuple(const.keys())


# ----------- draft
# @property
# def tides(self):
#     if not hasattr(self, '_tides'):
#         # get the first one you can find, since the Tides object is a
#         # singleton.
#         tides = None
#         for boundary in self.hgrid.boundaries.open.itertuples():
#             if boundary.iettype is not None:
#                 if hasattr(boundary.iettype, "tides"):
#                     tides = boundary.iettype.tides
#                     break
#                 elif boundary.ifltype is not None:
#                     if hasattr(boundary.ifltype, "tides"):
#                         tides = boundary.ifltype.tides
#                         break
#         self._tides = tides
#     return self._tides

# @property
# def tracers(self) -> List[Dict[Any, Union[bctides.itrtype.Itrtype, None]]]:
#     # if not hasattr(self, '_tracers'):
#     #     # tracers: List[Dict[Any, Union[itrtype.Itrtype, None]]] = []
#     #     boundary_data = {}
#     #     for boundary in self.hgrid.boundaries.open.itertuples():
#     #         itrtypes = boundary.itrtype
#     #         if itrtypes is None:
#     #             tracers.append({})
#     #         for tracer in boundary.itr
#     #             tracers.append()
#     #         tracer.setdefault(

#     #             )

#     # _itrtype = boundary.itrtype
#     # return self._tracers
#     # TODO: Cheating for now...
#     return []

# @property
# def ntip(self):
#     if self.tides is None:
#         return 0
#     return len(self.tides.get_active_potential_constituents())

# @property
# def nbfr(self):
#     if self.tides is None:
#         return 0
#     return self.tides.nbfr

# @property
# def Z0(self):
#     if hasattr(self.tides, '_Z0'):
#         return self.tides._Z0

# @Z0.setter
# def Z0(self, Z0):
#     self.tides.add_Z0(Z0)

# @property
# def cutoff_depth(self):
#     return self._cutoff_depth

# @cutoff_depth.setter
# def cutoff_depth(self, cutoff_depth: float):
#     self._cutoff_depth = float(cutoff_depth)

# @property
# def subtidal_database(self):
#     return self._subtidal_database

# @subtidal_database.setter
# def subtidal_database(self, subtidal_database: SubTidalDatabase):
#     if subtidal_database is not None:
#         # self._subtidal_database = Tides(subtidal_database=subtidal_database)
#     else:
#         self._subtidal_database = None

# @property
# def elevation(self):
#     return self._elevation

# @elevation.setter
# def elevation(self, elevation):
#     if elevation is not None:
#         assert isinstance(elevation, iettype.Iettype)
#     self._elevation = elevation

# @property
# def velocity(self):
#     return self._velocity

# @velocity.setter
# def velocity(self, velocity):
#     if velocity is not None:
#         assert isinstance(velocity, ifltype.Ifltype)
#     self._velocity = velocity

# @property
# def temperature(self):
#     return self._temperature

# @temperature.setter
# def temperature(self, temperature: Union[itetype.Itetype, None]):
#     if temperature is not None:
#         assert isinstance(temperature, itetype.Itetype)
#     self._temperature = temperature

# @property
# def salinity(self):
#     return self._salinity

# @salinity.setter
# def salinity(self, salinity: Union[isatype.Isatype, None]):
#     if salinity is not None:
#         assert isinstance(salinity, isatype.Isatype)
#     self._salinity = salinity


# class HgridDescriptor:

#     def __set__(self, obj, val: Hgrid):
#         if not isinstance(val, Hgrid):
#             raise TypeError(
#                 f'Argument hgrid must be of type {Hgrid}, not type '
#                 f'{type(val)}.')
#         obj.__dict__['hgrid'] = val

#     def __get__(self, obj, val):
#         return obj.__dict__['hgrid']


# class StartDateDescriptor:

#     def __set__(self, obj, val: datetime):
#         if not isinstance(val, datetime):
#             raise TypeError(
#                     f'Argument start_date must be of type {datetime}, '
#                     f'not type {type(val)}.')
#         if datetime_is_naive(val):
#             val = pytz.timezone('UTC').localize(val)
#         obj.__dict__['start_date'] = val

#     def __get__(self, obj, val):
#         return obj.__dict__['start_date']


# class RndayDescriptor:

#     def __set__(self, obj, val: Union[int, float, timedelta]):
#         if not isinstance(val, (int, float, timedelta)):
#             raise TypeError(
#                 f'Argument rnday must be of type {int}, {float} or '
#                 f'{timedelta}, not type {type(val)}.')
#         if not isinstance(val, timedelta):
#             val = timedelta(days=val)
#         obj.__dict__['rnday'] = val

#     def __get__(self, obj, val) -> timedelta:
#         return obj.__dict__['rnday']
