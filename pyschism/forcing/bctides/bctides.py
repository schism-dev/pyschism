from datetime import datetime, timedelta
from functools import cached_property, lru_cache
import pathlib
from typing import Dict, Union
import logging

from ordered_set import OrderedSet
import numpy as np

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
                    if not isinstance(ibctype, (self.bctype, type(None))):
                        raise TypeError(
                            f"Argument {ibctype} for boundary {bnd_id} must be of type {self.bctype} "
                            f" or None, not type {type(ibctype)}."
                        )
                    obj.gdf.at[np.where(obj.gdf["id"] == bnd_id)[0][0], self.bctype.__name__.lower()] = ibctype
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
    run_days = dates.RunDays()

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
        elev2D: Union[bool, str] = False,
        uv3D: Union[bool, str] = False,
        tem3D: Union[bool, str] = False,
        sal3D: Union[bool, str] = False,
        overwrite: bool = False,
        parallel_download=False,
        progress_bar=True,
    ):
        if start_date is not None:
            self.start_date = start_date
        if end_date is not None:
            self.run_days = end_date
        # self.tidal_database.write(path, )
        output_directory = pathlib.Path(output_directory)
        output_directory.mkdir(exist_ok=overwrite, parents=True)
        bctides = output_directory / "bctides.in" if bctides is True else bctides
        if bctides.exists() and not overwrite:
            raise IOError("path exists and overwrite is False")
        with open(bctides, "w") as f:
            f.write(str(self))
        # write nudge
        for bctype, tracer in {"itetype": "TEM", "isatype": "SAL"}.items():
            for boundary in self.gdf.itertuples():
                data_source = getattr(boundary, bctype)
                if data_source is not None:
                    import importlib
                    if hasattr(data_source, 'rlmax'):
                        # This generates all the nudges and writes the nudge files.
                        nudgemod = importlib.import_module('pyschism.forcing.bctides.nudge')
                        nudgeclass = getattr(nudgemod, f'{tracer}_Nudge')
                        _tracerfile = locals()[f'{tracer.lower()}3D']
                        if _tracerfile is False:
                            continue
                        elif _tracerfile is True:
                            _tracerfile = output_directory / f'{tracer}_nudge.gr3'
                        nudgeclass(self, data_source, rlmax=data_source.rlmax, rnu_day=data_source.rnu_day
                                   ).write(_tracerfile, overwrite=overwrite)
                        break

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
                return bctype.forcing_digit
            return "0"

        line = [
            f"{len(boundary.indexes)}",
            *[str(digit) for digit in map(get_focing_digit, bctypes)],
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
    def rnday(self):
        return self.run_days

    @cached_property
    def gdf(self):
        gdf = self.hgrid.boundaries.open.copy()
        gdf["iettype"] = None
        gdf["ifltype"] = None
        gdf["isatype"] = None
        gdf["itetype"] = None
        gdf["itrtype"] = None
        return gdf

    @property
    def ntip(self):
        return len(self.tides.get_active_potential_constituents())

    @property
    def nbfr(self):
        return len(self.tides.get_active_forcing_constituents())

    @cached_property
    def tides(self):
        return TidalConstituentCombiner(self.gdf)


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
        active_constituents = OrderedSet()
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
        active_constituents = OrderedSet()
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

    @lru_cache
    def get_active_constituents(self):
        return list(
            OrderedSet(
                [
                    *self.get_active_potential_constituents(),
                    *self.get_active_forcing_constituents(),
                ]
            )
        )


def ad_hoc_test():
    from datetime import datetime
    import logging

    from pyschism.mesh import Hgrid
    from pyschism.forcing.bctides import Bctides, iettype, ifltype, isatype, itetype

    # setup logging
    logging.basicConfig(
        format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
        force=True,
    )
    logging.getLogger("pyschism").setLevel(logging.DEBUG)

    startdate = datetime(2018, 8, 17)
    print(startdate)
    rnday = 61
    hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")

    # Bctides
    iet3 = iettype.Iettype3(constituents='major', database='tpxo')
    iet4 = iettype.Iettype4()
    iet5 = iettype.Iettype5(iettype3=iet3, iettype4=iet4)
    ifl3 = ifltype.Ifltype3(constituents='major', database='tpxo')
    ifl4 = ifltype.Ifltype4()
    ifl5 = ifltype.Ifltype5(ifltype3=ifl3, ifltype4=ifl4)
    isa3 = isatype.Isatype4()
    # ite3 = itetype.Itetype4()
    bctides = Bctides(hgrid, iettype={'1': iet5}, ifltype={'1': ifl5},
                      isatype=isa3,
                      itetype={'1': itetype.Itetype2(10, 1)}
                      )
    bctides.write(
            './',
            startdate,
            rnday,
            bctides=True,
            elev2D=False,
            uv3D=False,
            tem3D=False,
            sal3D=False,
            overwrite=True
            )


if __name__ == "__main__":
    ad_hoc_test()
