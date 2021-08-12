from datetime import timedelta
# from enum import Enum
import pathlib
import os
from typing import Union

import f90nml

from pyschism.enums import Stratification
from pyschism.param.schism_init import GitParamTemplate


# class IbcType(Enum):
#     BAROCLINIC = 0
#     BAROTROPIC = 1

#     @classmethod
#     def _missing_(self, name):
#         raise ValueError(f'{name} is not a valid integer for ibc. '
#                          'Valid integers are 0 or 1.')


class CORE:
    """Provides error checking implementation for CORE group"""

    mandatory = ["ipre", "ibc", "ibtp", "rnday", "dt", "nspool", "ihfskip"]

    def __init__(
        self,
        ipre: int = 0,
        ibc: int = 0,
        ibtp: int = 0,
        rnday: Union[float, timedelta] = 0.0,
        dt: Union[float, timedelta] = 150.0,
        nspool: Union[int, float, timedelta] = None,
        ihfskip: Union[int, timedelta] = None,
        template: Union[str, os.PathLike, bool] = None,
    ):
        self.ipre = ipre
        self.ibc = ibc
        self.ibtp = ibtp
        self.rnday = rnday
        self.dt = dt
        self.nspool = nspool
        self.ihfskip = ihfskip
        self.template = template

    def __str__(self):
        data = []
        for key, default in self.defaults.items():
            if hasattr(self, key):
                current = getattr(self, key)
            else:
                continue
            if key in self.mandatory:
                data.append(f"  {key}={str(current)}")
            elif default != current:
                data.append(f"  {key}={str(current)}")
        data = "\n".join(data)
        return f"&CORE\n{data}\n/"

    def to_dict(self):
        output = {}
        for key, default in self.template.items():
            if hasattr(self, f"{key}"):
                current = getattr(self, f"{key}")
            else:
                current = None
            if key in self.mandatory:
                output[key] = current
            elif default != current:
                output[key] = current
        return output

    @property
    def ipre(self) -> int:
        return self._ipre

    @ipre.setter
    def ipre(self, ipre: int):
        if ipre not in [0, 1]:
            raise ValueError("Argument to ipre attribute must be 0 or 1")
        self._ipre = ipre

    @property
    def ibc(self):
        return self._ibc

    @ibc.setter
    def ibc(self, ibc: Union[Stratification, int, str]):

        if isinstance(ibc, str):
            ibc = Stratification[ibc.upper()].value

        if isinstance(ibc, Stratification):
            ibc = ibc.value

        if isinstance(ibc, int):
            if ibc not in [0, 1]:
                raise ValueError(
                    "Argument to attribute ibc must be of type "
                    f"{Stratification} or an 0, 1 integer or a string "
                    "'barotropic', 'baroclinic', not type "
                    f"{type(ibc)}."
                )

        self._ibc = ibc

    @property
    def ibtp(self):
        return self._ibtp

    @ibtp.setter
    def ibtp(self, ibtp: int):

        if ibtp not in [0, 1]:
            raise TypeError(
                "Argument to attribute ibtp must be 0 or 1, not " f"{ibtp}."
            )
        if ibtp == 1 and self.ibc == 0:
            raise ValueError("ibtp cannot be set to 1 because ibc is equal to " "zero.")
        self._ibtp = ibtp

    @property
    def rnday(self) -> Union[float, None]:
        return self._rnday

    @rnday.setter
    def rnday(self, rnday: Union[float, timedelta, None]):
        if rnday is not None:
            if not isinstance(rnday, timedelta):
                rnday = timedelta(days=float(rnday))
            self._rnday = rnday / timedelta(days=1)
        else:
            self._rnday = None

    @property
    def dt(self) -> Union[float, None]:
        return self._dt

    @dt.setter
    def dt(self, dt: Union[float, timedelta, None]):
        if dt is None:
            dt = timedelta(seconds=150.0)

        if not isinstance(dt, timedelta):
            dt = timedelta(seconds=float(dt))

        self._dt = dt.total_seconds()

    @property
    def nspool(self):
        return self._nspool

    @nspool.setter
    def nspool(self, nspool: Union[int, float, timedelta, None]):
        if nspool is None and self.rnday is not None:
            nspool = int(round(self.rnday / self.dt))
        if isinstance(nspool, timedelta):
            nspool = int(round(nspool.total_seconds() / self.dt))
        if isinstance(nspool, float):
            nspool = int(round(timedelta(hours=nspool).total_seconds() / self.dt))
        if isinstance(nspool, (int, float)):
            if nspool < 0:
                raise ValueError("nspool must be positive.")
        if nspool is not None:
            self._nspool = int(nspool)
        else:
            self._nspool = None

    @property
    def ihfskip(self):
        return self._ihfskip

    @ihfskip.setter
    def ihfskip(self, ihfskip: Union[int, timedelta, None]):

        if not isinstance(ihfskip, (int, timedelta, type(None))):
            raise TypeError("Argument ihfskip must be int, timedelta or None.")

        if ihfskip is None and self.rnday is not None:
            ihfskip = int(round(timedelta(days=self.rnday).total_seconds() / self.dt))

        if isinstance(ihfskip, timedelta):
            ihfskip = int(round(ihfskip.total_seconds() / self.dt))

        if isinstance(self.nspool, int):
            if self.nspool > 0:
                if not (ihfskip / self.nspool).is_integer():
                    raise ValueError(
                        "ihfskip/nspool must be an integer but got "
                        "ihfskip/nspool="
                        f"{ihfskip}/{self.nspool}={ihfskip/self.nspool}"
                    )

        self._ihfskip = ihfskip

    @property
    def template(self):
        return self._template

    @template.setter
    def template(self, template: Union[str, os.PathLike, None]):
        if template is True:
            template = GitParamTemplate().path
        elif template is False or template is None:
            template = None
        else:
            template = pathlib.Path(template)
        self._template = template

    @property
    def defaults(self):
        if self.template is None:
            return GitParamTemplate().core
        return f90nml.read(self.template)['core']
