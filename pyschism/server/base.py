from typing import List, Union

from psutil import cpu_count  # type: ignore[import]


class ServerConfig:

    def __init__(
            self,
            nproc: int = None,
            # symlink_outputs: str = None,
            schism_binary: str = None,
            mpi_launcher: str = None,
            modules: Union[str, List[str]] = None,
            modulepath=None,
            modules_init=None
    ):
        self.nproc = nproc
        # self.symlink_outputs = symlink_outputs
        self.schism_binary = schism_binary
        self.mpi_launcher = mpi_launcher
        self.modules = modules
        self.modulepath = modulepath
        self.modules_init = modules_init

    def __str__(self):
        f = []
        f.append(self.MPI_LAUNCHER)
        f.append(self.SCHISM_BINARY)
        # if self.symlink_outputs is not None:
        #     f.append(self.SYMLINK_OUTPUTS)
        return "\n".join(f)

    @property
    def mpi_launcher(self):
        return self.__mpi_launcher

    @mpi_launcher.setter
    def mpi_launcher(self, mpi_launcher):
        self.__mpi_launcher = "mpiexec -n" if mpi_launcher is None \
            else mpi_launcher

    @property
    def nproc(self):
        return self.__nproc

    @nproc.setter
    def nproc(self, nproc):
        self.__nproc = cpu_count(logical=False) if nproc is None else nproc

    @property
    def SCHISM_BINARY(self):
        return f'SCHISM_BINARY={"pschism_TVD-VL" if self.schism_binary is None else self.schism_binary}'

    # @property
    # def SYMLINK_OUTPUTS(self):
    #     f = "SYMLINK_OUTPUTS="
    #     if self.symlink_outputs is not None:
    #         f += f"{self.symlink_outputs}"
    #     return f

    @property
    def MPI_LAUNCHER(self):
        f = "MPI_LAUNCHER="
        if self.mpi_launcher is not None:
            f += f'mpiexec -n {self.nproc}'
        return f
