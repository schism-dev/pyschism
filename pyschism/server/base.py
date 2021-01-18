from psutil import cpu_count  # type: ignore[import]


class ServerConfig:

    def __init__(
            self,
            nproc: int = None,
            symlink_outputs: str = None,
            schism_binary: str = None,
            mpi_launcher: str = None,
    ):
        self.nproc = nproc
        self.symlink_outputs = symlink_outputs
        self.schism_binary = schism_binary
        self.mpi_launcher = mpi_launcher

    def __str__(self):
        f = []
        f.append(self.MPI_LAUNCHER)
        f.append(self.SCHISM_BINARY)
        if self.symlink_outputs is not None:
            f.append(self.SYMLINK_OUTPUTS)
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
        f = "SCHISM_BINARY="
        if self.schism_binary is None:
            f += "pschism_TVD-VL"
        return f

    @property
    def SYMLINK_OUTPUTS(self):
        f = "SYMLINK_OUTPUTS="
        if self.symlink_outputs is not None:
            f += f"{self.symlink_outputs}"
        return f

    @property
    def MPI_LAUNCHER(self):
        f = "MPI_LAUNCHER="
        if self.mpi_launcher is not None:
            f += f'mpiexec -n {self.nproc}'
        return f
