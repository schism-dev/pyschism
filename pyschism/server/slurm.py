from datetime import timedelta
from typing import List, Union
import uuid


from pyschism.server.base import ServerConfig


# class SlurmJobList:

#     def __init__(self):
#         self.jobs = []

#     def __get__(self, obj, val):
#         return self.jobs


# class SlurmQueue:
#     jobs = SlurmJobList()

#     def __call__(self):
#         for job in self.jobs:
#             if job["depends_on"]


#     def add(self, callable, job_id=None, depends_on=None):
#         self.jobs.append({
#             id(callable) if job_id is None else job_id: callable,
#             "depends_on": depends_on})


class SlurmConfig(ServerConfig):
    """Configuration object for SLURM-enabled servers"""

    def __init__(
        self,
        account: str = None,
        ntasks: int = None,
        partition: str = None,
        walltime: timedelta = None,
        filename: str = None,
        run_directory: str = None,
        run_name: str = None,
        mail_type: str = None,
        mail_user: str = None,
        log_filename: str = None,
        modules: List[str] = None,
        modulepath=None,
        modules_init=None,
        schism_binary: str = None,
        extra_commands: List[str] = None,
        launcher: str = None,
        nodes: int = None,
        # symlink_outputs: str = None,
        # mpi_launcher: str = None,
    ):
        """
        Instantiate a new Slurm shell script (`*.job`).

        :param account: Slurm account name
        :param ntasks: number of total tasks for Slurm to run
        :param run_name: Slurm run name
        :param partition: partition to run on
        :param walltime: time delta
        :param driver_script_filename: file path to the driver shell script
        :param run_directory: directory to run in
        :param mail_type: email type
        :param mail_user: email address
        :param log_filename: file path to output log file
        :param modules: list of file paths to modules to load
        :param path_prefix: file path to prepend to the PATH
        :param extra_commands: list of extra shell commands to insert into script
        :param launcher: command to start processes on target system (`srun`, `ibrun`, etc.)
        :param nodes: number of total nodes
        :param schism_binary: path to schism binary to use
        """
        self.account = account
        self.nproc = ntasks
        self.run_name = run_name
        self.partition = partition
        self.walltime = walltime
        self.filename = filename
        self.run_directory = run_directory
        self.mail_type = mail_type
        self.mail_user = mail_user
        self.log_filename = log_filename
        self.modules = modules
        self.schism_binary = schism_binary
        self.extra_commands = extra_commands
        # self.launcher = launcher
        self.nodes = nodes
        # self.symlink_outputs = symlink_outputs
        self.mpi_launcher = launcher
        self.modulepath = modulepath
        self.modules_init = modules_init

    def __str__(self):
        f = [
            self.MPI_LAUNCHER,
            self.SCHISM_BINARY,
            # self.SYMLINK_OUTPUTS,
            self.SLURM_NTASKS,
            # self.SYMLINK_OUTPUTS,
            self.SLURM_ACCOUNT,
            self.SLURM_RUN_NAME,
            self.SLURM_PARTITION,
            self.SLURM_WALLTIME,
            self.SLURM_RUN_DIRECTORY,
            self.SLURM_MAIL_TYPE,
            self.SLURM_MAIL_USER,
            self.SLURM_LOG_FILE,
            self.SLURM_JOB_FILE,
        ]

        # if self.modules is not None:
        #     f.extend([
        #         '',
        #         f'module load {" ".join(module for module in self._modules)}'
        #         ])

        # if self._path_prefix is not None:
        #     f += f'\n' f'PATH={self._path_prefix}:$PATH\n'

        # if self._extra_commands is not None:
        #     f += '\n'
        #     for command in self._extra_commands:
        #         f += f'{command}\n'
        # self.modules = modules
        # self.path_prefix = path_prefix
        # self.extra_commands = extra_commands
        # self.launcher = launcher
        # self.nodes = nodes

        # f.extend([
        #     f if ,
        #     f,
        #     f"SLURM_JOB_FILE:={self.log_filename}",
        #     self.slurm,
        #     ])
        # f.append(self.slurm)
        return "\n".join(f)

    @property
    def walltime(self):
        if isinstance(self._walltime, timedelta):
            hours, remainder = divmod(self._walltime, timedelta(hours=1))
            minutes, remainder = divmod(remainder, timedelta(minutes=1))
            seconds = round(remainder / timedelta(seconds=1))
            return f"{hours:02}:{minutes:02}:{seconds:02}"
        return self._walltime

    @walltime.setter
    def walltime(self, walltime: Union[timedelta, str, None]):
        assert isinstance(walltime, (timedelta, str, type(None)))
        self._walltime = walltime

    @property
    def filename(self):
        return self.__filename

    @filename.setter
    def filename(self, filename):
        if filename is None:
            filename = "slurm.job"
        self.__filename = filename

    @property
    def run_name(self):
        return self.__run_name

    @run_name.setter
    def run_name(self, run_name):
        if run_name is None:
            run_name = uuid.uuid4().hex
        self.__run_name = run_name

    @property
    def run_directory(self):
        return self.__run_directory

    @run_directory.setter
    def run_directory(self, run_directory):
        if run_directory is None:
            run_directory = "."
        self.__run_directory = run_directory

    @property
    def log_filename(self):
        return self.__log_filename

    @log_filename.setter
    def log_filename(self, log_filename):
        if log_filename is None:
            log_filename = r"slurm.log"
        self.__log_filename = log_filename

    @property
    def mpi_launcher(self):
        return self.__mpi_launcher

    @mpi_launcher.setter
    def mpi_launcher(self, mpi_launcher):
        self.__mpi_launcher = "srun" if mpi_launcher is None else mpi_launcher

    @property
    def SLURM_NTASKS(self):
        f = "SLURM_NTASKS="
        if self.nproc is not None:
            f += f"{self.nproc}"
        return f

    @property
    def SLURM_ACCOUNT(self):
        f = "SLURM_ACCOUNT="
        if self.account is not None:
            f += f"{self.account}"
        return f

    @property
    def SLURM_RUN_NAME(self):
        f = "SLURM_RUN_NAME="
        if self.run_name is not None:
            f += f"{self.run_name}"
        return f

    @property
    def SLURM_PARTITION(self):
        f = "SLURM_PARTITION="
        if self.partition is not None:
            f += f"{self.partition}"
        return f

    @property
    def SLURM_WALLTIME(self):
        f = "SLURM_WALLTIME="
        if self.walltime is not None:
            f += f"{self.walltime}"
        return f

    @property
    def SLURM_RUN_DIRECTORY(self):
        f = "SLURM_RUN_DIRECTORY="
        if self.run_directory is not None:
            f += f"{self.run_directory}"
        return f

    @property
    def SLURM_MAIL_TYPE(self):
        f = "SLURM_MAIL_TYPE="
        if self.mail_type is not None:
            f += f"{self.mail_type}"
        return f

    @property
    def SLURM_MAIL_USER(self):
        f = "SLURM_MAIL_USER="
        if self.mail_user is not None:
            f += f"{self.mail_user}"
        return f

    @property
    def SLURM_LOG_FILE(self):
        f = "SLURM_LOG_FILE="
        if self.log_filename is not None:
            f += f"{self.log_filename}"
        return f

    @property
    def SLURM_JOB_FILE(self):
        f = "SLURM_JOB_FILE="
        if self.filename is not None:
            f += f"{self.filename}"
        return f

    @property
    def MPI_LAUNCHER(self):
        f = "MPI_LAUNCHER="
        if self.mpi_launcher is not None:
            f += f"{self.mpi_launcher}"
        return f
