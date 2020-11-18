from datetime import timedelta
from typing import List
import uuid


from pyschism.server.base import ServerConfig


class SlurmConfig(ServerConfig):
    """
    Object instance of a Slurm shell script (`*.job`).
    """

    def __init__(
            self,
            account: str,
            ntasks: int,
            partition: str,
            walltime: timedelta,
            filename: str = 'slurm.job',
            run_directory: str = '.',
            run_name: str = None,
            mail_type: str = None,
            mail_user: str = None,
            log_filename: str = None,
            modules: List[str] = None,
            path_prefix: str = None,
            extra_commands: List[str] = None,
            launcher: str = 'srun',
            nodes: int = None
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
        """
        self._account = account
        self._slurm_ntasks = ntasks
        self._run_name = run_name
        self._partition = partition
        self._walltime = walltime
        self._filename = filename
        self._run_directory = run_directory
        self._mail_type = mail_type
        self._mail_user = mail_user
        self._log_filename = log_filename
        self._modules = modules
        self._path_prefix = path_prefix
        self._extra_commands = extra_commands
        self._launcher = launcher
        self._nodes = nodes

    def __str__(self):
        f = ["#!/bin/bash --login",
             f'#SBATCH -D {self._run_directory}',
             f'#SBATCH -J {self._run_name}']

        if self._account is not None:
            f.append(f'#SBATCH -A {self._account}')
        if self._mail_type is not None:
            f.append(f'#SBATCH --mail-type={self._mail_type}')
        if self._mail_user is not None:
            f.append(f'#SBATCH --mail-user={self._mail_user}')
        if self._log_filename is not None:
            f.append(f'#SBATCH --output={self._log_filename}')

        f.append(f'#SBATCH -n {self._slurm_ntasks}')
        if self._nodes is not None:
            f.append(f'#SBATCH -N {self._nodes}')

        f.extend([f'#SBATCH --time={self._walltime}',
                  f'#SBATCH --partition={self._partition}',
                  'ulimit -s unlimited',
                  'set -e'])

        if len(self._modules) > 0:
            f.append(
                f'module load {" ".join(module for module in self._modules)}')

        if self._path_prefix is not None:
            f.append(f'PATH={self._path_prefix}:$PATH')

        if self._extra_commands is not None:
            for command in self._extra_commands:
                f.append(f'{command}')

        return '\n'.join(f)

    @property
    def nprocs(self):
        return self._slurm_ntasks

    @property
    def _walltime(self):
        return self.__walltime

    @_walltime.setter
    def _walltime(self, walltime):
        hours, remainder = divmod(walltime, timedelta(hours=1))
        minutes, remainder = divmod(remainder, timedelta(minutes=1))
        seconds = round(remainder / timedelta(seconds=1))
        self.__walltime = f'{hours:02}:{minutes:02}:{seconds:02}'

    @property
    def _filename(self):
        return self.__filename

    @_filename.setter
    def _filename(self, filename):
        if filename is None:
            filename = 'slurm.job'
        self.__filename = filename

    @property
    def _run_name(self):
        return self.__run_name

    @_run_name.setter
    def _run_name(self, run_name):
        if run_name is None:
            run_name = uuid.uuid4().hex
        self.__run_name = run_name

    @property
    def _run_directory(self):
        return self.__run_directory

    @_run_directory.setter
    def _run_directory(self, run_directory):
        if run_directory is None:
            run_directory = '.'
        self.__run_directory = run_directory

    @property
    def _log_filename(self):
        return self.__log_filename

    @_log_filename.setter
    def _log_filename(self, log_filename):
        if log_filename is None:
            log_filename = "slurm.log"
        self.__log_filename = log_filename
