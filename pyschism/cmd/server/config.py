import logging
import paramiko
import uuid
import tempfile
import pathlib
import os
from stat import S_ISDIR
# import getpass


class ServerConfig:
    """
    This class is used for configuring the server
    """

    def __init__(
        self,
        hostname,
        nprocs,
        wdir=None,
        binaries_prefix=None,
        port=22,
        username=None,
        password=None,
        pkey=None,
        writer_procs=None,
        source_script=None,
        additional_mpi_options=None,
        keep_wdir=False
    ):
        self._hostname = hostname
        self._nprocs = nprocs
        self._wdir = wdir
        self._binaries_prefix = binaries_prefix
        self._port = port
        self._username = username
        self._password = password
        self._pkey = pkey
        self._writer_procs = writer_procs
        self._source_script = source_script
        self._additional_mpi_options = additional_mpi_options
        self._keep_wdir = keep_wdir

    def run(
        self,
        driver,
        outdir,
        overwrite=False,
        coldstart=True,
        hotstart=True
    ):
        """
        puts driver outputs on outdir using a remote server for compute
        """
        outdir = pathlib.Path(outdir)
        if not outdir.exists():
            msg = f"{outdir} exists and overwrite is not enabled."
            raise IOError(msg)
        self.ssh.exec_command(f'mkdir -p {self.wdir}')
        self._deploy_files_to_server(driver)
        self._run_coldstart(driver)
        self._cleanup_rundir('coldstart')
        self._run_hotstart(driver)
        self._cleanup_rundir('hotstart')
        self._retrieve_files(outdir)
        if not self.keep_wdir:
            self.ssh.exec_command(f'rm -rf {self.wdir}')

    def _deploy_files_to_server(self, driver):
        outdir = tempfile.TemporaryDirectory()
        driver.dump(outdir.name)
        for item in pathlib.Path(outdir.name).glob('**/*'):
            self.sftp.put(item.absolute(), f'{self.wdir}/{item.name}')

    def _run_coldstart(self, driver):
        self._run_adcprep_command('coldstart')
        self._run_padcirc_command('coldstart', driver)

    def _run_hotstart(self, driver):
        self._run_adcprep_command('hotstart')
        self._run_padcirc_command('hotstart', driver)

    def _run_adcprep_command(self, runtype):
        cmd = f'rm -rf {self.wdir}/{runtype}; '
        cmd += f'mkdir -p {self.wdir}/{runtype}; '
        cmd += f"cd {self.wdir}/{runtype}; "
        cmd += f'ln -sf ../fort.14; '
        cmd += f'ln -sf ../fort.13; '
        cmd += f'ln -sf ../fort.15.{runtype} ./fort.15; '
        if runtype == 'hotstart':
            cmd += f'ln -sf ../coldstart/fort.67.nc ./fort.67.nc; '
        # if self.libraries_path:
        #     cmd += f'export LD_LIBRARY_PATH={self.libraries_path}:'
        #     cmd += "$LD_LIBRARY_PATH && "
        if self.source_script:
            cmd += f'source {self.source_script} && '
        cmd += f'{self.adcprep_binary} --np {self.nprocs} --partmesh && '
        cmd += f'{self.adcprep_binary} --np {self.nprocs} --prepall'
        stdin, stdout, stderr = self.ssh.exec_command(cmd)
        while True:
            out = stdout.readline()
            if not out:
                break
            print(out, end='')
        lines = stderr.readlines()
        if len(lines) > 0:
            msg = "\n"
            msg += "".join(lines)
            raise Exception(msg)

    def _run_padcirc_command(self, runtype, driver):
        cmd = ''
        if self.nprocs > 1:
            # if self.libraries_path:
            #     cmd += f'export LD_LIBRARY_PATH={self.libraries_path}:'
            #     cmd += "$LD_LIBRARY_PATH && "
            if self.source_script:
                cmd += f'source {self.source_script} && '
            cmd += f'mpiexec -n {self.nprocs} '
            if self.additional_mpi_options:
                mpi_opts = self.additional_mpi_options.strip("'\"")
                cmd += f'{mpi_opts} '
            cmd += f"--wdir {self.wdir}/{runtype} "
        cmd += f"{self.padcirc_binary}"
        self.logger.info(cmd)
        stdin, stdout, stderr = self.ssh.exec_command(cmd)
        while True:
            out = stdout.readline()
            if not out:
                break
            print(out, end='')
        lines = stderr.readlines()
        msg = "** ERROR: Elevation.gt.ErrorElev, ADCIRC stopping. **"
        if msg in "".join(lines):
            self.logger.warning(msg)
            driver._handle_blowup(lines)
        # filter IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
        msg = "Note: The following floating-point exceptions are signalling:"
        lines = [line for line in lines if msg not in line]
        if len(lines) > 0:
            if msg not in "".join(lines):
                msg = "\n"
                msg += "".join(lines)
                raise Exception(msg)
            else:
                raise Exception(msg)

    def _retrieve_files(self, outdir):
        self.sftp.chdir(str(self.wdir))
        for i, walker in enumerate(self._sftp_walk(str(self.wdir))):
            if i == 0:
                parent = walker[0]
            rdir = pathlib.Path(walker[0])
            for file in walker[1]:
                rfile = rdir / file
                rsubdir = str(rdir).split(parent)[1].strip('/')
                ldir = outdir / rsubdir
                if not ldir.exists():
                    ldir.mkdir()
                self.sftp.get(
                    str(rfile),
                    str(ldir / file))

    def _sftp_walk(self, remotepath):
        """
        https://techtalkontv.wordpress.com/2016/11/05/python-pramiko-sftp-copydownload-all-files-in-a-folder-recursively-from-remote-server/
        """
        path = remotepath
        files = []
        folders = []
        for f in self.sftp.listdir_attr(remotepath):
            if S_ISDIR(f.st_mode):
                folders.append(f.filename)
            else:
                files.append(f.filename)
        if files:
            yield path, files
        for folder in folders:
            new_path = os.path.join(remotepath, folder)
            for x in self._sftp_walk(new_path):
                yield x

    def _cleanup_rundir(self, runtype):
        self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/PE*')
        self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/partmesh.txt')
        self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/fort.13')
        self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/fort.14')
        self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/fort.15')
        if runtype == 'coldstart':
            self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/fort.68.nc')

    @property
    def hostname(self):
        return self._hostname

    @property
    def nprocs(self):
        return self._nprocs

    @property
    def wdir(self):
        return self._wdir

    @property
    def binaries_prefix(self):
        return self._binaries_prefix

    @property
    def port(self):
        return self._port

    @property
    def username(self):
        return self._username

    @property
    def password(self):
        return self._password

    @property
    def pkey(self):
        return self._pkey

    @property
    def writer_procs(self):
        return self._writer_procs

    @property
    def source_script(self):
        return self._source_script

    @property
    def additional_mpi_options(self):
        return self._additional_mpi_options

    @property
    def keep_wdir(self):
        return self._keep_wdir

    @property
    def logger(self):
        try:
            return self.__logger
        except AttributeError:
            self.__logger = logging.getLogger(
                __name__ + '.' + self.__class__.__name__)
            return self.__logger

    @property
    def ssh(self):
        try:
            return self.__ssh
        except AttributeError:
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            kwargs = {
                'hostname': self.hostname,
                'port': self.port,
                'username': self.username,
                'password': self.password
            }
            if self.pkey:
                kwargs.update({
                    'pkey': paramiko.RSAKey.from_private_key_file(
                        self.pkey)})
            # try:
            ssh.connect(**kwargs)
            # except paramiko.ssh_exception.SSHException:
            #     def auth_handler(title, _, fields):
            #         if len(fields) > 1:
            #             raise paramikoSSHException("Expecting one field only.")
            #         return [password]

            #     transport = ssh.get_transport()
            #     transport.auth_interactive(
            #         self.username,
            #         auth_handler)
            self.__ssh = ssh
            return self.__ssh

    @property
    def sftp(self):
        try:
            return self.__sftp
        except AttributeError:
            self.__sftp = self.ssh.open_sftp()
            return self.__sftp

    @property
    def padcirc_binary(self):
        if self.binaries_prefix:
            return self.binaries_prefix.absolute() / 'padcirc'
        else:
            return 'padcirc'

    @property
    def adcprep_binary(self):
        if self.binaries_prefix:
            return self.binaries_prefix.absolute() / 'adcprep'
        else:
            return 'adcprep'

    @property
    def _hostname(self):
        return self.__hostname

    @property
    def _nprocs(self):
        return self.__nprocs

    @property
    def _wdir(self):
        return self.__wdir

    @property
    def _binaries_prefix(self):
        return self.__binaries_prefix

    @property
    def _port(self):
        return self.__port

    @property
    def _username(self):
        return self.__username

    @property
    def _password(self):
        return self.__password

    @property
    def _pkey(self):
        return self.__pkey

    @property
    def _writer_procs(self):
        return self.__writer_procs

    @property
    def _source_script(self):
        return self.__source_script

    @property
    def _additional_mpi_options(self):
        return self.__additional_mpi_options

    @property
    def _keep_wdir(self):
        return self.__keep_wdir

    @_hostname.setter
    def _hostname(self, hostname):
        self.__hostname = hostname

    @_nprocs.setter
    def _nprocs(self, nprocs):
        self.__nprocs = nprocs

    @_wdir.setter
    def _wdir(self, wdir):
        if not wdir:
            wdir = f'/tmp/{uuid.uuid4().hex[:8]}'
        self.__wdir = pathlib.Path(wdir)

    @_binaries_prefix.setter
    def _binaries_prefix(self, binaries_prefix):
        if binaries_prefix:
            binaries_prefix = pathlib.Path(binaries_prefix)
        self.__binaries_prefix = binaries_prefix

    @_port.setter
    def _port(self, port):
        self.__port = port

    @_username.setter
    def _username(self, username):
        self.__username = username

    @_password.setter
    def _password(self, password):
        self.__password = password

    @_pkey.setter
    def _pkey(self, pkey):
        self.__pkey = pkey

    @_writer_procs.setter
    def _writer_procs(self, writer_procs):
        self.__writer_procs = writer_procs

    @_source_script.setter
    def _source_script(self, source_script):
        self.__source_script = source_script

    @_additional_mpi_options.setter
    def _additional_mpi_options(self, additional_mpi_options):
        self.__additional_mpi_options = additional_mpi_options

    @_keep_wdir.setter
    def _keep_wdir(self, keep_wdir):
        self.__keep_wdir = keep_wdir

















        # self._AdcircRun = AdcircRun
        # self._run_directory = run_directory
        # self._binaries_prefix = binaries_prefix
        # self._torque = torque
        # self._hostfile = hostfile
        # self._module_list = module_list
        # self._hostname = hostname
        # self._username = username
        # self._port = port
        # self._IdentityFile = IdentityFile
        # self._nprocs = nprocs
        # self._FORT14_PATH = FORT14_PATH
        # self._FORT13_PATH = FORT13_PATH
        # self._WRITER_PROCS = WRITER_PROCS
        # self._source_script = source_script
        # self.__init_nprocs()
        # self.__init_FORT14_PATH()
        # self.__init_FORT13_PATH()
        # self.__init_ServerConfiguration()

    # def deploy(run=True, monitor=True):
    #     self.__init_ssh_client()
    #     # self._ssh_client.

    # def dump(self, path, filename='ADCIRC_RUN.sh'):
    #     with open(path+'/'+filename, 'w') as f:
    #         for line in self._ServerConfiguration:
    #             f.write(line)

    # def printf(self):
    #     for line in self._ServerConfiguration:
    #         print(line, end="")

    # def __init_nprocs(self):
    #     if self._PBS is not None:
    #         self._nprocs = self._PBS._nprocs

    # def __init_FORT14_PATH(self):
    #     if self._FORT14_PATH is None:
    #         self._FORT14_PATH = '{}/fort.14'.format(self._run_directory)

    # def __init_FORT13_PATH(self):
    #     if self._FORT13_PATH is None:
    #         self._FORT13_PATH = '{}/fort.13'.format(self._run_directory)

    # def __init_ssh_client(self):
    #     self._ssh_client = paramiko.SSHClient()
    #     self._ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    #     try:
    #         self._ssh_client.connect(
    #             hostname=self._hostname,
    #             username=self._username,
    #             password=self._password,
    #             key_filename=self._IdentityFile
    #         )
    #     except paramiko.ssh_exception.SSHException:
    #         self._transport = self._ssh_client.get_transport()
    #         self._transport.auth_interactive(self._username, self.__interactive_auth_handler)

    # def __write_adcirc_run(self):
    #     self._ServerConfiguration.append('$ADCIRC_BINARIES_prefix/adcprep --np $RUN_PROCS --partmesh\n')
    #     self._ServerConfiguration.append('$ADCIRC_BINARIES_prefix/adcprep --np $RUN_PROCS --prepall\n')
    #     string='mpirun -np $NPROCS '
    #     if self._hostfile is not None:
    #         string+='--hostfile {} '.format(self._hostfile)
    #     if 300<=np.abs(self._AdcircRun.NWS)<=350:
    #         string+='$ADCIRC_BINARIES_prefix/padcswan -W $WRITER_PROCS\n'
    #     else:
    #         string+='$ADCIRC_BINARIES_prefix/padcirc -W $WRITER_PROCS\n'
    #     self._ServerConfiguration.append(string)
    #     self._ServerConfiguration.append('\n')

    # def __init_ServerConfiguration(self):
    #     self._ServerConfiguration=list()
    #     if self._PBS is not None:
    #         for line in self._PBS._PBS:
    #             self._ServerConfiguration.append(line)
    #     else:
    #         self._ServerConfiguration.append("#!/bin/bash\n")
    #     self._ServerConfiguration.append('\n')
    #     self._ServerConfiguration.append('set -e\n')
    #     self._ServerConfiguration.append('\n')
    #     self._ServerConfiguration.append('# Set path to ADCIRC binaries:\n')
    #     self._ServerConfiguration.append('ADCIRC_BINARIES_prefix={}\n'.format(self._ADCIRC_BINARIES_prefix))
    #     self._ServerConfiguration.append('\n')
    #     self._ServerConfiguration.append('# Set PATH to the location of the mesh file and fort.13:\n')
    #     self._ServerConfiguration.append('FORT14_PATH={}\n'.format(self._FORT14_PATH))
    #     self._ServerConfiguration.append('FORT13_PATH={}\n'.format(self._FORT13_PATH))
    #     self._ServerConfiguration.append('\n')
    #     if self._module_list is not None:
    #         self._ServerConfiguration.append('# Load HPC modules required to run ADCIRC\n')
    #         for module in self._module_list:
    #             self._ServerConfiguration.append('module load {}\n'.format(module))
    #         self._ServerConfiguration.append('\n')
    #     if self._source_script is not None:
    #         self._ServerConfiguration.append('source {}\n'.format(self._source_script))
    #         self._ServerConfiguration.append('\n')
    #     self._ServerConfiguration.append('#--------------------------------------------#\n')
    #     self._ServerConfiguration.append('#-------------Begin ADCIRC run---------------#\n')
    #     self._ServerConfiguration.append('#--------------------------------------------#\n')
    #     self._ServerConfiguration.append('\n')
    #     self._ServerConfiguration.append('RUNDIR={}\n'.format(self._run_directory))
    #     self._ServerConfiguration.append('WRITER_PROCS={}\n'.format(self._WRITER_PROCS))
    #     if self._PBS is not None:
    #         self._ServerConfiguration.append('RUN_PROCS="$(($PBS_NP-$WRITER_PROCS))"\n')
    #     elif self._nprocs is None:
    #         self._ServerConfiguration.append('RUN_PROCS="$(($(nproc)-$WRITER_PROCS))"\n')
    #     else:
    #         self._ServerConfiguration.append('RUN_PROCS="$(({}-$WRITER_PROCS))"\n'.format(self._nprocs))
    #     self._ServerConfiguration.append('NPROCS="$(($RUN_PROCS+$WRITER_PROCS))"\n')
    #     self._ServerConfiguration.append('COLDSTARTDIR=$RUNDIR/coldstart\n')
    #     self._ServerConfiguration.append('HOTSTARTDIR=$RUNDIR/hotstart\n')
    #     self._ServerConfiguration.append('OUTPUTDIR=$RUNDIR/outputs\n')
    #     self._ServerConfiguration.append('\n')
    #     self._ServerConfiguration.append('#----------------run coldstart---------------#\n')
    #     # self._ServerConfiguration.append('rm -rf $COLDSTARTDIR\n')
    #     self._ServerConfiguration.append('mkdir -p $COLDSTARTDIR\n')
    #     self._ServerConfiguration.append('cd $COLDSTARTDIR\n')
    #     self._ServerConfiguration.append('ln -sf $FORT14_PATH $COLDSTARTDIR/fort.14\n')
    #     self._ServerConfiguration.append('ln -sf $FORT13_PATH $COLDSTARTDIR/fort.13\n')
    #     self._ServerConfiguration.append('ln -sf $RUNDIR/fort.15.coldstart $COLDSTARTDIR/fort.15\n')
    #     self.__write_adcirc_run()
    #     self._ServerConfiguration.append('#----------------run hotstart----------------#\n')
    #     # self._ServerConfiguration.append('rm -rf $HOTSTARTDIR\n')
    #     self._ServerConfiguration.append('mkdir -p $HOTSTARTDIR\n')
    #     self._ServerConfiguration.append('cd $HOTSTARTDIR\n')
    #     self._ServerConfiguration.append('ln -sf $COLDSTARTDIR/fort.67.nc $HOTSTARTDIR/fort.67.nc\n')
    #     self._ServerConfiguration.append('ln -sf $RUNDIR/fort.15.hotstart $HOTSTARTDIR/fort.15\n')
    #     self._ServerConfiguration.append('ln -sf $FORT14_PATH $HOTSTARTDIR/fort.14\n')
    #     self._ServerConfiguration.append('ln -sf $FORT13_PATH $HOTSTARTDIR/fort.13\n')
    #     if hasattr(self._AdcircRun, 'fort22'):
    #         self._ServerConfiguration.append('ln -sf $RUNDIR/fort.22.best_track $HOTSTART/fort.22.best_track\n')
    #         self._ServerConfiguration.append('$ADCIRC_BINARIES_prefix/aswip -n 20 -m 4 -z 2 -w fort.22.best_track\n')
    #         self._ServerConfiguration.append('ln -sf $HOTSTARTDIR/NWS_20_fort.22 $HOTSTARTDIR/fort.22\n')
    #     self.__write_adcirc_run()
    #     self._ServerConfiguration.append('#----------------cleanup-------------------#\n')
    #     # self._ServerConfiguration.append('rm -rf $OUTPUTDIR\n')
    #     self._ServerConfiguration.append('mkdir -p $OUTPUTDIR\n')
    #     self._ServerConfiguration.append('ln -f $HOTSTARTDIR/*.nc $OUTPUTDIR/\n')
    #     self._ServerConfiguration.append('\n')

    # def __interactive_auth_handler(title, instructions, prompt_list):
    #     if prompt_list:
    #         return [getpass.getpass(prompt_list[0][0])]
    #     return []