import pathlib
from datetime import datetime, timedelta
from pyschism.mesh import Mesh
from pyschism.driver import SchismRun
from pyschism.forcing import TidalForcing
from pyschism.cmd.server import ServerConfig, SlurmConfig


class SchismBaseCommand:

    def __init__(self, args):
        self._args = args

    def run(self):

        # dump and exit if generate only
        if self.args.generate_only:
            self.driver.dump(
                self.args.output_directory,
                overwrite=self.args.overwrite
                )
            return

        outputs = self.driver.run(
            outdir=self.output_directory,
            nproc=self.args.nproc,
            overwrite=self.args.overwrite,
            # coldstart=self.coldstart,
            # hotstart=self.hotstart,
            server_config=self.server_config,
        )
        self._output_collection = outputs
        # outputs.maxele.make_plot(show=True)

    def _enable_outputs(self, driver):
        self._init_output_stations(driver)
        self._enable_output(driver, 'elevation', 'surface')

    def _enable_output(self, driver, name, _type):
        fs = getattr(self.args, f"{name}_{_type}_sampling_frequency")
        if fs is not None:
            fs = timedelta(minutes=fs)
        ha = getattr(self.args, f"{name}_{_type}_harmonic_analysis")
        getattr(driver, f"set_{name}_{_type}_output")(
            sampling_frequency=fs,
            netcdf=self.args.netcdf,
            harmonic_analysis=ha
            )

    def _init_output_stations(self, driver):
        pass

    @property
    def driver(self):
        try:
            return self.__driver
        except AttributeError:
            driver = SchismRun(
                self.mesh,
                self.start_date,
                self.end_date,
                self.spinup_time,
                self.tidal_forcing,
                self.wind_forcing,
                self.wave_forcing,
                # self.args.netcdf,
                )
            # self._enable_outputs(driver)
            self.__driver = driver
            return self.__driver

    @property
    def start_date(self):
        try:
            return self.__start_date
        except AttributeError:
            self.__start_date = datetime.strptime(
                self.args.start_date,
                "%Y-%m-%dT%H:%M"
                )
            return self.__start_date

    @property
    def end_date(self):
        try:
            return self.__end_date
        except AttributeError:
            self.__end_date = datetime.strptime(
                self.args.end_date,
                "%Y-%m-%dT%H:%M"
                )
            return self.__end_date

    @property
    def spinup_time(self):
        try:
            return self.__spinup_time
        except AttributeError:
            self.__spinup_time = timedelta(days=self.args.spinup_days)
            return self.__spinup_time

    @property
    def args(self):
        return self._args

    @property
    def mesh(self):
        return self._mesh

    @property
    def tidal_forcing(self):
        try:
            return self.__tidal_forcing
        except AttributeError:
            tidal_forcing = TidalForcing()
            for constituent in self.constituents:
                tidal_forcing.use_constituent(constituent)
            self.__tidal_forcing = tidal_forcing
            return self.__tidal_forcing

    @property
    def wind_forcing(self):
        return None

    @property
    def wave_forcing(self):
        return None

    @property
    def output_directory(self):
        if self.args.output_directory is not None:
            return pathlib.Path(self.args.output_directory).absolute()

    @property
    def constituents(self):
        try:
            return self.__constituents
        except AttributeError:
            # might be better to get these from TidalForcing()
            _major = ('Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2')
            _all = (*_major, 'Mm', 'Mf', 'M4', 'MN4', 'MS4', '2N2', 'S1')
            if ('all' in self.args.constituents
                    and len(self.args.constituents) > 1):
                msg = 'When using all, must only pass one'
                raise IOError(msg)

            elif ('major' in self.args.constituents
                    and len(self.args.constituents) > 1):
                msg = 'When using major, must only pass one'
                raise IOError(msg)
            if 'all' in self.args.constituents:
                constituents = _all
            elif 'major' in self.args.constituents:
                constituents = _major
            else:
                constituents = self.args.constituents
            self.__constituents = constituents
            return self.__constituents

    @property
    def server_config(self):
        if self.args.hostname:
            if (not self.args.use_slurm or
                    not self.args.use_torque or
                    not self.args.use_pbs):
                server_config = ServerConfig(
                    hostname=self.args.hostname,
                    nprocs=self.args.nproc,
                    wdir=self.args.wdir,
                    binaries_prefix=self.args.binaries_prefix,
                    source_script=self.args.source_script,
                    additional_mpi_options=self.args.additional_mpi_options,
                    )

            elif self.args.use_slurm:
                raise NotImplementedError
                server_config = SlurmConfig(
                    )

            elif self.args.use_torque or self.args.use_pbs:
                raise NotImplementedError

        else:
            server_config = None

        self.__server_config = server_config
        return self.__server_config

    @property
    def _args(self):
        return self.__args

    @property
    def _mesh(self):
        try:
            return self.__mesh
        except AttributeError:
            self.__mesh = Mesh.open(
                hgrid=self.args.hgrid,
                vgrid=self.args.vgrid,
                crs=self.args.crs
                )
            return self.__mesh

    @_args.setter
    def _args(self, args):
        self.__args = args
