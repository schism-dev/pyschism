from pyschism.cmd.server import ServerConfig


class SlurmConfig(ServerConfig):

    def __init__(
        self,
        hostname,
        numprocs,
        wdir,
        account,
        walltime,
        modules=[],
        **kwargs
    ):
        super().__init__(hostname, numprocs, wdir, **kwargs)
        self._account = account
        self._walltime = walltime
        self._modules = modules

    def _deploy_files_to_server(self):
        # avoids checking wdir on server.
        raise NotImplementedError
