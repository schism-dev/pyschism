import pathlib
import os
from typing import Union

from pyschism.server.base import ServerConfig


class ServerConfigDescriptor:

    def __set__(self, obj, server_config: Union[ServerConfig, None]):
        if not isinstance(server_config, (ServerConfig, type(None))):
            raise TypeError("Argument server_config mut be of type "
                            f"{ServerConfig} or None, not type "
                            f"{type(server_config)}.")
        self._server_config = server_config

    def __get__(self, obj, val):
        return self._server_config


class BashDriver:

    _server_config = ServerConfigDescriptor()

    def __init__(self, server_config: ServerConfig = None):
        self._server_config = server_config

    def __str__(self):
        f = []
        if self._server_config is not None:
            f.append(str(self._server_config))
        else:
            f.append("#!/bin/bash")

        return "\n".join(f)

    def write(self, path: Union[str, os.PathLike], overwrite: bool = False):
        path = pathlib.Path(path)
        if path.exists() and not overwrite:
            raise FileExistsError(
                f"File {str(path)} exists and overwrite is False.")
        with open(path, 'w') as fp:
            fp.write(str(self))

    @property
    def server_config(self):
        return self._server_config
