import os
import pathlib

from appdirs import user_data_dir  # type: ignore[import]
from daemons.prefab import run  # type: ignore[import]
from flask import Flask
import waitress

VIEWERD_DAEMON_LOCK_DIR = pathlib.Path(user_data_dir()) / 'pyschism'

app = Flask(__name__)


@app.route('/')
def home():
    return 'Could we bootstrap React here?'


class Viewerd(run.RunDaemon):

    def __init__(self, args):
        self._args = args
        # TODO: Put a system-lock in the pidfile
        super().__init__(pidfile=str(VIEWERD_DAEMON_LOCK_DIR / 'index.lock'))
        getattr(self, args.action)()

    def run(self):
        host = os.getenv('VIEWERD_HOSTNAME', 'http://127.0.0.1').split('/')[-1]
        port = os.getenv('VIEWERD_PORT', 5000)
        if self._args.deploy:
            waitress.serve(app, host=host, port=port)
        else:
            app.run(host=host, port=port)
