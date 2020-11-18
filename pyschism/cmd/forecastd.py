from argparse import Namespace
import atexit
# from datetime import timezone, timedelta
import hashlib
import os
import pathlib
from typing import List

from appdirs import user_data_dir  # type: ignore[import]
from apscheduler.schedulers.blocking import BlockingScheduler  # type: ignore[import]  # noqa: E501
from apscheduler.jobstores.sqlalchemy import SQLAlchemyJobStore  # type: ignore[import]  # noqa: E501
from apscheduler.job import Job  # type: ignore[import]
from daemons.prefab import run  # type: ignore[import]
import pytz
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import Column, String

from pyschism.forcing.tides import Tides
from pyschism.logger import get_logger

LOGGER = get_logger(__file__)

FORECASTD_DATA_DIR = pathlib.Path(user_data_dir()) / 'pyschism/forecastd'
FORECASTD_DATA_DIR.mkdir(exist_ok=True, parents=True)
SQLITE_DATABASE = FORECASTD_DATA_DIR / 'scheduler.sqlite'

DEBUG_CLEAR_DEFAULT_DATABASE = os.getenv('DEBUG_CLEAR_DEFAULT_DATABASE', '')
if DEBUG_CLEAR_DEFAULT_DATABASE == '.TRUE.':
    os.remove(SQLITE_DATABASE)

SQLITE_ENGINE = create_engine(f'sqlite:///{str(SQLITE_DATABASE)}', echo=False)


Base = declarative_base()


class HgridPointer(Base):
    __tablename__ = 'hgrid_pointers'
    id = Column(String, primary_key=True, nullable=False)
    md5 = Column(String, unique=True)
    tag = Column(String)
    path = Column(String)
    vgrid = relationship("VgridPointer")


class VgridPointer(Base):
    __tablename__ = 'vgrid_pointers'
    id = Column(String, primary_key=True, nullable=False)
    md5 = Column(String, unique=True)
    tag = Column(String)
    path = Column(String)


Base.metadata.create_all(SQLITE_ENGINE)


SCHEDULER = BlockingScheduler(
            jobstores={'default': SQLAlchemyJobStore(engine=SQLITE_ENGINE)},
            logger=LOGGER,
            timezone=pytz.utc,
            # jobstore_retry_interval=,
            # job_defaults=,
            # executors=,
        )


def get_forecast(*args, **kwargs):
    raise NotImplementedError('get_forecast()')


class Forecastd(run.RunDaemon):

    def __init__(self, args: Namespace):
        self.args = args
        # TODO: Put a system-lock in the pidfile
        super().__init__(pidfile=str(FORECASTD_DATA_DIR / 'index.lock'))

    def run(self):
        atexit.register(lambda: SCHEDULER.shutdown())
        SCHEDULER.start()

    def add(self):

        job_id = self._get_job_id()
        print(job_id)

        pid = self.pid
        if pid is not None:
            # if scheduler is running, stop it.
            self.stop()

        SCHEDULER.add_job(
            Job(
                get_forecast,
                # trigger=None,
                # args=None,
                # kwargs=None,
                id=job_id,
                # name=None,
                # misfire_grace_time=undefined,
                # coalesce=undefined,
                # max_instances=undefined,
                # next_run_time=undefined,
                # jobstore='default',
                # executor='default',
                # replace_existing=False,
                # **trigger_args
                )
            )

        if pid is not None:
            self.start()  # fork and exit parent

    def _get_job_id(self):

        # constituents
        job_id: List = []

        # compute hash for mesh file

        # hgrid_id =

        # compute ID for tidal constituent requests
        if self.args.constituents is not None:
            job_id.append(get_tidal_constituents_id(self.args.constituents))

        return hashlib.md5(''.join(job_id).encode('utf-8')).hexdigest()


def get_tidal_constituents_id(constituents):
    if 'all' in constituents:
        return ''.join(Tides.constituents)
    elif 'major' in constituents:
        return ''.join(Tides.major_constituents)
    else:
        f = []
        for constituent in Tides.constituents:
            if constituent in constituents:
                f.append(constituent)
        return ''.join(f)


def get_mesh_id(path):
    hashlib.md5(''.join(job_id).encode('utf-8')).hexdigest()