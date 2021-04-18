


# from argparse import Namespace
# import json
# import logging
# import pathlib

# from pyschism import dates
# from pyschism.mesh import Hgrid, Vgrid, Fgrid

# logger = logging.getLogger(__name__)

# CONFIG_FILE_NAME = 'config.json'
# STATIC_DIRECTORY = 'static'
# FORECAST_DIRECTORY = 'forecast'


# class ForecastInit:

#     def __init__(self, args: Namespace):
#         self.args = args
#         self.save_user_arguments()

#     def save_user_arguments(self):
#         logger.info(
#             f"Writing configuration file to path {self.config_file}")
#         with open(self.config_file, 'w') as fp:
#             json.dump(self.args.__dict__, fp, indent=4)

#     @property
#     def project_directory(self):
#         if not hasattr(self, '_project_directory'):
#             self._project_directory = pathlib.Path(self.args.project_directory)
#             self._project_directory.mkdir(exist_ok=self.args.overwrite)
#         return self._project_directory

#     @property
#     def config_file(self):
#         return self.project_directory / CONFIG_FILE_NAME

#     @property
#     def static_directory(self):
#         if not hasattr(self, '_static_directory'):
#             self._static_directory = self.project_directory / STATIC_DIRECTORY
#             self._static_directory.mkdir(exist_ok=self.args.overwrite)
#         return self._static_directory

#     @property
#     def forecast_directory(self):
#         return self.project_directory / FORECAST_DIRECTORY

#     @property
#     def hgrid_path(self):
#         if not hasattr(self, '_hgrid_path'):
#             self._hgrid_path = self.static_directory / 'hgrid.gr3'
#             if not self._hgrid_path.exists() or self.args.overwrite is True:
#                 logger.info(
#                     f'Writing initial hgrid file to {self._hgrid_path}.')
#                 Hgrid.open(self.args.hgrid, crs=self.args.hgrid_crs).write(
#                         self._hgrid_path,
#                         overwrite=self.args.overwrite)
#         return self._hgrid_path

#     @property
#     def hgrid(self):
#         if not hasattr(self, '_hgrid'):
#             self._hgrid = Hgrid.open(self.hgrid_path, crs=self.args.hgrid_crs)
#         return self._hgrid

    # @property
    # def current_forecast_directory(self):
    #     return self.forecast_directory / dates.nearest_cycle_date()

    # @property
    # def start_date(self):


    

















# from argparse import Namespace
# from datetime import datetime, timedelta, timezone
# import json
# import logging
# import os
# import pathlib
# import shutil
# import subprocess

# import numpy as np  # type: ignore[import]
# import pytz

# from pyschism.domain import ModelDomain
# from pyschism.driver import ModelDriver
# from pyschism.forcing import Tides
# from pyschism.mesh.vgrid import Vgrid
# from pyschism.server import SlurmConfig

# # from pyschism.logger import get_logger


# class ProjectDirectory:

#     def __get__(self, obj, val):
#         project_directory = obj.__dict__.get('project_directory')
#         if project_directory is None:
#             project_directory = pathlib.Path(obj.args.project_directory)
#             project_directory.mkdir(parents=True, exist_ok=True)
#             obj.__dict__['project_directory'] = project_directory
#         return project_directory


# class ConfigFile:

#     def __get__(self, obj, val):
#         config_file = obj.__dict__.get('config_file')
#         if config_file is None:
#             config_file = obj.project_directory / 'config.json'
#             obj.__dict__['config_file'] = config_file
#         return config_file


# class ColdstartDirectory:

#     def __get__(self, obj, val):
#         coldstart_directory = obj.__dict__.get('coldstart_directory')
#         if coldstart_directory is None:
#             timestamp = str(obj.target_datetime).replace(' ', 'T')
#             coldstart_directory = obj.project_directory / \
#                 'coldstart' / f'{timestamp}'
#             coldstart_directory.mkdir(
#                 parents=True, exist_ok=True)
#             obj.__dict__['coldstart_directory'] = coldstart_directory
#         return coldstart_directory


# class ColdstartDomain:

#     def __get__(self, obj, val):
#         model_domain = obj.__dict__.get('coldstart_domain')
#         if model_domain is None:
#             model_domain = ModelDomain.open(
#                 obj.hgrid_path,
#                 obj.fgrid_path,
#                 # TODO: Vgrid.open() raises NotImplementedError
#                 # obj.vgrid_path,
#                 hgrid_crs=obj.args.hgrid_crs,
#                 fgrid_crs=obj.args.fgrid_crs
#             )
#             if obj.tides is not None:
#                 model_domain.add_boundary_condition(obj.tides)
#         return model_domain


# class ColdstartDriver:

#     def __get__(self, obj, val):
#         coldstart = obj.__dict__.get('coldstart_driver')
#         if coldstart is None:
#             start_date = obj.target_datetime - timedelta(
#                 days=obj.args.spinup_days)
#             rnday = obj.target_datetime - start_date
#             coldstart = ModelDriver(
#                     model_domain=obj.coldstart_domain,
#                     dt=obj.args.timestep,
#                     rnday=rnday,
#                     dramp=(2./3.)*rnday,
#                     start_date=start_date,
#                     nhot_write=True,
#                     # ibc=obj.ibc,
#                     # drampbc=obj.drampbc,
#                     server_config=obj.server_config,
#                     )
#             obj.__dict__['coldstart_driver'] = coldstart
#         return coldstart


# class HgridPath:

#     def __get__(self, obj, val):
#         hgrid_path = obj.__dict__.get('hgrid_path')
#         if hgrid_path is None:
#             hgrid_path = obj.static_files_directory / 'hgrid.gr3'
#             if hgrid_path.is_file() and obj.args.overwrite is True:
#                 hgrid_path.unlink()
#             if not hgrid_path.is_file():
#                 shutil.copy2(obj.args.hgrid, hgrid_path,
#                              follow_symlinks=True)
#             obj.__dict__['hgrid_path'] = hgrid_path
#         return hgrid_path


# class VgridPath:

#     def __get__(self, obj, val):
#         vgrid_path = obj.__dict__.get('vgrid_path')
#         if vgrid_path is None:
#             vgrid_path = obj.static_files_directory / 'vgrid.in'
#             if vgrid_path.is_file():
#                 if obj.args.overwrite is True:
#                     vgrid_path.unlink()
#             if obj.args.vgrid is None:
#                 if not vgrid_path.exists():
#                     Vgrid().write(vgrid_path)
#             else:
#                 shutil.copy2(
#                     obj.args.vgrid, vgrid_path, follow_symlinks=True)
#             obj.__dict__['vgrid_path'] = vgrid_path
#         return vgrid_path


# class FgridPath:

#     def __get__(self, obj, val):
#         fgrid_path = obj.__dict__.get('fgrid_path')
#         if fgrid_path is None:
#             if obj.args.fgrid_type == 'auto':
#                 fgrid = pathlib.Path(obj.args.fgrid)
#                 fgrid_path = obj.static_files_directory / fgrid.name
#             else:
#                 fgrid_path = obj.static_files_directory / \
#                     (obj.args.fgrid_type + '.gr3')
#             if obj.args.overwrite is True:
#                 if fgrid_path.is_file():
#                     os.remove(fgrid_path)
#             shutil.copy2(obj.args.fgrid, fgrid_path, follow_symlinks=True)
#             obj.__dict__['fgrid_path'] = fgrid_path
#         return fgrid_path


# class StaticFilesDirectory:

#     def __get__(self, obj, val):
#         static_files = obj.__dict__.get('static_files')
#         if static_files is None:
#             static_files = obj.project_directory / 'static'
#             static_files.mkdir(exist_ok=True)
#             obj.__dict__['static_files'] = static_files
#         return static_files


# class TidesDescriptor:

#     def __get__(self, obj, val):
#         tides = obj.__dict__.get('tides')
#         if tides is None:
#             if not obj.args.all_constituents \
#                     and not obj.args.major_constituents \
#                     and not obj.args.constituents:
#                 return
#             else:
#                 tides = Tides(
#                     database=obj.args.tidal_database,
#                     velocity=obj.args.bnd_vel
#                 )
#                 if obj.args.all_constituents:
#                     tides.use_all()
#                 if obj.args.major_constituents:
#                     tides.use_major()
#                 if obj.args.constituents:
#                     for constituent in obj.args.constituents:
#                         tides.use_constituent(constituent)
#                 if obj.args.Z0 is not None:
#                     tides.add_Z0(obj.args.Z0)
#             obj.__dict__['tides'] = tides
#         return tides


# class ServerConfigDescriptor:

#     def __get__(self, obj, val):
#         server_config = obj.__dict__.get("server_config")
#         if server_config is None:
#             if obj.args.server_config == "slurm":
#                 kwargs = {
#                     "account": obj.args.account,
#                     "ntasks": obj.args.nproc,
#                     "partition": obj.args.partition,
#                     "walltime": obj.args.walltime,
#                     "mail_type": obj.args.mail_type,
#                     "mail_user": obj.args.mail_user,
#                     "log_filename": obj.args.log_filename,
#                     "modules": obj.args.modules,
#                     # "schism_binary": obj.args.schism_binary,
#                     "extra_commands": obj.args.extra_commands,
#                     "launcher": obj.args.slurm_launcher,
#                     "nodes": obj.args.slurm_nodes
#                 }
#                 # if obj.args.slurm_filename is not None:
#                 #     kwargs.update({"filename": obj.args.ntasks})
#                 # if obj.args.slurm_rundir is not None:
#                 #     kwargs.update({"run_directory": obj.args.slurm_rundir})
#                 # if obj.args.run_name is not None:
#                 #     kwargs.update({"run_name": obj.args.run_name})
#                 server_config = SlurmConfig(**kwargs)
#                 obj.__dict__['server_config'] = server_config
#         return server_config


# class ForecastInit:

#     project_directory = ProjectDirectory()
#     config_file = ConfigFile()
#     coldstart_directory = ColdstartDirectory()
#     coldstart_domain = ColdstartDomain()
#     coldstart_driver = ColdstartDriver()
#     static_files_directory = StaticFilesDirectory()
#     hgrid_path = HgridPath()
#     vgrid_path = VgridPath()
#     fgrid_path = FgridPath()
#     tides = TidesDescriptor()
#     server_config = ServerConfigDescriptor()

#     def __init__(self, args: Namespace):

#         """
#         TODO: Write a file that inidicates this project is currently being
#         accessed by a process. Write the pid of this process to the file.
#         May be implemented by a "check_project()" method.
#         """
#         self._args = args
#         self._write_config_file()
#         self._symlink_files(self.coldstart_directory,
#                             self.coldstart_domain.ics)
#         self.logger.info('Writting coldstart files to disk...')
#         self.coldstart_driver.write(
#             self.coldstart_directory,
#             hgrid=False,
#             vgrid=False,
#             fgrid=False,
#             wind_rot=False,
#             overwrite=self.args.overwrite,
#             use_param_template=self.args.use_param_template
#         )

#         if self.args.skip_run is False:
#             # release memory before launching SCHISM.
#             self.logger.info('Releasing memory before calling SCHISM...')
#             for item in list(self.__dict__.keys()):
#                 if not item.startswith('_'):
#                     del self.__dict__[item]
#             self.logger.info('Calling SCHISM using make.')
#             subprocess.check_call(
#                 ["make", "run"],
#                 cwd=self.coldstart_directory
#             )
#         # else:
#             # self.logger.info("Skipping coldstart run.")

#         # self.logger.info("Finished coldstart sequence.")

#     @property
#     def args(self):
#         return self._args

#     def _write_config_file(self):
#         if self.config_file.is_file() and self.args.overwrite is False:
#             raise IOError(
#                 f'The given directory {str(self.project_directory)} has '
#                 'been initialized previously. Please use\npyschism '
#                 'forecast --overwrite init [...]\nto allow overwrite of '
#                 'previous initialization options.')
#         self.logger.info(
#             f"Writting configuration file to path {self.config_file}")
#         with open(self.config_file, 'w') as fp:
#             json.dump(self.args.__dict__, fp, indent=4)

#     def _symlink_files(self, target_directory, ics):
#         self.logger.info(
#             f"Establishing symlinks to target_directory: {target_directory}")
#         hgrid_lnk = target_directory / 'hgrid.gr3'
#         vgrid_lnk = target_directory / 'vgrid.in'
#         fgrid_lnk = target_directory / f'{self.fgrid_path.name}'
#         if self.args.overwrite is True:
#             if hgrid_lnk.exists():
#                 hgrid_lnk.unlink()
#             if vgrid_lnk.exists():
#                 vgrid_lnk.unlink()
#             if fgrid_lnk.exists():
#                 fgrid_lnk.unlink()
#         os.symlink(
#             os.path.relpath(
#                 self.hgrid_path, target_directory), hgrid_lnk)
#         os.symlink(
#             os.path.relpath(
#                 self.vgrid_path, target_directory), vgrid_lnk)
#         os.symlink(
#             os.path.relpath(
#                 self.fgrid_path, target_directory), fgrid_lnk)
#         if ics == 2:
#             self._symlink_hgridll(target_directory)

#     def _symlink_hgridll(self, target_directory):
#         hgridll_lnk = target_directory / 'hgrid.ll'
#         if hgridll_lnk.exists():
#             hgridll_lnk.unlink()
#         if self.args.overwrite is True:
#             if self.args.overwrite is True:
#                 if hgridll_lnk.exists():
#                     hgridll_lnk.unlink()
#         os.symlink(os.path.relpath(
#             self.hgrid_path, target_directory), hgridll_lnk)

#     @property
#     def target_datetime(self):
#         if not hasattr(self, '_target_datetime'):
#             now_utc = datetime.now(timezone.utc)
#             nearest_cycle = int(6 * np.floor(now_utc.hour/6))
#             self._target_datetime = datetime(
#                 now_utc.year, now_utc.month, now_utc.day, nearest_cycle,
#                 tzinfo=timezone.utc)
#             self.logger.info(
#                 f'Target datetime is: {str(self._target_datetime)}.')
#         return self._target_datetime

#     @property
#     def logger(self):
#         if not hasattr(self, '_logger'):
#             self._logger = logging.getLogger(f"{self.__class__.__name__}")
#         return self._logger
