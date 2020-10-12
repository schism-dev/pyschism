#! /usr/bin/env python
from .cmd import plotting


def plot():
    plotting.main()

    



# import pathlib
# import argparse
# import json
# from functools import lru_cache
# import os
# import logging
# from pyschism import forcing, Mesh, Param
# from pyschism.driver.driver import SchismRun


# class PySchism:

#     def __init__(self, args):
#         self.args = args
#         self._init_logger()
#         self._init_forcings()

#     def main(self):
#         return self.driver.run()

#     @property
#     @lru_cache(maxsize=None)
#     def driver(self):
#         return SchismRun(self.mesh, self.param)

#     @property
#     @lru_cache(maxsize=None)
#     def mesh(self):
#         return Mesh.open(
#             self.hgrid,
#             self.vgrid,
#             self.fgrid,
#             crs=self.mesh_crs
#             )

#     @property
#     @lru_cache(maxsize=None)
#     def hgrid(self):
#         return pathlib.Path(
#             os.path.expandvars(self.mesh_opts["hgrid"])).resolve()

#     @property
#     @lru_cache(maxsize=None)
#     def vgrid(self):
#         vgrid = self.mesh_opts.get("vgrid", None)
#         if vgrid is not None:
#             vgrid = pathlib.Path(os.path.expandvars(vgrid)).resolve()
#         return vgrid

#     @property
#     def fgrid(self):
#         fgrid = self.mesh_opts.get("fgrid", None)
#         if fgrid is not None:
#             fgrid = pathlib.Path(os.path.expandvars(fgrid)).resolve()
#         return fgrid

#     @property
#     @lru_cache(maxsize=None)
#     def mesh_crs(self):
#         crs = self.mesh_opts.get('crs', None)
#         if crs is None:
#             msg = "Must specify 'crs' on 'mesh' section of the template."
#             raise AttributeError(msg)
#         return crs

#     @property
#     @lru_cache(maxsize=None)
#     def config(self):
#         src = pathlib.Path(self.args.config_file).resolve()
#         with open(src, 'r') as json_file:
#             config = json.load(json_file)
#         return config

#     @property
#     @lru_cache(maxsize=None)
#     def mesh_opts(self):
#         return self.config["mesh"]

#     @property
#     def meta(self):
#         return self.config.get("meta", {})

#     @property
#     def overwrite(self):
#         return self.meta.get("overwrite", False)

#     @property
#     def project_name(self):
#         return self.meta.get("project_name", False)

#     @property
#     @lru_cache(maxsize=None)
#     def output_directory(self):
#         outdir = self.meta.get("output_directory", "'.'")
#         outdir = pathlib.Path(os.path.expandvars(outdir)).resolve()
#         outdir.mkdir(parents=True, exist_ok=self.overwrite)
#         return outdir

#     @property
#     def server(self):
#         return self.config.get("server", False)

#     @property
#     @lru_cache(maxsize=None)
#     def param(self):
#         param = Param()

#         return param

#     @property
#     def hostname(self):
#         if self.server:
#             return self.server.get("hostname", "localhost")
#         return False

#     def _init_forcings(self):
#         forcings = self.config['forcing'].get("forcing", {})
#         self._init_elev_forcing(forcings)

#     def _init_elev_forcing(self, forcings):
#         if "tidal_constituents" in forcings:
#             tides = forcing.Tides()
#             if "all" in forcings["tidal_constituents"]:
#                 tides.use_all(**forcings["tidal_constituents"]["all"])
#             elif "major" in forcings["tidal_constituents"]:
#                 tides.use_major(
#                     tides.use_major(**forcings["tidal_constituents"]["major"]))
#             else:
#                 for const, kwargs in forcings["tidal_constituents"].items():
#                     tides.use_constituent(
#                         const,
#                         potential=kwargs.get("potential"),
#                         forcing=kwargs.get("forcing"),
#                         )
#             for data in forcings["tidal_constituents"].values():
#                 bnd_ids = data.get("boundary_id", [])
#                 if len(bnd_ids) == 0:
#                     self.mesh.add_forcing(tides)
#                 else:
#                     for bnd_id in bnd_ids:
#                         self.mesh.add_forcing(tides, id=bnd_id)

#     def _init_logger(self):
#         # log_level = {
#         #     "info": logging.INFO,
#         #     "debug": logging.DEBUG,
#         #     "warning": logging.WARNING,
#         # }[self.config["meta"].get("log_level", "info").lower()]
#         # self.logger = logging.getLogger(__name__)
#         # self.logger.setLevel(log_level)
#         # from geomesh import logger
#         # logger.setLevel(log_level)
#         logging.basicConfig(level={
#             "info": logging.INFO,
#             "debug": logging.DEBUG,
#             "warning": logging.WARNING,
#         }[self.config["meta"].get("log_level", "info").lower()])
#         logging.getLogger('matplotlib').setLevel(logging.WARNING)
#         logging.getLogger('fiona').setLevel(logging.WARNING)
#         logging.getLogger('rasterio').setLevel(logging.WARNING)
#         self.logger = logging.getLogger(__name__)


# def config_template():
#     config = dict()
#     config['meta'] = {
#         "_comments": "<> denotes optional parameters.",
#         "project_name": "<can be anything>",
#         "output_directory": "can contain environment variables",
#         "overwrite": "<bool>",
#         "log_level": "<warning | info | debug>",
#         "binary": "pschism_TVD-VL"
#         }

#     config['mesh'] = {
#         "hgrid": "Path to hgrid.",
#         "vgrid": "<Path to vgrid>",
#         "friction": {
#             "ftype": "manning | drag | rough",
#             "fgrid": "path to fgrid | float (constant)",
#         },

#     }

#     config["forcing"] = {
#             "tidal_constituents": {
#                 "<constituent | major | all >": {
#                     "boundary_id (optional)": [
#                         "list of boundary ids on which to apply forcing."
#                     ],
#                     "potential (required)": "<bool, use as potential>",
#                     "forcing (required)": "<bool, use as boundary forcing>"
#                 }
#             },
#             "winds": {
#                 "best_track": "<ATCF best track ID>"
#                 },
#             "waves": {
#                 "properties": "<To be discussed. Disabled for now.>"
#             },
#             "temperature": {
#                 "properties": "<To be discussed. Disabled for now.>"
#             },
#             "salinity": {
#                 "properties": "<To be discussed. Disabled for now.>"
#             },
#             "tracers": {
#                 "properties": "<To be discussed. Disabled for now.>"
#             }
#         }
#     config["start_date"] = r'%Y-%m-%dT%H:%M'
#     config["utc_start"] = '<float, defaults to 0.0>'
#     config["run_days"] = "float"
#     config["timestep"] = "float"
#     config["spinup_days"] = "float"
#     config["date_formatting"] = r'<%Y-%m-%dT%H:%M>'
#     config["outputs"] = {
#         "nspool": "int | float (seconds)",
#         "ihfskip": "<int>",
#         "variables": ["names of requested outputs"]

#     }

#     config['server'] = {
#         "hostname": "localhost",
#         "port": 22,
#     }
#     return config


# def parse_args():
#     parser = argparse.ArgumentParser()
#     parser.add_argument(
#         "config_file",
#         help="Path to configuration file.",
#         nargs='?')
#     return parser.parse_args()


# def main():
#     args = parse_args()
#     if args.config_file is None:
#         print(json.dumps(config_template(), indent=4))
#         exit()
#     exit(PySchism(args).main())


# if __name__ == '__main__':
#     main()
