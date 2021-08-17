from abc import ABC
from datetime import datetime, timedelta
import os
import pathlib
from typing import Dict, Union


# import cf
from netCDF4 import Dataset
import numpy as np
import xarray

from pyschism.enums import (
    IofWetdryVariables,
    IofZcorVariables,
    IofHydroVariables,
    IofDvdVariables,
    IofWwmVariables,
    IofGenVariables,
    IofAgeVariables,
    IofSedVariables,
    IofEcoVariables,
    IofIcmVariables,
    IofCosVariables,
    IofFibVariables,
    IofSed2dVariables,
    IofMarshVariables,
    IofIceVariables,
    IofAnaVariables,
)
from pyschism.mesh.base import Gr3


def get_stack_id_by_datetime(dt: datetime, flattened_timevector, stacks) -> str:
    output_stack = None
    for stack_id, stack_data in stacks.items():
        if dt in stack_data['timevector']:
            output_stack = stack_id
    if output_stack is None:
        raise ValueError(
            'The are no output slices corresponding to requested datetime:'
            f' {dt}. Available datetimes are '
            f'{flattened_timevector}.')
    return output_stack


# def aggregate_by_datetime(shape, filenames, n_local_to_global, hgrid, name, dt, stacks, flattened_timevector):
#     values = np.full(shape, np.nan)
#     stack_id = get_stack_id_by_datetime(dt, flattened_timevector, stacks)
#     local_time_index = get_local_time_index(dt, stacks, flattened_timevector)
#     for file in filenames:
#         cpu_id = file.name.split('_')[1]
#         _n_local_to_global = n_local_to_global[cpu_id]
#         local_node_ids = list(_n_local_to_global.keys())
#         global_node_ids = list(
#             map(lambda x: _n_local_to_global[x], local_node_ids))
#         idxs = list(map(hgrid.nodes.get_index_by_id, global_node_ids))
#         if file.name.endswith(f'{stack_id}.nc'):
#             values[idxs] = Dataset(file)[name][local_time_index, :]
#     return values


class OutputVariableCombiner(ABC):

    def __init__(self, parent, name, start_step=0):
        self.parent = parent
        self.name = name
        self.current_step = start_step

    @property
    def flattened_timevector(self):
        return self.parent.flattened_timevector

    @property
    def shape(self):
        return self.ncvar.shape[1:]

    @property
    def filenames(self):
        return self.parent.filenames

    @property
    def hgrid(self):
        return self.parent.hgrid

    @property
    def n_local_to_global(self):
        return self.parent.n_local_to_global

    @property
    def s_local_to_global(self):
        return self.parent.s_local_to_global

    @property
    def e_local_to_global(self):
        return self.parent.e_local_to_global

    @property
    def stacks(self):
        return self.parent.stacks

    def animation(
            self,
            save=False,
            fps=3,
            start_frame=0,
            end_frame=-1,
            figsize=None,
            wireframe=False,
            cmap='jet',
            levels=256,
            show=False,
            xmin=None,
            xmax=None,
            ymin=None,
            ymax=None,
            vmin=None,
            vmax=None,
    ):
        from matplotlib.animation import FuncAnimation
        import matplotlib.pyplot as plt
        self.current_step = start_frame
        fig = plt.figure(figsize)
        ax = fig.add_subplot(111)

        plt.tight_layout(pad=2)

        triangulation = self.hgrid.triangulation
        xmin = np.min(self.hgrid.x) if xmin is None else xmin
        xmax = np.max(self.hgrid.x) if xmax is None else xmax
        ymin = np.min(self.hgrid.y) if ymin is None else ymin
        ymax = np.max(self.hgrid.y) if ymax is None else ymax
        vmin = np.min(self.values) if vmin is None else vmin
        vmax = np.max(self.values) if vmax is None else vmax

        triangulation.set_mask(self.hgrid.elements.get_triangulation_mask(self.parent.wetdry_elem.values))

        def animate(index):
            self.current_step = index
            _ax = fig.get_axes()
            ax.clear()
            if len(_ax) > 1:
                cax = _ax[1]
                cax.cla()
            else:
                cax = None

            triangulation.set_mask(self.hgrid.elements.get_triangulation_mask(self.parent.wetdry_elem.values))

            if wireframe:
                ax.triplot(triangulation, color='k', linewidth=0.7)

            ax.tricontourf(
                triangulation,
                self.values,
                cmap=cmap,
                levels=levels,
                vmin=vmin,
                vmax=vmax
                )

            ax.set_ylim(ymin, ymax, auto=True)
            ax.set_xlim(xmin, xmax, auto=True)

            ax.set_xlabel('Longitude (°E)')
            ax.set_ylabel('Latitude (°N)')

            ax.set_title(self.flattened_timevector[index].strftime('%b %d, %Y %H:%M'))
            m = plt.cm.ScalarMappable(cmap=cmap)
            m.set_array(self.values[index])
            m.set_clim(vmin, vmax)
            cbar = fig.colorbar(
                m, cax=cax, format='%.1f',
                boundaries=np.linspace(vmin, vmax, levels)
            )
            ax.axis('scaled')
            # cbar = fig.colorbar(_ax)
            # cbar.ax.set_ylabel(f'{variable} [{unit}]', rotation=90)

        end_frame = end_frame % len(self.flattened_timevector) \
            if end_frame < 0 else end_frame
        start_frame = start_frame % len(self.flattened_timevector) \
            if start_frame < 0 else start_frame
        frames = range(start_frame, end_frame)
        anim = FuncAnimation(
            fig,
            animate,
            frames,
            blit=False
            )

        if save not in [False, None]:
            anim.save(
                pathlib.Path(save),
                writer='imagemagick',
                fps=fps
            )

        if show:
            plt.show()

        return anim

    def aggregate(self, dt: Union[datetime, timedelta, int]):
        if isinstance(dt, datetime):
            return self.aggregate_by_datetime(dt)
        # elif isinstance(dt, timedelta):
        #     return self.aggregate_by_timedelta(dt)
        elif isinstance(dt, int):
            return self.aggregate_by_index(dt)
        else:
            raise TypeError('Argument dt must be datetime, timedelta or int '
                            f'not type {type(dt)}')

    # def aggregate_parallel(self, nprocs=16):
    #     from time import time
    #     from multiprocessing import Pool
    #     start = time()
    #     with Pool(processes=nprocs) as p:
    #         # shape, filenames, n_local_to_global, hgrid, name, dt, stacks, flattened_timevector
    #         res = p.starmap(
    #             aggregate_by_datetime,
    #             [
    #                 (self.shape, self.filenames, self.n_local_to_global,
    #                  self.hgrid, self.name, dt, self.stacks,
    #                  self.flattened_timevector)
    #                 for dt in self.flattened_timevector
    #             ]
    #         )
    #     p.join()
    #     print(f'parallel took {time()-start}')
    #     start = time()
    #     for dt in self.flattened_timevector:
    #         aggregate_by_datetime(
    #             self.shape, self.filenames, self.n_local_to_global, self.hgrid,
    #             self.name, dt, self.stacks,
    #             self.flattened_timevector)
    #     print(f'serial took {time()-start}')
    #     # breakpoint()

    def aggregate_by_index(self, index):
        return self.aggregate_by_datetime(self.flattened_timevector[index])

    def aggregate_by_datetime(self, dt):

        stack_id = self.get_stack_id_by_datetime(dt)
        filenames = [file for file in self.filenames if str(file).endswith(f'{stack_id}.nc')]
        ncvar = Dataset(filenames[0])[self.name]
        shape = []
        if ncvar.dimensions[1] == 'nSCHISM_hgrid_node':
            shape.append(self.hgrid.nodes.values.shape[0])
        elif ncvar.dimensions[1] == 'nSCHISM_hgrid_face':
            shape.append(self.hgrid.elements.index.shape[0])
        else:
            raise Exception(f'Unhandled dimension1 {self.ncvar.dimensions[1]}')

        if len(ncvar.dimensions[1:]) > 1:
            raise Exception(f'Unhandled number of dimensions {self.ncvar.dimensions}')
        
        values = np.full(tuple(shape), np.nan)

        for i, file in enumerate(filenames):
            ncvar = Dataset(filenames[i])[self.name]
            cpu_id = file.name.split('_')[1]
            if ncvar.dimensions[1] == 'nSCHISM_hgrid_node':
                n_local_to_global = self.n_local_to_global[cpu_id]
                local_node_ids = list(n_local_to_global.keys())
                global_node_ids = list(
                    map(lambda x: n_local_to_global[x], local_node_ids))
                idxs = list(map(self.hgrid.nodes.get_index_by_id, global_node_ids))
            elif ncvar.dimensions[1] == 'nSCHISM_hgrid_face':
                e_local_to_global = self.e_local_to_global[cpu_id]
                local_elem_ids = list(e_local_to_global.keys())
                global_elem_ids = list(
                    map(lambda x: e_local_to_global[x], local_elem_ids))
                idxs = list(map(self.hgrid.elements.get_index_by_id, global_elem_ids))
            local_time_index = self.get_local_time_index(dt)
            values[idxs] = ncvar[local_time_index, :]

        return values

    def get_stack_id_by_datetime(self, dt: datetime) -> str:
        output_stack = None
        for stack_id, stack_data in self.stacks.items():
            if dt in stack_data['timevector']:
                output_stack = stack_id
        if output_stack is None:
            raise ValueError(
                'The are no output slices corresponding to requested datetime:'
                f' {dt}. Available datetimes are '
                f'{self.flattened_timevector}.')
        return output_stack

    def get_local_time_index(self, dt: datetime) -> str:
        for stack_data in self.stacks.values():
            if dt in stack_data['timevector']:
                return stack_data['timevector'].index(dt)
        raise ValueError(
                'The are no output slices corresponding to requested datetime:'
                f' {dt}. Available datetimes are '
                f'{self.flattened_timevector}.')

    def __iter__(self):
        return self

    def __next__(self):
        self.current_step += 1
        return self.values

    @property
    def values(self):
        return self.aggregate(self.current_step)

    @property
    def current_step(self):
        return self._current_step

    @current_step.setter
    def current_step(self, current_step: int):
        if current_step is None:
            current_step = self.timesteps[0]
        self._current_step = current_step

    # @property
    # def ncvar(self):
    #     return self.stacks

    # @property
    # def shape(self):
    #     if not hasattr(self, '_shape'):
    #         shape = [self.hgrid.values.shape[0]]
    #         if self.rank is not None:
    #             shape.append(self.rank)
    #         if self.nvrt is not None:
    #             shape.append(self.nvrt)
    #         self._shape = tuple(shape)
    #     return self._shape


class OutputsCollector:

    surface_output_vars = [
        IofWetdryVariables,
        IofZcorVariables,
        IofHydroVariables,
        IofDvdVariables,
        IofWwmVariables,
        IofGenVariables,
        IofAgeVariables,
        IofSedVariables,
        IofEcoVariables,
        IofIcmVariables,
        IofCosVariables,
        IofFibVariables,
        IofSed2dVariables,
        IofMarshVariables,
        IofIceVariables,
        IofAnaVariables,
    ]

    def __init__(self, path: Union[str, os.PathLike]):
        self.path = pathlib.Path(path)

        if not self.path.exists():
            raise ValueError(f'Specified directory {self.path} does not exist!')

        self.filenames = sorted(
            self.path.glob(r'schout_[0-9][0-9][0-9][0-9]_*.nc'))
        # self.fields = cf.read(self.filenames)
        nodes = {}
        elements = {}
        for file in sorted(self.path.glob(r'local_to_global_[0-9][0-9][0-9][0-9]')):
            with open(file) as f:
                ns_global, ne_global, np_global, nvrt, nproc, ntracers, \
                    T, S, GEN, AGE, SED3D, EcoSim, ICM, CoSINE, Feco, \
                    TIMOR, FABM, DVD = f.readline().split()
                f.readline()
                # elements
                ne_local = int(f.readline())
                e_local_to_global = {}
                for i in range(ne_local):
                    local_element_id, global_element_id = f.readline().split()
                    e_local_to_global[local_element_id] = global_element_id
                # points
                np_local = int(f.readline())
                n_local_to_global = {}
                for i in range(np_local):
                    local_node_id, global_node_id = f.readline().split()
                    n_local_to_global[local_node_id] = global_node_id
                # sides
                ns_local = int(f.readline())
                s_local_to_global = {}
                for i in range(ns_local):
                    local_side_id, global_side_id = f.readline().split()
                    s_local_to_global[local_side_id] = global_side_id
                f.readline()  # Header:
                line = f.readline().split()
                if len(line) != 5:
                    line.extend(f.readline().split())
                self.start_year, self.start_month, self.start_day, \
                    self.start_hour, self.utc_start = line
                # nrec, dtout, nspool, nvrt, kz, h0, h_s, h_c, theta_b, \
                #     theta_f, ics = f.readline().split()
                # f.readline()  # (ztot(k),k=1,kz-1),(sigma(k),k=1,nvrt-kz+1)
                # f.readline()  # repeat ne_local, np_local
                _ne_local = None
                _np_local = None
                while _ne_local != ne_local and _np_local != np_local:
                    line = f.readline().split()
                    _np_local = int(float(line[0]))
                    _ne_local = int(float(line[1]))
                for i in range(np_local):
                    x, y, z, flag = map(float, f.readline().split())
                    nodes.setdefault(
                        n_local_to_global[str(i+1)], ((x, y), -z))
                for i in range(ne_local):
                    eids = f.readline().split()[1:]
                    elements.setdefault(
                        e_local_to_global[str(i+1)],
                        list(map(lambda x: n_local_to_global[x], eids)))
            nproc_id = file.name.split('local_to_global_')[-1]
            self.n_local_to_global.setdefault(nproc_id, n_local_to_global)
            self.e_local_to_global.setdefault(nproc_id, e_local_to_global)
            self.s_local_to_global.setdefault(nproc_id, s_local_to_global)
        # sort (not necessary but done for consistency checks)
        nodes = {str(i+1): nodes[str(i+1)] for i in range(len(nodes))}
        elements = {str(i+1): elements[str(i+1)] for i in range(len(elements))}
        self.hgrid = Gr3(nodes=nodes, elements=elements, crs='epsg:4326')
        for vargroup in self.surface_output_vars:
            for vartype in vargroup:
                nc = Dataset(self.filenames[0])
                if vartype.name in nc.variables:
                    setattr(
                        self,
                        vartype.name,
                        OutputVariableCombiner(
                            self,
                            vartype.name,
                            start_step=0
                        ))

    @property
    def n_local_to_global(self):
        if not hasattr(self, '_n_local_to_global'):
            self._n_local_to_global = {}
        return self._n_local_to_global

    @property
    def s_local_to_global(self):
        if not hasattr(self, '_s_local_to_global'):
            self._s_local_to_global = {}
        return self._s_local_to_global

    @property
    def e_local_to_global(self):
        if not hasattr(self, '_e_local_to_global'):
            self._e_local_to_global = {}
        return self._e_local_to_global

    @property
    def stacks(self) -> Dict:
        if not hasattr(self, "_stacks"):
            self._stacks = {}
            for file in self.filenames:
                stack_id = file.name.split("_")[2].split('.')[0]
                nc = Dataset(file)
                self._stacks.setdefault(stack_id, {})
                self._stacks[stack_id].setdefault(
                    'timevector',
                    [self.start_date + timedelta(seconds=float(x))
                     for x in nc['time'][:]]
                    )
        return self._stacks

    @property
    def flattened_timevector(self):
        if not hasattr(self, "_flattened_timevector"):
            flattened_timevector = []
            for stack_id, stack_data in self.stacks.items():
                flattened_timevector.extend(stack_data['timevector'])
            self._flattened_timevector = flattened_timevector
        return self._flattened_timevector

    @property
    def start_date(self):
        if not hasattr(self, "_start_date"):
            start_hour = int(float(self.start_hour))
            self._start_date = datetime(
                int(self.start_year),
                int(self.start_month),
                int(self.start_day),
                start_hour, (start_hour*60) % 60,
                (start_hour*3600) % 60)
            self._start_date -= timedelta(hours=float(self.utc_start))
        return self._start_date


class CombinedOutputs:

    def __init__(self, resource):
        self.dataset = xarray.open_mfdataset(resource)
