import os
from time import time
import pathlib
from typing import Dict, Union
import glob

import numpy as np
import xarray as xr

from pyschism.mesh.base import Gr3

def combine(var, shape, l2g, name):
    values = np.full(tuple(shape), np.nan)

    local_ids = list(l2g.keys())
    for i, data in enumerate(var):
        cpu_id = str(i).zfill(6) 
        #print(cpu_id)
        n_local_to_global = l2g[cpu_id]
        #print(n_local_to_global[0])
        local_ids = list(n_local_to_global.keys())
        global_idxs = list(
            map(lambda x: int(n_local_to_global[x])-1, local_ids))
        #print(global_idxs[0])
        values[global_idxs] = data
    return values
    

def combine_(dst, var, shape, l2g, name):
    #if len(dst[0][var].shape) < 4:
    out = []
    for i in range(len(dst)):
        out.append(dst[i][var])
    r = combine(out, shape, l2g, name)
    return (xr.DataArray(r, dims = list(dst[0][var].dims), name = var))

class CombineOutputs:

    def __init__(self, path: Union[str, os.PathLike]):
        self.path = pathlib.Path(path)

        if not self.path.exists():
            raise ValueError(f'Directory {self.path} does not exist.')

        nodes = {}
        elements = {}
        for ifile in sorted(self.path.glob(r'local_to_global_[0-9][0-9][0-9][0-9][0-9][0-9]')):
            with open(ifile) as f:
                ns_global, ne_global, np_global, nvrt, nproc, ntracers, \
                    T, S, GEN, AGE, SED3D, EcoSim, ICM, CoSINE, Feco, \
                    TIMOR, FARM, DVD = f.readline().split()
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
                #old schism print to multiple lines not just one line
                if len(line) != 5:
                    line.extend(f.readline().split())
                self.start_year, self.start_month, self.start_day, \
                    self.start_hour, self.utc_start = line
                nrec, dtout, nspool, nvrt, kz, h0, h_s, h_c, theta_b, \
                    theta_f, ics = f.readline().split()
                #In the old version of schism, ztot was written to nvrt lines
                for i in np.arange(int(nvrt)):
                    f.readline()  # (ztot(k),k=1,kz-1),(sigma(k),k=1,nvrt-kz+1)
                f.readline()  # (ztot(k),k=1,kz-1),(sigma(k),k=1,nvrt-kz+1)
                
                #_ne_local = None
                #_np_local = None
                #while _ne_local != ne_local and _np_local != np_local:
                #    line = f.readline().split()
                #    _np_local = int(float(line[0])
                #    _ne_local = int(float(line[1])
                for i in range(np_local):
                    x, y, z, flag = map(float, f.readline().split())
                    nodes.setdefault(
                        n_local_to_global[str(i+1)], ((x, y), -z))
                for i in range(ne_local):
                    eids = f.readline().split()[1:]
                    elements.setdefault(
                        e_local_to_global[str(i+1)],
                        list(map(lambda x: n_local_to_global[x], eids)))
                nproc_id = ifile.name.split('local_to_global_')[-1]
                self.e_local_to_global.setdefault(nproc_id, e_local_to_global)
                self.n_local_to_global.setdefault(nproc_id, n_local_to_global)
                self.s_local_to_global.setdefault(nproc_id, s_local_to_global)
            #nodes = {str(i+1): nodes[str(i+1)] for i in range(len(nodes))}
            #elements = {str(i+1): elements[str(i+1)] for i in range(len(elements))}
        self.hgrid = Gr3(nodes=nodes, elements=elements, crs='epsg:4326')

    def hotstart(self, it=None):
        self.filenames = sorted(
            self.path.glob(r'hotstart_*_{}.nc'.format(it)))
        dst = []
        for i in range(len(self.filenames)):
            dst.append(xr.open_dataset(self.filenames[i]))
    
        #create dataset
        side = []
        node = []
        elem = []
        one = []
        #variables = ['eta2', 'su2', 'sv2']
        #variables = ['eta2', 'we', 'su2', 'tr_el', 'time', 'it', 'ifile', 'nsteps_from_cold']
        #for var in variables:
        for var in dst[0].variables:
            t0 = time()
            shape = []
            if 'nResident_elem' in dst[0][var].dims:
                shape.append(self.hgrid.elements.array.shape[0])
                if len(dst[0][var].shape) > 1:
                    for i in range(len(dst[0][var].shape)):
                        if i == 0:
                            continue
                        else:
                            shape.append(dst[0][var].shape[i])
                r = combine_(dst, var, shape, self.e_local_to_global, 'nResident_elem')
                elem.append(r)
            elif 'nResident_node' in dst[0][var].dims:
                shape.append(self.hgrid.nodes.values.shape[0])
                if len(dst[0][var].shape) > 1:
                    for i in range(len(dst[0][var].shape)):
                        if i == 0:
                            continue
                        else:
                            shape.append(dst[0][var].shape[i])
                r = combine_(dst, var, shape, self.n_local_to_global, 'nResident_node')
                node.append(r)
            elif 'nResident_side' in dst[0][var].dims:
                shape.append(self.hgrid.elements.sides.shape[0])
                if len(dst[0][var].shape) > 1:
                    for i in range(len(dst[0][var].shape)):
                        if i == 0:
                            continue
                        else:
                            shape.append(dst[0][var].shape[i])
                r = combine_(dst, var, shape, self.s_local_to_global, 'nResident_side')
                side.append(r)
            else:
                one.append(dst[0][var])
            print(f'It took {time()-t0} seconds to combine var {var} in file[{i}]')

        side = xr.merge(side).rename({'nResident_side': 'side'})
        elem = xr.merge(elem).rename({'nResident_elem': 'elem'})
        node = xr.merge(node).rename({'nResident_node': 'node'})
        one = xr.merge(one).rename({'one': 'one_new', 'it': 'iths'})

        xdat = xr.merge([side, elem, node, one])
        #xdat = xr.merge([node, one])
        hfile = 'hotstart_it={}.nc'.format(it)
        xdat.to_netcdf(f'./{hfile}')


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
