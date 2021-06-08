from abc import ABC, abstractmethod
from datetime import timedelta
import pathlib

from netCDF4 import Dataset


class MOD_3D(ABC):

    def __init__(self, hgrid, vgrid):
        self.hgrid = hgrid
        self.vgrid = vgrid

    def write(
            self,
            path,
            start_date,
            rnday,
            output_interval=timedelta(days=1),
            overwrite: bool = False,
            progress_bar: bool = True,
    ):  
        path = pathlib.Path(path)
        if path.exists() and overwrite is not True:
            raise IOError(f'File {path} exists and overwrite is not True.')

        # file_is_not_needed = True
        timevec = None
        for boundary in self.hgrid.boundaries.open.itertuples():
            obj = getattr(boundary, self.bctype)
            if obj is not None:
                bctype = getattr(obj, self.bctype)
                if bctype == 4:
                    # ds = getattr(obj.data_source, self.name)
                    datasets = obj.data_source.get_datasets(
                                start_date,
                                rnday,
                                output_interval
                                )
                    timevec = range(len(datasets))

        if timevec is None:
            return

        nOpenBndNodes = 0
        for boundary in self.hgrid.boundaries.open.itertuples():
            nOpenBndNodes += len(boundary.indexes)

        dst = Dataset(path, 'w', format='NETCDF4')
        # dimensions
        dst.createDimension('nOpenBndNodes', nOpenBndNodes)
        dst.createDimension('one', 1)
        dst.createDimension('time', None)
        dst.createDimension('nLevels', self.vgrid.nvrt)
        dst.createDimension('nComponents', self.nComponents)

        # variables
        dst.createVariable('time_step', 'f', ('one',))
        dst['time_step'][:] = int(output_interval.total_seconds())
        dst.createVariable('time', 'f', ('time',))
        dst['time'][:] = timevec
        dst.createVariable('time_series', 'f',
                           ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
        offset = 0
        for boundary in self.hgrid.boundaries.open.itertuples():
            obj = getattr(boundary, self.bctype)
            if obj is not None:
                bctype = getattr(obj, self.bctype)
                if bctype == 4:
                    obj.data_source.put_ncdata(
                        self.hgrid,
                        self.vgrid,
                        boundary,
                        dst,
                        start_date,
                        rnday,
                        overwrite=overwrite,
                        offset=offset,
                        output_interval=output_interval,
                        pixel_buffer=10,
                        progress_bar=progress_bar
                    )
            else:
                self.put_null_boundary_data(dst, len(boundary.indexes))
            offset += len(boundary.indexes)

    def put_null_boundary_data(self, dst, np):
        raise NotImplementedError('Must write null data.')

    @property
    @abstractmethod
    def nComponents(self):
        pass

    @property
    @abstractmethod
    def bctype(self):
        pass

    # @property
    # @abstractmethod
    # def name(self):
    #     pass


class TEM_3D(MOD_3D):

    @property
    def nComponents(self):
        return 1

    @property
    def bctype(self):
        return 'itetype'

    # @property
    # def name(self):
    #     return 'temperature'


class SAL_3D(MOD_3D):

    @property
    def nComponents(self):
        return 1

    @property
    def bctype(self):
        return 'isatype'

    # @property
    # def name(self):
    #     return 'salinity'
