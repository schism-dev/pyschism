from datetime import timedelta
import pathlib

# import geopandas as gpd
from netCDF4 import Dataset


class Elev2D:

    def __init__(self, hgrid):
        self.hgrid = hgrid

    def write(
            self,
            elev2D,
            start_date,
            rnday,
            output_interval=timedelta(days=1),
            overwrite: bool = False,
    ):
        elev2D = pathlib.Path(elev2D)
        if elev2D.exists() and overwrite is not True:
            raise IOError(f'File {elev2D} exists and overwrite is not True.')

        # file_is_not_needed = True
        timevec = None
        for boundary in self.hgrid.boundaries.open.itertuples():
            if boundary.iettype is not None:
                if boundary.iettype.iettype in [4, 5]:
                    ds = boundary.iettype.data_source
                    datasets = ds.get_datasets(
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

        dst = Dataset(elev2D, 'w', format='NETCDF4')
        # dimensions
        dst.createDimension('nOpenBndNodes', nOpenBndNodes)
        dst.createDimension('one', 1)
        dst.createDimension('time', None)
        dst.createDimension('nLevels', 1)
        dst.createDimension('nComponents', 1)

        # variables
        dst.createVariable('time_step', 'f', ('one',))
        dst['time_step'][:] = int(output_interval.total_seconds())
        dst.createVariable('time', 'f', ('time',))
        dst['time'][:] = timevec
        dst.createVariable('time_series', 'f',
                           ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
        offset = 0
        for boundary in self.hgrid.boundaries.open.itertuples():
            if boundary.iettype is not None:
                if boundary.iettype.iettype in [4, 5]:
                    boundary.iettype.data_source.put_ncdata(
                        boundary, dst, start_date, rnday, overwrite=overwrite,
                        offset=offset, output_interval=output_interval,
                        pixel_buffer=10)
            else:
                self.put_null_boundary_data(dst, len(boundary.indexes))
            offset += len(boundary.indexes)

    def put_null_boundary_data(self, dst, np):
        raise NotImplementedError('Must write null data.')
