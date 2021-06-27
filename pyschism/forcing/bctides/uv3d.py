from datetime import timedelta
import pathlib

# import geopandas as gpd
from netCDF4 import Dataset


class UV3D:

    def __init__(self, bctides):
        self.bctides = bctides

    def write(
            self,
            uv3d,
            start_date,
            rnday,
            output_interval=timedelta(days=1),
            overwrite: bool = False,
            progress_bar: bool = True,
    ):
        uv3d = pathlib.Path(uv3d)
        if uv3d.exists() and overwrite is not True:
            raise IOError(f'File {uv3d} exists and overwrite is not True.')

        # file_is_not_needed = True
        timevec = None
        for boundary in self.bctides.gdf.itertuples():
            if boundary.ifltype is not None:
                if boundary.ifltype.ifltype in [4, 5]:
                    ds = boundary.ifltype.data_component
                    datasets = ds.get_datasets(
                                start_date,
                                rnday,
                                output_interval
                                )
                    timevec = range(len(datasets))

        if timevec is None:
            return

        nOpenBndNodes = 0
        for boundary in self.bctides.gdf.itertuples():
            nOpenBndNodes += len(boundary.indexes)

        dst = Dataset(uv3d, 'w', format='NETCDF4')
        # dimensions
        dst.createDimension('nOpenBndNodes', nOpenBndNodes)
        dst.createDimension('one', 1)
        dst.createDimension('time', None)
        dst.createDimension('nLevels', self.bctides.vgrid.nvrt)
        dst.createDimension('nComponents', 2)

        # variables
        dst.createVariable('time_step', 'f', ('one',))
        dst['time_step'][:] = int(output_interval.total_seconds())
        dst.createVariable('time', 'f', ('time',))
        dst['time'][:] = timevec
        dst.createVariable('time_series', 'f',
                           ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
        offset = 0
        for boundary in self.bctides.gdf.itertuples():
            if boundary.ifltype is not None:
                if boundary.ifltype.ifltype in [4, 5]:
                    boundary.ifltype.data_component.put_boundary_ncdata(
                        self.bctides.hgrid,
                        self.bctides.vgrid,
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
