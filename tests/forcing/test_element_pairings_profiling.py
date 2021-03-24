#! /usr/bin/env python
from datetime import datetime
import logging
from multiprocessing import Pool, cpu_count
import os
import pathlib
from time import time

from netCDF4 import Dataset
# import numpy as np

from pyschism.mesh import Hgrid
from pyschism.forcing.hydrology import Hydrology
from pyschism.forcing.hydrology.nwm import (
    NWMElementPairings,
    AWSDataInventory,
    streamflow_lookup,
    localize_datetime,
    pivot_time
)


if __name__ == '__main__':

    logging.getLogger().setLevel(logging.INFO)

    logging.info('Read hgrid.')
    start = time()
    hgrid = Hgrid.open(
        os.getenv('NWM_TEST_MESH'),
        # '/ddnas/jreniel/ADCIRC/NetCDF_shinnecock_inlet/fort.14',
        crs='EPSG:4326'
    )
    logging.info(f'Read hgrid took {time() - start}.')

    pairings = NWMElementPairings(hgrid)

    start = time()
    inventory = AWSDataInventory(rnday=0.25)
    logging.info(f'Download AWS data took {time() - start}.')

    nprocs = cpu_count()
    source_indexes, sinks_indexes = inventory.get_nc_pairing_indexes(pairings)
    if nprocs > 1:
        logging.info('Running parallelized streamflow_lookup.')
        start = time()
        with Pool(processes=cpu_count()) as pool:
            sources = pool.starmap(
                streamflow_lookup,
                [(file, source_indexes) for file in inventory.files])
            sinks = pool.starmap(
                streamflow_lookup,
                [(file, sinks_indexes) for file in inventory.files])
        pool.join()
        logging.info(f'streamflow_lookup in parallel took {time() - start}.')

    else:
        logging.info('streamflow_lookup serial test.')
        start = time()
        sources = [streamflow_lookup(file, source_indexes)
                   for file in inventory.files]
        sinks = [streamflow_lookup(file, sinks_indexes)
                 for file in inventory.files]
        logging.info(f'streamflow_lookup in serial took {time() - start}.')

    # Pass NWM data to Hydrology class.
    logging.info('Generating per-element hydrologic timeseries...')
    start = time()
    hydro = Hydrology(pivot_time(), inventory.rnday)
    for i, file in enumerate(inventory.files):
        nc = Dataset(file)
        _time = localize_datetime(datetime.strptime(
            nc.model_output_valid_time,
            "%Y-%m-%d_%H:%M:%S"))
        for j, element_id in enumerate(pairings.sources.keys()):
            hydro.add_data(_time, element_id, sources[i][j], -9999, 0.)
        for k, element_id in enumerate(pairings.sinks.keys()):
            hydro.add_data(_time, element_id, -sinks[i][k])
    logging.info(
        f'Generating per-element hydrologic timeseries took {time() - start}.')

    # aggregate timeseries
    aggregation_radius = 4000.
    logging.info('Aggregating hydrology timeseries/elements using a radius of '
                 f'{aggregation_radius} meters.')
    start = time()
    hydro.aggregate_by_radius(hgrid, aggregation_radius)
    logging.info(f'Aggregating NWM elements took {time() - start}.')

    outdir = pathlib.Path('hydro_test')
    outdir.mkdir(exist_ok=True)
    hydro.write(outdir, overwrite=True)

    # verification plot:
    # gather extreme values
    # import numpy as np
    # source_max = {element_id: -float('inf') for element_id
    #               in hydro.sources.elements}
    # for element_data in hydro.sources.data.values():
    #     for element_id, data in element_data.items():
    #         source_max[element_id] = np.max(
    #             [source_max[element_id], data['flow']])

    # sink_max = {element_id: float('inf') for element_id
    #             in hydro.sinks.elements}
    # for element_data in hydro.sinks.data.values():
    #     for element_id, data in element_data.items():
    #         sink_max[element_id] = np.min([sink_max[element_id], data['flow']])

    # # create gdf based on extreme values
    # import geopandas as gpd
    # data = []
    # for i, element_id in enumerate(source_max):
    #     element_index = hgrid.elements.get_index_by_id(element_id)
    #     data.append({
    #         'geometry': hgrid.elements.gdf.loc[element_index].geometry,
    #         'max_flow': source_max[element_id]
    #     })
    # source_gdf = gpd.GeoDataFrame(data)
    # source_gdf.to_file('hydro_test/sources.shp', driver='ESRI Shapefile')

    # data = []
    # for i, element_id in enumerate(sink_max):
    #     element_index = hgrid.elements.get_index_by_id(element_id)
    #     data.append({
    #         'geometry': hgrid.elements.gdf.loc[element_index].geometry,
    #         'max_flow': sink_max[element_id]
    #     })
    # sink_gdf = gpd.GeoDataFrame(data)
    # sink_gdf.to_file('hydro_test/sinks.shp', driver='ESRI Shapefile')

    # ax = hgrid.elements.gdf.plot(facecolor='none', edgecolor='k')
    # source_gdf.plot(ax=ax, color='r', alpha=0.5)
    # sink_gdf.plot(ax=ax, color='b', alpha=0.5)
    # import matplotlib.pyplot as plt; plt.show()