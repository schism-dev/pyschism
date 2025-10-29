"""
Usage:
- change parameters in the main block at the bottom and run this script to generate sflux files from ERA5 data.
- if the download fails because the requested data is too large, increase n_chunks to split the task into smaller parts.
- check the docstring of gen_sflux_era5() for parameter details.
"""

import os
import gc
import psutil
from datetime import datetime, timedelta
from time import time
from pathlib import Path

from pyschism.forcing.nws.nws2.era5 import ERA5

from pyschism.mesh.hgrid import Hgrid
from matplotlib.transforms import Bbox

import re


def _replace_number(fname, n):
    # replace ".0001.nc" → ".1.nc" (no leading zeros)
    return re.sub(r'\.(\d{4})\.nc$', f'.{n}.nc', fname)


def _split_time(start, end, n):
    """
    Split the time period from start to end into n parts.
    """
    total_days = (end - start).days
    base = total_days // n
    remainder = total_days % n

    parts = []
    current = start
    for i in range(n):
        # Distribute the remainder among the first few parts
        length = base + (1 if i < remainder else 0)
        parts.append((current, length))
        current += timedelta(days=length)
    return parts


def gen_sflux_era5(
    startdate=datetime(2020, 12, 1), rnday=60, n_chunks=1,
    bbox=None, hgrid: Hgrid = None,
):
    """
    Generate sflux files from ERA5 data for a given period by splitting the task into n_chunks.
    Individual sflux files are symlinked into a single sflux directory.

    Parameters:
    - startdate: datetime, start date of the data
    - rnday: int, number of days to generate data for
    - n_chunks: int, number of chunks to split the data generation into
    - bbox: list of float, [[min_lon, min_lat], [max_lon, max_lat]], optional bounding box for the data
    - hgrid: Hgrid object, optional hgrid to derive bounding box from

    Either bbox or hgrid_object must be provided to define the spatial extent.
    If both are provided, bbox takes precedence.

    Outputs:
    - sflux/sflux_inputs.txt: configuration file for SCHISM
    - sflux/sflux_*.nc: symlinked sflux files for air, rad, and prc variables
    - sflux_era5_*/: directories containing the original sflux files for each chunk

    Note: the outdir is set to the current directory './', setting it to another path may require
      modifying deeper level code in ERA5DataInventory, which strips the path to only keep the last part.
    """

    t0 = time()

    # sanitize inputs
    if bbox is None and hgrid is None:
        raise ValueError("Either bbox or hgrid_object must be provided.")
    if bbox is None:
        bbox = hgrid.get_bbox('EPSG:4326', output_type='bbox')
        print(f"Derived bbox from hgrid: {bbox}")
    else:
        bbox = Bbox(bbox)
        print(f"Using provided bbox: {bbox}")

    outdir = Path('./')

    if any(outdir.glob('sflux*')):
        raise FileExistsError(
            f"sflux files or related directories already exist in {outdir}. "
            "Please remove them before running this script.")

    # split time periods and generate sflux files
    periods = _split_time(start=startdate, end=(startdate + timedelta(days=rnday)), n=n_chunks)
    process = psutil.Process(os.getpid())

    for i, [period_start, period_days] in enumerate(periods):
        print(f'\nGenerating sflux data from {period_start} for {period_days} days...')
        this_out_dir = outdir / f'sflux_era5_{i}'
        os.makedirs(this_out_dir, exist_ok=True)
        er = ERA5()
        if period_days < 1:
            continue  # skip zero-day periods
        er.write(
            outdir=this_out_dir,
            start_date=period_start, rnday=period_days-1,
            air=True, rad=True, prc=True,
            bbox=bbox, output_interval=1,  # raw data is hourly, so output_interval=1 means keeping all hours
            overwrite=True, tmpdir=this_out_dir
        )  # default tmpdir is system's tmp (/tmp), it is possbile /tmp has no enough space for large dataset.
        del er
        gc.collect()
        mem = process.memory_info().rss / 1024**2
        print(f"[Iter {i}] Memory usage after GC: {mem:.1f} MB")

    # consolidate sflux files
    # write sflux_inputs.txt
    os.makedirs(outdir / 'sflux', exist_ok=True)
    with open(outdir / 'sflux' / 'sflux_inputs.txt', 'w') as f:
        f.write("&sflux_inputs\n/\n")
    # link individual sflux periods to sflux directory
    file_id = 1
    for i, [_, _] in enumerate(periods):
        for var in ['air', 'rad', 'prc']:
            nc_files = sorted(Path(outdir / f'sflux_era5_{i}').glob(f'sflux*{var}*.nc'))
            for k, nc_file in enumerate(nc_files):
                src = f'../sflux_era5_{i}/{nc_file.name}'
                dst = outdir / 'sflux' / _replace_number(f'sflux_{var}.{file_id:04d}.nc', k + file_id)
                print(f"sym-linking: src: {src}, dest: {dst}")
                os.symlink(src, dst)
        file_id += len(nc_files)

    print(f'\nIt took {(time()-t0)/60} minutes to generate {rnday} days')


if __name__ == "__main__":
    hgrid = Hgrid.open('./hgrid.gr3', crs='EPSG:4326')

    # example using hgrid to define bbox
    gen_sflux_era5(
        startdate=datetime(2020, 12, 1), rnday=10, n_chunks=2,
        hgrid=hgrid,
    )

    # example using explicit bbox
    gen_sflux_era5(
        startdate=datetime(2020, 12, 1), rnday=3, n_chunks=1,
        bbox=[[-76.0, 34.0], [-70.0, 40.0]],
    )
