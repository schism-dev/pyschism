from datetime import date
import logging
import os
import pathlib
from typing import Union

from pyschism import Mesh
from pyschism.enums import ForecastCycle
from pyschism.logger import get_logger


def generate_forecast(
    date: date,
    cycle: Union[str, ForecastCycle],
    output_directory: os.PathLike,
    hgrid: os.PathLike,
    vgrid: os.PathLike = None,
    fgrid: os.PathLike = None,
    logger: Union[logging.Logger, dict] = None,
):
    """Generate input files for a given forcast cycle.

    This function is the main entrypoint to generate forecast cycles. Mainly,
    it is meant to be called by :class:`pyschism.cmd.forecastd.Forecastd` and
    by the forecast cli command. It can also be used standalone on a script.
    This function

    """

    if logger is None:
        logger = get_logger(__name__)

    elif isinstance(logger, dict):
        logger = get_logger(**logger)

    assert isinstance(logger, logging.Logger)

    outdir = pathlib.Path(output_directory)

    try:
        outdir.mkdir()
    except FileExistsError as err:
        logging.error(f'{err.message}')
        raise err

    logger.info(f'Output files will be written to {str(outdir)}')

    mesh = Mesh.open(
        hgrid=hgrid,
        vgrid=vgrid,
        fgrid=fgrid
    )
