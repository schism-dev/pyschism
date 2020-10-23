import pathlib
from typing import Union, List

import cf  # type: ignore[import]

from pyschism.forcing.atmosphere.nws2 import NWS2


def load_sflux(path: Union[str, pathlib.Path]) -> NWS2:
    """Factory method to load a sflux directory.

    Loads the data available from a local sflux directory.

    Args:
        path: Path to the local sflux directory.

    Returns:
        :class:`pyschism.forcing.NWS2`
    """
    if not isinstance(path, pathlib.Path):
        path = pathlib.Path(path)
    if not path.is_dir():
        raise OSError(f'Path provided {str(path)} is not a directory.')
    return NWS2(cf.read(path.glob('*.nc')))


def load_datasets(datasets: List[Union[str, pathlib.Path]],
                  ignore_read_error: bool = True) -> NWS2:
    """Factory functions to load datasets.

    Loads the data available from a list of paths. It will read the entire
    contents of a given directory, but it will not recurse.

    Args:
        datasets: list of paths to files or directories containing files.
    Returns:
        :class:`pyschism.forcing.NWS2`
    """
    for i, dataset in enumerate(datasets):
        if isinstance(dataset, str):
            datasets[i] = pathlib.Path(dataset)
    return NWS2(cf.read(datasets, ignore_read_error=ignore_read_error))


def load_fields(fields: cf.FieldList) -> NWS2:
    """Factory function to generate a :class:`pyschism.forcing.NWS2`


    Args:
        fields: Instance of :class:`cf.FieldList`
    Returns:
        :class:`pyschism.forcing.NWS2`
    """
    if not isinstance(fields, cf.FieldList):
        raise TypeError(f'Argument must be of type {cf.FieldList}')
    return NWS2(fields)


def fetch_storm_meta():
    pass
# def fetch_storm_meta(storm, clip: Union[Polygon, MultiPolygon]) -> namedtuple:
    
#     return namedtuple()()
