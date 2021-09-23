import pathlib

from pyschism.forcing.source_sink.nwm import NationalWaterModel, logger as nwm_logger, NWMElementPairings
from pyschism.forcing.source_sink.base import logger as source_sink_base_logger
from pyschism.mesh import Hgrid


for each file:
    my_idxs.append()

pool.starmap (
    func,
    [(nc, my_idxs[i]) for enumerate(each_file)]
)

def test_nwm():

    hgrid = Hgrid.open('https://raw.githubusercontent.com/geomesh/test-data/main/NWM/hgrid.ll')
    nwm = NationalWaterModel()
    cached_pairings = pathlib.Path('pairings.txt')
    nwm.fetch_data(
        hgrid,
        end_date=5.,
        pairings=None if cached_pairings.is_file() is False else NWMElementPairings.from_file(cached_pairings)
        )
    nwm.pairings.make_plot()
    if cached_pairings.is_file() is False:
        nwm.pairings.save('pairings.txt')


if __name__ == "__main__":
    import logging

    logging.basicConfig(
        format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
        force=True,
    )
    logging.captureWarnings(True)

    log_level = logging.DEBUG

    nwm_logger.setLevel(log_level)
    source_sink_base_logger.setLevel(log_level)

    test_nwm()
