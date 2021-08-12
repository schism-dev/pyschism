from pyschism.forcing.nwm import NationalWaterModel, logger as nwm_logger
from pyschism.forcing.base import logger as source_sink_base_logger
from pyschism.mesh import Hgrid


def test_nwm():

    nwm = NationalWaterModel()
    nwm.fetch_data(
        Hgrid.open('https://raw.githubusercontent.com/geomesh/test-data/main/NWM/hgrid.ll')
        )


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
