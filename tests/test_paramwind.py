from datetime import timedelta, datetime
from pathlib import Path

from pyschism.mesh import Hgrid
from pyschism.driver import ModelConfig
from pyschism.forcing.nws import BestTrackForcing


def test_paramwind():

    hgrid = Hgrid.open(
        'https://raw.githubusercontent.com/geomesh/test-data/main/NWM/hgrid.ll'
    )

    meteo = BestTrackForcing(storm='Florence2018')

    config = ModelConfig(
        hgrid=hgrid,
        vgrid=None,
        fgrid=None,
        iettype=None,
        ifltype=None,
        nws=meteo,
    )

    driver = config.coldstart(
        start_date=datetime(2018, 9, 8),
        end_date=datetime(2018, 9, 18),
        timestep=timedelta(seconds=150),
        nspool=24,
    )

    driver.write(Path(__file__).parent / 'paramwind', overwrite=True)

if __name__ == '__main__':
    test_paramwind()
