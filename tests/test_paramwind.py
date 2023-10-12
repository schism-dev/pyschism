import tempfile
from datetime import timedelta, datetime
from pathlib import Path
import unittest

from stormevents.nhc import VortexTrack

from pyschism.mesh import Hgrid
from pyschism.driver import ModelConfig
from pyschism.forcing.nws import BestTrackForcing

class ModelConfigurationTestCase(unittest.TestCase):
    def test_paramwind_from_stormname(self):

        hgrid = Hgrid.open(
            'https://raw.githubusercontent.com/geomesh/test-data/main/NWM/hgrid.ll'
        )

        meteo = BestTrackForcing(storm='Florence2018')

        config = ModelConfig(
            hgrid=hgrid,
            vgrid=None,
            fgrid=None,
            flags=[[3, 3, 0, 0]],
            nws=meteo,
        )

        driver = config.coldstart(
            start_date=datetime(2018, 9, 8),
            end_date=datetime(2018, 9, 18),
            timestep=timedelta(seconds=150),
            nspool=24,
        )

        with tempfile.TemporaryDirectory() as dn:
            tmpdir = Path(dn)
            driver.write(tmpdir / 'paramwind', overwrite=True)

    def test_paramwind_from_file(self):
        hgrid = Hgrid.open(
            'https://raw.githubusercontent.com/geomesh/test-data/main/NWM/hgrid.ll'
        )

        track = VortexTrack.from_storm_name('Florence', 2018)


        with tempfile.TemporaryDirectory() as dn:
            tmpdir = Path(dn)
            track.to_file(tmpdir / 'track.dat')

            meteo = BestTrackForcing.from_nhc_bdeck(nhc_bdeck=tmpdir / 'track.dat')
            config = ModelConfig(
                hgrid=hgrid,
                vgrid=None,
                fgrid=None,
                flags=[[3, 3, 0, 0]],
                nws=meteo,
            )

            driver = config.coldstart(
                start_date=datetime(2018, 9, 8),
                end_date=datetime(2018, 9, 18),
                timestep=timedelta(seconds=150),
                nspool=24,
            )

            driver.write(tmpdir / 'paramwind', overwrite=True)

if __name__ == '__main__':
    unittest.main()
