#! /usr/bin/env python
import unittest
import os


class OutputsApiTestCase(unittest.TestCase):

    def test_outputs_api(self):
        test_data = f'{os.getenv("ICOGS3D_HGRID")}'
        # test_data = "/sciclone/pscr/lcui01/ICOGS3D_dev/outputs_run3c"
        outputs = OutputCollection(test_data)
        # anim = outputs.elev.animation(vmin=-3, vmax=3, fps=6)
        from time import time
        start_total = time()
        print('Begin time aggregation.')
        for dt in outputs.elev.flattened_timevector:
            start_local = time()
            outputs.elev.aggregate(dt)
            print(f'took {time()-start_local} to aggregate elev for time {dt}.')
        print(f'Total aggregation time for variable elev took {time()-start_total}')


if __name__ == '__main__':
    unittest.main()
