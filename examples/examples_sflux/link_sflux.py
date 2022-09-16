from datetime import datetime, timedelta
import os

import numpy as np

if __name__ == "__main__":
    '''
    This script is used to link sflux
    srcdir: data directory, should have subfolder named as "20220201", "20220202"
    dstdir: path of run/sflux directory
    '''

    srcdir = '/sciclone/schism10/lcui01/schism20/ICOGS/ICOGS3D/Forecast/sflux'
    dstdir = '/sciclone/schism10/lcui01/schism20/ICOGS/ICOGS3D/Forecast/run/sflux'

    vars=["air", "prc", "rad"]

    startdate = datetime(2022, 2, 1)
    rnday = timedelta(days=5)

    timevector = np.arange(startdate, rnday, timedelta(days=1)).astype(datetime)
    print(timevector)

    for i, date in enumerate(timevector):
        srcgfs = f"{srcdir}/{date.strftime('%Y%m%d')}/gfs_{date.strftime('%Y%m%d%H')}.nc"
        srchrrr = f"{srcdir}/{date.strftime('%Y%m%d')}/hrrr_{date.strftime('%Y%m%d%H')}.nc"

        for var in vars:

            dstgfs = f"{dstdir}/sflux_{var}_1.{str(i+1).zfill(4)}.nc"
            os.symlink(srcgfs, dstgfs)

            dsthrrr = f"{dstdir}/sflux_{var}_2.{str(i+1).zfill(4)}.nc"
            os.symlink(srchrrr, dsthrrr)
