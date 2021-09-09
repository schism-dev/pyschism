from pyschism.mesh.hgrid import Hgrid
from pyschism.mesh.gridgr3 import Shapiro

hgrid=Hgrid.open('/sciclone/schism10/whuang07/NWM/Case1/RUN07b_ZG/hgrid.gr3', crs='epsg:4326')
shapiro=Shapiro.from_binary(hgrid)
regions=['coastal_0.2.cpp.reg', 'coastal_0.5_1.cpp.reg', 'coastal_0.5_2.cpp.reg']
values=[0.2, 0.5, 0.5]
flags=[0, 0, 0]
depth1=-99999
for reg, value, flag in zip(regions, values, flags):
    shapiro.modify_by_region(hgrid, f'Shapiro_python/{reg}', value, depth1, flag)

shapiro.write('shaprio.gr3', overwrite=True)
