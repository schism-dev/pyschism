from pyschism.mesh import Hgrid
from pyschism.mesh.gridgr3 import Shapiro

if __name__ == '__main__':
    hgrid = Hgrid.open('/sciclone/schism10/whuang07/NWM/Case1/RUN07b_ZG/hgrid.gr3', crs='EPSG:4326')

    shapiro_max = 0.5
    threshold_slope = 0.5
    depths = [-99999, 20, 50]  # tweaks in shallow waters
    shapiro_vals1 = [0.2, 0.2, 0.05]  # tweaks in shallow waters
    regdir = '/sciclone/data10/lcui01/schism/src/Utility/Pre-Processing/NWM/Shapiro_python'
    regions = [
        f'{regdir}/coastal_0.2.cpp.reg',
        f'{regdir}/coastal_0.5_1.cpp.reg',
        f'{regdir}/coastal_0.5_2.cpp.reg'
    ]  # tweaks in regions, the order matters
    shapiro_vals2 = [0.2, 0.5, 0.5]  # tweaks in regions, the order matters
    i_set_add_s = [0, 0, 0]

    shapiro = Shapiro.slope_filter(
        hgrid, shapiro_vals1, depths, shapiro_max, threshold_slope,
        regions, shapiro_vals2, i_set_add_s, lonc=77.07, latc=24.0)
    shapiro.write('shapiro_pyschism.gr3', overwrite=True)
