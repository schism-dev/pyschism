#! /usr/bin/env python
import unittest
import tempfile
import pathlib
import numpy as np
from pyschism.mesh import Gmesh


class GmeshTestCase(unittest.TestCase):

    def setUp(self):

        self.coords = {
            123964: (205774.81078231844, 4443752.534623082),
            132659: (206063.5749980295, 4443842.0242193565),
            132660: (205961.9161064203, 4443966.412854593),
            138088: (178078.4624008443, 4377433.342070606),
            138115: (178312.1315102501, 4377490.163494714),
            138116: (178273.18682817701, 4377348.881943488),
            139510: (199462.6522372386, 4436931.032656117),
            140443: (203187.56183790273, 4446165.092384715),
            140461: (205864.1247334437, 4445323.356508309),
            140462: (206119.8985982882, 4445353.242776221),
            141993: (201570.57090928522, 4444661.975748284),
            150761: (199459.23523418795, 4436762.253779167),
            150762: (199623.8859668278, 4436775.07013118),
            151406: (202372.86186895028, 4445742.149428694),
            151407: (202617.69929010718, 4445729.383170995),
            151410: (203247.0970142425, 4445945.949797375),
            151411: (203457.01591852625, 4446031.372848534),
            151426: (206024.5692980432, 4445141.813600075),
            152483: (201388.43549408024, 4444690.197771754),
            152484: (201486.6926371171, 4444556.231103622),
            152535: (199050.21183348674, 4441946.515479384),
            152536: (198972.56145616085, 4442047.694862202),
            159513: (202505.38260327603, 4445570.2157231495),
            160372: (199161.30104333046, 4442076.46013068),
            170399: (190302.56402584296, 4436544.191002601),
            170400: (190470.4690470995, 4436491.426104744),
            173152: (190365.53451060448, 4436388.863005645),
            174367: (211320.6036801472, 4378103.511599774),
            174368: (211279.74506958117, 4378246.304886084),
            175817: (211454.28782265307, 4378197.400590049),
            179719: (152032.81836110973, 4408228.504523149),
            179720: (152099.88111381052, 4408254.124981896),
            179721: (152166.94284954405, 4408279.733195366),
            179722: (152200.33388469755, 4408182.899772616),
            179723: (152233.72288592035, 4408086.066349866),
            179724: (152267.11188714317, 4407989.232927115),
            179725: (152189.0170960274, 4408068.984188688),
            179726: (152110.91823705027, 4408148.744355918),
            179727: (151799.47773146332, 4408468.438477338),
            179728: (151854.9084362822, 4408489.1419008365),
            179729: (151911.0510167384, 4408508.383683562),
            179730: (151967.83631919592, 4408526.235070765),
            179731: (152025.18502036916, 4408542.765081277),
            179732: (151955.0388292032, 4408308.486218568),
            179733: (151995.8699817119, 4408323.924173029),
            179734: (152036.9391040723, 4408338.871203209),
            179735: (152078.2207721591, 4408353.358478906),
            179736: (152119.69159577374, 4408367.4138302915),
            179737: (151877.25828033325, 4408388.457895126),
            179738: (151925.38920899705, 4408406.533036932),
            179739: (151973.996077372, 4408423.627443385),
            179740: (152023.0280371945, 4408439.80234087),
            179741: (152072.43830807143, 4408455.084446354),
            }

        self.triangles = {
            323078: [151410, 151411, 140443],
            323079: [151407, 151406, 159513],
            323080: [151426, 140462, 140461],
            323081: [152484, 141993, 152483],
            323082: [123964, 132659, 132660],
            323083: [152535, 160372, 152536],
            323084: [150762, 139510, 150761],
            323085: [173152, 170400, 170399],
            323086: [174368, 174367, 175817],
            323087: [138088, 138116, 138115],
            323088: [179719, 179726, 179720],
            323089: [179722, 179721, 179720],
            323090: [179720, 179726, 179722],
            323091: [179725, 179724, 179723],
            323092: [179723, 179722, 179725],
            323093: [179725, 179722, 179726],
            323094: [179719, 179733, 179732],
            323097: [179721, 179736, 179735]
          }

        self.quads = {
            323095: [179719, 179720, 179734, 179733],
            323096: [179720, 179721, 179735, 179734],
            323098: [179732, 179733, 179738, 179737],
            323099: [179733, 179734, 179739, 179738],
            323100: [179734, 179735, 179740, 179739],
            323101: [179735, 179736, 179741, 179740],
            323102: [179737, 179738, 179728, 179727],
            323103: [179738, 179739, 179729, 179728],
            323104: [179739, 179740, 179730, 179729],
            323105: [179740, 179741, 179731, 179730],
            }

    def test_init(self):
        gmsh = Gmesh(self.coords, self.triangles, self.quads)
        self.assertIsInstance(gmsh, Gmesh)

    def test_transform_to(self):
        gmsh = Gmesh(self.coords, self.triangles, self.quads, crs=3395)
        gmsh.transform_to(4326)
        self.assertIsInstance(gmsh, Gmesh)

    def test_xy(self):
        gmsh = Gmesh(self.coords, self.triangles, self.quads, crs=3395)
        values = [[x, y] for x, y in self.coords.values()]
        self.assertSequenceEqual(gmsh.xy.tolist(), values)

    def test_default_description(self):
        gmsh = Gmesh(self.coords, self.triangles, self.quads, crs=3395)
        gmsh.description = 'test'
        self.assertEqual(gmsh.description, 'test')

    def test_open_gr3(self):
        nodes = {
            '1': ((0., 0.), -99999.),
            '2': ((.5, 0.), -99999.),
            '3': ((1., 0.), -99999.),
            '4': ((1., 1.), -99999.),
            '5': ((0., 1.), -99999.),
            '6': ((.5, 1.5), -99999.),
            '7': ((.33, .33), -99999.),
            '8': ((.66, .33), -99999.),
            '9': ((.5, .66), -99999.),
            '10': ((-1., 1.), -99999.),
            '11': ((-1., 0.), -99999.),
            }
        elements = {
            '1': ['5', '7', '9'],
            '2': ['1', '2', '7'],
            '3': ['2', '3', '8'],
            '4': ['8', '7', '2'],
            '5': ['3', '4', '8'],
            '6': ['4', '9', '8'],
            '7': ['4', '6', '5'],
            '8': ['5', '10', '11', '1'],
            '9': ['9', '4', '5'],
            '10': ['5', '1', '7']
            }

        boundaries = dict()

        boundaries[None] = {  # "open" boundaries
                0: {'indexes': ['10', '11', '1', '2']},
                1: {'indexes': ['2', '3', '4']}
        }

        boundaries[0] = {  # "land" boundaries
            0: {'indexes': ['4', '6']},
            1: {'indexes': ['6',  '5', '10']}
        }

        boundaries[1] = { # "interior" boundary
            0: {'indexes': ['7', '8', '9', '7']}
        }
        f = "test\n"
        f += f'{len(elements):d} '
        f += f'{len(nodes):d}\n'
        for id, ((x, y), z) in nodes.items():
            f += f"{id} "
            f += f"{x} "
            f += f"{y} "
            f += f"{z}\n"
        for id, geom in elements.items():
            f += f"{id} "
            f += f"{len(geom)} "
            for idx in geom:
                f += f"{idx} "
            f += f"\n"
        if None in boundaries:
            f += f"{len(boundaries[None]):d} "
            f += "! total number of ocean boundaries\n"
            # count total number of ocean boundaries
            _sum = 0
            for bnd in boundaries[None].values():
                _sum += len(bnd['indexes'])
            f += f"{int(_sum):d} ! total number of ocean boundary nodes\n"
            # write ocean boundary indexes
            for i, boundary in boundaries[None].items():
                f += f"{len(boundary['indexes']):d}"
                f += f" ! number of nodes for ocean_boundary_{i}\n"
                for idx in boundary['indexes']:
                    f += f"{idx}\n"
        else:
            f += "0 ! total number of ocean boundaries\n"
            f += "0 ! total number of ocean boundary nodes\n"
        # remaining boundaries
        _cnt = 0
        for key in boundaries:
            if key is not None:
                for bnd in boundaries[key]:
                    _cnt += 1
        f += f"{_cnt:d}  ! total number of non-ocean boundaries\n"
        # count remaining boundary nodes
        _cnt = 0
        for ibtype in boundaries:
            if ibtype is not None:
                for bnd in boundaries[ibtype].values():
                    _cnt += np.asarray(bnd['indexes']).size
        f += f"{_cnt:d} ! Total number of non-ocean boundary nodes\n"
        # all additional boundaries
        for ibtype, boundaries in boundaries.items():
            if ibtype is None:
                continue
            for id, boundary in boundaries.items():
                f += f"{len(boundary['indexes']):d} "
                f += f"{ibtype} "
                f += f"! boundary {ibtype}:{id}\n"
                for idx in boundary['indexes']:
                    f += f"{idx}\n"
        tmpdir = tempfile.TemporaryDirectory()
        gr3 = pathlib.Path(tmpdir.name) / 'gr3.gr3'
        with open(gr3.absolute(), 'w') as h:
            h.write(f)
        msh = Gmesh.open_gr3(gr3.absolute())
        self.assertIsInstance(msh, Gmesh)

    def test_add_existing_boundary_type_raises(self):
        msh = Gmesh(self.coords, self.triangles)
        self.assertRaises(Exception, msh.add_boundary_type, None)

    def test_boundary_type(self):
        msh = Gmesh(self.coords, self.triangles)
        msh.delete_boundary_type(None)

    def test_set_boundary_data_raises_bad_indexes(self):
        msh = Gmesh(self.coords, self.triangles)
        data = ['10', '11', '1', '2']
        self.assertRaises(
            AssertionError,
            msh.set_boundary_data,
            None, 0, data)

    def test_set_boundary_data(self):
        msh = Gmesh(self.coords, self.triangles)
        data = [139510, 140443, 140461, 140462, 141993, 150761]
        msh.set_boundary_data(None, 0, data)

    def test_delete_boundary_data(self):
        msh = Gmesh(self.coords, self.triangles)
        data = [139510, 140443, 140461, 140462, 141993, 150761]
        msh.set_boundary_data(None, 0, data)
        msh.delete_boundary_data(None, 0)

if __name__ == '__main__':
    unittest.main()
