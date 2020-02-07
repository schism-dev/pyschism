from pyschism.mesh.mesh import Mesh


class Bctides:

    def __init__(self, mesh):
        self.mesh = mesh

    def write(self, path, overwrite=False):
        print(self.bctides)
        raise NotImplementedError

    @property
    def mesh(self):
        return self.__mesh

    @mesh.setter
    def mesh(self, mesh):
        assert isinstance(mesh, Mesh)
        self.__mesh = mesh

    @property
    def bctides(self):
        print(self.mesh)
        # f = f"{self.tidal_forcing.forcing_start_date}\n"
        # f += f"{self.tidal_forcing.ntip} "
        # f += f"{self.tidal_forcing.cutoff_depth}\n"
        # active = self.tidal_forcing.get_active_potential_constituents()
        # for constituent in active:
        #     forcing = self.tidal_forcing(constituent)
        #     f += f'{constituent} \n'
        #     f += f'{forcing[0]:G} '
        #     f += f"{forcing[1]:G} "
        #     f += f'{forcing[2]:G} '
        #     f += f'{forcing[3]:G} '
        #     f += f'{forcing[4]:G}'
        #     f += '\n'
        # f += f'{self.tidal_forcing.nbfr:d}\n'
        # active = self.tidal_forcing.get_active_forcing_constituents()
        # for constituent in active:
        #     f += f'{constituent} \n'
        #     f += f"{forcing[2]:G} "
        #     f += f'{forcing[3]:G} '
        #     f += f'{forcing[4]:G}'
        #     f += '\n'
        # f += f"{len(self.mesh.ocean_boundaries)}\n"
        # for ftype, boundary in self.mesh.open_boundaries.items():
        #     f += f"{len(boundary)} "
        #     f += 
        # # for constituent in self.get_active_constituents():
        # #     f += f'{constituent}\n'
        # #     for boundary in self.mesh.ocean_boundaries:
        # #         vertices = self.mesh.hgrid.get_xy(crs='EPSG:4326')[boundary, :]
        # #         amp, phase = self.tpxo(constituent, vertices)
        # #         for i in range(len(vertices)):
        # #             f += f'{amp[i]:.8e} {phase[i]:.8e}\n'
        # return f
