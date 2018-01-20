# reading spectra file
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')

class SpectraReader:
    def __init__(self):
        self._n_voxels = None
        self._spectra_matrix = None
        self._n_dims = None

    @property
    def spectra_matrix(self):
        return self._spectra_matrix

    @property
    def n_voxels(self):
        return self._n_voxels

    @property
    def n_dims(self):
        return self._n_dims

    def read(self, filename):
        # todo add value checking
        dfe = pd.read_csv(filename, comment="#", header=None)
        spectra_values = dfe.values[:, 2:]
        voxel_idx = np.unique(spectra_values[:, 0])
        fluence = spectra_values[:, 1]
        self._n_dims = sum(spectra_values[:, 0] < 1)
        self._n_voxels = len(voxel_idx)
        # test reshape index
        fluence_map = fluence.reshape(self.n_voxels, self.n_dims, order="F")
        self._spectra_matrix = fluence_map
        #
        return fluence_map

    def plot_spectra(self, at, voxel_width, energy_bins):
        # plotting spectra
        Nz = self.n_voxels
        z = np.arange(Nz)
        # extent x,y,z

        zz = (-Nz + 1 + 2 * z) * voxel_width / 2
        # first dose value at voxel_z_Width / 2
        zz += zz.max() + voxel_width / 2

        # nearest index at position in mm
        idx_at = np.argmin(abs(zz - at))
        plt.plot(energy_bins, self.spectra_matrix[idx_at], label=str(zz[idx_at]) + " mm")
        plt.title("Energy spectra")
        plt.xlabel("Energy [MeV]")
        plt.legend()


if __name__ == '__main__':
    spectra_file = '/home/victor/Dropbox/victorgabr-geant4/build-SimulationBase-GEANT4_Clang-Release/PromptGammaSpectra.csv'
    spectra_reader = SpectraReader()
    spectra_reader.read(spectra_file)
    # plot proton spectra at
    pos = 50  # mm
    voxel_width = 2  # mm
    # energy in keV
    energy_bins = np.linspace(0, 10000, spectra_reader.n_dims)
    spectra_reader.plot_spectra(pos, voxel_width, energy_bins)
    plt.show()

    plt.style.use("ggplot")
    # test gamma histogram
    plt.plot(energy_bins / 1000, spectra_reader.spectra_matrix[0])
    plt.title("Prompt γ-ray spectra around 175 MeV proton Bragg peak in water ")
    plt.xlabel("Energy [MeV]")
    plt.ylabel("γ-ray count [a.u]")
    plt.show()
