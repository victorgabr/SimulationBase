import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# plt.style.use('ggplot')
class Geant4Reader:
    def __init__(self):
        self.size = np.zeros(3)
        self.voxels = np.zeros(3)
        self.voxel_width = np.zeros(3)
        self.xx = []
        self.yy = []
        self.zz = []
        self.matrix3D = []

    def read(self, filename):
        self.read_header(filename)
        self.read_matrix(filename)
        self.calc_coordinates()

    def read_header(self, filename):
        # GETTING scoring sizes
        with open(filename) as f:
            txt = f.readlines()
            self.size = np.asarray(eval(txt[2].split(":")[-1]))
            self.voxels = np.asarray(eval(txt[3].split(":")[-1]))
            # TODO add unit specification
        # TODO add unit tests
        # GEANT4 specifies volumes in half sizes
        self.voxel_width = self.size / self.voxels * 2

    def read_matrix(self, filename):
        # reading output file
        df = pd.read_csv(filename, skiprows=5, header=None)
        vals = df.values
        tmp = np.zeros((self.voxels[0], self.voxels[1], self.voxels[2]))
        for row in vals:
            r, c, d = row[:3].astype(int)
            tmp[r, c, d] = row[3]

        self.matrix3D = tmp

        return tmp

    def calc_coordinates(self):
        # calculating coordinates
        x = np.arange(self.voxels[0])
        y = np.arange(self.voxels[1])
        z = np.arange(self.voxels[2])
        # extent x,y,z
        self.xx = (-self.voxels[0] + 1 + 2 * x) * self.voxel_width[0] / 2
        self.yy = (-self.voxels[1] + 1 + 2 * y) * self.voxel_width[1] / 2
        self.zz = (-self.voxels[2] + 1 + 2 * z) * self.voxel_width[2] / 2


# plotting values
if __name__ == '__main__':
    # reader instance
    kerma_reader = Geant4Reader()
    arq = arq = '/home/victor/Dropbox/victorgabr-geant4/build-BrachySourceKerma-GEANT4_Clang-Debug/KermaDeposition.csv'
    kerma_reader.read(arq)

    dose_map = kerma_reader.matrix3D
    xx = kerma_reader.xx
    yy = kerma_reader.yy
    zz = kerma_reader.zz
    numberOfVoxel_x = kerma_reader.voxels[0]
    numberOfVoxel_y = kerma_reader.voxels[1]

    plt.set_cmap("nipy_spectral")
    # XY plane - axial
    z_idx = 0
    im = dose_map[:, :, z_idx]
    plt.imshow(im, interpolation="bicubic", extent=[xx.min(), xx.max(), yy.min(), yy.max()])
    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
    plt.title("Depth z %s mm" % str(zz[z_idx]))

    # plotting values
    # xy plane
    plt.figure()
    y_idx = int(numberOfVoxel_y / 2)
    x_proj = dose_map.sum(axis=0)

    plt.imshow(x_proj, interpolation="bicubic", extent=[zz.min(), zz.max(), xx.min(), xx.max()])
    plt.xlabel("z [mm]")
    plt.ylabel("x [mm]")
    plt.title("X-axis projection")

    # plotting values
    # zx plane
    plt.figure()
    x_idx = int(numberOfVoxel_x / 2)
    x_proj = dose_map.sum(axis=1)
    plt.imshow(x_proj, interpolation="bicubic", extent=[zz.min(), zz.max(), yy.min(), yy.max()])
    plt.xlabel("z [mm]")
    plt.ylabel("y [mm]")
    # plt.title("x  %s mm" % str(xx[x_idx]))
    plt.title("Y-axis projection")

    # # plot bragg curve
    # bragg_curve = dose_map.sum(axis=1).sum(axis=0)
    # bragg_curve /= bragg_curve.max()
    # bragg_curve *= 100
    # plt.figure()
    # plt.plot(zz, bragg_curve)
    # plt.xlabel("Z-Depth [mm]")
    # plt.ylabel("PDD(%)")
    # plt.title("Bragg curve")
    # plt.show()
