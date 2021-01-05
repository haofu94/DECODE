import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import sys
sys.path.append("/Users/haofu/GalaxyProjects/Discrete-STEEL/D-STEEL/dsteel/")
from DM_to_SM import SMHM_Grylls_param #(DM, Params, z, scatter)

SMHM_params = [11.92, 0.58, 0.032, -0.014, 1.64, -0.69, 0.53, 0.03] #PyMorph

plt.ion()

z_array = np.arange(0., 6., 0.1)

Mstar = np.arange(6, 12.5, 0.05)

Mhalo_interp = np.arange(9, 16, 0.01)

matrix_out = np.zeros((Mstar.size + 1, z_array.size + 1))

matrix_out[0,-1] = -np.inf
for j in range(z_array.size):
    matrix_out[0,j] = z_array[j]
for i in range(Mstar.size):
    matrix_out[i+1,-1] = Mstar[i]

for j,z in enumerate(z_array):

    interpolator = interp1d(SMHM_Grylls_param(Mhalo_interp, SMHM_params, z, 0.), Mhalo_interp, fill_value="extrapolate")
    Mhalo = interpolator(Mstar)

    for i in range(Mstar.size):

        matrix_out[i+1, j] = Mhalo[i]

np.savetxt("SMHM_Grylls_extrapolated.txt", matrix_out)
