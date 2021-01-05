import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

plt.ion()


#z_array = np.loadtxt("SMHM_Tomczak_test.txt", skiprows=1, max_rows=1)[:-1]
#matrix = np.loadtxt("SMHM_Tomczak_test.txt", skiprows=2)

#z_array = np.loadtxt("SMHM_null_sigma.txt", skiprows=1, max_rows=1)[:-1]
#matrix = np.loadtxt("SMHM_null_sigma.txt", skiprows=2)

z_array = np.loadtxt("SMHM_constant_sigma.txt", skiprows=1, max_rows=1)[:-1]
matrix = np.loadtxt("SMHM_constant_sigma.txt", skiprows=2)


new_z = np.arange(2.1, 6., 0.1)

Mstar = matrix[:,-1]

matrix_out = np.zeros((Mstar.size + 1, z_array.size + new_z.size + 1))

matrix_out[0,-1] = -np.inf
for j in range(z_array.size):
    matrix_out[0,j] = z_array[j]
    for i in range(Mstar.size):
        matrix_out[i+1,j] = matrix[i,j]
for j in range(new_z.size):
    matrix_out[0,z_array.size+j] = new_z[j]
for i in range(Mstar.size):
    matrix_out[i+1,-1] = Mstar[i]

for j,z in enumerate(new_z):

    Mhalo = np.array([interp1d(z_array, matrix[i,:-1], fill_value="extrapolate")(z) for i in range(Mstar.size)])

    idx = z_array.size + j

    for i in range(Mstar.size):

        matrix_out[i+1, idx] = Mhalo[i]

#np.savetxt("SMHM_Tomczak_extrapolated.txt", matrix_out)

#np.savetxt("SMHM_null_sigma_extrapolated.txt", matrix_out)

np.savetxt("SMHM_constant_sigma_extrapolated.txt", matrix_out)
