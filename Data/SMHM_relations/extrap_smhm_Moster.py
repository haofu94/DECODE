import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

plt.ion()



def SMHM_Moster(DM, z):

    # MCMC FIT
    """M = 11.339 + 0.692 * z / (z+1.)
    e = 0.005 + 0.689 * z / (z+1.)
    b = 3.344 - 2.079 * z / (z+1.)
    g = 0.966"""

    # ALL CENTRALS
    """if 0 <= z < 0.3:
        M = 11.8; e = 0.14; b = 1.75; g = 0.57
    elif 0.3 <= z < 0.75:
        M = 11.85; e = 0.16; b = 1.7; g = 0.58
    elif 0.75 <= z < 1.5:
        M = 11.95; e = 0.18; b = 1.6; g = 0.6
    elif 1.5 <= z < 3.:
        M = 12.; e = 0.18; b = 1.55; g = 0.62
    elif 3. <= z < 6.:
        M = 12.05; e = 0.19; b = 1.5; g = 0.64
    elif 6 <= z:
        M = 12.1; e = 0.24; b = 1.3; g = 0.64"""

    # ALL GALAXIES
    if 0 <= z < 0.3:
        M = 11.78; e = 0.15; b = 1.78; g = 0.57
    elif 0.3 <= z < 0.75:
        M = 11.86; e = 0.18; b = 1.67; g = 0.58
    elif 0.75 <= z < 1.5:
        M = 11.98; e = 0.19; b = 1.53; g = 0.59
    elif 1.5 <= z < 3.:
        M = 11.99; e = 0.19; b = 1.46; g = 0.59
    elif 3. <= z < 6.:
        M = 12.07; e = 0.2; b = 1.36; g = 0.6
    elif 6 <= z:
        M = 12.1; e = 0.24; b = 1.3; g = 0.6

    SM = np.log10( 2. * 10**DM * e / ( np.power(10**DM/10**M, -b) + np.power(10**DM/10**M, g) ) )
    #SM = 2. * e / ( np.power(10**DM/10**M, -b) + np.power(10**DM/10**M, g) )

    return SM

"""
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

    interpolator = interp1d(SMHM_Moster(Mhalo_interp, z), Mhalo_interp, fill_value="extrapolate")
    Mhalo = interpolator(Mstar)

    for i in range(Mstar.size):

        matrix_out[i+1, j] = Mhalo[i]

np.savetxt("SMHM_Moster_2018_extrapolated.txt", matrix_out)
"""



Mpeak = np.linspace(10, 14.5, 1000)
Mstar = SMHM_Moster(Mpeak, 1.1)

plt.figure(1); plt.clf()
plt.plot(Mpeak, Mstar)
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.xlabel(r"$\log (M_{\rm peak} / M_\odot)$", fontsize=15)
plt.ylabel(r"$\log (M_* / M_\odot)$", fontsize=15)
plt.rc('font', family="serif")
#plt.yscale("log")
