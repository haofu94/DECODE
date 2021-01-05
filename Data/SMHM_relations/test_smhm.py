import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("/Users/haofu/GalaxyProjects/Discrete-STEEL/D-STEEL/dsteel/")
from DM_to_SM import SMHM_Moster

SMHM_params = [11.92, 0.58, 0.032, -0.014, 1.64, -0.69, 0.53, 0.03, 0.] #PyMorph

plt.ion()

mrange = np.linspace(10, 16, 1000)


#smhm = np.loadtxt("SMHM_null_sigma_extrapolated.txt", skiprows=2)
smhm = np.loadtxt("SMHM_Tomczak_extrapolated.txt", skiprows=2)

z = np.arange(0, 6, 0.1)

plt.figure(1); plt.clf()


#for i in np.arange(0, 128, 10):
#    plt.plot(z, smhm[i,:-1])


for i in np.arange(0, 50, 5):
    plt.plot(smhm[:,i], smhm[:,-1], label=r"{:.2f}".format(z[i]))
    plt.plot(mrange, SMHM_Moster(mrange, SMHM_params, z[i]), "-.")


plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.xlabel(r"$\log (M_{\rm h} / M_\odot)$", fontsize=15)
plt.ylabel(r"$\log (M_* / M_\odot)$", fontsize=15)
plt.rc('font', family="serif")
plt.legend(fontsize=12, frameon=False)
plt.xlim(10, 15.8); plt.ylim(6)
