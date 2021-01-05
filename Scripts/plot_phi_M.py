import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("../dsteel/")

from colossus.cosmology import cosmology
cosmology.setCosmology("planck18")

from analytic_SHMFs import *

plt.ion()


with open("default.conf", "r") as conf_file:
    lines = conf_file.readlines()
    for i in range(len(lines)):
        line = lines[i].split()
        if "path_to_folder" in line: path_to_folder = line[2]
        break
if path_to_folder[-1] != "/": path_to_folder+="/"


Params = [0.22, -0.91, 6., 3., 1.] #Jiang & van den Bosch 2016 Table A1, Unevolved, total

M_h = 14.
bin = 0.1
Msh_bins = np.arange(M_h-3., M_h, 0.1)
psi = 10**Msh_bins/10**M_h

analytical_SHMF = vdB_USHMF(Params, M_h, Msh_bins)

Params_1st_order = [0.13, -0.83, 1.33, -0.02, 5.67, 1.19]
analytical_SHMF_1st_order = vdB_USHMF_1st_order(Params_1st_order, M_h, Msh_bins)
analytical_SHMF_2nd_order = vdB_USHMF_ith_order(Params_1st_order, M_h, Msh_bins, 2)

data = np.loadtxt(path_to_folder+"data/unevolved-shmf_Mh_"+str(M_h)+".txt")


plt.figure(1); plt.clf()
plt.rc('font', family="serif")
plt.rc('text', usetex=True)
plt.plot(np.log10(psi), np.log10(analytical_SHMF), "-", color="black", alpha=0.6, label = r"analytical shmf (Jiang \& vdB 2016)")
plt.plot(np.log10(psi), np.log10(analytical_SHMF_1st_order), "-", color="red", label = "1st order")
plt.plot(np.log10(psi), np.log10(analytical_SHMF_2nd_order), "-", color="red", label = "2nd order")
plt.plot(np.log10(data[:,0]), np.log10(data[:,1]), "-.", color="C1", label=r"$\log M_h = $ {:.1f} (D-STEEL)".format(M_h))
plt.xlim(right=0); plt.ylim(-1,3.5)
plt.legend(frameon=False)
plt.xlabel(r"$\log (M_{\rm subhalo} / M_{\rm host})$", fontsize=15)
plt.ylabel(r"$\log \phi(M_{\rm subhalo}) \;\; [{\rm dex}^{-1}]$", fontsize=15)
plt.title("Unevolved subhalo mass function")
