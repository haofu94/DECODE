import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

plt.ion()


with open("default.conf", "r") as conf_file:
    lines = conf_file.readlines()
    for i in range(len(lines)):
        line = lines[i].split()
        if "path_to_folder" in line: path_to_folder = line[2]
        break
if path_to_folder[-1] != "/": path_to_folder+="/"


z = 0.1

data = np.loadtxt(path_to_folder+"data/SFR_Mstar_relation_for_z_{:.2f}.txt".format(z))

plt.figure(1); plt.clf()
plt.rc('font', family="serif")
plt.rc('text', usetex=True)
plt.plot(data[:,0], np.log10(data[:,1]), ".", label="z = {:.2f}".format(z))
plt.xlim(8,12); plt.ylim(-1,3)
plt.xlabel(r"$\log (M_*/M_\odot)$")
plt.ylabel(r"$\log {\rm SFR} [M_\odot {\rm yr}^{-1}]$")
plt.legend(frameon=False, fontsize=12)



"""
# Moster fit with scipy.optimize.curve_fit
# guess parameters (from Grylls 2019) for scipy.optimize.curve_fit
p0 = np.array([], dtype=float)
p0 = np.append(p0, np.power(10, 10.65 + 0.33*z_for_SFR - 0.08*z_for_SFR**2))
p0 = np.append(p0, np.power(10, 0.69 + 0.71*z_for_SFR - 0.088*z_for_SFR**2))
p0 = np.append(p0, 1. - 0.022*z_for_SFR + 0.009*z_for_SFR**2)
p0 = np.append(p0, 1.8 - 1.*z_for_SFR - 0.1*z_for_SFR**2)
plt.plot(stellar_masses_for_SFR, np.log10(SFR_Moster(stellar_masses_for_SFR, p0[0], p0[1], p0[2], p0[3])), "--", color=color_idx)
"""
