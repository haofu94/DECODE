import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("../dsteel/")

from colossus.cosmology import cosmology
cosmology.setCosmology("planck18")
Cosmo = cosmology.getCurrent()

plt.ion()


with open("default.conf", "r") as conf_file:
    lines = conf_file.readlines()
    for i in range(len(lines)):
        line = lines[i].split()
        if "path_to_folder" in line: path_to_folder = line[2]
        break
if path_to_folder[-1] != "/": path_to_folder+="/"


Mstar = 12.

data = np.loadtxt(path_to_folder+"data/SFH_for_Mstar_{:.2f}.txt".format(Mstar))

sdss = np.loadtxt("../Data/SDSS_data/galaxy_growth.dat")

umachine = np.loadtxt("../Data/UMachine_release-sfh_z0_z8_052913/Berhoozi_Mstar_tracks.txt")

z = data[:-1,0] + np.mean(np.diff(data[:,0]))/2.
tot = data[:,1]
mer = data[:,2]
sfh = data[:,3]

delta_tot = -np.diff(10**tot)
delta_mer = -np.diff(10**mer)
delta_sfh = -np.diff(10**sfh)

plt.figure(2); plt.clf()
plt.rc('font', family="serif")
plt.rc('text', usetex=True)

#plt.plot(z, delta_tot/delta_tot, "-", color="red", label = "Total growth")
plt.plot(z, delta_mer/delta_tot, "--", color="red", label = r"$\dot M$ of mergers accretion / $\dot M_{\rm tot}$")
plt.plot(z, delta_sfh/delta_tot, "-.", color="red", label = r"$\dot M$ of SFH accretion / $\dot M_{\rm tot}$")
plt.hlines(1, 0, 10, color="black", linewidth=1)




plt.xlim(0.1, 4)
#plt.ylim(0.1)
plt.xscale("log")
plt.yscale("log")
plt.xticks([0.1, 0.3, 1, 2, 3], ['0.1', '0.3', '1', '2', '3'], fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel(r"$z$", fontsize=15)
plt.ylabel(r"$(\dot M_* / M_\odot)$", fontsize=15)
plt.legend(frameon=False, fontsize=12)
plt.title(r"$\log (M_*/M_\odot)$ = {:.2f}".format(Mstar))
