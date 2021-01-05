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

plt.figure(1); plt.clf()
plt.rc('font', family="serif")
plt.rc('text', usetex=True)

#plt.plot(data[:,0], 10**data[:,1]/10**data[0,1], "-", color="C0", label = "Accreted mass fraction")

plt.plot(data[:,0], data[:,1], "-", color="red", label = "Total growth")
plt.plot(data[:,0], data[:,2], "--", color="red", label = "Accretion from mergers")
plt.plot(data[:,0], data[:,3], "-.", color="red", label = "Star formation history")

"""plt.plot(Cosmo.lookbackTime(data[:,0]), data[:,1], "-", color="C0", label = "Total growth")
plt.plot(Cosmo.lookbackTime(data[:,0]), data[:,2], "--", color="C0", label = "Accretion from mergers")
plt.plot(Cosmo.lookbackTime(data[:,0]), data[:,3], "-.", color="C0", label = "Star formation history")"""

"""Mstar = Mstar+0.1
idx = int(Mstar-9.+1.)
plt.plot(Cosmo.lookbackTime(sdss[:,0], inverse=True), sdss[:,idx], "^", color="blue")
plt.plot(Cosmo.lookbackTime(sdss[:,0], inverse=True), sdss[:,idx], ":", color="blue", label = "SDSS total growth")"""

idx = 3
#plt.plot(umachine[:,0], umachine[:,idx], "^", color="green")
#plt.plot(umachine[:,0], umachine[:,idx], "-", color="green", label = "U-Machine total growth")



plt.xlim(0.1, 4)
plt.ylim(Mstar-1.8)
plt.xscale("log")
plt.xticks([0.1, 0.3, 1, 2, 3], ['0.1', '0.3', '1', '2', '3'], fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel(r"$z$", fontsize=15)
plt.ylabel(r"$\log (M_*/M_\odot)$", fontsize=15)
plt.legend(frameon=False, fontsize=12)
plt.title(r"$\log (M_*/M_\odot)$ = {:.2f}".format(Mstar))
