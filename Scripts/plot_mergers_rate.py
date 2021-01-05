import numpy as np
import numpy.ma as ma
from scipy.integrate import trapz
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("../dsteel/")

from colossus.cosmology import cosmology
Cosmo = cosmology.setCosmology("planck18")
h = Cosmo.h

from halo_growth import Mass_acc_history_VDB_FS as VDB
#from analysis import *

plt.ion()


def merger_rate_FM(Mh, xi, z, z_bin , Params, major_merger_frac):
    """
    Fakhouri & Ma et al. 2010, Eq. (1) integrated over xi
    """
    # M [log10M_sun/h]
    a, b, g, n, A, xi0 = Params
    Mh = Mh - np.log10(h)
    idx = np.argsort(np.abs(xi-major_merger_frac))[0]
    result = A* np.power(10.**(Mh-12.), a) * trapz( np.power(xi[idx:], b) * np.exp(np.power(xi[idx:]/xi0, g)), xi[idx:]) *np.power(1+z, n)
    Result = np.zeros(z.size-1)
    for i in range(z.size-1):
        Result[i] = trapz(result[i:i+2], z[i:i+2]) / (Cosmo.age(z[i])-Cosmo.age(z[i]+z_bin))

    return Result #[N/Gyr]


with open("default.conf", "r") as conf_file:
    lines = conf_file.readlines()
    for i in range(len(lines)):
        line = lines[i].split()
        if "path_to_folder" in line: path_to_folder = line[2]
        break
if path_to_folder[-1] != "/": path_to_folder+="/"


M_h = 11.

z_bin = 0.1
z = np.arange(0, 3.+z_bin, z_bin)
major_merger_frac = 0.3
M_h_track = VDB(M_h, z, Cosmo.h, Cosmo.Om0)

Msh_bins = np.arange(9,M_h+0.05,0.01)
psi = 10**Msh_bins/10**M_h
z__bin = 0.01; z__ = np.arange(0,2.6, z__bin)
M_h_track_ = VDB(M_h, z__, Cosmo.h, Cosmo.Om0)
ana_merger_rate = merger_rate_FM(M_h_track_, psi, z__, z__bin, [0.133, -1.995, 0.263, 0.0993, 0.0104, 9.72e-3], major_merger_frac)

data = np.loadtxt(path_to_folder+"data/merger_rate_Mh_"+str(M_h)+".txt")

plt.figure(1); plt.clf()
plt.rc('font', family="serif")
plt.rc('text', usetex=True)
plt.plot(data[:,0], data[:,1], "-.", color="C0", label=r"$\log M_h = $ {:.1f} (D-STEEL)".format(M_h))
plt.plot(z__[0:-1], ana_merger_rate, "-", color="C0", label=r"Fakhouri \& Ma 2010")
plt.yscale("log")
plt.xlim(0, 2.5)
plt.ylim(1e-2,1e0)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.xlabel(r"$z$", fontsize=15)
plt.ylabel(r"$dN/dt \;\; [{\rm Gyr}^{-1}]$", fontsize=15)
plt.legend(frameon=False, fontsize=12, loc="lower right")
