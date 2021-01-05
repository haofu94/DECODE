import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

plt.ion()


def forward(x):
    return x**(1/2)

def inverse(x):
    return x**2


with open("default.conf", "r") as conf_file:
    lines = conf_file.readlines()
    for i in range(len(lines)):
        line = lines[i].split()
        if "path_to_folder" in line: path_to_folder = line[2]
        break
if path_to_folder[-1] != "/": path_to_folder+="/"


Mstar = 11

data = np.loadtxt(path_to_folder+"data/SFR_z_relation_for_Mstar_{:.2f}.txt".format(Mstar))

plt.figure(1); plt.clf()
plt.rc('font', family="serif")
plt.rc('text', usetex=True)

plt.plot(data[:,0], np.log10(data[:,2]), "-", label=r"$\log M_* = $ {:.1f}".format(Mstar))
plt.xlim(0,9)
plt.xlabel("z")
plt.xscale("symlog")
plt.xscale("function", functions=(forward, inverse))
plt.ylabel(r"$\log {\rm SFR} [M_\odot {\rm yr}^{-1}]$")
plt.xticks([0,1,2,3,4,5,6,7,8,9], ['0', '1', '2','3','4', '5', '6', '7', '8', '9'])

"""
plt.plot(data[:,0], data[:,1], "-", label=r"$\log M_* = $ {:.1f}".format(Mstar))
plt.xlim(0,12); plt.ylim(0)
plt.xlabel("lookback time [Gyr]")
plt.ylabel(r"${\rm SFR} [M_\odot {\rm yr}^{-1}]$")
"""

plt.xlabel(r"$\log (M_*/M_\odot)$")
plt.legend(frameon=False, fontsize=12)
