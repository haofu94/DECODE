import numpy as np
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

data = np.loadtxt(path_to_folder+"data/pdf_redshift_13.0.txt")

plt.figure(1); plt.clf()
plt.rc('font', family="serif")
plt.rc('text', usetex=True)
plt.plot(data[:,0], np.log10(data[:,1]))
plt.xlim(0,1)
#plt.ylim(-2,1)
plt.xlabel(r"$\log_{10} (z+1)$", fontsize=15)
plt.ylabel(r"$\log_{10} \phi(z) \;\; [{\rm dex}^{-1}]$", fontsize=15)
#plt.legend(frameon=False, fontsize=10)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.show()
