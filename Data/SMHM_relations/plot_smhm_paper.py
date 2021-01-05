import numpy as np
import matplotlib.pyplot as plt

plt.ion()

plt.rc('text', usetex=True)
plt.rc('font', family="serif")
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)


smhm_null_sigma = np.loadtxt("SMHM_null_sigma.txt", skiprows=2)
smhm_const_sigma = np.loadtxt("SMHM_constant_sigma_extrapolated.txt", skiprows=2)
smhm_Tomczak = np.loadtxt("SMHM_Tomczak_extrapolated.txt", skiprows=2)
smhm_var_sigma = np.loadtxt("SMHM_variable_sigma_extrapolated.txt", skiprows=2)
smhm_Tomczak_var_sigma = np.loadtxt("SMHM_Tomczak_variable_sigma_extrapolated.txt", skiprows=2)

z = np.arange(0, 6, 0.1)

fig, ax = plt.subplots(4, 1, figsize=(5.5,9.5), sharex=True, sharey=True, gridspec_kw={'hspace': 0, 'height_ratios': [1, 1, 1, 1]})

idx_colour = 0
colours = ["red", "blue", "green", "darkviolet"]
#colours = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"]
for i in np.arange(0, 15, 5):
    c = colours[idx_colour]
    ax[0].plot(smhm_const_sigma[:,i], smhm_const_sigma[:,-1], color=c, label=r"$z$ = {:.1f}".format(z[i]))
    ax[1].plot(smhm_Tomczak[:,i], smhm_Tomczak[:,-1], color=c, label=r"{:.1f}".format(z[i]))
    ax[2].plot(smhm_var_sigma[:,i], smhm_var_sigma[:,-1], color=c, label=r"{:.1f}".format(z[i]))
    ax[3].plot(smhm_Tomczak_var_sigma[:,i], smhm_Tomczak_var_sigma[:,-1], color=c, label=r"{:.1f}".format(z[i]))
    idx_colour += 1



ax[0].set_ylabel(r"$\log (M_* / M_\odot)$", fontsize=15)
ax[0].legend(fontsize=12, frameon=False, loc="upper left")
ax[0].set_xlim(10, 15.8)
ax[0].set_ylim(7, 12.3)
ax[0].text(14, 8, "Model 1", fontsize=15)

ax[1].set_ylabel(r"$\log (M_* / M_\odot)$", fontsize=15)
ax[1].text(14, 8, "Model 2", fontsize=15)

ax[2].set_ylabel(r"$\log (M_* / M_\odot)$", fontsize=15)
ax[2].text(14, 8, "Model 3", fontsize=15)

ax[3].set_xlabel(r"$\log (M_{\rm h} / M_\odot)$", fontsize=15)
ax[3].set_ylabel(r"$\log (M_* / M_\odot)$", fontsize=15)
ax[3].text(14, 8, "Model 4", fontsize=15)
