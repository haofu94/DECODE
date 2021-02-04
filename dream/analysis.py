"""
@ file analysis.py

Written by Hao fu

UNDER DEVELOPMENT
"""

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from tqdm import tqdm

import analytic_SHMFs
from halo_growth import Mass_acc_history_VDB_FS as VDB
from DM_to_SM import *


import warnings
warnings.filterwarnings("ignore")


import ctypes
from numpy.ctypeslib import ndpointer

analyze = ctypes.CDLL("dream/C_functions/analyze.so")

analyze.get_evolved_mergers_number.argtypes = [ctypes.c_double,
                                               ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                               ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                               ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                               ctypes.c_int]

analyze.get_evolved_mergers.argtypes = [ctypes.c_double,
                                        ctypes.c_double,
                                        ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                        ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                        ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                        ctypes.c_double,
                                        ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                        ctypes.c_double,
                                        ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                        ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                        ctypes.c_double,
                                        ctypes.c_int,
                                        ctypes.c_int]



from colossus.cosmology import cosmology
Cosmo = cosmology.getCurrent()
h = Cosmo.h


#plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)




def vdB_USHMF(Params, psi):
    """
    Jiang & van den Bosch 2016
    """
    gamma, alpha, beta, omega, a = Params
    Part1 = gamma*np.power(a*psi, alpha)
    Part2 = np.exp(-beta*np.power(a*psi, omega))
    dn_dlnX_arr = Part1*Part2
    dn_dlogX_arr = dn_dlnX_arr*np.log(10.)
    return dn_dlogX_arr #N dex-1


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


def SFR_Moster(mstar, M, N, a, b):
    """
    M = np.power(10, 10.65 + 0.33*z - 0.08*z**2)
    N = np.power(10, 0.69 + 0.71*z - 0.088*z**2)
    a = 1. - 0.022*z + 0.009*z**2
    b = 1.8 - 1.*z - 0.1*z**2
    """
    return 2*N / (np.power(10**mstar/M, -a) + np.power(10**mstar/M, b))


def f_loss(tau):
    # tau in [Gyr]
    return 0.05*np.log(1. + tau/1.4e-3) #Moster 2018


def forward(x):
    return x**(1/2)

def inverse(x):
    return x**2



def get_z_infall_distribution(z_infall,
                              data_folder,
                              file_name):

    logz_bin = 0.1
    logz_bins = np.arange(0., 1.+logz_bin, logz_bin)
    phi_logz = np.histogram(z_infall, bins = logz_bins)[0] / logz_bin
    logz_bins = logz_bins[:-1] + logz_bin/2.

    header = "1) 1st column: log(1+z)\n2) 2nd column: phi(z) [dex^-1]\n"

    np.savetxt(data_folder+file_name, np.column_stack((logz_bins, phi_logz)), header = header)

    return



def get_unevolved_shmf(mergers_array,z_at_merge,
                       M_h, Msh_bins, bin,
                       N_halos, idx_sub,
                       data_folder):

    header = "1) 1st column: Msubhalo/Mparent\n2) 2nd column: phi(Msubhalo) [dex^-1]\n"

    psi = 10**Msh_bins[:-1] / 10**M_h

    phi_M = np.histogram(mergers_array[idx_sub], bins=Msh_bins)[0] / bin / N_halos
    np.savetxt(data_folder+"unevolved-shmf_Mh_"+str(M_h)+".txt", np.column_stack((psi, phi_M)), header = header)

    mask = ma.array(z_at_merge, mask=z_at_merge<0.).recordmask
    mask = np.logical_and(mask, idx_sub)
    surv_phi_M = np.histogram(mergers_array[mask], bins=Msh_bins)[0] / bin / N_halos
    np.savetxt(data_folder+"unevolved-surviving-shmf_Mh_"+str(M_h)+".txt", np.column_stack((psi, surv_phi_M)), header = header)

    return



def get_evolved_shmf(z_,
                     mergers_array,
                     z_infall,
                     z_at_merge,
                     M_h, Msh_bins,
                     N_halos, idx_sub,
                     data_folder):

    z_arr = np.arange(0,20,0.1)
    Idx = np.argsort(np.abs(z_arr-z_))[0]
    M_h_track = VDB(M_h, z_arr, Cosmo.h, Cosmo.Om0)
    ev_M_h = M_h_track[Idx]
    Psi = 10**(Msh_bins) /10**ev_M_h
    fs = analytic_SHMFs.get_fs(M_h_track, z_arr)
    Params = [0.280*fs, -0.86, 50., 4., 1.] #Jiang & van den Bosch 2016 Table A1, Evolved
    ana_esHMF = vdB_USHMF(Params, Psi)

    analyze.get_evolved_mergers_number.restype = ctypes.c_int
    evolved_mergers_number = analyze.get_evolved_mergers_number(z_,
                                                                mergers_array[idx_sub],
                                                                z_at_merge[idx_sub],
                                                                z_infall[idx_sub],
                                                                mergers_array[idx_sub].size)

    analyze.get_evolved_mergers.restype = ndpointer(dtype=ctypes.c_double, shape = np.zeros(evolved_mergers_number).shape)
    #print(evolved_mergers_number)
    ev_mergers_array = analyze.get_evolved_mergers(M_h, z_,
                                                   mergers_array[idx_sub],
                                                   z_at_merge[idx_sub],
                                                   z_infall[idx_sub],
                                                   Cosmo.age(z_),
                                                   Cosmo.age(z_infall[idx_sub]),
                                                   Cosmo.h,
                                                   Cosmo.Om(z_infall[idx_sub]),
                                                   Cosmo.Hz(z_infall[idx_sub]),
                                                   Cosmo.H0,
                                                   mergers_array[idx_sub].size,
                                                   evolved_mergers_number)

    """
    ev_mergers_array = np.array([], dtype=float)
    for i,m in enumerate(mergers_array[idx_sub]):
        if z_at_merge[idx_sub][i]<z_ and z_<z_infall[idx_sub][i]:
            ev_mergers_array = np.append(ev_mergers_array, mass_at_t(m, M_h, z_, z_infall[idx_sub][i]))"""

    ev_hist = np.histogram(ev_mergers_array, bins=Msh_bins)
    ev_sHMF_from_hist = ev_hist[0]/0.1/N_halos

    #if output_verbose in [2,3]:
    if True:

        """if output_verbose in [2,3]:
            if use_latex:
                plt.rc('text', usetex=True)
                ampersand = "\&"
            else:
                ampersand = "&"
        plt.figure(3); plt.clf(); plt.figsize=((12,20))
        plt.plot(1,100,"-", color="black", label = "analytical shmf (Jiang "+ampersand+" vdB 2016)", alpha=0.8)
        plt.plot(np.log10(Psi), np.log10(ana_esHMF), "-", color="black", label=r"$\log M_h = $ {:.1f}".format(M_h), alpha=0.6)
        plt.plot(np.log10(Psi[0:-1]), np.log10(ev_sHMF_from_hist), ".", color="C1", label=r"$\log M_h = $ {:.1f} (from histogram)".format(M_h))
        plt.xlim(right=0); plt.ylim(-2,2)
        plt.xlabel(r"$\log_{10} (M_{\rm subhalo} / M_{\rm host})$", fontsize=15)
        plt.ylabel(r"$\log_{10} \phi(M_{\rm subhalo}) \;\; [{\rm dex}^{-1}]$", fontsize=15)
        plt.title("Evolved subhalo mass function")
        plt.rc('font', family=font_for_plot)
        plt.legend(frameon=want_frame)
        plt.savefig(plots_folder+"evolved-shmf_Mh_"+str(M_h)+".pdf", dpi=300)"""
        ########################
        fp = open(data_folder+"evolved-shmf_Mh_"+str(M_h)+".txt", "w")
        for i in range(ev_sHMF_from_hist.size):
            fp.write("{:.6f} {:.6f}\n".format(np.log10(Psi[i]), np.log10(ev_sHMF_from_hist[i])))
        fp.close()
        ########################

        """plt.figure(4); plt.clf(); plt.figsize=((12,20))
        plt.plot(np.log10(Psi[0:-1]), np.log10(ana_esHMF[0:-1]/ev_sHMF_from_hist), label=r"$\log M_h = $ {:.1f}".format(M_h))
        plt.xlabel(r"$\log (M_{\rm subhalo} / M_{\rm host})$", fontsize=15)
        plt.ylabel("log[ratio]")
        plt.ylim(-0.8,0.8)
        plt.rc('font', family=font_for_plot)
        plt.title("Ratio of evolved-shmf")
        plt.legend(frameon=want_frame)
        plt.savefig(plots_folder+"ratio_evolved-shmf_Mh_"+str(M_h)+".pdf", dpi=300)"""

    return



def get_mergers_rate(mergers_array,
                     z_infall,
                     M_h, N_halos, idx_sub,
                     data_folder):

    z_bin = 0.1
    z = np.arange(0, 3.2, z_bin)
    major_merger_frac = 0.3
    M_h_track = VDB(M_h, z, Cosmo.h, Cosmo.Om0)
    mergers_list_z = np.zeros(z.size)
    for j in range(mergers_array[idx_sub].size):
        for k, Z in enumerate(z):
            if Z <= z_infall[idx_sub][j] <Z+z_bin:
                if 10**mergers_array[idx_sub][j]/10**M_h_track[k] > major_merger_frac:
                    mergers_list_z[k] += 1
    merger_rate = mergers_list_z/(Cosmo.age(z)-Cosmo.age(z+z_bin)) / N_halos
    z = z + z_bin/2.
    header = "1) 1st column: redshift\n2) 2nd column: dN/dt [Gyr^-1]\n"
    np.savetxt(data_folder+"merger_rate_Mh_"+str(M_h)+".txt", np.column_stack((z, merger_rate)), header = header)

    return
