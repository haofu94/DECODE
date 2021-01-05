"""
@ file halo_growth.py

Written by:
Chris Marsden, Francesco Shankar, Hao Fu

The main goal of this module is to provide a quick
computation of the accretion history of dark matter haloes.
"""

import numpy as np
from scipy import interpolate


try:
    #Test if a cosmology exists yet
    from colossus.cosmology import cosmology
    Cosmo = cosmology.getCurrent()
except:
    #Cosmology doesn't exist, get it from file
    from input_parameters import read_cosmo_model
    cosmo_str = read_cosmo_model(sys.argv[2:])
    cosmology.setCosmology(cosmo_str)
    Cosmo = cosmology.getCurrent()



def fit_var_sem(Mlog, h, omega_m):
    """
    Direct fitting formula by van den Bosch [note that the masses are in units of Msun/h]
    """

    omega_b=Cosmo.Ob0
    sigma_8=Cosmo.sigma8
    spectral_index=Cosmo.ns

    Mh = np.linspace(17, 5)

    gammam = omega_m*h*np.exp(-omega_b*(1+np.sqrt(2.*h)/omega_m))

    c = 3.804e-4
    x = (c*gammam*((10**(Mh/3.))))/((omega_m)**(1./3.))
    g1 = 64.087 * ((1 + 1.074*(x**0.3) - 1.581*(x**0.4) + 0.954*(x**0.5) - 0.185*(x**(0.6)))**(-10))
    x = (32*gammam)
    g2 = 64.087 * ((1+1.074*(x**0.3) - 1.581*(x**0.4) + 0.954*(x**0.5) - 0.185*(x**0.6))**-10)
    f = (g1*g1)/(g2*g2)


    s = np.sqrt(f*sigma_8*sigma_8)
    sig=s*(10**((Mh-14.09)*(1-spectral_index)/9.2))/(1+(1.00-spectral_index)/9.2)

    sig2Mh = interpolate.interp1d(Mh, sig, fill_value="extrapolate")
    res = sig2Mh(Mlog)

    return res


def D_z_White(z, omega_m):
    """
    Mo & White 2002 the normalized Growth factor which is the one to be included
    in The P&S MF (see Eke, Coles & Frenk 1996).
    The growth factor of perturbations is however just given by gz!!!
    """
    omega_l = 1-omega_m
    omega_k = 0
    Ez = (omega_l+omega_k*(1+z)**2.+omega_m*(1+z)**3)**0.5
    omega_m_z=omega_m*(1+z)**3/Ez**2.
    omega_l_z=omega_l/Ez**2
    gz=2.5*omega_m_z*(omega_m_z**(4/7)-omega_l_z+(1+omega_m_z/2)*(1+omega_l_z/70))**(-1)

    gz0=2.5*omega_m*(omega_m**(4/7)-omega_l+(1+omega_m/2)*(1+omega_l/70))**(-1)

    dz=(gz/gz0)/(1+z)

    return dz


def Mass_acc_history_VDB_FS(M0, zc, h, omega_m):

    # Parameters for average
    apar1 = 3.2954
    apar2 = 0.1975
    apar3 = 0.7542
    apar4 = 0.0898
    apar5 = 0.4415

    z0=0.00

    logM0h = M0+np.log10(h)  # Msun/h

    var = fit_var_sem(logM0h, h, omega_m)
    s0 = var**2
    #print(var)

    dc0 = 1.686/D_z_White(0, omega_m)

    logMzc = np.linspace(M0-0.0001, M0-7., 1000) # Possible range of masses
    logpsiz = logMzc - M0 # In log, the rato of the possible range of masses to M0
    psiz = 10**logpsiz # The value un-logged
    Fpsi = apar1*(1-apar2*logpsiz)**apar3*(1-psiz**apar4)**apar5

    logMzch = logMzc + np.log10(h)
    varz = fit_var_sem(logMzch, h, omega_m)
    sz = varz**2
    Gz=0.57*(sz/s0)**0.19*(dc0/varz)**(-0.01)
    nz=len(zc)
    logMz=np.zeros(nz)
    logMz[0]=M0 # The first element is the last element

    # Loop over every redshift interval beyond the first
    for j in range(nz-1):
        dcz = 1.686/D_z_White(zc[j+1], omega_m)
        wz = (dcz-dc0)/np.sqrt(sz-s0)
        gam = 0.4
        wtz = wz*Gz**gam
        interpolator = interpolate.interp1d(Fpsi-wtz, logpsiz, fill_value="extrapolate")
        logpsi = interpolator(0)

        logMz[j+1] = logpsi + M0

    return logMz


def getMassEvolution(m0_range, z_range, h, omega_m):
    results = np.zeros((len(m0_range), len(z_range)))

    print("Results Shape", results.shape)
    results = Parallel(n_jobs = -1)(delayed(Mass_acc_history_VDB_FS)(m0, z_range, h, omega_m) for m0 in m0_range)

    results = np.array(results)

    print("Results Shape", results.shape)

    return results



def get_mean_halo_growths(mass_params, z_range):
    """
    Calculate mean halo growth history
    """

    halo_mass_range = np.arange(mass_params[0], mass_params[1], mass_params[2])
    N = halo_mass_range.size - 1

    mean_halo_tracks = {
        "mass_bin": [],
        "mass_track": []
    }

    for i in range(N):
        mean_halo_tracks["mass_bin"].append([halo_mass_range[i], halo_mass_range[i+1]])
        mean_mass = (halo_mass_range[i] + halo_mass_range[i+1]) /2.
        mean_halo_tracks["mass_track"].append(Mass_acc_history_VDB_FS(mean_mass, z_range, Cosmo.h, Cosmo.Ob0))

    return mean_halo_tracks
