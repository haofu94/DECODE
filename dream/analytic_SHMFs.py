"""
@ file analytic_SHMFs.py

This module provides the analytical form of different subhalo mass functions,
as well as other useful quantities concerning mergers.
"""

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from colossus.halo.mass_so import M_to_R
from colossus.cosmology import cosmology


Cosmo = cosmology.getCurrent()
h = Cosmo.h


def Delta_vir(z):
    """
    Bryan & Norman 1998, ApJ, 495, 80
    """
    x = Cosmo.Om(z) - 1
    return (18*np.power(np.pi, 2) + 82*x - 39*np.power(x, 2))


def tau_dyn(Redshift):
    """
    Jiang & van den Bosch 2016, MNRAS 458, 2870–2884

    Return
    DynamicalTime: dynamical timescale in [Gyr]
    """
    G_MPC_Msun_Gyr = -14.34693206746732 #np.log10(constants.G.value*np.power(constants.pc.value*10**6,-3)*(constants.M_sun.value)*np.power(3.15576*10**16, 2))
    DynamicalTime = 1.628*(h**-1)*np.power(Delta_vir(Redshift)/178, -0.5)*np.power(Cosmo.Hz(Redshift)/Cosmo.H0, -1)
    return DynamicalTime #Jiang & VDB 2016 #Gyr


def mass_at_t(m, M, z, z_inf):

  m = np.power(10., m)
  M = np.power(10., M)
  xi = 0.07
  A = 1.54
  #A = 0.86
  #xi = 0.36
  #A = 23.7
  dt = Cosmo.age(z) - Cosmo.age(z_inf)
  tau = tau_dyn(z_inf) / A

  return np.log10( m* np.power(1.+ xi * np.power(m/M, xi) *(dt/tau), -1./xi) )
  #return log10( m * exp(-dt/tau))

def f_function(x):
    return np.log(1.+x) - x / (1.+x)



def tau_merge(central_mass, sub_mass, fudge, z):
    """
    Merging time-scale, Boylan-Kolchin et al. 2008

    Input params
    central_mass: central halo mass [log10(Msun)]
    sub_mass: subhalo mass [log10(Msun)]

    Return
    Merging timescale in units of [Gyr]
    """

    # Note that this needs to be in h^-1
    A = 0.9; B = 1.0; C = 0.6; D = 0.1 #McCavana et al. 2012 values

    MassRatio = 10**(central_mass - sub_mass) #unitless

    VR = np.divide(M_to_R(10**central_mass, z, 'vir'), h*(10**3) ) # mega parsecs - Not sure what this is?
    Tdyn = tau_dyn(z) #Gyr

    #dynamical_friction /= 1. + z;
    """
    if Paramaters['Aldynamical_frictionamicalTime'] != False:
        dynamical_friction = dynamical_friction * Paramaters['Aldynamical_frictionamicalTime']
    if Paramaters['Aldynamical_frictionamicalTimeB']:
        dynamical_friction = dynamical_friction * (1/(1+Redshift))
    """

    #fudge = 0.5
    if fudge<0.:
        x = 1./MassRatio; m = 4.032021; a = 0.03724761
        fudge = m*(np.log10(x+a) -np.log10(a))
    Tdyn = fudge*Tdyn

    NormalRnd = 0.5 #We are considering only avarage halos on circular orbits

    OrbitalEnergy = (VR * NormalRnd**2.17) / (1 - np.sqrt( 1 - (NormalRnd**2) ) ) #energy

    Part1 = MassRatio**B
    Part2 = np.log( 1 + MassRatio )
    Part3 = np.exp( C * NormalRnd ) #exponent of orbital circularity
    Part4 = ( OrbitalEnergy/VR )**D #(OrbitalEnergy/virialradius)^D

    return Tdyn*A*(Part1/Part2)*Part3*Part4 #returns a timescale of infall [Gyr]


def get_fs(M,z):
    """
    Function that calculates fs for first order subhalos, according to Jiang & van den Bosch 2016 Eq. (19)
    """
    M = 10.**M
    M_half = M[0]/2.

    interpolator = interp1d(M, z)

    z_form = interpolator(M_half)
    z_arr = np.linspace(0, z_form, 10001)
    dt_dz = np.zeros(z_arr.size-1)

    for i in range(z_arr.size-1):
        dt_dz[i] = -(Cosmo.age(z_arr[i+1])-Cosmo.age(z_arr[i])) /(z_arr[i+1]-z_arr[i])

    z_arr = np.linspace(0, z_form, 1000)
    N_tau = quad(tau_dyn, 0., z_form)[0]
    fs = 0.325/np.power(N_tau, 0.6) - 0.075

    return fs


def vdB_USHMF(Params, M_host, M_subhalos_range):
    """
    unevolved SHMF, Jiang & van den Bosch 2016 Appendix A, Eq. (A1)

    Input params
    Params: list of shmf parameters
    M_host: mass of the parent halo in [log10(Msun)]
    M_subhalos_range: mass interval of subhalos in [log10(Msun)] (numpy.ndarray)

    Return
    dn_dlogX_arr: subhalo mass function in units of [N/dex] as function of psi
    """
    psi = 10.**M_subhalos_range/10.**M_host
    gamma, alpha, beta, omega, a = Params
    Part1 = gamma*np.power(a*psi, alpha)
    Part2 = np.exp(-beta*np.power(a*psi, omega))
    dn_dlnX_arr = Part1*Part2
    dn_dlogX_arr = dn_dlnX_arr*np.log(10.)
    return dn_dlogX_arr


def vdB_USHMF_double(Params, M_host, M_subhalos_range):
    psi = 10.**M_subhalos_range/10.**M_host
    gamma1, alpha1, gamma2, alpha2, beta, xi = Params
    Part1 = gamma1*np.power(psi, alpha1) + gamma2*np.power(psi, alpha2)
    Part2 = np.exp(-beta*np.power(psi, xi))
    dn_dlnX_arr = Part1*Part2
    dn_dlogX_arr = dn_dlnX_arr*np.log(10.)
    return dn_dlogX_arr


def vdB_USHMF_1st_order(Params, M_host, M_subhalos_range):
    """
    Jiang \& van den Bosch 2016
    """
    psi = 10.**M_subhalos_range/10.**M_host
    gamma1, alpha1, gamma2, alpha2, beta, omega = Params
    Part1 = gamma1*np.power(psi, alpha1)
    Part2 = gamma2*np.power(psi, alpha2)
    dn_dlnX_arr = (Part1 + Part2) * np.exp(-beta * np.power(psi, omega))
    dn_dlogX_arr = dn_dlnX_arr*np.log(10.)
    return dn_dlogX_arr #N dex-1



def vdB_USHMF_ith_order(Params, M_host, M_subhalos_range, order):

    m_for_integrate = np.arange(0, M_host, 0.01)


    phi_1 = []
    for m in m_for_integrate:
        phi_1.append(vdB_USHMF_1st_order(Params, m, M_subhalos_range))
    phi_1 = np.array(phi_1)

    phi_i_1_for_integrate = vdB_USHMF_1st_order(Params, M_host, m_for_integrate)
    phi_i_1 = np.zeros(M_subhalos_range.size)

    if order > 1:
        for i in range(order-1):
            for j in range(phi_i_1.size):
                phi_i_1[j] = trapz( phi_1[:,j] * phi_i_1_for_integrate * m_for_integrate / M_subhalos_range[j] , m_for_integrate)
            interpolator = interp1d(M_subhalos_range, phi_i_1, fill_value="extrapolate")
            phi_i_1_for_integrate = interpolator(m_for_integrate)

    return phi_i_1
