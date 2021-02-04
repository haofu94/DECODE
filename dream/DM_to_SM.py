"""
@ file DM_to_SM.py

The goal of this module is to assign galaxy stellar mass to dark matter haloes.
"""

import numpy as np
from scipy.interpolate import interp1d
#import sys
from bisect import bisect

"""
from input_parameters import read_cosmo_model
try:
    from colossus.cosmology import cosmology
    cosmology = cosmology.getCurrent()
except:
    cosmo_str = read_cosmo_model(sys.argv[2:])
    cosmology.setCosmology(cosmo_str)
    cosmology = cosmology.getCurrent()"""



def SMHM_Grylls_param(DM, Params, z, scatter):

    """
    Grylls parametrization for the SMHM relation
    Centered in z = 0.1
    """

    M10, M11, SHMnorm10, SHMnorm11, beta10, beta11, gamma10, gamma11 = Params
    zparameter = np.divide(z-0.1, z+1)
    M = M10 + M11*zparameter
    N = SHMnorm10 + SHMnorm11*zparameter
    b = beta10 + beta11*zparameter
    g = gamma10 + gamma11*zparameter

    SM = np.power(10, DM) * (2*N*np.power( (np.power(np.power(10,DM-M), -b) + np.power(np.power(10,DM-M), g)), -1))
    SM = np.log10(SM)

    if scatter !=0.:
        return np.random.normal(loc = SM, scale = scatter, size = np.shape(SM))
    else:
        return SM



def SMHM_Moster(DM, Params, z):

    """
    Moster SMHM relation
    """

    M10, M11, SHMnorm10, SHMnorm11, beta10, beta11, gamma10, gamma11, Scatter = Params
    zparameter = np.divide(z, z+1)
    M = M10 + M11*zparameter
    N = SHMnorm10 + SHMnorm11*zparameter
    b = beta10 + beta11*zparameter
    g = gamma10 + gamma11*zparameter

    SM = np.power(10, DM) * 2*N* np.power( (np.power(np.power(10,DM-M), -b) + np.power(np.power(10,DM-M), g)), -1)
    SM = np.log10(SM)

    if Scatter !=0.:
        return np.random.normal(loc = SM, scale = Scatter, size = np.shape(SM))
    else:
        return SM


def SMHM_Behroozi_2013(DarkMatter, z, Scatter):

    e = np.array([-1.777, -0.006, 0.000, -0.119])
    M = np.array([11.514,-1.793,-0.251])
    alpha = np.array([-1.412,0.731])
    delta = np.array([3.508,2.608,-0.043])
    gamma_b = np.array([0.361,1.391,0.279])
    ep = np.array([0.218,-0.023])

    a      = 1/(1+z)
    afac   = a-1
    v_     = np.exp(-4*np.power(a,2))
    M_     = M[0]       + (M[1]      *afac +M[2]*z      )*v_
    e_     = e[0]       + (e[1]      *afac +e[2]*z      )*v_ + e[3]*afac
    alpha_ = alpha[0]   + (alpha[1]  *afac              )*v_
    delta_ = delta[0]   + (delta[1]  *afac +delta[2]*z  )*v_
    gamma_ = gamma_b[0] + (gamma_b[1]*afac +gamma_b[2]*z)*v_
    ep_    = ep[0]      + (ep[1]*afac                   )


    e_ = np.power(10, e_)
    M_ = np.power(10, M_)

    def f(x, a = alpha_, d = delta_, g = gamma_):
        Part1 = np.log10(np.power(10, a*x) + 1)
        Part2 = np.power(np.log10(1+np.exp(x)),g)
        Part3 = 1 + np.exp(np.power(10, -x))
        return -Part1 + d*np.divide(Part2, Part3)

    Part1 = np.log10(e_*M_)
    Part2 = f( np.log10( np.divide(np.power(10,DarkMatter), M_) ) )
    Part3 = f(0)

    M_Star = Part1 + Part2 - Part3
    if Scatter!=0.:
        Scatter = np.random.normal(scale = Scatter, size = np.shape(M_Star))
        return M_Star + Scatter
    else:
        return M_Star


def SMHM_Behroozi(DM, z, scatter):

    #Behroozi et al. 2019, Obs. Q/SF Cen/Sat
    eps = np.array([-1.435, 1.831, 1.368, -0.217])
    M = np.array([12.035, 4.556, 4.417, -0.731])
    alpha = np.array([1.963, -2.316, -1.732, 0.178])
    beta = np.array([0.482, -0.841, -0.471])
    gamma = np.array([-1.034, -3.100, -1.055])
    d = 0.411

    #Behroozi et al. 2019, True Q/SF Cen/Sat
    """eps = np.array([-1.430, 1.796, 1.360, -0.216])
    M = np.array([12.040, 4.675, 4.513, -0.744])
    alpha = np.array([1.973, -2.353, 1.783, 0.186])
    beta = np.array([0.473, -0.884, -0.486])
    gamma = np.array([-1.088, -3.241, -1.079])
    d = 0.407"""

    a = 1./(1.+z)

    logM1 = M[0] + M[1]*(a-1.) - M[2]*np.log(a) + M[3]*z
    e = eps[0] + eps[1]*(a-1.) - eps[2]*np.log(a) + M[3]*z
    a = alpha[0] + alpha[1]*(a-1.) - alpha[2]*np.log(a) + alpha[3]*z
    b = beta[0] + beta[1]*(a-1.) + beta[2]*z
    g = np.power(10., gamma[0] + gamma[1]*(a-1.) + gamma[2]*z)

    x = DM - logM1

    SM = logM1 + e - np.log10(np.power(10., -a*x) + np.power(10., -b*x)) + g * np.exp(-0.5 * np.power(np.divide(x,d), 2))

    if scatter != 0.:
        SM += np.random.normal(0, scatter)

    return SM



def SMHM_matrix(DM, matrix, z_array, z, scatter):

    Mstar = matrix[:,-1]

    Mhalo = np.array([interp1d(z_array, matrix[i,:-1], fill_value="extrapolate")(z) for i in range(Mstar.size)])

    SM = interp1d(Mhalo, Mstar, fill_value="extrapolate")(DM) + np.random.normal(0, scatter)

    return SM



def SMHM_matrix_v2(DM, matrix, z_array, z, scatter):

    Mstar = matrix[:,-1]

    idx = bisect(z_array, z)

    Mhalo = matrix[:,idx]

    SM = interp1d(Mhalo, Mstar, fill_value="extrapolate")(DM) + np.random.normal(0, scatter)

    return SM
