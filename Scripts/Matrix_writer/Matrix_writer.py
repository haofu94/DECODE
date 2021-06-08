"""
Written by Max Dickson
"""

import numpy as np
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from colossus.halo import mass_adv
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
from scipy.special import erfc
import sys
import matplotlib.pyplot as plt

plt.ion()

cosmo = cosmology.setCosmology("planck18")

import warnings
warnings.filterwarnings("ignore")

from time import time


def HMF(z , bins = 0.05 ,cum_var = 1.0):

    #Produces a Halo Mass Function (HMF) using a Tinker et al. (2008) HMF and a correction for subhalos from Behroozi et al. (2013)

    M = np.arange(9. , 15.5+bins ,bins) #Sets HMF mass range

    if z<=1.:
        phi = mass_function.massFunction((10.0**M)*cosmo.h, z , mdef = 'vir' , model = 'tinker08' , q_out="dndlnM") * np.log(10) * cosmo.h**3.0 * cum_var #Produces the Tinker et al. (2008) HMF
    else:
        phi = mass_function.massFunction((10.0**M)*cosmo.h, z , mdef = '200m' , model = 'tinker08' , q_out="dndlnM") * np.log(10) * cosmo.h**3.0 * cum_var
        M = np.log10(mass_adv.changeMassDefinitionCModel(10.**M, z, '200m','vir'))[0]

    a = 1./(1.+z) #Constants for Behroozi et al. (2013) Appendix G, eqs. G6, G7 and G8 subhalos correction
    C = np.power(10., -2.415 + 11.68*a - 28.88*a**2 + 29.33*a**3 - 10.56*a**4)
    logMcutoff = 10.94 + 8.34*a - 0.36*a**2 - 5.08*a**3 + 0.75*a**4
    correction = phi * C * (logMcutoff - M)

    return M , np.log10(phi + correction)

def derivative(x,y):

    #Returns the derivative of any input function

    func = interp1d(x,y,fill_value="extrapolate")
    dx = 0.1
    x_calc = np.arange(x[0],x[-1]+dx,dx)
    y_calc = func(x_calc)
    dydx = np.diff(y_calc)/np.diff(x_calc)

    dydx = np.append(dydx, dydx[-1]) #Preventing boundary abnormalities in the returned function
    dydx_low = np.mean(dydx[:10])
    dydx[0] = dydx_low
    dydx[1] = dydx_low
    dydx_high = np.mean(dydx[-10:])

    return interp1d(x_calc,dydx,fill_value=(dydx_low,dydx_high),bounds_error = False)

def integrals(SMF , HMF , sig):

    #Returns the integrals of the Stellar Mass Function (SMF) and HMF as calculated in Aversa et al. (2015) eq. 37

    M_s , phi_s = SMF[0] , SMF[1]
    M_h , phi_h = HMF[0] , HMF[1]
    phi_s , phi_h = 10.0**phi_s , 10.0**phi_h

    I_phi_s = np.flip(cumtrapz(np.flip(phi_s), M_s))

    I_phi_h = np.array([])
    for m in M_h:
        I = np.trapz(phi_h*0.5*erfc((m-M_h)/(np.sqrt(2)*sig)) , M_h)
        I_phi_h = np.append(I_phi_h,I)

    M_s , M_h = M_s[:-1]+0.025 , M_h+0.025

    return I_phi_s , I_phi_h , M_s , M_h

def SMFfromSMHM(M_s , M_h , sig , z , M_hforsigma = None):

    #Reconstructs the SMF from the Stellar Mass Halo Mass (SMHM) relationship

    bins = 0.1
    volume = (200*cosmo.h)**3

    M_hmf , phi = HMF(z , bins , cum_var = volume)
    M_hmf = M_hmf #Cutting the low mass end of the HMF for increased speed
    phi = 10**phi

    cum_phi = np.cumsum(phi)
    max_number = np.floor(np.trapz(phi, M_hmf)) #np.floor(np.max(cum_phi))
    #plt.plot(M_hmf, np.log10(phi))
    #print(max_number, np.floor(np.max(cum_phi)))
    #sys.exit()
    #if (np.random.uniform(0,1) > np.max(cum_phi)-max_number): #Calculating number of halos to compute
    if (np.random.uniform(0,1) > np.trapz(phi, M_hmf)-max_number): #Calculating number of halos to compute
        max_number += 1

    int_cum_phi = interp1d(cum_phi, M_hmf)
    range_numbers = np.random.uniform(np.min(cum_phi), np.max(cum_phi), int(max_number))
    halo_masses = int_cum_phi(range_numbers)

    if M_hforsigma is not None:
        #interpolating the values of sigma over the halo masses if M_hforsigma is given
        sig = interp1d(M_hforsigma,sig,fill_value="extrapolate")(halo_masses)

    M_smf = np.arange(6., 12.5 , 0.1) #SMF histogram bins
    SMHMinterp = interp1d(M_h , M_s , fill_value="extrapolate")
    stellar_masses = SMHMinterp(halo_masses) + np.random.normal(0., sig, halo_masses.size) #Calculating stellar masses using SMHM with scatter
    phi_smf = np.histogram(stellar_masses , bins = M_smf)[0]/0.1/volume #SMF histogram

    return M_smf[:-1]+0.05 , np.log10(phi_smf)

def writematrix(SMF , sigma , z_range , title , M_hforsigma = None):

    #Writes a SMHM matrix to produce a SMHM for any redshift within z_range

    M_s = SMF[0]
    abun = np.loadtxt("SMHM_data/Baldry + Bernardi SMHM.txt") #Reads SMHM produced from abundance matching at z = 0.1
    matrix = np.empty((0,M_s.size-1))

    if type(sigma) == float:
        #Adaptation for the code to use constant and variable values of sigma
        sigma = np.repeat(sigma,z_range.size)

    for i in range(z_range.size):
        #Calculates SMHM for each redshift in z_range

        z = z_range[i]
        M_h , phi_h = HMF(z)#0.1

        if M_hforsigma is not None:
            #interpolating the values of sigma over the halo masses if M_hforsigma is given
            if sigma.shape[1] == 1:
                sig = interp1d(M_hforsigma,sigma[:,0],fill_value="extrapolate")(M_h)
            else:
                sig = interp1d(M_hforsigma,sigma[:,i],fill_value="extrapolate")(M_h)
        else:
            sig = sigma[i]

        if SMF.shape[0] == 2:
            phi_s = SMF[1]
        else:
            phi_s = SMF[i+1]

        deri = derivative(abun[:,0],abun[:,1])(M_h)
        n=0
        e=1.

        while n < 3:
            #Calls the integrals function and matches integral values of SMF to interpolated values in the integral of HMF, iterating multiple times
            I_phi_s , I_phi_h , M_s_temp , M_h_temp = integrals(np.array([M_s , phi_s]) , np.array([M_h , phi_h]) , sig/deri)
            int_I_phi_h = interp1d(I_phi_h , M_h_temp , fill_value="extrapolate")
            M_h_match = np.array([])
            for m in range(M_s_temp.size):
                M_h_match = np.append(M_h_match , int_I_phi_h(I_phi_s[m]))

            time_ini = time()
            if M_hforsigma is not None:
                M_s_iter , phi_s_iter = SMFfromSMHM(M_s_temp , M_h_match , sig , z , M_h)
            else:
                M_s_iter , phi_s_iter = SMFfromSMHM(M_s_temp , M_h_match , sig , z) #Reconstructing SMF
            time_fin = time()
            print("time in seconds: ", time_fin-time_ini)
            int_phi_s_iter = interp1d(M_s_iter , phi_s_iter , fill_value="extrapolate")
            #e_temp = max((phi_s - int_phi_s_iter(M_s))/phi_s) #Calculating relative error between reconstructed SMF and the input SMF
            e_array = np.abs((phi_s - int_phi_s_iter(M_s))/phi_s)
            mask = np.where(np.isfinite(e_array))
            e_temp = max(e_array[mask]) #Calculating relative error between reconstructed SMF and the input SMF
            #print(M_s[np.argmax(e_array[mask])], e_temp)
            #plt.figure(1); plt.clf()
            #plt.plot(M_s, phi_s); plt.plot(M_s, int_phi_s_iter(M_s))
            #print(M_s)
            #print(e_array)

            if e_temp < e:
                #Only accepts iterations that decrease the relative error
                e = e_temp
                deri = derivative(M_h_match,M_s_temp)(M_h)
            n += 1
        matrix = np.append(matrix,[M_h_match], axis=0) #Appending each halo mass array from the SMHM at each redshift to the matrix array

        percent = float(i+1)*100./float(z_range.size)
        sys.stdout.write("\r{per:.2f}%".format(per = percent))
        sys.stdout.flush()

    matrix = np.append(matrix,[M_s_temp], axis=0) # Appending the stellar mass array to the end of the matrix array

    with open("SMHM_data/"+title+".txt" ,"w+") as file:
        #Writes the matrix to a file
        np.savetxt(file,np.array([M_s.size+1,z_range.size+1])[None],fmt="%i") #Writing the rows and columns of the file
        np.savetxt(file,np.append(z_range,np.inf)[None],fmt = "%f") #Writing the redshift of each SMHM
        np.savetxt(file,matrix.T) #Writing the SMHM matrix


if __name__ == "__main__":

    SMF = np.transpose(np.loadtxt("baldry_bernardi_SMF.txt"))
    for i in range(SMF.shape[0]-1):
        SMF[i+1,:] = np.log10(SMF[i+1,:])

    sigma = 0.15
    z_range = np.arange(0, 6., 0.1)
    #z_range = np.array([6.])
    title = "SMHM_out"

    writematrix(SMF, sigma , z_range , title)
