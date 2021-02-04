import numpy as np
import matplotlib.pyplot as plt

plt.ion()

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



Mpeak = np.linspace(10, 15, 1000)

Mstar = SMHM_Behroozi(Mpeak, 0.1, 0)


plt.figure(1); plt.clf()
plt.plot(Mpeak, Mstar)
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.xlabel(r"$\log (M_{\rm peak} / M_\odot)$", fontsize=15)
plt.ylabel(r"$\log (M_* / M_\odot)$", fontsize=15)
plt.rc('font', family="serif")
