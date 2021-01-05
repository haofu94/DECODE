import numpy as np
from scipy.interpolate import interp1d
from ctypes import *
#import ctypes
from numpy.ctypeslib import ndpointer


class stellar_mass_halo_mass(Structure):

    _fields_ = [("is_analytical", c_int),
                ("model", c_char_p),
                ("file", c_char_p),
                ("Mhalo", POINTER(c_double)),
                ("Mstar", POINTER(c_double)),
                #("Mhalo", np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C")),
                #("Mstar", np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C")),
                ("length", c_int),
                ("constant_scatter", c_int),
                ("scatter", c_double)]

    def __init__(self, is_analytical, model, file, Mhalo, Mstar, length, constant_scatter, scatter):

        #print(Mhalo); print(Mstar); print(Mstar.size)
        self.is_analytical = is_analytical
        self.model = model.encode('utf-8')
        self.file = file.encode('utf-8')
        Mhalo = Mhalo.astype(np.float64); Mstar = Mstar.astype(np.float64)
        self.Mhalo = Mhalo.ctypes.data_as(POINTER(c_double))
        self.Mstar = Mstar.ctypes.data_as(POINTER(c_double))
        #self.Mhalo = Mhalo
        #self.Mhalo = Mstar
        self.length = length
        self.constant_scatter = constant_scatter
        self.scatter = scatter




class SMHM_matrix(Structure):

    _fiels = [("rows", c_int),
              ("cols", c_int),
              ("matrix", POINTER(c_double)),
              ("redshift", POINTER(c_double)),
              ("Mstar", POINTER(c_double))]

    def __init__(self, rows, cols, matrix, redshift, Mstar):

        self.rows = rows
        self.cols = cols
        self.matrix = matrix.ctypes.data_as(POINTER(c_double))
        self.redshift = redshift.ctypes.data_as(POINTER(c_double))
        self.Mstar = Mstar.ctypes.data_as(POINTER(c_double))




def SMHM_read_matrix(filename):

    rows = int(np.loadtxt(filename, usecols=0, max_rows=1))
    cols = int(np.loadtxt(filename, usecols=1, max_rows=1))
    print(rows, cols)
    redshift = np.loadtxt(filename, usecols=range(0,cols-1), skiprows=1, max_rows=1)
    Mstar = np.loadtxt(filename, usecols=-1, skiprows=2)
    matrix_2D = np.loadtxt(filename, usecols=range(0,cols-1), skiprows=2)
    matrix = np.zeros((rows-1)*(cols-1))
    for i in range(rows-1):
        for j in range(cols-1):
            matrix[i*(cols-1)+j] = matrix_2D[i][j]

    smhm_data = [rows, cols, matrix, redshift, Mstar]
    smhm_data = SMHM_matrix(*smhm_data)

    return smhm_data



def SMHM_for_z(z, matrix):

    z_range = np.arange(0., 1.6, 0.1)
    #Mstar = matrix[:,-1]
    Mhalo = np.array([], dtype=float)

    for i in range(matrix[:,0].size):
        Mhalo = np.append(Mhalo, interp1d(z_range, matrix[i,:-1])(z))

    #return Mhalo, Mstar
    return Mhalo



def get_SMHM_numerical(SMHM_is_analytical, SMHM_file, SMHM_model, constant_scatter, scatter, z):

    if constant_scatter:
        constant_scatter = 1
    elif not constant_scatter:
        constant_scatter = 0

    if SMHM_is_analytical:
        is_analytical = 1
        Mhalo = np.array([], dtype=float)
        Mstar = np.array([], dtype=float)
    elif not SMHM_is_analytical:
        is_analytical = 0
        SMHM_matrix = np.loadtxt(SMHM_file, skiprows=2)
        Mstar = np.loadtxt(SMHM_file, usecols=-1, skiprows=2) #SMHM_matrix[:,-1]
        #Mhalo = SMHM_for_z(z, SMHM_matrix)
        z_range = np.arange(0., 1.6, 0.1)
        mhalo = []
        for i in range(SMHM_matrix[:,0].size):
            mhalo.append(interp1d(z_range, SMHM_matrix[i,:-1])(z))
        Mhalo = np.array(mhalo.copy())

    #print(Mhalo); print(Mstar)

    SMHM = [is_analytical, SMHM_model, SMHM_file, Mhalo, Mstar, Mhalo.size, constant_scatter, scatter]
    SMHM = stellar_mass_halo_mass(*SMHM)

    return SMHM
