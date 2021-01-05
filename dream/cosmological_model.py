from ctypes import *


class cosmological_parameters(Structure):

    _fields_ = [("Om0", c_double),
                ("Ob0", c_double),
                ("sigma8", c_double),
                ("ns", c_double),
                ("h", c_double),
                ("H0", c_double),
                ("Om", POINTER(c_double)),
                ("Ob", POINTER(c_double)),
                ("Hz", POINTER(c_double)),
                ("z", POINTER(c_double)),
                ("length", c_int)]

    def __init__(self, Om0, Ob0, sigma8, ns, h, H0, Om, Ob, Hz, z, length):

        self.Om0 = Om0
        self.Ob0 = Ob0
        self.sigma8 = sigma8
        self.ns = ns
        self.h = h
        self.H0 = H0
        self.Om = Om.ctypes.data_as(POINTER(c_double))
        self.Ob = Ob.ctypes.data_as(POINTER(c_double))
        self.Hz = Hz.ctypes.data_as(POINTER(c_double))
        self.z = z.ctypes.data_as(POINTER(c_double))
        self.length = length


class cosmological_time(Structure):

    _fields_ = [("redshift", POINTER(c_double)), ("lookback_time",POINTER(c_double)), ("age",POINTER(c_double)), ("length",c_int)]

    def __init__(self, redshift, lookback_time, age, length):

        self.redshift = redshift.ctypes.data_as(POINTER(c_double))
        self.lookback_time = lookback_time.ctypes.data_as(POINTER(c_double))
        self.age = age.ctypes.data_as(POINTER(c_double))
        self.length = length
