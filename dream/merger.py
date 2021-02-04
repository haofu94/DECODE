import numpy as np
from ctypes import *


class mergers_parameters(Structure):

    _fields_ = [("id", c_int),
                ("halo_mass_at_z0", c_double),
                ("halo_mass_at_z", c_double),
                ("halo_mass_at_z_plus_dz", c_double),
                ("redshift", c_double),
                ("redshift_bin", c_double),
                ("redshift_max", c_double),
                ("subhalo_mass_range", POINTER(c_double)),
                ("subhalo_mass_bin", c_double),
                ("length", c_int),
                ("type_orbital_circularity", c_int),
                ("orbital_circularity", c_double),
                ("fudge", c_double),
                ("max_order", c_int)]

    def __init__(self, id, halo_mass_at_z0, halo_mass_at_z, halo_mass_at_z_plus_dz,
                 redshift, redshift_bin, redshift_max, sub_mass_params, type_orbital_circularity, orbital_circularity, fudge, max_order):

        subhalo_mass_range = np.arange(sub_mass_params[0], sub_mass_params[1] + \
                                sub_mass_params[2], sub_mass_params[2]) #[log(M/Msun)]

        #print(subhalo_mass_range)

        self.id = id
        self.halo_mass_at_z0 = halo_mass_at_z0
        self.halo_mass_at_z = halo_mass_at_z
        self.halo_mass_at_z_plus_dz = halo_mass_at_z_plus_dz
        self.redshift = redshift
        self.redshift_bin = redshift_bin
        self.redshift_max = redshift_max
        self.subhalo_mass_range = subhalo_mass_range.ctypes.data_as(POINTER(c_double))
        self.subhalo_mass_bin = sub_mass_params[2]
        self.length = subhalo_mass_range.size
        self.type_orbital_circularity = type_orbital_circularity
        self.orbital_circularity = orbital_circularity
        self.fudge = fudge
        self.max_order = max_order


class subhalo_mass_functions(Structure):

    _fields_ = [("total", POINTER(c_double)),
                ("first", POINTER(c_double)),
                ("second", POINTER(c_double)),
                ("third", POINTER(c_double)),
                ("fourth", POINTER(c_double)),
                ("fifth", POINTER(c_double)),
                ("psi", POINTER(c_double)),
                ("length_psi", c_int)]

    def __init__(self, total, first, second, third, fourth, fifth, psi, length_psi):

        self.total = total.ctypes.data_as(POINTER(c_double))
        self.first = first.ctypes.data_as(POINTER(c_double))
        self.second = second.ctypes.data_as(POINTER(c_double))
        self.third = third.ctypes.data_as(POINTER(c_double))
        self.fourth = fourth.ctypes.data_as(POINTER(c_double))
        self.fifth = fifth.ctypes.data_as(POINTER(c_double))
        self.psi = psi.ctypes.data_as(POINTER(c_double))
        self.length_psi = length_psi


class DM_catalogue(Structure):

    _fields_ = [("len_parents", c_int),
                ("id_parents", POINTER(c_int)),
                ("mass_parents", POINTER(c_double)),
                ("len_mergers", c_int),
                ("id_mergers", POINTER(c_int)),
                ("order_mergers", POINTER(c_int)),
                ("mass_mergers", POINTER(c_double)),
                ("z_infall", POINTER(c_double)),
                ("tau_merge", POINTER(c_double)),
                ("z_at_merge", POINTER(c_double))]

    def __init__(self, len_parents, id_parents, mass_parents, len_mergers, id_mergers, order_mergers, mass_mergers, z_infall, tau_merge, z_at_merge):

        self.len_parents = len_parents
        self.id_parents = id_parents.ctypes.data_as(POINTER(c_int))
        self.mass_parents = mass_parents.ctypes.data_as(POINTER(c_double))
        self.len_mergers = len_mergers
        self.id_mergers = id_mergers.ctypes.data_as(POINTER(c_int))
        self.order_mergers = order_mergers.ctypes.data_as(POINTER(c_int))
        self.mass_mergers = mass_mergers.ctypes.data_as(POINTER(c_double))
        self.z_infall = z_infall.ctypes.data_as(POINTER(c_double))
        self.tau_merge = tau_merge.ctypes.data_as(POINTER(c_double))
        self.z_at_merge = z_at_merge.ctypes.data_as(POINTER(c_double))
