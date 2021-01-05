"""
@ file DREAM_galaxies.py

Written by:
Hao Fu
"""

import numpy as np
import sys
import os
from tqdm import tqdm


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


from cosmological_model import *
from merger import *
from other_functions import *
from stellar_mass_halo_mass_relation import *

from ctypes import *
from numpy.ctypeslib import ndpointer


central_galaxies = CDLL("dream/C_functions/baryonic_matter/central_galaxies.so")

central_galaxies.get_central_galaxies.argtypes = [c_char_p, c_char_p, c_char_p, c_int, POINTER(stellar_mass_halo_mass)]#c_char_p]
                                                  #ndpointer(np.float64, flags="C_CONTIGUOUS")]

satellite_galaxies = CDLL("dream/C_functions/baryonic_matter/satellite_galaxies.so")

satellite_galaxies.get_satellite_galaxies_at_z.argtypes = [c_char_p, c_char_p, c_char_p, c_int, c_int,
                                                           #c_char_p,
                                                           POINTER(stellar_mass_halo_mass),
                                                           c_double, c_int, c_int, c_int,
                                                           POINTER(mergers_parameters),
                                                           POINTER(cosmological_time),
                                                           POINTER(cosmological_parameters)]


def DREAM_galaxies(input_params_run):

    logfile_name = input_params_run.output_folder + "output.log"
    file_name = input_params_run.output_folder + "data/output_parents.txt"
    len_centrals = int(np.loadtxt(logfile_name, usecols=1, max_rows=1))

    if os.path.isfile(input_params_run.output_folder + "data/output_centrals.txt"):
        os.remove(input_params_run.output_folder + "data/output_centrals.txt")


    print_centrals_header(input_params_run.output_folder)

    #SMHM_params = np.array([11.92, 0.58, 0.032, -0.014, 1.64, -0.69, 0.53, 0.03, 0.]) #PyMorph
    #SMHM = get_SMHM_numerical(SMHM_model = "Fu_Dickson_2020", scatter = 0., z = 0.1)

    #print(logfile_name)
    #smhm_data = SMHM_read_matrix("Data/SMHM_relations/SMHM_constant_sigma.txt")

    central_galaxies.get_central_galaxies(create_string_buffer(logfile_name.encode('utf-8')),
                                          create_string_buffer(file_name.encode('utf-8')),
                                          create_string_buffer(input_params_run.output_folder.encode('utf-8')),
                                          len_centrals,
                                          input_params_run.SMHM)
                                          #create_string_buffer(input_params_run.SMHM_model.encode('utf-8')))
                                          #SMHM_params)

    file_name = input_params_run.output_folder + "data/output_mergers.txt"
    len_mergers = int(np.loadtxt(logfile_name, usecols=1, skiprows=3, max_rows=1))

    if os.path.isfile(input_params_run.output_folder + "data/output_satellites.txt"):
        os.remove(input_params_run.output_folder + "data/output_satellites.txt")

    print_satellites_header(input_params_run.output_folder, input_params_run.satellites_redshift)

    type_orbital_circularity = input_params_run.merging_timescale_params[0]
    orbital_circularity = input_params_run.merging_timescale_params[1]
    fudge = input_params_run.merging_timescale_params[2]
    mergers_params = [0, 0., -1., -1., -1., input_params_run.z_range[1] - input_params_run.z_range[0], \
                        input_params_run.z_range[-1], input_params_run.sub_mass_params, type_orbital_circularity, orbital_circularity, fudge, input_params_run.max_order]
    mergers_params = mergers_parameters(*mergers_params)

    cosmo_params = [Cosmo.Om0, Cosmo.Ob0, Cosmo.sigma8, Cosmo.ns, Cosmo.h, Cosmo.H0,
                    Cosmo.Om(input_params_run.z_range), Cosmo.Ob(input_params_run.z_range), Cosmo.Hz(input_params_run.z_range), input_params_run.z_range, input_params_run.z_range.size]
    cosmo_params = cosmological_parameters(*cosmo_params)

    input_params_run.want_z_at_merge = 1
    z_for_interp = np.linspace(0,20,10000)
    t_for_interp = Cosmo.lookbackTime(z_for_interp)
    age_for_interp = Cosmo.age(z_for_interp)
    cosmo_time = [z_for_interp, t_for_interp, age_for_interp, z_for_interp.size]
    cosmo_time = cosmological_time(*cosmo_time)
    #print("ciao 3")

    satellite_galaxies.get_satellite_galaxies_at_z(create_string_buffer(logfile_name.encode('utf-8')),
                                                   create_string_buffer(file_name.encode('utf-8')),
                                                   create_string_buffer(input_params_run.output_folder.encode('utf-8')),
                                                   len_mergers,
                                                   len_centrals,
                                                   #create_string_buffer(input_params_run.SMHM_model.encode('utf-8')),
                                                   input_params_run.SMHM,
                                                   input_params_run.satellites_redshift,
                                                   input_params_run.ignore_high_orders,
                                                   input_params_run.include_quenching,
                                                   input_params_run.include_stripping,
                                                   mergers_params,
                                                   cosmo_time,
                                                   cosmo_params)

    #print("\nritornato in python")

    return
