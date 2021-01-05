"""
@ file info.py

Written by Hao fu

UNDER DEVELOPMENT
"""


import numpy as np
#import sys
import os
from tqdm import tqdm
import analysis
import other_functions as others

from cosmological_model import *
from merger import *
from stellar_mass_halo_mass_relation import *

from colossus.cosmology import cosmology
Cosmo = cosmology.getCurrent()
h = Cosmo.h


import ctypes
from numpy.ctypeslib import ndpointer

analyze = ctypes.CDLL("dream/C_functions/analyze.so")

analyze.get_mask_in_range_.argtypes = [ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                      ctypes.c_double,
                                      ctypes.c_double,
                                      ctypes.c_int]

analyze.array1_elements_in_array2_.argtypes = [ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                              ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                              ctypes.c_int,
                                              ctypes.c_int]

star_formation = ctypes.CDLL("dream/C_functions/star_formation/star_formation.so")
star_formation.compute_star_formation.argtypes = [POINTER(stellar_mass_halo_mass),
                                                  #ctypes.c_char_p,
                                                  #ctypes.c_double,
                                                  ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                                  ctypes.c_int,
                                                  ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                                  ctypes.c_int,
                                                  ctypes.c_char_p,
                                                  POINTER(DM_catalogue),
                                                  POINTER(cosmological_parameters),
                                                  POINTER(cosmological_time)]



def analyse_mergers(input_params_info):

    data_folder = input_params_info.folder_name+"data/"

    id_halo, halo_catalog, id_mergers, mergers_array, z_infall, merging_timescale, z_at_merge, orders = others.read_data(input_params_info.folder_name, input_params_info.exist_z_at_merge)

    if input_params_info.want_unevolved_shmf or input_params_info.want_evolved_shmf:
        analyze.get_mask_in_range_.restype = ndpointer(dtype=ctypes.c_bool, shape = id_halo.shape)
        analyze.array1_elements_in_array2_.restype = ndpointer(dtype=ctypes.c_bool, shape = id_mergers.shape)

    # Z INFALL DISTRIBUTIONS TOTAL
    if input_params_info.want_z_infall_pdf:
        analysis.get_z_infall_distribution(z_infall,
                                           data_folder,
                                           "pdf_redshift_total.txt")

    if input_params_info.want_unevolved_shmf or input_params_info.want_evolved_shmf or input_params_info.want_mergers_rate:

        print("\nCalculating mass function and merger rates...")
        #for M_h in halo_masses:
        for q in tqdm(range(len(input_params_info.halo_masses))):

            M_h = input_params_info.halo_masses[q]
            bin = 0.2
            sub_mass_params = (10., M_h+bin, bin)
            Msh_bins = np.arange(sub_mass_params[0], sub_mass_params[1], sub_mass_params[2])

            """
            # OLD METHOLOGY, KEEP IT HERE, NEED TO FURTHER TEST THE SPEED
            idx_1 = ma.array(halo_catalog, mask = halo_catalog >= M_h).recordmask
            idx_2 = ma.array(halo_catalog, mask = halo_catalog < M_h+0.1).recordmask
            idx_halo = np.logical_and(idx_1, idx_2)

            idx_sub = np.array([], dtype=bool)
            for i, id in enumerate(id_mergers):
                idx_sub = np.append(idx_sub, id in id_halo[idx_halo])"""

            idx_halo = analyze.get_mask_in_range_(halo_catalog,
                                                 M_h, M_h+bin,
                                                 halo_catalog.size)

            idx_sub = analyze.array1_elements_in_array2_(id_mergers,
                                                        id_halo[idx_halo],
                                                        id_mergers.size,
                                                        id_halo.size)

            N_halos = len(halo_catalog[idx_halo])

            #UN-EVOLVED SHMF
            if input_params_info.want_unevolved_shmf:
                analysis.get_unevolved_shmf(mergers_array, z_at_merge,
                                            M_h, Msh_bins, bin,
                                            N_halos, idx_sub,
                                            data_folder)


            # Z INFALL DISTRIBUTIONS
            #order_mask = ma.array(orders, mask=orders>=1).recordmask
            #mask = np.logical_and(order_mask, idx_sub)
            analysis.get_z_infall_distribution(z_infall[idx_sub],
                                               data_folder,
                                               "pdf_redshift_{:.1f}.txt".format(M_h))

            print("ciao\n")

            #EVOLVED SHMF
            if input_params_info.want_evolved_shmf:
                z_ = 0.
                analysis.get_evolved_shmf(z_,
                                          mergers_array,
                                          z_infall,
                                          z_at_merge,
                                          M_h, Msh_bins,
                                          N_halos, idx_sub,
                                          data_folder)

            #MERGERS RATE
            if input_params_info.want_mergers_rate:
                analysis.get_mergers_rate(mergers_array,
                                          z_infall,
                                          M_h, N_halos, idx_sub,
                                          data_folder)


    #STAR FORMATION RATE
    if input_params_info.compute_SF:

        z_for_interp = np.linspace(0,20,10000)
        t_for_interp = Cosmo.lookbackTime(z_for_interp)
        age_for_interp = Cosmo.age(z_for_interp)
        z_range = np.array([])
        cosmo_params = [Cosmo.Om0, Cosmo.Ob0, Cosmo.sigma8, Cosmo.ns, Cosmo.h, Cosmo.H0,
                        Cosmo.Om(z_range), Cosmo.Ob(z_range), Cosmo.Hz(z_range), z_range, z_range.size]
        cosmo_params = cosmological_parameters(*cosmo_params)
        cosmo_time = [z_for_interp, t_for_interp, age_for_interp, z_for_interp.size]
        cosmo_time = cosmological_time(*cosmo_time)

        DM_catalog = [id_halo.size, id_halo, halo_catalog, id_mergers.size, id_mergers, np.array([]), mergers_array, z_infall, np.array([]), z_at_merge]
        DM_catalog = DM_catalogue(*DM_catalog)

        star_formation.compute_star_formation(input_params_info.SMHM,
                                              #create_string_buffer(input_params_info.SMHM_model.encode('utf-8')),
                                              #input_params_info.SMHM.scatter,
                                              np.array(input_params_info.stellar_masses_for_SF),
                                              len(input_params_info.stellar_masses_for_SF),
                                              np.array(input_params_info.redshifts_for_SFR),
                                              len(input_params_info.redshifts_for_SFR),
                                              create_string_buffer(data_folder.encode('utf-8')),
                                              DM_catalog,
                                              cosmo_params,
                                              cosmo_time)

    return
