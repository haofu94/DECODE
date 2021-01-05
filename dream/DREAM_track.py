"""
@ file DREAM_track.py

Written by:
Hao Fu, Chris Marsden
"""

import numpy as np
import sys
import os
from tqdm import tqdm
from colossus.halo.mass_so import M_to_R


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


import halo_growth
import analytic_SHMFs
import other_functions as others

from cosmological_model import *
from merger import *

from ctypes import *
from numpy.ctypeslib import ndpointer

centrals = CDLL("dream/C_functions/dark_matter/centrals.so")
centrals.get_halo_index.argtypes = [ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                    c_int,
                                    ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                    c_int]
centrals.get_halo_IDs.argtypes = [c_int, c_int, c_int]

mergers = CDLL("dream/C_functions/dark_matter/mergers.so")
mergers.generate_mergers.argtypes = [POINTER(mergers_parameters),
                                     ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                     ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                     c_int,
                                     c_int,
                                     POINTER(cosmological_parameters),
                                     #ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                     c_char_p,
                                     c_int,
                                     POINTER(cosmological_time),
                                     POINTER(subhalo_mass_functions)]



def DREAM_track(halo_catalog,
                 input_params_run):
    """
    DREAM Tracks

    This function generates all the (currently implemented) data for DREAM, for all central halos.
    This function can be parallelized.

    Args:
        halo_catalog (list): The mass of the centrals at z = 0, in [log10 Msun/h]
        z_range (array/numpy array/list): The redshifts at which to generate the information [dimensionless]
        subhalo_range (array/numpy array/list): The bins upon which the subhahalo masses will be sampled.
        output_folder
        want_z_at_merge: state where want to get the redshift at merging
    Returns:
        none
    """


    type_orbital_circularity = input_params_run.merging_timescale_params[0]
    orbital_circularity = input_params_run.merging_timescale_params[1]
    fudge = input_params_run.merging_timescale_params[2]

    print("Calculating accretion tracks...")

    if input_params_run.use_mean_track:

        min_mass = 9.; max_mass = 17.
        mass_range = np.arange(min_mass, max_mass+0.05, 0.05)
        sample_tracks = []

        for i,m in enumerate(mass_range):
            sample_tracks.append(halo_growth.Mass_acc_history_VDB_FS(m, input_params_run.z_range, Cosmo.h, Cosmo.Om0))
        sample_tracks = np.array(sample_tracks)

        centrals.get_halo_index.restype = ndpointer(dtype=c_int, shape = np.array(halo_catalog).shape)
        halo_index = centrals.get_halo_index(np.array(halo_catalog),
                                             len(halo_catalog),
                                             mass_range,
                                             mass_range.size)

    elif not use_mean_track:

        accretion_tracks = []
        for i in tqdm(range(len(halo_catalog))):
            i_th_track = halo_growth.Mass_acc_history_VDB_FS(halo_catalog[i], input_params_run.z_range, Cosmo.h, Cosmo.Om0)
            accretion_tracks.append(i_th_track)

    print("Accretion tracks calculated")


    # Generate random ID
    print("\nGenerating halo IDs...")
    id_list = []
    Num_haloes = len(halo_catalog)
    Num_digits = int(np.floor(np.log10(Num_haloes)))+2
    low = np.power(10, Num_digits-1); high = np.power(10, Num_digits) - 1

    centrals.get_halo_IDs.restype = ndpointer(dtype=c_int, shape = np.array(halo_catalog).shape)
    id_list = centrals.get_halo_IDs(Num_haloes, low, high)
    print("Halo IDs generated")


    if os.path.isfile(input_params_run.output_folder + "data/output_mergers.txt"):
        os.remove(input_params_run.output_folder + "data/output_mergers.txt")
    if os.path.isfile(input_params_run.output_folder + "data/output_parents.txt"):
        os.remove(input_params_run.output_folder + "data/output_parents.txt")

    # print parents to file and create header for mergers file
    others.print_parents(input_params_run.output_folder, id_list, halo_catalog)
    others.print_mergers_header(input_params_run.output_folder, input_params_run.want_z_at_merge)

    if input_params_run.use_merger_tree == True:
        input_params_run.use_merger_tree = 1
    elif input_params_run.use_merger_tree == False:
        input_params_run.use_merger_tree = 0

    if input_params_run.want_z_at_merge == "yes":
        input_params_run.want_z_at_merge = 1
        z_for_interp = np.linspace(0,20,10000)
        t_for_interp = Cosmo.lookbackTime(z_for_interp)
        age_for_interp = Cosmo.age(z_for_interp)
    elif input_params_run.want_z_at_merge == "no":
        input_params_run.want_z_at_merge = 0
        z_for_interp = np.array([])
        t_for_interp = np.array([])
        age_for_interp = np.array([])

    mergers_params = [0, 0., -1., -1., -1., input_params_run.z_range[1] - input_params_run.z_range[0], \
                        input_params_run.z_range[-1], input_params_run.sub_mass_params, type_orbital_circularity, orbital_circularity, fudge, input_params_run.max_order]
    mergers_params = mergers_parameters(*mergers_params)

    cosmo_params = [Cosmo.Om0, Cosmo.Ob0, Cosmo.sigma8, Cosmo.ns, Cosmo.h, Cosmo.H0,
                    Cosmo.Om(input_params_run.z_range), Cosmo.Ob(input_params_run.z_range), Cosmo.Hz(input_params_run.z_range), input_params_run.z_range, input_params_run.z_range.size]
    cosmo_params = cosmological_parameters(*cosmo_params)

    cosmo_time = [z_for_interp, t_for_interp, age_for_interp, z_for_interp.size]
    cosmo_time = cosmological_time(*cosmo_time)

    ##################################
    # get sub halo mass functions
    Params = [0.22, -0.91, 6., 3., 1.] #Jiang \& van den Bosch 2016 Table A1, Unevolved, total
    Params_1st_order = [0.13, -0.83, 1.33, -0.02, 5.67, 1.19]
    psi = 10**np.arange(8,16,0.1)/10**14
    shmf_all = analytic_SHMFs.vdB_USHMF(Params, 14, np.arange(8,16,0.1))
    shmf_1st = analytic_SHMFs.vdB_USHMF_1st_order(Params_1st_order, 14, np.arange(8,16,0.1))
    shmf_2nd = analytic_SHMFs.vdB_USHMF_ith_order(Params_1st_order, 14, np.arange(8,16,0.1), 2)
    shmf_3rd = analytic_SHMFs.vdB_USHMF_ith_order(Params_1st_order, 14, np.arange(8,16,0.1), 3)
    shmf_4th = analytic_SHMFs.vdB_USHMF_ith_order(Params_1st_order, 14, np.arange(8,16,0.1), 4)
    shmf_5th = analytic_SHMFs.vdB_USHMF_ith_order(Params_1st_order, 14, np.arange(8,16,0.1), 5)
    #shmf_all = shmf_1st + shmf_2nd + shmf_3rd
    SHMF = [shmf_all, shmf_1st, shmf_2nd, shmf_3rd, shmf_4th, shmf_5th, psi, psi.size]
    SHMF = subhalo_mass_functions(*SHMF)
    ##################################

    print("\nGenerating mergers...")
    mergers.generate_mergers.restype = c_int
    Num_mergers = 0
    for i in tqdm(range(len(halo_catalog))):

        if input_params_run.use_mean_track:
            assert min_mass<=halo_catalog[i]<=max_mass, "Halo mass out of range. Min: log(M/Msun) = 9. Max: log(M/Msun) = 17."
            track = sample_tracks[halo_index[i]]
        elif not input_params_run.use_mean_track:
            track = accretion_tracks[i]

        mergers_params.id = id_list[i]
        mergers_params.halo_mass_at_z0 = track[0]

        N = mergers.generate_mergers(mergers_params,
                                     track,
                                     input_params_run.z_range,
                                     track.size,
                                     input_params_run.use_merger_tree,
                                     cosmo_params,
                                     #M_to_R(10**track, z_range, 'vir'),
                                     create_string_buffer(input_params_run.output_folder.encode('utf-8')),
                                     input_params_run.want_z_at_merge,
                                     cosmo_time,
                                     SHMF)

        Num_mergers += N

    others.print_log_file(input_params_run.output_folder, len(halo_catalog), Num_mergers, input_params_run.want_z_at_merge)
    print("\nData stored in:\n{}".format(os.getcwd() + "/" + input_params_run.output_folder + "data/"))

    return
