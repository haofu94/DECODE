import numpy as np
import sys
from colossus.cosmology import cosmology

Cosmo = cosmology.getCurrent()
h = Cosmo.h


import ctypes
from numpy.ctypeslib import ndpointer

num_c = ctypes.CDLL("dream/C_functions/numerical/num_c.so")

num_c.read_data.argtypes = [ctypes.c_int,
                            ctypes.c_char_p,
                            ctypes.c_char_p,
                            ctypes.c_int]

mergers = ctypes.CDLL("dream/C_functions/dark_matter/mergers.so")

mergers.get_redshift_at_merging.argtypes = [ctypes.c_double,
                                  ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                  ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                  ndpointer(np.float64, flags="C_CONTIGUOUS"),
                                  ctypes.c_int,
                                  ctypes.c_int]


def initialize(cosmo_str, args):

    print("\nRunning DREAM (DiscRete sEmi-empiricAl Model)\n\nCosmological model set to: {}\n".format(cosmo_str))

    return args[1]


def print_command_error():

    ErrMsg = "\n  /!\    Sub-command not recognized :(\n / ! \   See README.md for help.\n"
    sys.exit(ErrMsg)


def print_parents(output_folder, id_list, halo_catalog):

    fp = open(output_folder+"data/output_parents.txt", "w")
    fp.write("# Catalogue of the parent halos.\n# Colums:\n# 1) ID number identifying each halo\n# 2) halo mass in units of log10[M/Msun]\n\n")
    for i, M in enumerate(halo_catalog):
        fp.write("{:d} {:.6f}\n".format(id_list[i], M))
    fp.close()

    return


def print_parents_with_mah(output_folder, id_list, accretion_tracks, z):

    fp = open(output_folder+"data/output_parents.txt", "w")
    fp.write("# Catalogue of the parent halos.\n")
    fp.write("# Colums and rows:\n")
    fp.write("# 1st col: ID number identifying each halo\n")
    fp.write("# 2nd col: halo mass in units of log10[M/Msun]\n")
    fp.write("# Following columns: halo mass accretion histories\n")
    fp.write("# 1st row: redshift for halo mass accretion histories\n")

    for i in range(z.size-1):
        fp.write("{:.6f} ".format(z[i]))
    fp.write("{:.6f}".format(z[-1]))

    for i in range(id_list.size):
        fp.write("\n{:d}".format(id_list[i]))
        for j in range(accretion_tracks[i].size):
            fp.write(" {:.6f}".format(accretion_tracks[i][j]))

    fp.close()

    return


def print_mergers_header(output_folder, want_z_at_merge):

    fp = open(output_folder+"data/output_mergers.txt", "w")

    if want_z_at_merge == "yes":
        fp.write("# Mergers information.\n# Column:\n# 1) ID number identifying the corresponding parent halo\n# 2) ID identifying the (i-1)-th order progenitor\n# 3) ID identifying the subhalo itself\n# 4) order of the subhalo\n# 5) subhalo mass in units of log10[M/Msun]\n# 6) redshift at first accretion\n# 7) merging timescale in units of [Gyr]\n# 8) redshift at full merging, set to -1 if time at merging is later than today\n\n")
    elif want_z_at_merge == "no":
        fp.write("# Mergers information.\n# Column:\n# 1) ID number identifying the corresponding parent halo\n# 2) ID identifying the (i-1)-th order progenitor\n# 3) ID identifying the subhalo itself\n# 4) order of the subhalo\n# 5) subhalo mass in units of log10[M/Msun]\n# 6) redshift at first accretion\n# 7) merging timescale in units of [Gyr]\n\n")

    fp.close()

    return


def print_centrals_header(output_folder):

    fp = open(output_folder+"data/output_centrals.txt", "w")

    fp.write("# Catalogue of the central galaxies.\n# Colums:\n# 1) ID number identifying each galaxy\n# 2) stellar mass in units of log10[M/Msun]\n\n")

    fp.close()

    return


def print_satellites_header(output_folder, satellites_redshift):

    fp = open(output_folder+"data/output_satellites.txt", "w")

    fp.write("# Catalogue of the satellite galaxies at redshift z = {:.2f}.\n# Colums:\n# 1) ID number identifying the corresponding central galaxy\n# 2) order of the satellite\n# 3) stellar mass in units of log10[M/Msun]\n\n".format(satellites_redshift))

    fp.close()

    return


def print_log_file(output_folder, parents_len, mergers_len, want_z_at_merge):

    fp = open(output_folder+"output.log", "w")

    fp.write("# Number of rows and columns of parents output file.\n{:d} {:d} {:d}\n".format(0, parents_len, 2))

    if want_z_at_merge == 1:
        fp.write("# Number of rows and columns of mergers output file.\n{:d} {:d} {:d}\n".format(1, mergers_len, 8))
    elif want_z_at_merge == 0:
        fp.write("# Number of rows and columns of mergers output file.\n{:d} {:d} {:d}\n".format(1, mergers_len, 7))

    fp.close()

    return


def print_time(secs):

    print("\nProcess completed correctly")

    if secs < 60.:
        print("Time elapsed: {} sec\n".format(np.int(secs)))

    elif secs >= 60.:
        mins = int(secs/60.)
        secs = np.around(secs-60.*mins)
        print("Time elapsed: {} min {} sec\n".format(mins, np.int(secs)))

    return



def read_data(folder_name, exist_z_at_merge):

    print("\nReading data from file...")

    logfile_name = folder_name+"output.log"

    file_name = folder_name+"data/output_parents.txt"
    len = int(np.loadtxt(logfile_name, usecols=1, max_rows=1))
    halo_type = 0

    num_c.read_data.restype = ndpointer(dtype=ctypes.c_double, shape = np.zeros(len).shape)
    id_halo = num_c.read_data(halo_type,
                              ctypes.create_string_buffer(logfile_name.encode('utf-8')),
                              ctypes.create_string_buffer(file_name.encode('utf-8')),
                              0)
    id_halo = np.int32(id_halo)

    halo_catalog = num_c.read_data(halo_type,
                                   ctypes.create_string_buffer(logfile_name.encode('utf-8')),
                                   ctypes.create_string_buffer(file_name.encode('utf-8')),
                                   1)

    file_name = folder_name+"data/output_mergers.txt"
    len = int(np.loadtxt(logfile_name, usecols=1, skiprows=3, max_rows=1))
    halo_type = 1

    num_c.read_data.restype = ndpointer(dtype=ctypes.c_double, shape = np.zeros(len).shape)
    id_mergers = num_c.read_data(halo_type,
                                 ctypes.create_string_buffer(logfile_name.encode('utf-8')),
                                 ctypes.create_string_buffer(file_name.encode('utf-8')),
                                 0)
    id_mergers = np.int32(id_mergers)

    orders = num_c.read_data(halo_type,
                             ctypes.create_string_buffer(logfile_name.encode('utf-8')),
                             ctypes.create_string_buffer(file_name.encode('utf-8')),
                             1)

    mergers_array = num_c.read_data(halo_type,
                                    ctypes.create_string_buffer(logfile_name.encode('utf-8')),
                                    ctypes.create_string_buffer(file_name.encode('utf-8')),
                                    2)

    z_infall = num_c.read_data(halo_type,
                               ctypes.create_string_buffer(logfile_name.encode('utf-8')),
                               ctypes.create_string_buffer(file_name.encode('utf-8')),
                               3)

    merging_timescale = num_c.read_data(halo_type,
                                        ctypes.create_string_buffer(logfile_name.encode('utf-8')),
                                        ctypes.create_string_buffer(file_name.encode('utf-8')),
                                        4)


    if exist_z_at_merge == "yes":
        z_at_merge = num_c.read_data(halo_type,
                                     ctypes.create_string_buffer(logfile_name.encode('utf-8')),
                                     ctypes.create_string_buffer(file_name.encode('utf-8')),
                                     5)

    elif exist_z_at_merge == "no":
        age_at_merge = Cosmo.age(z_infall) + merging_timescale
        age_today = Cosmo.age(0)

        z_for_interp = np.linspace(0,20,10000)
        age_for_interp = Cosmo.age(z_for_interp)
        mergers.get_redshift_at_merging.restype = ndpointer(dtype=ctypes.c_double, shape = age_at_merge.shape)
        z_at_merge = mergers.get_redshift_at_merging(age_today,
                                            age_at_merge,
                                            z_for_interp,
                                            age_for_interp,
                                            z_for_interp.size,
                                            age_at_merge.size)

    return id_halo, halo_catalog, id_mergers, mergers_array, z_infall, merging_timescale, z_at_merge, orders
