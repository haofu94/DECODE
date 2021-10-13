"""
@ input_parameters.py

Written by Hao Fu

The goal of this module is to read the input parameters from file.
"""

import numpy as np
import getopt
import sys
import os
import ast


from stellar_mass_halo_mass_relation import *


class input_parameters_run:

    """
    Return parameters:
    @ halo_mass_function -> halo mass function (string)
    @ cube_side -> box size [Mpc/h]
    @ z_range -> redshift array
    @ host_mass_params -> minimum, maximum and bin width for halo mass range log[Msun] (tuple)
    @ use_mean_track -> state if use mean track for each halo mass bin (bool)
    @ sub_mass_params -> minimum, maximum and bin width for subhalo mass range log[Msun] (tuple)
    @ output_folder -> output folder path (string)
    @ max_order -> maximum desired subhalo order (int)
    @ use_merger_tree -> use merger tree (bool)
    @ merging_timescale_params -> merging timescale parameters (tuple)
    @ want_z_at_merge -> state if save redshift at merging (string)
    """

    def __init__(self, halo_mass_function, cube_side, z_range, host_mass_params, mass_definition, use_mean_track, path_to_diffmah, mah_type, sub_mass_res, output_folder, max_order, use_merger_tree, merging_timescale_params, want_z_at_merge, want_galaxies, exist_DM, ignore_high_orders, SMHM, include_sat_SF, include_mass_loss, include_quenching, include_stripping, satellites_redshift):

        self.halo_mass_function = halo_mass_function
        self.cube_side = cube_side
        self.z_range = z_range
        self.host_mass_params = host_mass_params
        self.mass_definition = mass_definition
        self.use_mean_track = use_mean_track
        self.path_to_diffmah = path_to_diffmah
        self.mah_type = mah_type
        self.sub_mass_res = sub_mass_res
        self.output_folder = output_folder
        self.max_order = max_order
        self.use_merger_tree = use_merger_tree
        self.merging_timescale_params = merging_timescale_params
        self.want_z_at_merge = want_z_at_merge
        self.want_galaxies = want_galaxies
        self.exist_DM = exist_DM
        self.ignore_high_orders = ignore_high_orders
        self.SMHM = SMHM
        self.include_sat_SF = include_sat_SF
        self.include_mass_loss = include_mass_loss
        self.include_quenching = include_quenching
        self.include_stripping = include_stripping
        self.satellites_redshift = satellites_redshift



class input_parameters_info:

    def __init__(self, folder_name, exist_z_at_merge, want_z_infall_pdf, want_unevolved_shmf, want_evolved_shmf, want_mergers_rate, halo_masses, halo_mass_bin, SMHM, compute_SF, redshifts_for_SFR, stellar_masses_for_SF):

        self.folder_name = folder_name
        self.exist_z_at_merge = exist_z_at_merge
        self.want_z_infall_pdf = want_z_infall_pdf
        self.want_unevolved_shmf = want_unevolved_shmf
        self.want_evolved_shmf = want_evolved_shmf
        self.want_mergers_rate = want_mergers_rate
        self.halo_masses = halo_masses
        self.halo_mass_bin = halo_mass_bin
        self.SMHM = SMHM
        self.compute_SF = compute_SF
        self.redshifts_for_SFR = redshifts_for_SFR
        self.stellar_masses_for_SF = stellar_masses_for_SF



def test_latex():
    use_latex = True

    try:
        import matplotlib.pyplot as plt
        plt.rc('text', usetex=True)
    except:
        ErrMsg = "\n  /!\    Warning: LaTex not installed :(\n / ! \   Not using LaTex for plots text.\n"
        use_latex = False

    return use_latex


def initialize_opts(opts):

    for i in range(len(opts)):
        if opts[i][0] == "-o":
            idx = i
            break

    new_opts = []
    new_opts.append(opts[idx])

    for i in range(len(opts)):
        if i != idx:
            new_opts.append(opts[i])

    return new_opts


def test_output_folder(output_folder):

    try:
        if output_folder[-1] != "/":
            output_folder += "/"

        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)

        if not os.path.isdir(output_folder+"data/"):
            os.mkdir(output_folder+"data/")

    except:
        ErrMessage = "  /!\    Error in the command line:\n / ! \   output folder not given\n"
        sys.exit(ErrMessage)

    return output_folder


def get_save_parameters(input_lines):

    for i in range(len(input_lines)):
        line = input_lines[i].split()
        if len(line)!=0 and line[0][0]!="#":
            if "save_parameters" in line: save_parameters = line[2]
    try:
        if save_parameters == "yes":
            save_parameters = True
        elif save_parameters == "no":
            save_parameters = False
    except:
        save_parameters = True

    return save_parameters



def read_cosmo_model(argv):
    try:
        opts, args = getopt.getopt(argv, "p:o:", ["ifile=", "ofile"])
    except getopt.GetoptError:
        ErrMessage = "Error in reading input file."
        sys.exit(ErrMessage)
    for opt, arg in opts:
        if "-p" in opt:
            with open(arg) as input_file:
                input_lines = input_file.readlines()
            for i in range(len(input_lines)):
                line = input_lines[i].split()
                if len(line)!=0 and line[0][0]!="#":
                    if "cosmological_model" in line: cosmo_str = line[2]
    try:
        cosmo_str
    except:
        cosmo_str = "planck18" #Default
        print("\n  /!\    Attention: cosmological model not give in file.\n / ! \   Will be set to default (Planck 2018).")

    return cosmo_str



def read_parameters_run(argv):

    try:
        opts, args = getopt.getopt(argv, "p:o:", ["ifile=", "ofile"])
    except getopt.GetoptError:
        ErrMessage = "Error in reading input file."
        sys.exit(ErrMessage)

    opts = initialize_opts(opts)

    for opt, arg in opts:

        if "-o" in opt:
            output_folder = arg
            output_folder = test_output_folder(output_folder)

        if "-p" in opt and output_folder != "":
            print("Reading parameters from file:\n{}\n".format(os.getcwd()+"/"+arg))
            with open(arg) as input_file:
                input_lines = input_file.readlines()

            save_parameters = get_save_parameters(input_lines)

            if save_parameters:
                fp = open(output_folder+"log_param.ini", "w")

            for i in range(len(input_lines)):

                line = input_lines[i].split()

                if len(line)!=0 and line[0][0]!="#":

                    if "halo_mass_function" in line: halo_mass_function = line[2]
                    elif "cube_side" in line: cube_side = float(line[2])
                    elif "z_min" in line: z_min = float(line[2])
                    elif "z_max" in line: z_max = float(line[2])
                    #elif "z_bin" in line: z_bin = float(line[2])
                    elif "M0_host_min" in line: M0_host_min = float(line[2])
                    elif "M0_host_max" in line: M0_host_max = float(line[2])
                    elif "M0_host_bin" in line: M0_host_bin = float(line[2])
                    elif "mass_definition" in line: mass_definition = line[2]
                    elif "use_mean_track" in line: use_mean_track = line[2]
                    elif "path_to_diffmah" in line: path_to_diffmah = line[2]
                    elif "mah_type" in line: mah_type = line[2]
                    elif "M_sub_res" in line: M_sub_res = float(line[2])
                    elif "max_order" in line: max_order = int(line[2])
                    elif "use_merger_tree" in line: use_merger_tree = line[2]
                    elif "type_orbital_circularity" in line: type_orbital_circularity = line[2]
                    elif "orbital_circularity" in line: orbital_circularity = float(line[2])
                    elif "fudge" in line: fudge = float(line[2])
                    elif "want_z_at_merge" in line: want_z_at_merge = line[2]
                    elif "want_galaxies" in line: want_galaxies = line[2]
                    elif "exist_DM" in line: exist_DM = line[2]
                    elif "satellites_redshift" in line: satellites_redshift = float(line[2])
                    elif "ignore_high_orders" in line: ignore_high_orders = line[2]
                    elif "SMHM_file" in line: SMHM_file = line[2]
                    elif "SMHM_model" in line: SMHM_model = line[2]
                    elif "scatter" in line: scatter = float(line[2])
                    elif "include_sat_SF" in line: include_sat_SF = line[2]
                    elif "include_mass_loss" in line: include_mass_loss = line[2]
                    elif "include_quenching" in line: include_quenching = line[2]
                    elif "include_stripping" in line: include_stripping = line[2]

                    if save_parameters:
                        fp.write("{0} = {1}\n".format(line[0], line[2]))

            if save_parameters:
                fp.close()

    try:
        halo_mass_function
    except:
        halo_mass_function = "tinker08"

    try:
        cube_side
    except:
        ErrMessage = "  /!\    Error in the input parameter file:\n / ! \   cube side length not given\n"
        sys.exit(ErrMessage)

    try:
        z_bin = 0.1 #default
        z_range = np.arange(z_min, z_max, z_bin)
    except:
        ErrMessage = "  /!\    Error in the input parameter file:\n / ! \   redshift parameters missing\n"
        sys.exit(ErrMessage)

    try:
        host_mass_params = (M0_host_min, M0_host_max, M0_host_bin)
    except:
        ErrMessage = "  /!\    Error in the input parameter file:\n / ! \   mass parameters missing\n"
        sys.exit(ErrMessage)

    try:
        M_sub_res
    except:
        M_sub_res = 1.e-3
        print("  /!\    Subhalo mass resolution not given.\n / ! \   Set to default = 1e-3.\n")

    try:
        mass_definition
    except:
        mass_definition = "vir"

    try:
        if use_mean_track == "yes":
            use_mean_track = True
        elif use_mean_track == "no":
            use_mean_track = False
    except:
        use_mean_track = True

    if not use_mean_track:
        try:
            path_to_diffmah
        except:
            ErrMessage = "  /!\    Error in the input parameter file:\n / ! \   path to diffmah missing\n"
            sys.exit(ErrMessage)
        if path_to_diffmah[-1] != "/":
            path_to_diffmah += "/"
    else:
        path_to_diffmah = ""

    if not use_mean_track:
        try:
            mah_type
            if mah_type == "None":
                mah_type = None
        except:
            mah_type = None

    try:
        max_order
    except:
        max_order = 3

    try:
        if use_merger_tree == "yes":
            use_merger_tree = True
        elif use_merger_tree == "no":
            use_merger_tree = False
    except:
        use_merger_tree = False

    try:
        type_orbital_circularity
    except:
        type_orbital_circularity = "constant"
        orbital_circularity = 0.5
    if type_orbital_circularity == "constant":
        type_orbital_circularity = 0
    elif type_orbital_circularity == "gaussian":
        type_orbital_circularity = 1
        orbital_circularity = 0.
    else:
        ErrMessage = "  /!\    Error in the input parameter file:\n / ! \   non existing orbital circularity type\n"
        sys.exit(ErrMessage)

    try:
        fudge
    except:
        fudge = 0.5

    merging_timescale_params = (type_orbital_circularity, orbital_circularity, fudge)

    try:
        want_z_at_merge
    except:
        want_z_at_merge = "no"

    try:
        if want_galaxies == "yes":
            want_galaxies = True
        elif want_galaxies == "no":
            want_galaxies = False
    except:
        want_galaxies = False

    try:
        if exist_DM == "yes":
            exist_DM = True
        elif exist_DM == "no":
            exist_DM = False
    except:
        exist_DM = True

    try:
        if ignore_high_orders == "yes":
            ignore_high_orders = 1
        elif ignore_high_orders == "no":
            ignore_high_orders = 0
    except:
        ignore_high_orders = 1

    try:
        SMHM_file
    except:
        ErrMessage = "  /!\    Error in the input parameter file:\n / ! \   SMHM file name missing\n"
        sys.exit(ErrMessage)

    try:
        scatter
        constant_scatter = True
    except:
        scatter = 0.
        constant_scatter = True

    try:
        if include_sat_SF == "yes":
            include_sat_SF = 1
        elif include_sat_SF == "no":
            include_sat_SF = 0
    except:
        include_sat_SF = 0

    try:
        if include_mass_loss == "yes":
            include_mass_loss = 1
        elif include_mass_loss == "no":
            include_mass_loss = 0
    except:
        include_mass_loss = 0

    try:
        if include_quenching == "yes":
            include_quenching = 1
        elif include_quenching == "no":
            include_quenching = 0
    except:
        include_quenching = 0

    try:
        if include_stripping == "yes":
            include_stripping = 1
        elif include_stripping == "no":
            include_stripping = 0
    except:
        include_stripping = 0

    try:
        satellites_redshift
    except:
        ErrMessage = "  /!\    Error in the input parameter file:\n / ! \   satellites redshift missing\n"
        sys.exit(ErrMessage)

    SMHM = get_SMHM_numerical(SMHM_file, constant_scatter, scatter, z = satellites_redshift)

    input_params_run = input_parameters_run(halo_mass_function, cube_side, z_range, host_mass_params, mass_definition, use_mean_track, path_to_diffmah, mah_type, M_sub_res, output_folder, max_order, use_merger_tree, merging_timescale_params, want_z_at_merge, want_galaxies, exist_DM, ignore_high_orders, SMHM, include_sat_SF, include_mass_loss, include_quenching, include_stripping, satellites_redshift)

    return input_params_run




#---------------------------------------------------------------

def read_parameters_info(argv):

    """
    Return parameters for info
    """

    try:
        opts, args = getopt.getopt(argv, "p:o:", ["ifile=", "ofile"])
    except getopt.GetoptError:
        ErrMessage = "Error in reading input file."
        sys.exit(ErrMessage)

    for opt, arg in opts:

        if "-p" in opt:
            print("Reading parameters from file:\n{}\n".format(os.getcwd()+"/"+arg))
            with open(arg) as input_file:
                input_lines = input_file.readlines()
            for i in range(len(input_lines)):
                line = input_lines[i].split()
                if len(line)!=0 and line[0][0]!="#":
                    if "exist_z_at_merge" in line: exist_z_at_merge = line[2]
                    elif "want_z_infall_pdf" in line: want_z_infall_pdf = line[2]
                    elif "want_unevolved_shmf" in line: want_unevolved_shmf = line[2]
                    elif "want_evolved_shmf" in line: want_evolved_shmf = line[2]
                    elif "want_mergers_rate" in line: want_mergers_rate = line[2]
                    elif "halo_masses" in line:
                        Line = line[2]
                        for i in range(3, len(line)):
                            Line = Line + line[i]
                        halo_masses = ast.literal_eval(Line)
                    elif "halo_mass_bin" in line: halo_mass_bin = float(line[2])
                    elif "SMHM_file" in line: SMHM_file = line[2]
                    elif "SMHM_model" in line: SMHM_model = line[2]
                    elif "scatter" in line: scatter = float(line[2])
                    elif "compute_SF" in line: compute_SF = line[2]
                    elif "redshifts_for_SFR" in line:
                        Line = line[2]
                        for i in range(3, len(line)):
                            Line = Line + line[i]
                        redshifts_for_SFR = ast.literal_eval(Line)
                    elif "stellar_masses_for_SF" in line:
                        Line = line[2]
                        for i in range(3, len(line)):
                            Line = Line + line[i]
                        stellar_masses_for_SF = ast.literal_eval(Line)

    # Check if parameters are given and set defaults

    try:
        exist_z_at_merge
    except:
        exist_z_at_merge = "no"

    try:
        if want_z_infall_pdf == "yes":
            want_z_infall_pdf = True
        elif want_z_infall_pdf == "no":
            want_z_infall_pdf = False
    except:
        want_z_infall_pdf = True

    try:
        if want_unevolved_shmf == "yes":
            want_unevolved_shmf = True
        elif want_unevolved_shmf == "no":
            want_unevolved_shmf = False
    except:
        want_unevolved_shmf = False

    try:
        if want_evolved_shmf == "yes":
            want_evolved_shmf = True
        elif want_evolved_shmf == "no":
            want_evolved_shmf = False
    except:
        want_evolved_shmf = False

    try:
        if want_mergers_rate == "yes":
            want_mergers_rate = True
        elif want_mergers_rate == "no":
            want_mergers_rate = False
    except:
        want_mergers_rate = False

    try:
        halo_masses
        halo_mass_bin
    except:
        ErrMessage = "  /!\    Error in the input parameter file:\n / ! \   mass parameters missing\n"
        sys.exit(ErrMessage)

    try:
        SMHM_file
    except:
        ErrMessage = "  /!\    Error in the input parameter file:\n / ! \   SMHM file name missing\n"
        sys.exit(ErrMessage)

    try:
        scatter
        constant_scatter = True
    except:
        scatter = 0.

    try:
        if compute_SF == "yes":
            compute_SF = True
        elif compute_SF == "no":
            compute_SF = False
    except:
        compute_SF = True

    try:
        redshifts_for_SFR
    except:
        redshifts_for_SFR = [0, 1, 2, 3]

    try:
        stellar_masses_for_SF
    except:
        stellar_masses_for_SF = [11.]

    folder_name = argv[-1]
    if folder_name[-1] != "/":
        folder_name += "/"
    print("Analysing data from folder: {}".format(folder_name+"data/"))

    SMHM = get_SMHM_numerical(SMHM_file, constant_scatter, scatter, z = 0.1) # TO DO: EXTRAPOLATE OUT TO EVERY z

    input_params_info = input_parameters_info(folder_name, exist_z_at_merge, want_z_infall_pdf, want_unevolved_shmf, want_evolved_shmf, want_mergers_rate, halo_masses, halo_mass_bin, SMHM, compute_SF, redshifts_for_SFR, stellar_masses_for_SF)

    return input_params_info
