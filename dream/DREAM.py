"""
DREAM (DiscRete statistical sEmi-empiricAl Model)

Authors:
Hao Fu <h.fu@soton.ac.uk>
Chris Marsden <c.marsden@soton.ac.uk>
Francesco Shankar <f.shankar@soton.ac.uk>
Max Dickson <md2g17@soton.ac.uk>
"""


# Generic Imports
import sys
from time import time

from colossus.cosmology import cosmology
from input_parameters import *
cosmo_str = read_cosmo_model(sys.argv[2:])
cosmology.setCosmology(cosmo_str)
Cosmo = cosmology.getCurrent()

# Local
import semi_analytic_catalog
from DREAM_track import *
from DREAM_galaxies import *
import info
import other_functions as others



if __name__ == "__main__":

    time_ini = time()

    command = others.initialize(cosmo_str, sys.argv)

    if command == "run":

        input_params_run = read_parameters_run(sys.argv[2:])

        if (input_params_run.want_galaxies and not input_params_run.exist_DM) or not input_params_run.want_galaxies:

            # Generate the correct distribution of host central halos at z=0
            halo_catalog = semi_analytic_catalog.generate_parents_catalogue(input_params_run, Cosmo.h)

            # Generate the mergers
            DREAM_track(halo_catalog, input_params_run)

        if input_params_run.want_galaxies:

            DREAM_galaxies(input_params_run)

    elif command == "info":

        input_params_info = read_parameters_info(sys.argv[2:])

        info.analyse_mergers(input_params_info)

    else:

        others.print_command_error()


    time_fin = time(); others.print_time(time_fin-time_ini)
