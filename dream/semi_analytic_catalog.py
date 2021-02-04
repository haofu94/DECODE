"""
@ file semi_analytic_catalog.py

Written by Chris Marsden, Hao Fu
"""

import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
from colossus.lss import mass_function


def generate_parents_catalogue(input_params_run, h):

    """
    Function to generate the semi analytic halo catalogue (without coordinates) for galaxy testing

    :param catalogue_volume: float, cosmological volume within which to generate the catalog. [Mpc/h]^3
    :param mass_params: tuple, (mass low, mass high, spacing). log[Msun]
    :param z: float, redshift.
    :param h: float, reduced hubble constant.
    :return array, of halo masses. log[Msun]
    """

    print("Generating catalogue for a volume of ({:.2f} Mpc/h)^3\n".format(input_params_run.cube_side))

    catalogue_volume = input_params_run.cube_side**3

    # Get the bin width and generate the bins.
    bin_width = input_params_run.host_mass_params[2]
    mass_range = 10 ** np.arange(input_params_run.host_mass_params[0], input_params_run.host_mass_params[1] + bin_width, bin_width) #log[M/Msun]

    # Generate the mass function itself - this is from the colossus toolbox
    local_mass_function = mass_function.massFunction(mass_range*h, input_params_run.z_range[0], mdef = input_params_run.mass_definition, model = input_params_run.halo_mass_function, q_out='dndlnM') * np.log(10) # dn/dlog10M [dex^-1 (Mpc/h)^-3]

    # Calculate cumulative of the halo mass function
    cumulative_mass_function = np.cumsum(local_mass_function * bin_width * catalogue_volume)

    # Get the maximum cumulative number.
    max_number = np.floor(np.max(cumulative_mass_function))
    if (np.random.uniform(0,1) > np.max(cumulative_mass_function)-max_number):
        max_number += 1

    interpolator = interp1d(cumulative_mass_function, mass_range)
    range_numbers = np.random.uniform(np.min(cumulative_mass_function), np.max(cumulative_mass_function), int(max_number))
    mass_catalog = interpolator(range_numbers)

    print("Number of halos generated: {:d}\n".format(len(mass_catalog)))

    mass_catalog = np.log10(mass_catalog) #log[M/Msun]

    # MAGHEGGIO
    #mass_catalog = np.random.uniform(11.5, 16., 1500)
    """mass_catalog = np.random.uniform(11.5, 12., 500)
    mass_catalog = np.concatenate((mass_catalog, np.random.uniform(12., 13., 500)))
    mass_catalog = np.concatenate((mass_catalog, np.random.uniform(13., 14., 500)))
    mass_catalog = np.concatenate((mass_catalog, np.random.uniform(14., 15., 500)))
    mass_catalog = np.concatenate((mass_catalog, np.random.uniform(15., 16., 300)))"""

    return mass_catalog
