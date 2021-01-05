#ifndef error_messages
# define error_messages


#define _alloc_error_message_ "Error: array dynamical allocation failed!"
#define _dealloc_error_message_ "Error: array deallocation failed!"

#define _load_data_error_message_ "\nFailed to load data from file."

#define _mean_error_message_ "Failed to compute mean."

#define _mass_accretion_error_message_ "\nError in function get_subhaloes\nThe mass at z+dz is higher than that at z.\n"

#define _mass_function_error_message_ "\nError in calculating mass function.\n"

#define _haloes_number_error_message_ "\nError in getting total haloes number.\n"
#define _mass_from_pdf_error_message_ "\nError in getting mass from mass function.\n"
#define _assign_order_error_message_ "\nError in assigning subhalo order.\n"
#define _dynamical_timescale_error_message_ "\nError in computing dynamical timescale.\n"
#define _merging_timescale_error_message_ "\nError in calculating merging timescale.\n"
#define _quenching_timescale_error_message_ "\nError in calculating quenching timescale.\n"
#define _stripping_timescale_error_message_ "\nError in calculating stripping timescale.\n"
#define _infall_redshift_error_message_ "\nError in calculating infall redshift.\n"
#define _infall_redshift_from_pdf_error_message_ "\nError in calculating infall redshift from PDF.\n"
#define _redshift_pdf_error_message_ "\nError in computing infall redshift PDF.\n"
#define _delta_mass_func_error_message_ "\nError in calculating delta mass function.\n"
#define _redshift_merging_error_message_ "\nError in computing redshift at merging.\n"

#define _cumsum_error_message_ "\nError in calculating cumulative sum.\n"
#define _linear_space_error_message "\nFailed to allocate linear space.\n"
#define _arange_error_message "\nFailed to allocate arange.\n"

#define _get_SMHM_params_error_message_ "\nFailed to get SMHM parameters from model name.\n"
#define _SMHM_error_message_ "\nFailed to assign stellar mass via the SMHM relation.\n"
#define _SMHM_read_matrix_error_message_ "\nFailed to read SMHM matrix from file.\n"
#define _get_subhaloes_error_message_ "\nFailed to get subhaloes!\n"

#define _get_mask_error_message_ "\nFailed to get mask.\n"
#define _append_error_message_ "\nFailed to append object to array.\n"

#define _mstar_at_z_error_message_ "\nFailed to compute evolved stellar mass."
#define _mass_acc_error_message_ "\nFailed to computed halo accretion track.\n"

#define _Mhalo_track_to_Mstar_track_error_message_ "\nFailed to convert Mhalo track to Mstar track.\n"
#define _merged_stellar_mass_error_message_ "\nFailed to compute stellar mass from mergers.\n"
#define _compute_SFR_error_message_ "\nFailed to compute SFR.\n"


#endif
