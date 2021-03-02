/** @ file mergers.c
  *
  * Written by Hao Fu
  *
  * The main goal of this module is to generate a dark matter subhalo catalogue.
  **/


#include "../include/dream.h"

#include "centrals.c"


/**
  * Compute the total unevolved SHMF, Jiang & van den Bosch 2016 Appendix A, Eq. (A1)
  *
  * Input:
  * @ Params -> list of SHMF parameters
  * @ M_host -> mass of the parent halo [log10(M/Msun)]
  * @ M_subhalos_range -> mass interval of subhalos [log10(M/Msun)]
  * @ len -> length of dn_dlogX_arr array
  * @ dn_dlogX_arr -> pointer to the subhalo mass function array
  **/
int vdB_USHMF(double *Params,
              double M_host,
              double *M_subhalos_range,
              int len,
              double **dn_dlnX_arr){

  int i;
  double psi;
  double gamma = Params[0];
  double alpha = Params[1];
  double beta = Params[2];
  double omega = Params[3];
  double a = Params[4];
  double Part1, Part2;

  for (i=0; i<len; i++){
    psi = pow(10., *(M_subhalos_range+i)) / pow(10., M_host);
    Part1 = gamma*pow(a*psi, alpha);
    Part2 = exp(-beta*pow(a*psi, omega));
    *(*dn_dlnX_arr+i) = Part1*Part2;
    *(*dn_dlnX_arr+i) = *(*dn_dlnX_arr+i) * log(10.);
  }

  return _success_;
}



/**
  * Compute the i-th order unevolved SHMF, Jiang & van den Bosch 2016 Eq. (14)
  *
  * Input:
  * @ Params -> list of SHMF parameters
  * @ M_host -> mass of the parent halo [log10(M/Msun)]
  * @ M_subhalos_range -> mass interval of subhalos [log10(M/Msun)]
  * @ len -> length of dn_dlogX_arr array
  * @ dn_dlogX_arr -> pointer to the subhalo mass function array
  **/
int vdB_USHMF_1st_order(double *Params,
                        double M_host,
                        double *M_subhalos_range,
                        int len,
                        double **dn_dlnX_arr){

  int i;
  double psi;
  double gamma1 = Params[0];
  double alpha1 = Params[1];
  double gamma2 = Params[2];
  double alpha2 = Params[3];
  double beta = Params[4];
  double omega = Params[5];
  double Part1, Part2;
  for (i=0; i<len; i++){
    psi = pow(10., *(M_subhalos_range+i)) / pow(10., M_host);
    Part1 = gamma1*pow(psi, alpha1);
    Part2 = gamma2*pow(psi, alpha2);
    *(*dn_dlnX_arr+i) = (Part1 + Part2) * exp(-beta * pow(psi, omega));
    *(*dn_dlnX_arr+i) = *(*dn_dlnX_arr+i) * log(10.);
  }

  return _success_;
}



/**
  * Compute the delta-SHMF between z and z+dz [N]
  *
  * Input:
  * @ length -> length of the mass function array
  * @ mass_function_at_z -> mass function at redshift z
  * @ mass_function_at_z_plus_dz -> mass function at redshift z+dz
  * @ mass_bin -> mass bin width [log(M/Msun)]
  * @ delta_mass_function -> pointer to the delta_mass_function array
  **/
int compute_delta_mass_function(int length,
                                double **mass_function_at_z,
                                double **mass_function_at_z_plus_dz,
                                double mass_bin,
                                double **delta_mass_function){

  int i;

  for (i=0; i<length; i++){

    *(*delta_mass_function+i) = ( (*(*mass_function_at_z+i)) - \
                          (*(*mass_function_at_z_plus_dz+i)) ) * mass_bin; //[N]

  }

  return _success_;

}



/**
  * Check that M(z) >= M(z+dz)
  * mass_at_z -> mass at redshift z [log(M/Msun)]
  * mass_at_z_plus_dz -> mass at redshift z+dz [log(M/Msun)]
  **/
int check_accretion(double mass_at_z, double mass_at_z_plus_dz){

  if (mass_at_z < mass_at_z_plus_dz) {

    return _failure_;

  } else {

    return _success_;

  }
}


/**
  * Compute the probability distribution of the infall redshifts.
  * Infall redshift stays for the redshift at which the j-th order subhalo
  * first falls into the (j-1)-th order "parent" subhalo.
  * The PDF depends on the subhalo order and mass.
  * Parameters fitted by running merger trees via DREAM [H. Fu].
  *
  * Input:
  * @ redshift_pdf_range -> redshift range array for PDF
  * @ length -> length of redshift_pdf_range and redshift_pdf arrays
  * @ subhalo_mass -> mass of the subhalo [log(M/Msun)]
  * @ subhalo_order -> order of the subhalo
  * @ redshift_pdf -> pointer to the redshift PDF array
  **/
int compute_analytical_redshift_PDF(double *redshift_pdf_range,
                                    int length,
                                    double subhalo_mass,
                                    int subhalo_order,
                                    double **redshift_pdf) {

  int i;
  double A, beta, delta, gamma, alpha;

  if ((subhalo_order == 2) && (subhalo_mass < 11.)) {

    A = 17.59386; beta = 1.25453; delta = 13.06571; gamma = 2.27080; alpha = 1.39963;

  } else if ((subhalo_order == 2) && (11. <= subhalo_mass) && (subhalo_mass < 12.)) {

    A = 23.86118; beta = 1.64827; delta = 9.47616; gamma = 0.69930; alpha = 1.39298;

  } else if ((subhalo_order == 2) && (12. <= subhalo_mass)) {

    A = 22.32932; beta = 2.40317; delta = 2.97333; gamma = -4.17881; alpha = 1.19073;

  } else if ((subhalo_order >= 3) && (subhalo_mass < 11.)) {

    A = 23.26552; beta = 1.66035; delta = 11.89017; gamma = 1.24507; alpha = 2.32610;

  } else if ((subhalo_order >= 3) && (11. <= subhalo_mass)) {

    A = 24.47686; beta = 2.20970; delta = 4.40460; gamma = -2.77968; alpha = 2.06852;

  }

  /*else if (subhalo_order == 1) {
    A = 9.206457497136162; beta = 0.8316294547478578; delta = 13.180614924713375; gamma = 1.2384602297980127; alpha = 0.9695950516941695;
  }*/

   else {
    printf("\nError in subhalo order or mass!\n\n");
    exit(0);
  }

  for (i=0; i<length; i++) {

    *(*redshift_pdf+i) = A * pow(*(redshift_pdf_range+i), alpha) / \
                        (delta * exp(beta * (*(redshift_pdf_range+i))) - gamma);

  }

  return _success_;

}


/**
  * Generate infall redshift according to the PDF,
  * given the subhalo mass and order.
  *
  * Input:
  * @ subhalo_mass -> mass of the subhalo [log(M/Msun)]
  * @ subhalo_order -> order of the subhalo
  * @ redshift_lower_limit -> lower limit of the infall redshift range
  * @ redshift_upper_limit -> upper limit of the infall redshift range
  * @ infall_redshift -> pointer to the infall redshift variable
  **/
int get_infall_redshift_from_pdf(double subhalo_mass,
                                 int subhalo_order,
                                 double redshift_lower_limit,
                                 double redshift_upper_limit,
                                 double *infall_redshift) {

  double low, high;

  int length,i;

  double redshift_bin = 0.1;

  double *redshift_pdf_range, *redshift_pdf, *cumulative;

  length = (int) ((redshift_upper_limit - redshift_lower_limit) / redshift_bin);

  dream_call(double_malloc(length,
                            &redshift_pdf_range),
              _alloc_error_message_);

  dream_call(arange(redshift_lower_limit,
                     redshift_upper_limit,
                     redshift_bin,
                     length,
                     &redshift_pdf_range),
              _arange_error_message);

  dream_call(double_malloc(length,
                            &redshift_pdf),
              _alloc_error_message_);

  dream_call(compute_analytical_redshift_PDF(redshift_pdf_range,
                                              length,
                                              subhalo_mass,
                                              subhalo_order,
                                              &redshift_pdf),
              _redshift_pdf_error_message_);

  dream_call(double_calloc(length,
                            &cumulative),
              _alloc_error_message_);

  dream_call(cumsum(redshift_pdf,
                     length,
                     &cumulative),
              _cumsum_error_message_);

  low = min(cumulative, length);
  high = max(cumulative, length);

  *infall_redshift = linear_interp(get_random_uniform(low, high),
                                   cumulative, redshift_pdf_range, length);

  dream_call(dealloc(redshift_pdf_range),
              _dealloc_error_message_);

  dream_call(dealloc(redshift_pdf),
              _dealloc_error_message_);

  dream_call(dealloc(cumulative),
              _dealloc_error_message_);

  return _success_;
}



/**
  * Generate infall redshift for subhaloes of order higher than 1.
  *
  * Input:
  * subhalo_mass -> subhalo mass [log(M/Msun)]
  * subhalo_order -> subhalo order
  * redshift -> redshift
  * redshift_bin -> redshift bin
  * redshift_max -> maximum redshift
  * redshift_infall -> pointer to the infall redshift variable
  **/
int generate_infall_redshift_higher_order(double subhalo_mass,
                                          int subhalo_order,
                                          double redshift,
                                          double redshift_bin,
                                          double redshift_max,
                                          double *redshift_infall){

  dream_call(get_infall_redshift_from_pdf(subhalo_mass,
                                           subhalo_order,
                                           0., 20.,
                                           &(*redshift_infall)),
              _infall_redshift_from_pdf_error_message_);

  if (*redshift_infall < get_random_uniform(redshift, redshift+redshift_bin)) {

    *redshift_infall = -1.;
  }

  return _success_;
}



/**
  * Compute redshift when full merging happens.
  *
  * Input:
  * @ universe_age_today -> age of the Universe at z=0 [Gyr]
  * @ universe_age_at_merging -> age of the Universe at merging [Gyr]
  * @ redshift_for_interpolation -> redshift array of the age-redshift relation
  * @ universe_age_for_interpolation -> age array of the age-redshift relation [Gyr]
  * @ interp_length -> length of redshift_for_interpolation
                       and universe_age_for_interpolation arrays
  * @ length -> length of universe_age_at_merging array
  * Return:
  * £ redshift_at_merging -> redshift at merging
  **/
double *get_redshift_at_merging(double universe_age_today,
                                double *universe_age_at_merging,
                                double *redshift_for_interpolation,
                                double *universe_age_for_interpolation,
                                int interp_length,
                                int length) {

  // Note to myself:
  // TO BE DELETED, modify in other_functions.py

  int i;

  double *redshift_at_merging;

  dream_call(double_malloc(length,
                            &redshift_at_merging),
              _alloc_error_message_);

  for (i=0; i<length; i++) {
    if (*(universe_age_at_merging+i) > universe_age_today) {
      *(redshift_at_merging+i) = -1.;
    } else if (*(universe_age_at_merging+i) <= universe_age_today) {
      *(redshift_at_merging+i) = linear_interp(*(universe_age_at_merging+i), \
        universe_age_for_interpolation, redshift_for_interpolation, interp_length);
    }
  }

  return redshift_at_merging;

}



/**
  * Compute redshift at which full merging happens.
  *
  * Input:
  * @ redshift_infall -> infall redshift array
  * @ merging_timescale -> merging timescale array [Gyr]
  * @ length -> length of the previous arrays
  * @ cosmo_time -> cosmological time
  * @ z_at_merge -> pointer to the redshift at merging array
  **/
int compute_redshift_at_merging(double *redshift_infall,
                                double *merging_timescale,
                                int length,
                                cosmological_time *cosmo_time,
                                double **z_at_merge){

  int i;
  double age_at_merge, age_today;

  age_today = linear_interp(0., cosmo_time->redshift, cosmo_time->age, cosmo_time->length);

  for (i=0; i<length; i++){

    age_at_merge = linear_interp(*(redshift_infall+i), cosmo_time->redshift, cosmo_time->age, cosmo_time->length) + (*(merging_timescale+i));

    if (age_at_merge > age_today) {

      *(*z_at_merge+i) = -1.;

    } else if (age_at_merge <= age_today) {

      *(*z_at_merge+i) = linear_interp(age_at_merge, cosmo_time->age, cosmo_time->redshift, cosmo_time->length);

    }
  }

  return _success_;
}



/**
  * Assign order to subhalo.
  *
  * Input:
  * @ max_order -> maximum order required for subhaloes
  * @ subhalo_mass -> mass of the subhalo [log(M/Msun)]
  * @ halo_mass_at_z0 -> mass of the parent halo at z=0 [log(M/Msun)]
  * @ SHMF -> subhalo mass functions
  * @ order -> pointer to the order variable
  **/
int assign_subhalo_order(int max_order,
                         double subhalo_mass,
                         double halo_mass_at_z0,
                         subhalo_mass_functions *SHMF,
                         int *order){

  int i, index_psi;
  double this_psi, decide_order;
  double ith_order_cumulative_probability = 0.;

  *order = 0;

  decide_order = get_random_uniform(0., 1.);

  //this_psi = pow(10., subhalo_mass) / pow(10., halo_mass_at_z);
  this_psi = pow(10., subhalo_mass) / pow(10., halo_mass_at_z0);

  index_psi = find_index(SHMF->psi, this_psi, SHMF->length_psi);

  if ( (*(SHMF->first+index_psi)+ (*(SHMF->second+index_psi)) + (*(SHMF->third+index_psi))) > (*(SHMF->total+index_psi)) ) {

    *(SHMF->total+index_psi) = *(SHMF->first+index_psi)+ (*(SHMF->second+index_psi)) + (*(SHMF->third+index_psi));

  }

  for (i=0; i<max_order; i++) {

    if (i==0) {

      ith_order_cumulative_probability += *(SHMF->first+index_psi) / (*(SHMF->total+index_psi));

    } else if (i==1) {

      ith_order_cumulative_probability += *(SHMF->second+index_psi) / (*(SHMF->total+index_psi));

    } else if (i==2) {

      ith_order_cumulative_probability += *(SHMF->third+index_psi) / (*(SHMF->total+index_psi));

    } else if (i==3) {

      ith_order_cumulative_probability += *(SHMF->fourth+index_psi) / (*(SHMF->total+index_psi));

    } else if (i==4) {

      ith_order_cumulative_probability += *(SHMF->fifth+index_psi) / (*(SHMF->total+index_psi));

    }

    if (decide_order <= ith_order_cumulative_probability) {

      *order = i+1;

      break;
    }
  }

  return _success_;
}



/**
  * Bryan & Norman 1998, ApJ, 495, 80
  **/
double Delta_vir(double z,
                 double Om){

  double x;

  x = Om - 1.;

  return 18. * pow(pi, 2.) + 82. * x - 39. * pow(x, 2.); // for curvature = 0
  //return 18. * pow(pi, 2.) + 60. * x - 32. * pow(x, 2.); // for dark energy = 0
}



/**
  * Compute dynamical friction timescale. [Gyr]
  * Jiang & van den Bosch 2016, MNRAS 458, 2870–2884
  *
  * Input:
  * @ redshift -> redshift
  * @ cosmo_params -> cosmological parameters
  * @ tau_dyn -> pointer to the dynamical timescale pointer
  **/
int compute_dynamical_timescale(double redshift,
                                cosmological_parameters *cosmo_params,
                                double *tau_dyn){

  double Om, Hz;

  Om = linear_interp(redshift, cosmo_params->z, cosmo_params->Om, cosmo_params->length);
  Hz = linear_interp(redshift, cosmo_params->z, cosmo_params->Hz, cosmo_params->length);

  *tau_dyn = 1.628 / cosmo_params->h / sqrt(Delta_vir(redshift, Om)/178.) * (cosmo_params->H0 / Hz);

  return _success_;
}



/**
  * Jiang & van den Bosch 2016 Eq. (3)
  **/
double mass_at_t(double m,
                 double M,
                 double z,
                 double z_inf,
                 double age_at_z,
                 double age_at_z_inf,
                 cosmological_parameters *cosmo_params){

  m = pow(10., m);
  M = pow(10., M);
  double xi = 0.07;
  double A = 1.54;
  //double xi = 0.36;
  //double A = 23.7;
  double dt = age_at_z - age_at_z_inf;
  double tau;

  dream_call(compute_dynamical_timescale(z_inf,
                                          cosmo_params,
                                          &tau),
              _dynamical_timescale_error_message_);

  tau /= A;

  return log10( m* pow(1.+ xi*pow(m/M, xi) *(dt/tau), -1./xi) );
  //return log10( m * exp(-dt/tau));
}
double get_mass_at_t(double m,
                 double M,
                 double z,
                 double z_inf,
                 cosmological_time *cosmo_time,
                 cosmological_parameters *cosmo_params){

  m = pow(10., m);
  M = pow(10., M);
  double xi = 0.07;
  double A = 1.54;
  //double xi = 0.36;
  //double A = 23.7;
  double age_at_z = linear_interp(z, cosmo_time->redshift, \
              cosmo_time->age, cosmo_time->length);
  double age_at_z_inf = linear_interp(z_inf, cosmo_time->redshift, \
              cosmo_time->age, cosmo_time->length);
  double dt = age_at_z - age_at_z_inf;
  double tau;

  dream_call(compute_dynamical_timescale(z_inf,
                                          cosmo_params,
                                          &tau),
              _dynamical_timescale_error_message_);

  tau /= A;

  return log10( m* pow(1.+ xi*pow(m/M, xi) *(dt/tau), -1./xi) );
  //return log10( m * exp(-dt/tau));
}



/**
  * Compute merging time-scale [Gyr]
  * Boylan-Kolchin et al. 2008
  * McCavana et al. 2012
  * Grylls et al. 2018 Sec. 3.1.4 (fudge factor)
  *
  * Input:
  * @ central_mass -> central halo mass [log(M/Msun)]
  * @ sub_mass -> subhalo mass [log(M/Msun)]
  * @ z -> redshift
  * @ mergers_params -> mergers parameters
  * @ cosmo_params -> cosmological parameters
  * @ merging_timescale -> pointer to merging timescale variable
  **/
int compute_merging_timescale(double parent_mass,
                              double subhalo_mass,
                              double z,
                              mergers_parameters *mergers_params,
                              cosmological_parameters *cosmo_params,
                              double *merging_timescale){

  // Note that this needs to be in h^-1
  double A = 0.9;
  double B = 1.0;
  double C = 0.6;
  double D = 0.1;

  //Boylan-Kolchin et al. 2008
  //A = 0.216; B = 1.3; C = 1.9; D = 1.0;

  double mass_ratio, fudge, dynamical_friction, orbital_circularity, orbital_energy;

  mass_ratio = pow(10., parent_mass - subhalo_mass);

  dream_call(compute_dynamical_timescale(z,
                                          cosmo_params,
                                          &dynamical_friction),
              _dynamical_timescale_error_message_); //[Gyr]

  if (mergers_params->fudge<0.){

    // linear relation on mass ratio
    /*double fudge1 = 0.65; //alti psi
    double fudge2 = 1.; //bassi psi
    fudge = (mass_ratio - 1.) / (1. - 1000.) * (fudge1 - fudge2) + fudge1;*/
    fudge = 0.00035035 * (mass_ratio - 1.) + 0.65;

  } else if (mergers_params->fudge>=0.) {

    fudge = mergers_params->fudge;

  }

  dynamical_friction *= fudge;

  if (mergers_params->type_orbital_circularity == 0) {

    orbital_circularity = mergers_params->orbital_circularity;

  } else if (mergers_params->type_orbital_circularity == 1) {

    orbital_circularity = get_random_gaussian(0.5, 0.23, 0., 2.); //Shankar et al 2014

  }

  orbital_energy = pow(orbital_circularity, 2.17) / (1. - sqrt( 1. - pow(orbital_circularity, 2.) ) ); // energy

  *merging_timescale =  dynamical_friction * A * pow(mass_ratio, B) / \
    log( 1. + mass_ratio ) * exp(C * orbital_circularity) * pow(orbital_energy, D);

  return _success_;
}



/**
  * Compute the toal number of haloes given the mass function [N]
  *
  * Input:
  * mass_function -> mass function
  * length -> length of the mass function array
  * count -> pointer to the count variable
  **/
int compute_total_haloes_number(double *mass_function,
                                int length,
                                int *count){

  double * cumulative_mass_function;

  dream_call(double_calloc(length,
                            &cumulative_mass_function),
              _alloc_error_message_);

  dream_call(cumsum(mass_function,
                     length,
                     &cumulative_mass_function),
              _cumsum_error_message_); //Cumulative delta-SHMF

  *count = (int)max(cumulative_mass_function, length);

  //Fractional part
  //Random number between the highest integer and itself +1
  if (get_random_uniform((double)(*count), (double)(*count+1)) < max(cumulative_mass_function, length)) {

    *count += 1;
  }

  dream_call(dealloc(cumulative_mass_function),
              _dealloc_error_message_);

  return _success_;
}



/**
  * Generate a halo mass given the mass function.
  *
  * Input:
  * cumulative_mass_function -> cumulative mass function
  * mass_range -> mass range array
  * length -> length of the previous arrays
  * halo_mass -> pointer to the halo mass variable
  **/
int generate_halo_from_mass_function(double *cumulative_mass_function,
                                     double *mass_range,
                                     int length,
                                     double *halo_mass){

  double low, high;

  low = min(cumulative_mass_function, length);
  high = max(cumulative_mass_function, length);

  *halo_mass = linear_interp(get_random_uniform(low, high),
                             cumulative_mass_function, mass_range, length);

  return _success_;
}



/**
  * Generate subhaloes at given redshift for a given parent halo via merger tree.
  * Save the information to file in a subfolder named "data".
  *
  * Input:
  * @ mergers_params -> parameters for subhaloes catalogue
  * @ file_pointer -> pointer to the file to write
  * @ M_to_R -> see https://bdiemer.bitbucket.io/colossus/halo_mass_so.html
  * @ z_array -> redshift array
  * @ len_track -> length of the accretion track array
  * @ want_redshift_at_merging -> switch for choosing whether to calculate redshift at merging
  * @ cosmo_params -> cosmological parameters
  * @ cosmo_time -> cosmological time
  * Return:
  * £ count -> number of the generated subhaloes
  **/
int get_subhaloes_1st_order(mergers_parameters *mergers_params,
                            int order,
                            FILE *file_pointer,
                            //double *M_to_R,
                            double *z_array,
                            int len_track,
                            int want_redshift_at_merging,
                            cosmological_parameters *cosmo_params,
                            struct cosmological_time *cosmo_time,
                            int *final_count) {

  // This function has been used to run initially a merger tree.
  // Its purpose is beyond the statistical semi-empirical approach.
  // Do not use it, unless specific needs.

  // variables declaration
  int i, j, N;
  double *delta_shmf;
  double *shmf, *shmf2;
  double *cum_dshmf;
  double max_cumu;
  double frac_part;
  double sub_mass;
  double random;
  int int_part, count, count_higher_order;
  double low;
  double high;
  double redshift_infall;
  double merging_timescale;
  //double M_to_R_at_this_z;

  int len_merger_track;
  double *redshift_for_merger_track;
  double *this_merger_track;
  //double *M_to_R_for_merger_track;

  //Check that the mass of the halo at z is higher than at z+dz
  dream_call(check_accretion(mergers_params->halo_mass_at_z, mergers_params->halo_mass_at_z_plus_dz), _mass_accretion_error_message_);

  //Choose the sHMF, default: total
  //double Params[5] = {0.22, -0.91, 6., 3., 1.}; //Jiang & van den Bosch 2016 Table A1, Unevolved, total
  double Params[6] = {0.13, -0.83, 1.33, -0.02, 5.67, 1.19}; //Jiang & van den Bosch 2016 Eq. (14), Unevolved, 1st order

  dream_call(double_malloc(mergers_params->length, &shmf), _alloc_error_message_);
  dream_call(double_malloc(mergers_params->length, &shmf2), _alloc_error_message_);
  dream_call(double_malloc(mergers_params->length, &delta_shmf), _alloc_error_message_);

  dream_call(vdB_USHMF_1st_order(Params, mergers_params->halo_mass_at_z, mergers_params->subhalo_mass_range, mergers_params->length, &shmf), _mass_function_error_message_);
  dream_call(vdB_USHMF_1st_order(Params, mergers_params->halo_mass_at_z_plus_dz, mergers_params->subhalo_mass_range, mergers_params->length, &shmf2), _mass_function_error_message_);

  //delta sHMF at M_host(z=z)
  for (i=0; i<mergers_params->length; i++){
    *(delta_shmf+i) = ( (*(shmf+i)) - (*(shmf2+i)) ) * mergers_params->subhalo_mass_bin; //[N]
  }

  dream_call(double_calloc(mergers_params->length, &cum_dshmf), _alloc_error_message_);

  dream_call(cumsum(delta_shmf, mergers_params->length, &cum_dshmf), _cumsum_error_message_); //Cumulative dshmf

  dream_call(dealloc(shmf), _dealloc_error_message_);
  dream_call(dealloc(shmf2), _dealloc_error_message_);
  dream_call(dealloc(delta_shmf), _dealloc_error_message_);

  //count = compute_total_haloes_number(delta_shmf, mergers_params->length);
  dream_call(compute_total_haloes_number(delta_shmf, mergers_params->length, &count), _haloes_number_error_message_);

  count_higher_order = 0;

  //M_to_R_at_this_z = *(M_to_R+find_index(z_array, redshift, len_track));

  if (want_redshift_at_merging == _False_) {

    for (i=0; i<count; i++){

      dream_call(generate_halo_from_mass_function(cum_dshmf, mergers_params->subhalo_mass_range, mergers_params->length, &sub_mass), _mass_from_pdf_error_message_);
      redshift_infall = get_random_uniform(mergers_params->redshift, mergers_params->redshift + mergers_params->redshift_bin);
      //merging_timescale = ging_timescale(halo_mass_at_z, sub_mass, fudge, redshift, h, Om, Hz, H0, M_to_R_at_this_z);
      dream_call(compute_merging_timescale(mergers_params->halo_mass_at_z, sub_mass, mergers_params->redshift, mergers_params, cosmo_params, &merging_timescale), _merging_timescale_error_message_);

      fprintf(file_pointer, "%d %d %lf %lf %.5le\n", mergers_params->id, order, sub_mass, redshift_infall, merging_timescale);


      if ((order < mergers_params->max_order) && (redshift_infall < mergers_params->redshift_max)) {

        len_merger_track = (int)((mergers_params->redshift_max - redshift_infall) / mergers_params->redshift_bin);
        dream_call(double_calloc(len_merger_track, &redshift_for_merger_track), _alloc_error_message_);
        dream_call(double_calloc(len_merger_track, &this_merger_track), _alloc_error_message_);
        //M_to_R_for_merger_track = (double *)calloc(len_merger_track, sizeof(double));

        dream_call(arange(redshift_infall, mergers_params->redshift_max, mergers_params->redshift_bin, len_merger_track, &redshift_for_merger_track), _arange_error_message);

        this_merger_track = Mass_acc_history_VDB_(sub_mass, redshift_for_merger_track, cosmo_params->h, cosmo_params->Om0, cosmo_params->Ob0, cosmo_params->sigma8, cosmo_params->ns, len_merger_track);

        //for (j=0; j<len_merger_track; j++) {
        //  *(M_to_R_for_merger_track+len_merger_track-1+j) = *(M_to_R+len_merger_track-1+j);
        //}

        for (j=0; j<len_merger_track-1; j++){

          dream_call(get_subhaloes_1st_order(mergers_params,
                                      order+1,
                                      file_pointer,
                                      //M_to_R_for_merger_track,
                                      redshift_for_merger_track,
                                      len_merger_track,
                                      want_redshift_at_merging,
                                      cosmo_params,
                                      cosmo_time,
                                      &N), _get_subhaloes_error_message_);

          count_higher_order += N;

        }

        dream_call(dealloc(redshift_for_merger_track), _dealloc_error_message_);
        dream_call(dealloc(this_merger_track), _dealloc_error_message_);
        //free(M_to_R_for_merger_track);
        //M_to_R_for_merger_track = NULL;

      }

    }

  } else if (want_redshift_at_merging == _True_) {

    double redshift_at_merging, age_at_merge;
    double age_at_0 = linear_interp(0., cosmo_time->redshift, cosmo_time->age, cosmo_time->length);

    for (i=0; i<count; i++){

      dream_call(generate_halo_from_mass_function(cum_dshmf, mergers_params->subhalo_mass_range, mergers_params->length, &sub_mass), _mass_from_pdf_error_message_);
      redshift_infall = get_random_uniform(mergers_params->redshift, mergers_params->redshift + mergers_params->redshift_bin);
      //merging_timescale = compute_merging_timescale(halo_mass_at_z, sub_mass, fudge, redshift, h, Om, Hz, H0, M_to_R_at_this_z);
      dream_call(compute_merging_timescale(mergers_params->halo_mass_at_z, sub_mass, mergers_params->redshift, mergers_params, cosmo_params, &merging_timescale), _merging_timescale_error_message_);
      age_at_merge = linear_interp(redshift_infall, cosmo_time->redshift, cosmo_time->age, cosmo_time->length) + merging_timescale;

      if (age_at_merge > age_at_0) {
        redshift_at_merging = -1.;
      } else if (age_at_merge <= age_at_0) {
        redshift_at_merging = linear_interp(age_at_merge, cosmo_time->age, cosmo_time->redshift, cosmo_time->length);
      }

      fprintf(file_pointer, "%d %lf %lf %.5le %lf\n", mergers_params->id, sub_mass, redshift_infall, merging_timescale, redshift_at_merging);


      if ((order < mergers_params->max_order) && (redshift_infall < mergers_params->redshift_max)) {

        len_merger_track = (int)((mergers_params->redshift_max - redshift_infall) / mergers_params->redshift_bin);
        dream_call(double_calloc(len_merger_track, &redshift_for_merger_track), _alloc_error_message_);
        dream_call(double_calloc(len_merger_track, &this_merger_track), _alloc_error_message_);
        //M_to_R_for_merger_track = (double *)calloc(len_merger_track, sizeof(double));

        dream_call(arange(redshift_infall, mergers_params->redshift_max, mergers_params->redshift_bin, len_merger_track, &redshift_for_merger_track), _arange_error_message);

        this_merger_track = Mass_acc_history_VDB_(sub_mass, redshift_for_merger_track, cosmo_params->h, cosmo_params->Om0, cosmo_params->Ob0, cosmo_params->sigma8, cosmo_params->ns, len_merger_track);

        //for (j=0; j<len_merger_track; j++) {
        //  *(M_to_R_for_merger_track+len_merger_track-1+j) = *(M_to_R+len_merger_track-1+j);
        //}

        for (j=0; j<len_merger_track-1; j++){

          dream_call(get_subhaloes_1st_order(mergers_params,
                                      order+1,
                                      file_pointer,
                                      //M_to_R_for_merger_track,
                                      redshift_for_merger_track,
                                      len_merger_track,
                                      want_redshift_at_merging,
                                      cosmo_params,
                                      cosmo_time, &N), _get_subhaloes_error_message_);

          count_higher_order += N;

        }

        dream_call(dealloc(redshift_for_merger_track), _dealloc_error_message_);
        dream_call(dealloc(this_merger_track), _dealloc_error_message_);
        //free(M_to_R_for_merger_track);
        //M_to_R_for_merger_track = NULL;

      }

    }

  }

  count += count_higher_order;

  dream_call(dealloc(cum_dshmf), _dealloc_error_message_);

  //return count;
  *final_count = count;

  return _success_;

}



/**
  * Generate subhaloes at given redshift for a given parent halo statistically.
  * Save the information to file in a subfolder named "data".
  *
  * Input:
  * @ mergers_params -> parameters for subhaloes catalogue
  * @ file_pointer -> pointer to the file to write
  * @ want_redshift_at_merging -> switch for choosing whether to calculate redshift at merging
  * @ SHMF -> subhalo mass functions
  * @ cosmo_params -> cosmological parameters
  * @ cosmo_time -> cosmological time
  * * final_count -> pointer to the number of the generated subhaloes
  **/
int get_subhaloes(mergers_parameters *mergers_params,
                  FILE *file_pointer,
                  int want_redshift_at_merging,
                  subhalo_mass_functions *SHMF,
                  cosmological_parameters *cosmo_params,
                  cosmological_time *cosmo_time,
                  int *final_count) {

  int i, j;
  int order;
  int int_part, count;
  double *delta_shmf;
  double *shmf_at_z, *shmf_at_z_plus_dz;
  double *cum_dshmf, *cum_shmf_1st;
  double max_cumu, frac_part;
  double sub_mass, random, its_parent_mass;
  double redshift_infall, merging_timescale;
  double redshift_infall_lower_limit;
  double redshift_at_merging, age_at_merge, age_at_0;

  dream_call(check_accretion(mergers_params->halo_mass_at_z,
                              mergers_params->halo_mass_at_z_plus_dz),
              _mass_accretion_error_message_);

  dream_call(double_malloc(mergers_params->length,
                            &shmf_at_z),
              _alloc_error_message_);

  dream_call(double_malloc(mergers_params->length,
                            &shmf_at_z_plus_dz),
              _alloc_error_message_);

  dream_call(double_malloc(mergers_params->length,
                            &delta_shmf),
              _alloc_error_message_);

  dream_call(vdB_USHMF(params_jiang_vdb_total_unevolved,
                        mergers_params->halo_mass_at_z,
                        mergers_params->subhalo_mass_range,
                        mergers_params->length,
                        &shmf_at_z),
              _mass_function_error_message_);

  dream_call(vdB_USHMF(params_jiang_vdb_total_unevolved,
                        mergers_params->halo_mass_at_z_plus_dz,
                        mergers_params->subhalo_mass_range,
                        mergers_params->length,
                        &shmf_at_z_plus_dz),
              _mass_function_error_message_);

  //delta SHMF at M_host(z=z) [N]
  dream_call(compute_delta_mass_function(mergers_params->length,
                                          &shmf_at_z,
                                          &shmf_at_z_plus_dz,
                                          mergers_params->subhalo_mass_bin,
                                          &delta_shmf),
              _delta_mass_func_error_message_);

  dream_call(double_calloc(mergers_params->length,
                            &cum_dshmf),
              _alloc_error_message_);

  dream_call(cumsum(delta_shmf,
                     mergers_params->length,
                     &cum_dshmf),
              _cumsum_error_message_); //Cumulative delta-SHMF

  dream_call(compute_total_haloes_number(delta_shmf,
                                          mergers_params->length,
                                          &count),
              _haloes_number_error_message_);

  dream_call(double_calloc(mergers_params->length,
                            &cum_shmf_1st),
              _alloc_error_message_);

  dream_call(cumsum(SHMF->first,
                     mergers_params->length,
                     &cum_shmf_1st),
              _cumsum_error_message_); //Cumulative 1st order SHMF

  if (want_redshift_at_merging == _True_) {

    age_at_0 = linear_interp(0., cosmo_time->redshift, cosmo_time->age, cosmo_time->length);
  }

  int count_unwanted_orders = 0;




  /*****************************************/
  double try_sub_mass[count];
  double mass_sum;
  int times = 0, times_threshold=20;
  do {
    mass_sum = 0.;
    for (i=0; i<count; i++){
      dream_call(generate_halo_from_mass_function(cum_dshmf,
                                                   mergers_params->subhalo_mass_range,
                                                   mergers_params->length,
                                                   &try_sub_mass[i]),
                  _mass_from_pdf_error_message_);
      mass_sum += pow(10., try_sub_mass[i]);
    }
    mass_sum = log10(mass_sum);
    times += 1;
  } while ((mass_sum>mergers_params->halo_mass_at_z) && (times < times_threshold));

  if (times >= times_threshold){
    int loop_to_be_repeated;
    double available_mass;
    do {
      loop_to_be_repeated = _False_;
      available_mass = pow(10., mergers_params->halo_mass_at_z);
      for (i=0; i<count; i++){
        //printf("starting %d %d %lf %lf\n", i, count, log10(available_mass), mergers_params->subhalo_mass_range[0]);
        if (available_mass<pow(10., mergers_params->subhalo_mass_range[0])){
          loop_to_be_repeated = _True_; break;
        }
        do {
          dream_call(generate_halo_from_mass_function(cum_dshmf,
                                                       mergers_params->subhalo_mass_range,
                                                       mergers_params->length,
                                                       &try_sub_mass[i]),
                      _mass_from_pdf_error_message_);
        } while (pow(10.,try_sub_mass[i]) >= available_mass);
        available_mass -= pow(10., try_sub_mass[i]);
        //printf("finished %d %lf %lf\n", i, log10(available_mass), try_sub_mass[i]);
      }
    } while (loop_to_be_repeated==_True_);
  }
  /*****************************************/

  for (i=0; i<count; i++){

    /*dream_call(generate_halo_from_mass_function(cum_dshmf,
                                                 mergers_params->subhalo_mass_range,
                                                 mergers_params->length,
                                                 &sub_mass),
                _mass_from_pdf_error_message_);*/
    sub_mass = try_sub_mass[i]; /******************************/

    dream_call(assign_subhalo_order(mergers_params->max_order,
                                     sub_mass,
                                     mergers_params->halo_mass_at_z0,
                                     SHMF, &order),
                _assign_order_error_message_);

    if (order == 1) {

      redshift_infall = get_random_uniform(mergers_params->redshift, \
                       mergers_params->redshift + mergers_params->redshift_bin);

      dream_call(compute_merging_timescale(mergers_params->halo_mass_at_z,
                                            sub_mass,
                                            redshift_infall,
                                            mergers_params,
                                            cosmo_params,
                                            &merging_timescale),
                  _merging_timescale_error_message_);

    } else if (order>1) {

      dream_call(generate_infall_redshift_higher_order(sub_mass,
                                                        order,
                                                        mergers_params->redshift,
                                                        mergers_params->redshift_bin,
                                                        mergers_params->redshift_max,
                                                        &redshift_infall),
                  _infall_redshift_error_message_);

      dream_call(generate_halo_from_mass_function(cum_shmf_1st,
                                                   mergers_params->subhalo_mass_range,
                                                   mergers_params->length,
                                                   &its_parent_mass),
                  _mass_from_pdf_error_message_);

      dream_call(compute_merging_timescale(its_parent_mass,
                                            sub_mass,
                                            redshift_infall,
                                            mergers_params,
                                            cosmo_params,
                                            &merging_timescale),
                  _merging_timescale_error_message_);

    }

    if ((order > 0) && (redshift_infall < mergers_params->redshift_max) && (redshift_infall >= 0.)){

      if (want_redshift_at_merging == _False_) {

        fprintf(file_pointer, "%d %d %lf %lf %.5le\n", mergers_params->id, \
                           order, sub_mass, redshift_infall, merging_timescale);

      } else if (want_redshift_at_merging == 1) {

        age_at_merge = linear_interp(redshift_infall, cosmo_time->redshift, \
                       cosmo_time->age, cosmo_time->length) + merging_timescale;

        if (age_at_merge > age_at_0) {

          redshift_at_merging = -1.;

        } else if (age_at_merge <= age_at_0) {

          redshift_at_merging = linear_interp(age_at_merge, \
                     cosmo_time->age, cosmo_time->redshift, cosmo_time->length);
        }

        fprintf(file_pointer, "%d %d %lf %lf %.5le %lf\n", mergers_params->id, \
          order, sub_mass, redshift_infall, merging_timescale, redshift_at_merging);

      }

    } else {
      count_unwanted_orders += 1;
    }
  }

  dream_call(dealloc(shmf_at_z),
              _dealloc_error_message_);

  dream_call(dealloc(shmf_at_z_plus_dz),
              _dealloc_error_message_);

  dream_call(dealloc(delta_shmf),
              _dealloc_error_message_);

  dream_call(dealloc(cum_dshmf),
              _dealloc_error_message_);

  dream_call(dealloc(cum_shmf_1st),
              _dealloc_error_message_);

  *final_count = count - count_unwanted_orders;

  return _success_;

}




/**
  * Generate subhaloes for a given parent halo statistically.
  * Save the information to file in a subfolder named "data".
  *
  * Input:
  * @ mergers_params -> parameters for subhaloes catalogue
  * @ track -> accretion track of the parent halo [log(M/Msun)]
  * @ z_range -> redshift array of the accretion track
  * @ len_track -> length of the track array
  * @ use_merger_tree -> switch for choosing whether to use merger tree
  * @ cosmo_params -> cosmological parameters
  * @ output_folder -> name of the output folder
  * @ want_redshift_at_merging -> switch for choosing whether to calculate redshift at merging
  * @ cosmo_time -> cosmological time
  * @ SHMF -> subhalo mass functions
  * Return:
  * £ Num_mergers -> number of the generated subhaloes
  **/
int generate_mergers(mergers_parameters *mergers_params,
                     double *track,
                     double *z_range,
                     int len_track,
                     int use_merger_tree,
                     cosmological_parameters *cosmo_params,
                     char *output_folder,
                     int want_redshift_at_merging,
                     cosmological_time *cosmo_time,
                     subhalo_mass_functions *SHMF) {


  int i, N, Num_mergers;
  FILE *file_pointer;

  char filename[24] = "data/output_mergers.txt";
  const size_t len1 = strlen(output_folder);
  const size_t len2 = strlen(filename);
  char *output_filename;

  dream_call(char_malloc(len1 + len2 + 1,
                          &output_filename),
              _alloc_error_message_);

  strcpy(output_filename, output_folder);
  strcpy(output_filename + len1, filename);

  file_pointer = fopen(output_filename, "a");

  Num_mergers = 0;

  for (i=0; i<len_track-1; i++){

    mergers_params->halo_mass_at_z = *(track+i);
    mergers_params->halo_mass_at_z_plus_dz = *(track+i+1);
    mergers_params->redshift = *(z_range+i);

    N = 0;

    if (use_merger_tree == 1) {

      dream_call(get_subhaloes_1st_order(mergers_params,
                                          1,
                                          file_pointer,
                                          //M_to_R,
                                          z_range,
                                          len_track,
                                          want_redshift_at_merging,
                                          cosmo_params,
                                          cosmo_time,
                                          &N),
                  _get_subhaloes_error_message_);

      Num_mergers += N;

    } else if (use_merger_tree == 0) {

      dream_call(get_subhaloes(mergers_params,
                                file_pointer,
                                want_redshift_at_merging,
                                SHMF,
                                cosmo_params,
                                cosmo_time,
                                &N),
                  _get_subhaloes_error_message_);

      Num_mergers += N;

    }

  }

  fclose(file_pointer);

  dream_call(dealloc(output_filename),
              _dealloc_error_message_);

  return Num_mergers;
}
