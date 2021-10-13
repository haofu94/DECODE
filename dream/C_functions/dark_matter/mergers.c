/** @ file mergers.c
  *
  * Written by Hao Fu
  *
  * The main goal of this module is to generate a dark matter subhalo catalogue.
  **/


#ifndef dm_mergers
# define dm_mergers


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

    //A = 17.59386; beta = 1.25453; delta = 13.06571; gamma = 2.27080; alpha = 1.39963;
    A = 10.05658; beta = 1.06991; delta = 13.78276; gamma = 1.85808; alpha = 1.76572;

  } else if ((subhalo_order == 2) && (11. <= subhalo_mass) && (subhalo_mass < 12.)) {

    //A = 23.86118; beta = 1.64827; delta = 9.47616; gamma = 0.69930; alpha = 1.39298;
    A = 24.92779; beta = 1.76825; delta = 8.86702; gamma = 0.42264; alpha = 1.93852;

  } else if ((subhalo_order == 2) && (12. <= subhalo_mass)) {

    //A = 22.32932; beta = 2.40317; delta = 2.97333; gamma = -4.17881; alpha = 1.19073;
    A = 24.80613; beta = 2.67003; delta = 2.30625; gamma = -2.59638; alpha = 1.85692;

  } else if ((subhalo_order >= 3) && (subhalo_mass < 11.)) {

    //A = 23.26552; beta = 1.66035; delta = 11.89017; gamma = 1.24507; alpha = 2.32610;
    A = 8.37148; beta = 1.49595; delta = 13.49528; gamma = 2.99628; alpha = 3.28954;

  } else if ((subhalo_order >= 3) && (11. <= subhalo_mass)) {

    //A = 24.47686; beta = 2.20970; delta = 4.40460; gamma = -2.77968; alpha = 2.06852;
    A = 23.71313; beta = 2.54297; delta = 3.40928; gamma = -3.56339; alpha = 3.12188;

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

  //*order = 1; /////////////////// DELETE THIS

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

  //printf("%lf %lf %lf\n", redshift, Om, Hz);

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

  //printf("%lf %lf %lf\n", parent_mass, subhalo_mass, mergers_params->fudge);

  dream_call(compute_dynamical_timescale(z,
                                          cosmo_params,
                                          &dynamical_friction),
              _dynamical_timescale_error_message_); //[Gyr]
              //printf("%lf %lf\n", z, dynamical_friction);

  if (mergers_params->fudge<0.){

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
  * mass_function -> mass function (MUST be in units of [dex^-1]
  * mass_range -> mass range in units of log10(M/Msun)
  * length -> length of the mass function array
  * count -> pointer to the count variable
  **/
int compute_total_haloes_number(double *mass_function,
                                double *mass_range,
                                int length,
                                int *count){

  int i;
  double sum;

  //printf("length = %d\n", length);
  /*for (i=0; i<length; i++){
    printf("%lf %lf\n", mass_range[i], mass_function[i]);
  }*/

  sum = integrate_trapz(mass_function, mass_range, length);

  *count = (int)sum;

  //Fractional part
  //Random number between the highest integer and itself +1
  if (get_random_uniform((double)(*count), (double)(*count+1)) < sum) {

    *count += 1;
  }

  //printf("number = %d\n", *count);
  //exit(0);

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
                            double main_progenitor,
                            double parent_halo_mass,
                            double *z_track,
                            double *track,
                            int len_track,
                            int order,
                            FILE *file_pointer,
                            int want_redshift_at_merging,
                            subhalo_mass_functions *SHMF,
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
  double *cum_shmf;
  double max_cumu;
  double frac_part;
  double sub_mass;
  double random;
  int int_part, count, count_higher_order;
  double low;
  double high;
  double redshift_infall;
  double merging_timescale;

  int len_merger_track;
  double *redshift_for_merger_track;
  double *this_merger_track;

  double Params[6] = {0.13, -0.83, 1.33, -0.02, 5.67, 1.19}; //Jiang & van den Bosch 2016 Eq. (14), Unevolved, 1st order

  double res = -4.;
  int len_range = (int)((parent_halo_mass - main_progenitor - res) / 0.1);

  double *subhalo_mass_range;

  dream_call(double_malloc(len_range,
                           &subhalo_mass_range),
             _alloc_error_message_);

  dream_call(arange(main_progenitor+res,
                    parent_halo_mass,
                    0.1, len_range,
                    &subhalo_mass_range),
             _arange_error_message);

  double *shmf_1st;
  dream_call(double_malloc(len_range, &shmf_1st),
             _alloc_error_message_);

  { double this_psi; int index_psi;
    for(i=0; i<len_range; i++){
      this_psi = pow(10., subhalo_mass_range[i]) / pow(10., track[0]);
      index_psi = find_index(SHMF->psi, this_psi, SHMF->length_psi);
      shmf_1st[i] = *(SHMF->first+index_psi);
    }
  }

  dream_call(double_calloc(len_range, &cum_shmf),
             _alloc_error_message_);

  dream_call(cumsum(shmf_1st, len_range, &cum_shmf),
             _cumsum_error_message_); //Cumulative dshmf

  dream_call(compute_total_haloes_number(shmf_1st, subhalo_mass_range,
                                         len_range, &count),
             _haloes_number_error_message_);

  count_higher_order = 0;

  if (want_redshift_at_merging == _False_) {

    for (i=0; i<count; i++){

      dream_call(generate_halo_from_mass_function(cum_shmf, subhalo_mass_range,
                                                  len_range, &sub_mass),
                 _mass_from_pdf_error_message_);

      double dphi, dz;
      double z_array[len_track-1], PDFz[len_track-1];
      double this_psi, this_psi_1;
      for (j=0; j<len_track-1; j++){
        this_psi = pow(10.,sub_mass - track[j]);
        this_psi_1 = pow(10.,sub_mass - track[j+1]);
        dphi = linear_interp(this_psi, SHMF->psi, SHMF->first, SHMF->length_psi) - \
                linear_interp(this_psi_1 , SHMF->psi, SHMF->first, SHMF->length_psi);
        dz = *(z_track+j+1) - (*(z_track+j));
        if (this_psi < 1.) { *(PDFz+j) = dphi / dz; } else { PDFz[j] = 0.; }
        *(z_array+j) = 0.5 * ( *(z_track+j) + (*(z_track+j+1)) );
      }
      redshift_infall = get_random_from_distribution(z_array, PDFz, len_track-1);

      //fprintf(file_pointer, "%d %d %lf %lf %.5le\n", mergers_params->id, order, sub_mass, redshift_infall, merging_timescale);
      fprintf(file_pointer, "%d %d %.4lf %.4lf\n", mergers_params->id, order, sub_mass, redshift_infall);


      if ((order < mergers_params->max_order) && (redshift_infall < mergers_params->redshift_max)) {

        len_merger_track = (int)((mergers_params->redshift_max - redshift_infall) / mergers_params->redshift_bin);

        if ( ( sub_mass > (main_progenitor+res) ) && (len_merger_track > 0) ) {

          dream_call(double_calloc(len_merger_track,
                                   &redshift_for_merger_track),
                     _alloc_error_message_);

          dream_call(double_calloc(len_merger_track,
                                   &this_merger_track),
                     _alloc_error_message_);

          dream_call(arange(redshift_infall,
                            mergers_params->redshift_max,
                            mergers_params->redshift_bin,
                            len_merger_track,
                            &redshift_for_merger_track),
                     _arange_error_message);

          dream_call(Mass_acc_history_VDB_(sub_mass,
                                           redshift_for_merger_track,
                                           cosmo_params->h,
                                           cosmo_params->Om0,
                                           cosmo_params->Ob0,
                                           cosmo_params->sigma8,
                                           cosmo_params->ns,
                                           len_merger_track,
                                           &this_merger_track),
                     "error\n");

          dream_call(get_subhaloes_1st_order(mergers_params,
                                             main_progenitor,
                                             sub_mass,
                                             redshift_for_merger_track,
                                             this_merger_track,
                                             len_merger_track,
                                             order+1,
                                             file_pointer,
                                             want_redshift_at_merging,
                                             SHMF, cosmo_params, cosmo_time, &N),
                     _get_subhaloes_error_message_);

          count_higher_order += N;

          dream_call(dealloc(redshift_for_merger_track),
                     _dealloc_error_message_);

          dream_call(dealloc(this_merger_track),
                     _dealloc_error_message_);

        }
      }

    }

  }

  count += count_higher_order;

  dream_call(dealloc(cum_shmf),
             _dealloc_error_message_);

  dream_call(dealloc(shmf_1st),
             _dealloc_error_message_);

  dream_call(dealloc(subhalo_mass_range),
             _dealloc_error_message_);

  *final_count = count;

  return _success_;

}



/**
  * Generate subhaloes at given redshift for a given parent halo statistically.
  * Save the information to file in a subfolder named "data".
  *
  * Input:
  * @ mergers_params -> parameters for subhaloes catalogue
  * @ halo_accretion -> dark matter halo accretion track
  * @ file_pointer -> pointer to the file to write
  * @ want_redshift_at_merging -> switch for choosing whether to calculate redshift at merging
  * @ SHMF -> subhalo mass functions
  * @ cosmo_params -> cosmological parameters
  * @ cosmo_time -> cosmological time
  * * final_count -> pointer to the number of the generated subhaloes
  **/
int get_subhaloes(mergers_parameters *mergers_params,
                  DM_halo_accretion *halo_accretion,
                  FILE *file_pointer,
                  int want_redshift_at_merging,
                  subhalo_mass_functions *SHMF,
                  cosmological_parameters *cosmo_params,
                  cosmological_time *cosmo_time,
                  int *final_count) {

  int i, j;
  int order;
  int int_part, count;
  double *subhalo_mass_function_total;
  double *cum_shmf;
  double *cum_shmf_1st;
  double max_cumu, frac_part;
  double random, its_parent_mass;
  double redshift_infall_lower_limit;
  double redshift_at_merging, age_at_merge, age_at_0;

  int *ID_first_orders;
  int *all_orders;
  double *all_masses, *all_zinfalls, *all_tau_merges;


  do { //repeat loop until the generated number of subhaloes is positive

    dream_call(double_malloc(mergers_params->length,
                             &subhalo_mass_function_total),
               _alloc_error_message_);

    { double this_psi; int index_psi;

      for(i=0; i<mergers_params->length; i++){

        this_psi = pow(10., mergers_params->subhalo_mass_range[i]) / pow(10., mergers_params->halo_mass_at_z0);
        index_psi = find_index(SHMF->psi, this_psi, SHMF->length_psi);

        *(subhalo_mass_function_total+i) = *(SHMF->total+index_psi);
      }
    }

    dream_call(compute_total_haloes_number(subhalo_mass_function_total,
                                           mergers_params->subhalo_mass_range,
                                           mergers_params->length,
                                           &count),
               _haloes_number_error_message_);

     dream_call(double_calloc(mergers_params->length,
                              &cum_shmf),
                _alloc_error_message_);

     dream_call(cumsum(subhalo_mass_function_total,
                       mergers_params->length,
                       &cum_shmf),
                _cumsum_error_message_); //Cumulative SHMF

    dream_call(double_calloc(SHMF->length_psi,
                             &cum_shmf_1st),
               _alloc_error_message_);

    dream_call(cumsum(SHMF->first,
                      SHMF->length_psi,
                      &cum_shmf_1st),
               _cumsum_error_message_); //Cumulative 1st order SHMF

    if (want_redshift_at_merging == _True_) {

      age_at_0 = linear_interp(0., cosmo_time->redshift, cosmo_time->age, cosmo_time->length);
    }

    int count_unwanted_orders = 0;

    int count_first_orders = 0;

    dream_call(int_calloc(count_first_orders,
                          &ID_first_orders),
               _alloc_error_message_);

    dream_call(int_calloc(count,
                          &all_orders),
               _alloc_error_message_);

    dream_call(double_calloc(count,
                             &all_masses),
               _alloc_error_message_);

    dream_call(double_calloc(count,
                             &all_zinfalls),
               _alloc_error_message_);

    dream_call(double_calloc(count,
                             &all_tau_merges),
               _alloc_error_message_);


    for (i=0; i<count; i++){

      dream_call(generate_halo_from_mass_function(cum_shmf,
                                                  mergers_params->subhalo_mass_range,
                                                  mergers_params->length,
                                                  &all_masses[i]),
                 _mass_from_pdf_error_message_);

      dream_call(assign_subhalo_order(mergers_params->max_order,
                                      all_masses[i],
                                      mergers_params->halo_mass_at_z0,
                                      SHMF, &all_orders[i]),
                 _assign_order_error_message_);

      if (all_orders[i] == 1) {

        count_first_orders += 1;
        ID_first_orders = realloc(ID_first_orders, count_first_orders*sizeof(int));
        ID_first_orders[count_first_orders-1] = get_int_random_uniform(100, 999);

        double dphi, dz;
        double z_array[halo_accretion->length-1], PDFz[halo_accretion->length-1];
        double this_psi, this_psi_1;

        for (j=0; j<halo_accretion->length-1; j++){

          this_psi = pow(10., all_masses[i] - halo_accretion->track[j]);
          this_psi_1 = pow(10., all_masses[i] - halo_accretion->track[j+1]);
          dphi = linear_interp(this_psi, SHMF->psi, SHMF->first, SHMF->length_psi) - \
                  linear_interp(this_psi_1 , SHMF->psi, SHMF->first, SHMF->length_psi);
          dz = *(halo_accretion->redshift+j+1) - (*(halo_accretion->redshift+j));

          if (this_psi < 1.) { *(PDFz+j) = dphi / dz; } else { PDFz[j] = 0.; }
          *(z_array+j) = 0.5 * ( *(halo_accretion->redshift+j) + (*(halo_accretion->redshift+j+1)) );

        }

        all_zinfalls[i] = get_random_from_distribution(z_array, PDFz, halo_accretion->length-1);

        dream_call(compute_merging_timescale(linear_interp(all_zinfalls[i],
                                                           halo_accretion->redshift,
                                                           halo_accretion->track,
                                                           halo_accretion->length),
                                             all_masses[i],
                                             all_zinfalls[i],
                                             mergers_params,
                                             cosmo_params,
                                             &all_tau_merges[i]),
                   _merging_timescale_error_message_);

      } else if (all_orders[i]>1) {

        dream_call(generate_infall_redshift_higher_order(all_masses[i],
                                                         all_orders[i],
                                                         mergers_params->redshift,
                                                         mergers_params->redshift_bin,
                                                         mergers_params->redshift_max,
                                                         &all_zinfalls[i]),
                   _infall_redshift_error_message_);

        do {

          dream_call(generate_halo_from_mass_function(cum_shmf_1st,
                                                      mergers_params->subhalo_mass_range,
                                                      mergers_params->length,
                                                      &its_parent_mass),
                     _mass_from_pdf_error_message_);

        } while (its_parent_mass <= all_masses[i]);

        dream_call(compute_merging_timescale(its_parent_mass,
                                             all_masses[i],
                                             all_zinfalls[i],
                                             mergers_params,
                                             cosmo_params,
                                             &all_tau_merges[i]),
                   _merging_timescale_error_message_);

      }

    }

    int num_first_orders = count_first_orders;
    count_first_orders = 0;

    for (i=0; i<count; i++){

      if ((all_orders[i] > 0) && (all_zinfalls[i] < mergers_params->redshift_max) && (all_zinfalls[i] >= 0.)){

        if (all_orders[i] == 1) {

          count_first_orders += 1;

          fprintf(file_pointer, "%d -1 %d %d %lf %lf %.5le\n", mergers_params->id, ID_first_orders[count_first_orders-1], \
                             all_orders[i], all_masses[i], all_zinfalls[i], all_tau_merges[i]);

        } else if ( (all_orders[i] > 1) && (num_first_orders > 0) ) {

          fprintf(file_pointer, "%d %d -1 %d %lf %lf %.5le\n", mergers_params->id, ID_first_orders[get_int_random_uniform(1,num_first_orders)-1], \
                             all_orders[i], all_masses[i], all_zinfalls[i], all_tau_merges[i]);

        } else {

          count_unwanted_orders += 1;

        }

      } else {
        count_unwanted_orders += 1;
      }

    }

    *final_count = count - count_unwanted_orders;

  } while (*final_count<0);

  dream_call(dealloc(subhalo_mass_function_total),
             _dealloc_error_message_);

  dream_call(dealloc(cum_shmf),
              _dealloc_error_message_);

  dream_call(dealloc(cum_shmf_1st),
              _dealloc_error_message_);

  dream_call(dealloc(all_orders),
             _dealloc_error_message_);

  dream_call(dealloc(all_masses),
             _dealloc_error_message_);

  dream_call(dealloc(all_zinfalls),
             _dealloc_error_message_);

  dream_call(dealloc(all_tau_merges),
             _dealloc_error_message_);

  dream_call(dealloc(ID_first_orders),
             _dealloc_error_message_);


  if (*final_count >= 0) {

    return _success_;

  } else {

    return _failure_;

  }

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
                     DM_halo_accretion *halo_accretion,
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

  if (use_merger_tree == 1) {

    dream_call(get_subhaloes_1st_order(mergers_params, halo_accretion->track[0], halo_accretion->track[0],
                                        halo_accretion->redshift, halo_accretion->track, halo_accretion->length,
                                        1, file_pointer, want_redshift_at_merging,
                                        SHMF, cosmo_params, cosmo_time,
                                        &Num_mergers),
                _get_subhaloes_error_message_);

  } else if (use_merger_tree == 0) {

    dream_call(get_subhaloes(mergers_params,
                             halo_accretion,
                             file_pointer,
                              want_redshift_at_merging,
                              SHMF,
                              cosmo_params,
                              cosmo_time,
                              &Num_mergers),
                _get_subhaloes_error_message_);

  }

  fclose(file_pointer);

  dream_call(dealloc(output_filename),
              _dealloc_error_message_);

  return Num_mergers;
}


#endif
