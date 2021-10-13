/** @ file star_formation.c
  *
  * Written by Hao Fu
  *
  * The main goal of this module is to compute the evolution of satellite galaxies.
  **/

#ifndef star_formation_satellites
# define star_formation_satellites


#include "../include/dream.h"

#include "../dark_matter/Mhalo_to_Mstar.c"
#include "../dark_matter/mergers.c"


#define pi M_PI



/**
  * Compute the derivative of f_loss, Moster et al. 2018
  *
  * Input:
  * tau -> timescale [yr]
  * Return:
  * derivative of f_loss [yr^-1]
**/
double fdot(double tau){

  return 0.05 / (tau + 1.4e6);
  //return 0.046 / (tau + 2.76e5);
}



/**
  * Compute the quenching timescale [Gyr]
  *
  * Input:
  * @ Mstar -> stellar mass [Msun]
  * @ tau_quenching -> pointer to the tau_quenching variable
**/
int compute_quenching_timescale(double Mstar,
                                double *tau_quenching){

  *tau_quenching = -0.5 * log10(Mstar) + 5.7;

  return _success_;
}



/**
  * Compute stripping timescale [Gyr]
  *
  * Input:
  * @ Mpar -> mass of the parent halo [log(M/Msun)]
  * @ Msub -> mass of the subhalo [log(M/Msun)]
  * @ tau_stripping -> pointer to the tau_stripping variable
**/
int compute_stripping_timescale(double Mpar, double Msub,
                                double *tau_stripping){
  Mpar = pow(10., Mpar);
  Msub = pow(10., Msub);

  *tau_stripping = (1.428/2.*pi) * (Mpar/Msub) / log(1. + Mpar/Msub);
  //*tau_stripping = (1.227/2.*pi) * (Mpar/Msub) / log(1. + Mpar/Msub);

  return _success_;
}



/**
  * Compute analytical star formation rate for satellites [Msun/yr]
  * Grylls et al. 2019, Paper I, Eq. (7)
  *
  * Input:
  * @ Mstar -> stellar mass [Msun]
  * @ z -> redshift
**/
int get_star_formation_rate_satellites(double Mstar, double z, double *SFR) {

  /*
    Tomczak et al. 2016
  */

  double s0 = 0.195 + 1.157*z - 0.143*z*z;
  double M0 = pow(10., 9.244 + 0.753*z - 0.090*z*z);
  double alpha = 1.118;
  //double alpha = 0.8;

  /*
    Leja et al. 2015
  */
  /*double s0 = 0.6 + 1.22*z - 0.2*z*z;
  double M0 = pow(10., 10.3 + 0.753*z - 0.15*z*z);
  double alpha = 1.3 - 0.1*z;*/

  *SFR = pow(10., s0 - log10(1. + pow(Mstar/M0, -alpha)) ) ;//+ get_random_gaussian(0., 0.3, -5., 5.);

  return _success_;
}


/*double retta(double x){
  return - (x - 12.) * (0.3) + 0.1;
}*/

/**
  * Compute the evolved satellite mass at redshift z [log(M/Msun)]
  *
  * Input:
  * @ SMHM_params -> parameters of the SMHM relation
  * @ progenitor_mass -> mass of the progenitor [log(M/Msun)]
  * @ subhalo_mass -> mass of the subhalo [log(M/Msun)]
  * @ z -> redshift at which to compute the evolved mass
  * @ z_infall -> infall redshift of the subhalo
  * @ include_stripping -> switch to include stripping in the evolution
  * @ include_quenching -> switch to include quenching in the evolution
  * @ cosmo_time -> cosmological time
  * @ mass_at_z -> pointer to the mass_at_z variable
**/

int get_evolved_satellite_mass(stellar_mass_halo_mass *SMHM,
                               SMHM_matrix *smhm_data,
                               double progenitor_mass,
                               double subhalo_mass,
                               double z,
                               double z_infall,
                               int include_sat_SF,
                               int include_mass_loss,
                               int include_stripping,
                               int include_quenching,
                               cosmological_parameters *cosmo_params,
                               cosmological_time *cosmo_time,
                               double *mass_at_z){

  double Mstar_satellite;

  dream_call(SMHM_numerical_interp(subhalo_mass, SMHM, smhm_data, z_infall, &Mstar_satellite),
             _SMHM_error_message_);

  //if (Mstar_satellite < 10.5){Mstar_satellite += 0.5;}
  //if (Mstar_satellite < 10.8){Mstar_satellite += retta(subhalo_mass);}
  // TNG
  /*if (progenitor_mass<13.) {double Params[9] = {11.12, 0., 0.108, 0., 27.5, 0., 15.6, 0., 0.}; Mstar_satellite = SMHM_Moster(subhalo_mass, Params, 0.);}
  if ( (13<=progenitor_mass) && (progenitor_mass<14.) ) {double Params[9] = {10.85, 0., 0.127, 0., 23.6, 0., 10.1, 0., 0.}; Mstar_satellite = SMHM_Moster(subhalo_mass, Params, 0.);}
  if ( (14<=progenitor_mass) && (progenitor_mass<14.5) ) {double Params[9] = {10.93, 0., 0.137, 0., 30.5, 0., 10.9, 0., 0.}; Mstar_satellite = SMHM_Moster(subhalo_mass, Params, 0.);}
  if (progenitor_mass>=14.5) {double Params[9] = {10.85, 0., 0.129, 0., 38.6, 0., 9.5, 0., 0.}; Mstar_satellite = SMHM_Moster(subhalo_mass, Params, 0.);}*/

  if (z_infall <= z) {

    *mass_at_z = log10(Mstar_satellite);

    return _success_;
  }

  int i, j, idx;

  double tau_strip;

  dream_call(compute_stripping_timescale(progenitor_mass,
                                         subhalo_mass,
                                         &tau_strip),
             _stripping_timescale_error_message_);

  double subhalo_mass_at_t, f_DM = 1., f_new;

  double t_infall = linear_interp(z_infall, cosmo_time->redshift, cosmo_time->age, cosmo_time-> length)*1.e9;
  double tau_q = 0.; // [Gyr]
  double tau_f;
  double fq;
  double formed_stars = 0.;

  //double lookbacktime_z = linear_interp(z, cosmo_time->redshift, cosmo_time->lookback_time, cosmo_time->length);
  //double lookback_time_infall = linear_interp(z_infall, cosmo_time->redshift, cosmo_time->lookback_time, cosmo_time->length);
  double cosmic_time_z = linear_interp(z, cosmo_time->redshift, cosmo_time->age, cosmo_time->length);
  double cosmic_time_infall = linear_interp(z_infall, cosmo_time->redshift, cosmo_time->age, cosmo_time->length);

  //double zbin = 0.1;
  double dt = 1.e-2;

  //int len = 20;
  int len = (int) ( (cosmic_time_z - cosmic_time_infall) / dt ) ;
  double *z_range, *lookback_time, *cosmic_time, *SFR, *SFH, *MLR;

  dream_call(double_malloc(len,
                           &z_range),
             _alloc_error_message_);

  dream_call(double_malloc(len,
                           &lookback_time),
             _alloc_error_message_);

  dream_call(double_malloc(len,
                           &cosmic_time),
             _alloc_error_message_);

  dream_call(double_calloc(len,
                           &SFR),
             _alloc_error_message_);

  dream_call(double_calloc(len,
                           &SFH),
             _alloc_error_message_);

  dream_call(double_calloc(len,
                           &MLR),
             _alloc_error_message_);

  /*dream_call(linear_space(z,
                          z_infall,
                          len,
                          &z_range),
             _linear_space_error_message);*/

  //printf("Debug flag 1\n");
  //printf("%lf %lf %lf %d\n", cosmic_time_infall, cosmic_time_z, dt, len);

  //dream_call(arange(lookbacktime_z, lookback_time_infall, dt, len, &lookback_time), _arange_error_message);
  dream_call(arange(cosmic_time_infall, cosmic_time_z, dt, len, &cosmic_time), _arange_error_message);

  for (i=0; i<len; i++){

    *(z_range+i) = linear_interp(*(cosmic_time+i), cosmo_time->age, cosmo_time->redshift, cosmo_time->length);

    *(lookback_time+i) = linear_interp(*(cosmic_time+i), cosmo_time->age, cosmo_time->lookback_time, cosmo_time->length);

  }

  double Mstar_ini = Mstar_satellite;

  Mstar_satellite = pow(10., Mstar_satellite);

  if (include_sat_SF == _True_){

    if (len > 0) {
      if (Mstar_satellite < 1.e+9){
        tau_q = *lookback_time - 2.;
      } else if (Mstar_satellite >= 1.e+9){
        tau_q = *lookback_time - 1.;
      }
    }

    if ( (include_quenching == _True_) && (tau_q > 0.) && (len>0) ) {

      dream_call(compute_quenching_timescale(Mstar_satellite,
                                             &tau_f),
                 _quenching_timescale_error_message_);

      for (i=0; i<len; i++) {

        if (*(lookback_time+i) > tau_q) {

          fq = 1.;

        } else {

          fq = exp( -(tau_q - (*(lookback_time+i)) ) / tau_f);

        }

        dream_call(get_star_formation_rate_satellites(Mstar_satellite, *(z_range+i), &(*(SFR+i))),
                   _compute_SFR_error_message_);

        if (include_mass_loss == _True_){

          for (j=0; j<=i; j++){

           *(MLR+i) += *(SFR+j) * (cosmic_time[j+1]-cosmic_time[j]) * 1.e+9 * fdot( (cosmic_time[j+1]-cosmic_time[j]) * 1.e+9 );

          }

          *(SFR+i) -= *(MLR+i);

        }

        Mstar_satellite += (*(SFR+i)) * fq * dt * 1.e+9;

        if (include_stripping == _True_){

          f_new = (*(SFR+i)) * dt / Mstar_satellite;
          f_DM = f_DM + f_new * f_DM;
          Mstar_satellite = Mstar_satellite * (1. - exp(-14.2 * f_DM) );

        }

      }

    } else {

      for (i=0; i<len; i++) {

        dream_call(get_star_formation_rate_satellites(Mstar_satellite, *(z_range+i), &(*(SFR+i))),
                   _compute_SFR_error_message_);

        if (include_mass_loss == _True_){

          for (j=0; j<=i; j++){

            *(MLR+i) += *(SFR+j) * (cosmic_time[j+1]-cosmic_time[j]) * 1.e+9 * fdot( (cosmic_time[j+1]-cosmic_time[j]) * 1.e+9 );

          }

        }

        *(SFR+i) -= *(MLR+i);

        Mstar_satellite += (*(SFR+i)) * dt * 1.e+9;

        if (include_stripping == _True_){

          f_new = (*(SFR+i)) * dt / Mstar_satellite;
          f_DM = f_DM + f_new * f_DM;
          Mstar_satellite = Mstar_satellite * (1. - exp(-14.2 * f_DM) );

        }

      }

    }

  } else if (include_sat_SF == _False_){

    if (include_stripping == _True_){

      subhalo_mass_at_t = mass_at_t(subhalo_mass, progenitor_mass, z, z_infall, cosmic_time_z, cosmic_time_infall, cosmo_params);

      f_DM = pow(10., subhalo_mass_at_t) / pow(10., subhalo_mass);

      Mstar_satellite = Mstar_satellite * (1. - exp(-14.2 * f_DM) );

    }

  }

  /*if (Mstar_ini < 10.) {
    if (Mstar_satellite / pow(10., Mstar_ini) > 2.) {
      printf("%lf %lf %lf\n", Mstar_ini, log10(Mstar_satellite),  Mstar_satellite / pow(10., Mstar_ini));
    }
  }*/

  dream_call(dealloc(z_range),
             _dealloc_error_message_);

  dream_call(dealloc(SFR),
             _dealloc_error_message_);

  dream_call(dealloc(SFH),
             _dealloc_error_message_);

  dream_call(dealloc(MLR),
             _dealloc_error_message_);

  dream_call(dealloc(lookback_time),
             _dealloc_error_message_);

  dream_call(dealloc(cosmic_time),
             _dealloc_error_message_);

  *mass_at_z = log10(Mstar_satellite);

  return _success_;
}



#endif
