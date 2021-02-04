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
}



/**
  * Compute the quenching timescale [Gyr]
  *
  * Input:
  * @ Mstar -> stellar mass [log(M/Msun)]
  * @ tau_quenching -> pointer to the tau_quenching variable
**/
int compute_quenching_timescale(double Mstar,
                                double *tau_quenching){

  *tau_quenching = -0.5 * Mstar + 5.7;

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
  * @ Mstar -> stellar mass [log(M/Msun)]
  * @ z -> redshift
**/
int get_star_formation_rate_satellites(double Mstar, double z, double *SFR) {

  double s0 = 0.6 + 1.22*z - 0.2*z*z;
  double M0 = pow(10., 10.3 + 0.753*z - 0.15*z*z);
  double alpha = 1.3 - 0.1*z;
  Mstar = pow(10., Mstar);

  *SFR = pow(10., s0 - log10(1. + pow(Mstar/M0, -alpha)) ) ;//+ get_random_gaussian(0., 0.3, -5., 5.);

  return _success_;
}



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
                               int include_stripping,
                               int include_quenching,
                               cosmological_time *cosmo_time,
                               double *mass_at_z){

  int i, j, idx;

  double tau_strip;

  dream_call(compute_stripping_timescale(progenitor_mass,
                                         subhalo_mass,
                                         &tau_strip),
             _stripping_timescale_error_message_);

  double eta_strip = 0.4;

  double t_infall = linear_interp(z_infall, cosmo_time->redshift, cosmo_time->age, cosmo_time-> length)*1.e9;
  double tau_q = 2.e9;
  double tau_f;

  double f2, f1;
  double dt = 1.e4;

  double zbin = 0.1;

  int len = 20;
  double *z_range, *lookback_time, *Mstar, *SFR, *SFH, *MLR;

  dream_call(double_malloc(len,
                           &z_range),
             _alloc_error_message_);

  dream_call(double_malloc(len,
                           &lookback_time),
             _alloc_error_message_);

  dream_call(double_calloc(len,
                           &Mstar),
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

  dream_call(linear_space(z,
                          z_infall,
                          len,
                          &z_range),
             _linear_space_error_message);

  for (i=0; i<len; i++){

    *(lookback_time+i) = linear_interp(0., cosmo_time->redshift, cosmo_time->age, cosmo_time-> length) - \
                linear_interp(*(z_range+i), cosmo_time->redshift, cosmo_time->age, cosmo_time-> length);

    *(lookback_time+i) *= 1.e9; //[yr]
  }

  dream_call(SMHM_numerical_interp(subhalo_mass, SMHM, smhm_data, z_infall, &(*(Mstar+len-1))),
             _SMHM_error_message_);

  for (i=len-2; i>=0; i--){

    dream_call(get_star_formation_rate_satellites(*(Mstar+i+1), *(z_range+i+1), &(*(SFR+i))),
               _compute_SFR_error_message_);

    //GAS LOSS
    for (j=len-2; j>=i; j--){

      *(MLR+i) += *(SFR+j) * (lookback_time[j+1]-lookback_time[j]) * fdot(lookback_time[j+1]-lookback_time[j]);

    }

    *(SFR+i) -= *(MLR+i);

    if ((include_quenching == _True_) && (*(lookback_time+i) < tau_q)){

      if (*(Mstar+i+1) < 9.) {

        dream_call(compute_quenching_timescale(*(Mstar+i),
                                               &tau_f),
                   _quenching_timescale_error_message_);

        tau_f *= 1.e9;

        *(SFR+i) = *(SFR+i) * exp(-(tau_q - (*(lookback_time+i))) / tau_f);
      }
    }

    *(SFH+i) = integrate_trapz(SFR, lookback_time, len); //[Msun]

    *(Mstar+i) = log10(pow(10., *(Mstar+len-1)) + (*(SFH+i)));

    if (include_stripping == _True_){
      if (pow(1. - eta_strip, tau_strip)>1.e-4){
        //printf("%lf %lf %.5e\n", 1.-eta_strip, tau_strip, pow(1. - eta_strip, tau_strip));
        *(Mstar+i) = log10(pow(10., *(Mstar+i)) * pow(1. - eta_strip, tau_strip));
      }
      //printf("%lf %lf\n", tau_strip, pow(10., *(Mstar+i)) * pow(1. - eta_strip, tau_strip));
    }

  }

  dream_call(dealloc(z_range),
             _dealloc_error_message_);

  dream_call(dealloc(Mstar),
             _dealloc_error_message_);

  dream_call(dealloc(SFR),
             _dealloc_error_message_);

  dream_call(dealloc(SFH),
             _dealloc_error_message_);

  dream_call(dealloc(MLR),
             _dealloc_error_message_);

  dream_call(dealloc(lookback_time),
             _dealloc_error_message_);

  *mass_at_z = *Mstar;

  return _success_;
}






//--------------------------------------------------------
// For using numerical SMHM from Python (IGNORE)
//--------------------------------------------------------

double get_evolved_satellite_mass_for_python(double *SMHM_params,
                                  double progenitor_mass,
                                  double subhalo_mass,
                                  double z,
                                  double z_infall,
                                  int include_stripping,
                                  int include_quenching,
                                  cosmological_time *cosmo_time,
                                  double *mhalo,
                                  double *mstar,
                                  int size){

  int i, j, idx;

  double tau_strip;

  dream_call(compute_stripping_timescale(progenitor_mass,
                                          subhalo_mass,
                                          &tau_strip),
              _stripping_timescale_error_message_);
  //printf("%lf\n", tau_strip);
  double eta_strip = 0.4;

  double t_infall = linear_interp(z_infall, cosmo_time->redshift, cosmo_time->age, cosmo_time-> length)*1.e9;
  double tau_q = 2.e9;
  double tau_f;

  double f2, f1;
  double dt = 1.e4;

  double zbin = 0.1;

  //int len = (int)((z_infall-z)/zbin);
  int len = 20;
  double *z_range, *lookback_time, *Mstar, *SFR, *SFH, *MLR;

  dream_call(double_malloc(len,
                            &z_range),
              _alloc_error_message_);

  dream_call(double_malloc(len,
                            &lookback_time),
              _alloc_error_message_);

  dream_call(double_calloc(len,
                            &Mstar),
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

  //z_range = arange(z, z_infall, zbin, len, z_range);
  dream_call(linear_space(z,
                           z_infall,
                           len,
                           &z_range),
              _linear_space_error_message);

  for (i=0; i<len; i++){

    *(lookback_time+i) = linear_interp(0., cosmo_time->redshift, cosmo_time->age, cosmo_time-> length) - \
                linear_interp(*(z_range+i), cosmo_time->redshift, cosmo_time->age, cosmo_time-> length);

    *(lookback_time+i) *= 1.e9; //[yr]
  }

  //*(Mstar+len-1) = SMHM_Moster(subhalo_mass, SMHM_params, z_infall); //[log(M/Msun)]
  *(Mstar+len-1) = linear_interp(subhalo_mass, mhalo, mstar, size); //[log(M/Msun)]
  for (i=len-2; i>=0; i--){

    //*(SFR+i) = get_star_formation_rate_satellites(*(Mstar+i+1), *(z_range+i+1));
    dream_call(get_star_formation_rate_satellites(*(Mstar+i+1), *(z_range+i+1), &(*(SFR+i))),
                _compute_SFR_error_message_);

    //GAS LOSS
    for (j=len-2; j>=i; j--){

      *(MLR+i) += *(SFR+j) * (lookback_time[j+1]-lookback_time[j]) * fdot(lookback_time[j+1]-lookback_time[j]);

    }

    *(SFR+i) -= *(MLR+i);

    if ((include_quenching == _True_) && (*(lookback_time+i) < tau_q)){

      if (*(Mstar+i+1) < 9.) {

        tau_f = -0.5 * (*(Mstar+i)) + 5.7; tau_f *= 1.e9;

        *(SFR+i) = *(SFR+i) * exp(-(tau_q - (*(lookback_time+i))) / tau_f);
      }
    }

    *(SFH+i) = integrate_trapz(SFR, lookback_time, len); //[Msun]

    *(Mstar+i) = log10(pow(10., *(Mstar+len-1)) + (*(SFH+i)));

    if (include_stripping == _True_){
      printf("%lf %lf %lf\n", 1.-eta_strip, tau_strip, pow(1. - eta_strip, tau_strip));
      if (pow(1. - eta_strip, tau_strip)>0.){
        *(Mstar+i) = log10(pow(10., *(Mstar+i)) * pow(1. - eta_strip, tau_strip));
      }
      //printf("%lf %lf\n", tau_strip, pow(10., *(Mstar+i)) * pow(1. - eta_strip, tau_strip));
    }

    //printf("%.3e\n", *(Mstar+i));

  }

  double m_at_z;

  if (len==0){
    m_at_z = linear_interp(subhalo_mass, mhalo, mstar, size);//SMHM_Moster(subhalo_mass, SMHM_params, z_infall);
  } else {
    m_at_z = *Mstar;
  }
  if (m_at_z < 1.){
  //if (SMHM_Moster(subhalo_mass, SMHM_params, z_infall)< 10.){
    //printf("%d\n", len);
    //printf("Evolved mass: %lf %lf\n", linear_interp(subhalo_mass, mhalo, mstar, size), m_at_z);
  }

  dream_call(dealloc(z_range),
              _dealloc_error_message_);

  dream_call(dealloc(Mstar),
              _dealloc_error_message_);

  dream_call(dealloc(SFR),
              _dealloc_error_message_);

  dream_call(dealloc(SFH),
              _dealloc_error_message_);

  dream_call(dealloc(MLR),
              _dealloc_error_message_);

  dream_call(dealloc(lookback_time),
              _dealloc_error_message_);

  return m_at_z;
}




#endif
