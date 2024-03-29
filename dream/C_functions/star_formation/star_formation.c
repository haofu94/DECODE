/** @ file star_formation.c
  *
  * Written by Hao Fu
  *
  * The main goal of this module is to compute the star formation quantities.
  **/

#ifndef star_formation
# define star_formation


#include <stdio.h>
#include <math.h>

#include "../include/dream.h"

#include "../dark_matter/Mhalo_to_Mstar.c"
#include "../dark_matter/centrals.c"
#include "star_formation_satellites.c"
#include "../analyze.c"



/**
  * Compute fraction of mass lost, Moster et al. 2018
  *
  * Input:
  * tau -> timescale [Gyr]
  * Return:
  * f_loss -> fraction of mass lost
**/
double f_loss(double tau){

  return 0.05 * log(1. + tau/1.4e-3);

}



/**
  * Compute accretion contribution from mergers [log(M/Msun)]
  *
  * Input:
  * @ reduced_catalog -> catalogue within desired mass range
  * @ SMHM_params -> parameters of the SMHM relation
  * @ scatter -> scatter [dex]
  * @ z -> redshift array
  * @ len_z -> length of the redshift array
  * @ merged_stellar_mass -> pointer to array of stellar mass from mergers
**/
int compute_accretion_from_mergers(DM_catalogue *reduced_catalog,
                                   //double *SMHM_params,
                                   stellar_mass_halo_mass *SMHM,
                                   SMHM_matrix *smhm_data,
                                   double scatter,
                                   double *z,
                                   int len_z,
                                   double **merged_stellar_mass, double mstar,
                                   cosmological_parameters *cosmo_params) {

  int i, j;

  double SMHM_params[8];

  double merged_stellar_mass_at_z;

  if (reduced_catalog->len_parents == 0){

    for (i=0; i<len_z; i++) {

      *(*merged_stellar_mass+i) = -INFINITY;

    }

    printf("\n--> Zero galaxies in mass bin log(Mstar/Msun) = %lf\n", mstar);

    return _success_;
  }

  for (i=0; i<len_z; i++) {

    *(*merged_stellar_mass+i) = 0.;

    for (j=0; j<reduced_catalog->len_mergers; j++) {

      if ( *(z+i) <= (*(reduced_catalog->z_at_merge+j)) ) {
      //if ( *(z+i) <= (*(reduced_catalog->z_infall+j)) ) {

        dream_call(SMHM_numerical_interp(*(reduced_catalog->mass_mergers+j),
                                            SMHM, smhm_data, *(reduced_catalog->z_infall+j), &merged_stellar_mass_at_z),
                    _SMHM_error_message_);

        *(*merged_stellar_mass+i) += pow(10., merged_stellar_mass_at_z)*0.61;

      }

    }

    //printf("%lf %lf\n", z[i], *(*merged_stellar_mass+i) / ((double)reduced_catalog->len_parents));
    *(*merged_stellar_mass+i) = log10(*(*merged_stellar_mass+i) / ((double)reduced_catalog->len_parents));
    //printf("%lf %lf\n", z[i], *(*merged_stellar_mass+i));
    //printf("%d %lf %lf\n", reduced_catalog->len_parents, z[i], *(*merged_stellar_mass+i));

  }


  // BELOW TEST FRACTION OF ELLIPTICALS

  /*
  double frac;
  int k, N = 100000;
  int ID[N];
  int counter;
  for (i=0; i<N; i++) {ID[i] = 0;}
  frac = 0.; counter = 0;
  int has_major_merger;
  double *Mhalo_track, *Mstar_track;
  double Mhalo;
  dream_call(SMHM_numerical_interp_inverse(mstar, SMHM, smhm_data, *(reduced_catalog->z_infall+j), &Mhalo),
              _SMHM_error_message_);
  dream_call(double_malloc(len_z, &Mhalo_track), _alloc_error_message_);
  dream_call(double_malloc(len_z, &Mstar_track), _alloc_error_message_);
  dream_call(Mass_acc_history_VDB(Mhalo, z, len_z, cosmo_params, &Mhalo_track), _mass_acc_error_message_);
  dream_call(Mhalo_track_to_Mstar_track(Mhalo_track, SMHM, smhm_data, z, len_z, &Mstar_track),
               _Mhalo_track_to_Mstar_track_error_message_);

  for (i=0; i<len_z-1; i++) {

    for (j=0; j<reduced_catalog->len_mergers; j++) {

      has_major_merger = 0;

      if ( *(z+i) <= (*(reduced_catalog->z_at_merge+j)) ) {

        dream_call(SMHM_numerical_interp(*(reduced_catalog->mass_mergers+j), SMHM, smhm_data, *(reduced_catalog->z_infall+j), &merged_stellar_mass_at_z),
                    _SMHM_error_message_);

        mstar = linear_interp(*(reduced_catalog->z_at_merge+j), z, Mstar_track, len_z);
        if (pow(10., merged_stellar_mass_at_z) / pow(10., mstar) > 0.25 ) { //MM

          has_major_merger += 1;
        }
      }

      if (i==0){
        if ( (has_major_merger>0.) && (j==0) ) {
          frac += 1.;
          ID[counter] = reduced_catalog->id_mergers[j];
          counter += 1;
        } else if ( (has_major_merger>0.) && (j>0) ) {
          int id_exist = _False_;
          for (k=0; k<counter; k++){
            if (ID[k] == reduced_catalog->id_mergers[j]) {id_exist = _True_; break;}
          }
          if (id_exist == _False_) {
            frac += 1.;
            ID[counter] = reduced_catalog->id_mergers[j];
            counter += 1;
          }
        }
      }

    }
  }*/

  printf("\n--> Computed mergers for log(Mstar/Msun) = %lf\n", mstar);
  printf(" -> number of centrals = %d\n", reduced_catalog->len_parents);
  /*printf("number of ellipticals = %d\n", (int)frac);

  frac /= (double)reduced_catalog->len_parents;
  printf("fraction of ellipticals = %lf\n", frac);

  dream_call(dealloc(Mhalo_track), _dealloc_error_message_);
  dream_call(dealloc(Mstar_track), _dealloc_error_message_);
  */

  return _success_;
}



/**
  * Computer star formation rate as function of time [Msun/yr]
  *
  * Input:
  * @ star_formation_history -> star formation history [Msun]
  * @ lookback_time -> lookback_time [Gyr]
  * @ len_sfr -> length of the star formation rate array
  * @ len_sfh -> length of the star formation history array
  * @ star_formation_rate -> pointer to the star formation rate array
**/
int compute_star_formation_rate(double *star_formation_history,
                                double *lookback_time,
                                int len_sfh,
                                double **star_formation_rate) {

  int i, j, k, idx;

  for (j=1; j<(len_sfh-2); j++) {

    idx = len_sfh-j-1;

    *(*star_formation_rate+idx) = *(star_formation_history+idx-1);

    for (k=idx+1; k<(len_sfh-1); k++) {

      *(*star_formation_rate+idx) -= *(*star_formation_rate+k) * \
                  ((*(lookback_time+k)) - (*(lookback_time+k-1)))*1.e+9 * \
                      (1. - f_loss((*(lookback_time+k)) - (*(lookback_time+idx-1))));

    }

    *(*star_formation_rate+idx) /= ((*(lookback_time+idx)) - (*(lookback_time+idx-1))) \
                  * 1.e+9 * (1. - f_loss((*(lookback_time+idx+1)) - (*(lookback_time+idx))));

  }

  /*double diff_sum_i;
  for (idx = len_sfh-2; idx>=0; idx--) {
    diff_sum_i = 0.;
    for (i=idx+1; i<len_sfh; i++){
      //diff_sum_i += *(*star_formation_rate+i) * (1. - f_loss(lookback_time[len_sfh-1] - lookback_time[i])) * (lookback_time[len_sfh-1] - lookback_time[i]) * 1e+9;
      diff_sum_i += *(*star_formation_rate+i) * (1. - f_loss(lookback_time[i] - lookback_time[idx])) * (lookback_time[i] - lookback_time[idx]) * 1e+9;
    }
    if (diff_sum_i <= 0.){ diff_sum_i = 1e-10;}
    //diff_sum_i = star_formation_history[idx+1];
    *(*star_formation_rate+idx) = (star_formation_history[idx] - diff_sum_i) / (1. - f_loss(lookback_time[len_sfh-1] - lookback_time[idx])) / (lookback_time[len_sfh-1] - lookback_time[idx]) / 1e+9;
    if (*(*star_formation_rate+idx) <= 0.){ *(*star_formation_rate+idx) = 1e-10;}
    printf("%lf %lf %lf\n", log10(star_formation_history[idx]), log10(diff_sum_i), log10(star_formation_history[idx] - diff_sum_i));
  }*/

  return _success_;
}




/**
  * Compute star formation histories and rates
  *
  * Input:
  * @ SMHM_params -> parameters of the SMHM relation
  * @ scatter -> scatter [dex]
  * @ stellar_masses_for_plot -> stellar masses for which SF is desired [log(M/Msun)]
  * @ len_stellar_masses_for_plot -> length of stellar_masses_for_plot array
  * @ redshifts_for_SFR -> redshifts for which SFR-Mstar relation is resired
  * @ len_redshifts_for_SFR -> length of redshifts_for_SFR array
  * @ data_folder -> path to data folder
  * @ DM_catalog -> dark matter catalogue
  * @ cosmo_params -> cosmological parameters
  * @ cosmo_time -> cosmological time
**/
void compute_star_formation(stellar_mass_halo_mass *SMHM,
                            double *stellar_masses_for_plot,
                            int len_stellar_masses_for_plot,
                            double *redshifts_for_SFR,
                            int len_redshifts_for_SFR,
                            char *data_folder,
                            DM_catalogue *DM_catalog,
                            cosmological_parameters *cosmo_params,
                            cosmological_time *cosmo_time) {

  int i, j, k;
  double SMHM_params[8];
  double m, z_bin, m_for_plot, z_, z_for_SFR;
  double Mhalo;
  FILE *file_pointer;

  SMHM_matrix *smhm_data;
  smhm_data = (SMHM_matrix *)malloc(sizeof(SMHM_matrix));

  SMHM->scatter = 0.; //force scatter to 0

  dream_call(SMHM_read_matrix(SMHM, &smhm_data),
             _SMHM_read_matrix_error_message_);

  double *central_gal_masses;
  double *z, *lookback_time, *Mhalo_track;
  double *stellar_masses_for_SFR;
  int len_z, len_stellar_masses_for_SFR;

  dream_call(double_malloc(DM_catalog->len_parents,
                           &central_gal_masses),
             _alloc_error_message_);

    for (i=0; i<DM_catalog->len_parents; i++){

      dream_call(SMHM_numerical_interp(*(DM_catalog->mass_parents+i),
                                        SMHM, smhm_data,
                                        0.1, &(*(central_gal_masses+i))),
                  _SMHM_error_message_);

    }

  //SMHM->scatter = 0.; //set scatter back to 0

  z_bin = 0.001;
  double z_MAX = 6.;
  len_z = (int)(z_MAX/z_bin);

  dream_call(double_malloc(len_z, &z),
             _alloc_error_message_);

  dream_call(arange(0., z_MAX, z_bin, len_z, &z),
             _arange_error_message);

  dream_call(double_malloc(len_z, &lookback_time),
             _alloc_error_message_);

  for(i=0; i<len_z; i++){
    *(lookback_time+i) = linear_interp(*(z+i), cosmo_time->redshift, cosmo_time->lookback_time, cosmo_time->length);
  }

  double stellar_mass_MIN = 10.; //default 6.
  double stellar_mass_MAX = 12.2; //default 12.5
  len_stellar_masses_for_SFR = (int)((stellar_mass_MAX - stellar_mass_MIN) / 0.1);

  dream_call(double_malloc(len_stellar_masses_for_SFR,
                            &stellar_masses_for_SFR),
              _alloc_error_message_);

  dream_call(arange(stellar_mass_MIN, stellar_mass_MAX, 0.1,
                    len_stellar_masses_for_SFR,
                    &stellar_masses_for_SFR),
              _arange_error_message);

  dream_call(double_malloc(len_z, &Mhalo_track),
             _alloc_error_message_);


  int **idx_cen = malloc(len_stellar_masses_for_SFR * sizeof(int *));
  int **idx_sat = malloc(len_stellar_masses_for_SFR * sizeof(int *));

  printf("\nInitializing star formation calculation...\n");

  for (j=0; j<len_stellar_masses_for_SFR; j++){

    dream_call(int_malloc(DM_catalog->len_parents, &(*(idx_cen+j))),
               _alloc_error_message_);

    dream_call(int_malloc(DM_catalog->len_mergers, &(*(idx_sat+j))),
               _alloc_error_message_);

    m = *(stellar_masses_for_SFR+j);

    dream_call(get_mask_in_range(central_gal_masses, m, m+0.1,
                                 DM_catalog->len_parents, &(*(idx_cen+j))),
               _get_mask_error_message_);

    int *id_halo_idx_cen;
    int len_id_halo_idx_cen = 0;

    dream_call(int_malloc(len_id_halo_idx_cen,
                          &id_halo_idx_cen),
               _alloc_error_message_);

    for (i=0; i<DM_catalog->len_parents; i++){

      if (*(*(idx_cen+j)+i) == _True_){

        dream_call(int_append(len_id_halo_idx_cen,
                              &id_halo_idx_cen,
                              *(DM_catalog->id_parents+i)),
                   _append_error_message_);

        len_id_halo_idx_cen += 1;

      }

    }

    dream_call(array1_elements_in_array2(DM_catalog->id_mergers,
                                         id_halo_idx_cen,
                                         DM_catalog->len_mergers,
                                         len_id_halo_idx_cen,
                                         &(*(idx_sat+j))),
               _get_mask_error_message_);

    dream_call(dealloc(id_halo_idx_cen),
               _dealloc_error_message_);
  }


  double **Mstar_track = malloc(len_stellar_masses_for_SFR * sizeof(double *));
  double **merged_stellar_mass = malloc(len_stellar_masses_for_SFR * sizeof(double *));
  double **stellar_mass = malloc(len_stellar_masses_for_SFR * sizeof(double *));
  double **star_formation_rate = malloc(len_stellar_masses_for_SFR * sizeof(double *));

  printf("\nLooping on stellar mass range...\n");

  for (i=0; i<len_stellar_masses_for_SFR; i++){

    m = *(stellar_masses_for_SFR+i);

    int *id_halo_catalog_idx_cen;
    double *halo_catalog_idx_cen;
    int len_halo_catalog_idx_cen = 0;

    dream_call(int_malloc(len_halo_catalog_idx_cen,
                           &id_halo_catalog_idx_cen),
               _alloc_error_message_);

    dream_call(double_malloc(len_halo_catalog_idx_cen,
                             &halo_catalog_idx_cen),
               _alloc_error_message_);

    for (j=0; j<DM_catalog->len_parents; j++){

      if (*(*(idx_cen+i)+j) == _True_){

        dream_call(int_append(len_halo_catalog_idx_cen,
                              &id_halo_catalog_idx_cen,
                              *(DM_catalog->id_parents+j)),
                   _append_error_message_);

        dream_call(double_append(len_halo_catalog_idx_cen,
                                 &halo_catalog_idx_cen,
                                 *(DM_catalog->mass_parents+j)),
                   _append_error_message_);

        len_halo_catalog_idx_cen += 1;

      }

    }

    Mhalo = mean(halo_catalog_idx_cen, len_halo_catalog_idx_cen); //log(M/Msun)

    dream_call(Mass_acc_history_VDB(Mhalo, z, len_z, cosmo_params, &Mhalo_track),
               _mass_acc_error_message_);

    dream_call(double_calloc(len_z, &(*(Mstar_track+i))),
               _alloc_error_message_);

    //SMHM->scatter = 0.;
    dream_call(Mhalo_track_to_Mstar_track(Mhalo_track,
                                          SMHM, smhm_data,
                                          z, len_z, &(*(Mstar_track+i))),
               _Mhalo_track_to_Mstar_track_error_message_);
               //int w; for (w=0; w<len_z; w++){printf("%lf %lf\n", Mhalo_track[w], Mstar_track[i][w]);}
               //SMHM->scatter = 0.11;

    int *id_mergers_idx_sat;
    double *mergers_array_idx_sat;
    double *z_infall_idx_sat;
    double *z_at_merge_idx_sat;
    int len_mergers_array_idx_sat = 0;

    dream_call(int_malloc(len_mergers_array_idx_sat,
                          &id_mergers_idx_sat),
               _alloc_error_message_);

    dream_call(double_malloc(len_mergers_array_idx_sat,
                             &mergers_array_idx_sat),
               _alloc_error_message_);

    dream_call(double_malloc(len_mergers_array_idx_sat,
                             &z_infall_idx_sat),
               _alloc_error_message_);

    dream_call(double_malloc(len_mergers_array_idx_sat,
                             &z_at_merge_idx_sat),
               _alloc_error_message_);

    for (j=0; j<DM_catalog->len_mergers; j++){

      if (*(*(idx_sat+i)+j) == _True_){

        dream_call(int_append(len_mergers_array_idx_sat,
                              &id_mergers_idx_sat,
                              *(DM_catalog->id_mergers+j)),
                   _append_error_message_);

        dream_call(double_append(len_mergers_array_idx_sat,
                                 &mergers_array_idx_sat,
                                 *(DM_catalog->mass_mergers+j)),
                   _append_error_message_);

        dream_call(double_append(len_mergers_array_idx_sat,
                                 &z_infall_idx_sat,
                                 *(DM_catalog->z_infall+j)),
                   _append_error_message_);

        dream_call(double_append(len_mergers_array_idx_sat,
                                 &z_at_merge_idx_sat,
                                 *(DM_catalog->z_at_merge+j)),
                   _append_error_message_);

        len_mergers_array_idx_sat += 1;

      }

    }

    dream_call(double_malloc(len_z, &(*(merged_stellar_mass+i))),
               _alloc_error_message_);

    DM_catalogue *reduced_catalog = malloc(sizeof(DM_catalogue));
    reduced_catalog->id_parents = id_halo_catalog_idx_cen;
    reduced_catalog->mass_parents = halo_catalog_idx_cen;
    reduced_catalog->len_parents = len_halo_catalog_idx_cen;
    reduced_catalog->len_mergers = len_mergers_array_idx_sat;
    reduced_catalog->id_mergers = id_mergers_idx_sat;
    reduced_catalog->mass_mergers = mergers_array_idx_sat;
    reduced_catalog->z_infall = z_infall_idx_sat;
    reduced_catalog->z_at_merge = z_at_merge_idx_sat;

    dream_call(compute_accretion_from_mergers(reduced_catalog, SMHM, smhm_data,
                                              SMHM->scatter, z, len_z,
                                              &(*(merged_stellar_mass+i)), m, cosmo_params),
               _merged_stellar_mass_error_message_);

    //printf(" -> Computed mergers for log(Mstar/Msun) = %lf\n\n", m);

    dream_call(double_calloc(len_z, &(*(stellar_mass+i))),
               _alloc_error_message_);

    for(j=0; j<len_z; j++){
      if ( (pow(10., *(*(Mstar_track+i)+j)) - pow(10., *(*(merged_stellar_mass+i)+j))) > 0. ) {
        *(*(stellar_mass+i)+j) = pow(10., *(*(Mstar_track+i)+j)) - pow(10., *(*(merged_stellar_mass+i)+j));
      }
      //printf("%lf %lf %lf %lf\n", z[j], *(*(Mstar_track+i)+j), *(*(merged_stellar_mass+i)+j), log10(*(*(stellar_mass+i)+j)));
      //printf("%lf %lf %lf %lf\n", z[j], Mstar_track[i][j], merged_stellar_mass[i][j], log10(stellar_mass[i][j]));
    }

    dream_call(double_calloc(len_z, &(*(star_formation_rate+i))),
               _alloc_error_message_);

    dream_call(compute_star_formation_rate(*(stellar_mass+i),
                                           lookback_time, len_z,
                                           &(*(star_formation_rate+i))),
               _compute_SFR_error_message_);

    for(k=0; k<len_stellar_masses_for_plot; k++){

      m_for_plot = *(stellar_masses_for_plot+k);

      if (fabs(m - m_for_plot) < 0.01){

        char m_str[16], file_name[256];

        sprintf(m_str, "%.2lf", m_for_plot);

        snprintf(file_name, sizeof file_name, "%s%s%s%s", data_folder, \
                                    "SFR_z_relation_for_Mstar_", m_str, ".txt");

        file_pointer = fopen(file_name, "w");

        fprintf(file_pointer, "# 1) first column: redshift\n# 2) second column: lookback time [Gyr]\n# 3) third column: SFR [Msun/yr]\n\n");

        for (j=0; j<len_z-1; j++){
          fprintf(file_pointer, "%.10le %.10le %.10le\n", \
                    *(z+j), *(lookback_time+j), *(*(star_formation_rate+i)+j));
        }

        fclose(file_pointer);
      }

    }

    dream_call(dealloc(id_halo_catalog_idx_cen),
               _dealloc_error_message_);

    dream_call(dealloc(halo_catalog_idx_cen),
               _dealloc_error_message_);

    dream_call(dealloc(id_mergers_idx_sat),
               _dealloc_error_message_);

    dream_call(dealloc(mergers_array_idx_sat),
               _dealloc_error_message_);

    dream_call(dealloc(z_infall_idx_sat),
               _dealloc_error_message_);

    dream_call(dealloc(z_at_merge_idx_sat),
               _dealloc_error_message_);

    dream_call(dealloc(reduced_catalog),
               _dealloc_error_message_);

  }

  double *DM_linspace, *SM_linspace;
  int len_linspace = 1000;

  dream_call(double_malloc(len_linspace, &DM_linspace),
             _alloc_error_message_);

  dream_call(double_malloc(len_linspace, &SM_linspace),
             _alloc_error_message_);

  dream_call(linear_space(10., 16., len_linspace, &DM_linspace),
             _linear_space_error_message);

  for(j=0; j<len_z; j++){

    z_ = *(z+j);

    for (k=0; k<len_redshifts_for_SFR; k++){

      z_for_SFR = *(redshifts_for_SFR+k);

      if (fabs(z_ - z_for_SFR) < 1.e-4){

        double *stellar_masses_to_plot;

        if (z_>0.) {

          double *halo_masses_for_SFR, *z_for_evolution, *this_halo_track;

          dream_call(double_malloc(len_stellar_masses_for_SFR,
                                   &stellar_masses_to_plot),
                     _alloc_error_message_);

          dream_call(double_malloc(len_stellar_masses_for_SFR,
                                   &halo_masses_for_SFR),
                     _alloc_error_message_);

          int len_z_for_evolution = (int)((z_+0.1)/0.05);

          dream_call(double_malloc(len_z_for_evolution,
                                   &z_for_evolution),
                     _alloc_error_message_);

          dream_call(double_malloc(len_z_for_evolution,
                                   &this_halo_track),
                     _alloc_error_message_);

          dream_call(arange(0., z_+0.1, 0.05,
                            len_z_for_evolution,
                            &z_for_evolution),
                     _arange_error_message);

          for (i=0; i<len_linspace; i++){

            dream_call(SMHM_numerical_interp(*(DM_linspace+i),
                                             SMHM, smhm_data,
                                             z_, &(*(SM_linspace+i))),
                       _SMHM_error_message_);

          }

          for (i=0; i<len_stellar_masses_for_SFR; i++){

            *(halo_masses_for_SFR+i) = linear_interp(*(stellar_masses_for_SFR+i), \
                                        SM_linspace, DM_linspace, len_linspace);

            dream_call(Mass_acc_history_VDB(*(halo_masses_for_SFR+i),
                                            z_for_evolution, len_z_for_evolution,
                                            cosmo_params, &this_halo_track),
                       _mass_acc_error_message_);

            dream_call(SMHM_numerical_interp(*(this_halo_track+len_z_for_evolution-1),
                                             SMHM, smhm_data,
                                             z_, &(*(stellar_masses_to_plot+i))),
                       _SMHM_error_message_);

          }

          dream_call(dealloc(halo_masses_for_SFR),
                     _dealloc_error_message_);

          dream_call(dealloc(z_for_evolution),
                     _dealloc_error_message_);

          dream_call(dealloc(this_halo_track),
                     _dealloc_error_message_);

        } else {

          stellar_masses_to_plot = stellar_masses_for_SFR;

        }

        char z_str[16], file_name[256];

        sprintf(z_str, "%.2lf", z_for_SFR);

        snprintf(file_name, sizeof file_name, "%s%s%s%s", data_folder, \
                                    "SFR_Mstar_relation_for_z_", z_str, ".txt");

        file_pointer = fopen(file_name, "w");

        fprintf(file_pointer, "# 1) first column: log(Mstar/Msun)\n# 2) second column: SFR [Msun/yr]\n\n");

        for (i=0; i<len_stellar_masses_for_SFR; i++){

          fprintf(file_pointer, "%.10le %.10le\n", *(stellar_masses_to_plot+i), *(*(star_formation_rate+i)+j));

        }

        fclose(file_pointer);

        if (z_>0.) {
          dream_call(dealloc(stellar_masses_to_plot),
                     _dealloc_error_message_);
        }

      }
    }
  }



  for (k=0; k<len_stellar_masses_for_plot; k++){

    for (i=0; i<len_stellar_masses_for_SFR; i++){

      if (fabs(stellar_masses_for_SFR[i] - stellar_masses_for_plot[k]) < 1.e-4){

        m_for_plot = *(stellar_masses_for_SFR+i);

        char m_str[16], file_name[256];

        sprintf(m_str, "%.2lf", m_for_plot);

        snprintf(file_name, sizeof file_name, "%s%s%s%s", data_folder, \
                                    "SFH_for_Mstar_", m_str, ".txt");

        file_pointer = fopen(file_name, "w");

        fprintf(file_pointer, "# 1) first column: redshift\n# 2) second column: total growth history [log(M/Msun)]\n# 3) third column: accretion from mergers [log(M/Msun)]\n# 4) fourth column: star formation history [log(M/Msun)]\n\n");

        for (j=0; j<len_z; j++){

          fprintf(file_pointer, "%.10le %.10le %.10le %.10le\n", *(z+j), *(*(Mstar_track+i)+j), *(*(merged_stellar_mass+i)+j), log10(*(*(stellar_mass+i)+j)) );
        }

        fclose(file_pointer);
      }
    }
  }



  dream_call(dealloc(central_gal_masses),
             _dealloc_error_message_);

  dream_call(dealloc(z),
             _dealloc_error_message_);

  dream_call(dealloc(lookback_time),
             _dealloc_error_message_);

  dream_call(dealloc(stellar_masses_for_SFR),
             _dealloc_error_message_);

  for (j=0; j<len_stellar_masses_for_SFR; j++){
    dream_call(dealloc(*(idx_cen+j)),
               _dealloc_error_message_);
    dream_call(dealloc(*(idx_sat+j)),
               _dealloc_error_message_);
    dream_call(dealloc(*(Mstar_track+j)),
               _dealloc_error_message_);
    dream_call(dealloc(*(merged_stellar_mass+j)),
               _dealloc_error_message_);
    dream_call(dealloc(*(stellar_mass+j)),
               _dealloc_error_message_);
    dream_call(dealloc(*(star_formation_rate+j)),
               _dealloc_error_message_);
  }
  dream_call(dealloc(idx_cen),
             _dealloc_error_message_);
  dream_call(dealloc(idx_sat),
             _dealloc_error_message_);
  dream_call(dealloc(Mstar_track),
             _dealloc_error_message_);
  dream_call(dealloc(merged_stellar_mass),
             _dealloc_error_message_);
  dream_call(dealloc(stellar_mass),
             _dealloc_error_message_);
  dream_call(dealloc(star_formation_rate),
             _dealloc_error_message_);

  dream_call(dealloc(Mhalo_track),
             _dealloc_error_message_);

  dream_call(dealloc(DM_linspace),
             _dealloc_error_message_);

  dream_call(dealloc(SM_linspace),
             _dealloc_error_message_);

  dream_call(dealloc(smhm_data->redshift),
             _dealloc_error_message_);
  dream_call(dealloc(smhm_data->matrix),
             _dealloc_error_message_);
  dream_call(dealloc(smhm_data->Mstar),
             _dealloc_error_message_);
  dream_call(dealloc(smhm_data),
             _dealloc_error_message_);

  printf("\nStar formation successfully computed.\n");

  return;
}
























//*************************** RUBBISH ******************************************


double f_dot_Leitner(double tau){

  //
  //  * Fraction of mass lost
  //  * Eq. (1) from
  //  * Leinet et al. 2011, Astrophysical Journal, 734:48
  //

  return 0.046 / (tau + 2.76e-4);
}



double f_dot_Moster(double tau){

  //
  //  * Fraction of mass lost
  //  * Moster et al. 2018
  //

  return 0.05 / (tau + 1.4e-3);
}


double f_Moster(double tau){
  return 0.05*log(1. + tau/1.4e-3);
}


double Mdot_Leitner(double Mstar,
                    double t,
                    double sfr,
                    double z,
                    double *t_array,
                    int idx,
                    int len){

  //
  //  * Star formation rate and mass loss rate
  //  * Eqs. (2) and (3) from
  //  * Leinet et al. 2011, Astrophysical Journal, 734:48
  //

  double A = 2.8;  //[Msun/yr]
  double a = 3.4;
  double b = -0.25;
  double z_break = 2.;
  double Mdot, psi, gmlr;
  int i;


  if (z < z_break) {
    gmlr = 0.;
    psi = A * pow(z+1., a) * pow(pow(10., Mstar) / pow(10., 10.75), 1.+b);

    for (i=idx; i<len; i++){
      gmlr += sfr * (f_dot_Leitner((*(t_array+i) - t)/1e+9) + f_dot_Leitner((*(t_array+i-1) - t)/1e+9) ) * 0.5 * (*(t_array+i) - (*(t_array+i-1))) / 1e+9;
    }
    Mdot = psi - gmlr;
  } else {
    gmlr = 0.;
    psi = A * pow(pow(10., Mstar) / pow(10., 10.75), 1.+b);
    for (i=idx; i<len; i++){
      gmlr += sfr * (f_dot_Leitner((*(t_array+i) - t)/1e+9) + f_dot_Leitner((*(t_array+i-1) - t)/1e+9) ) * 0.5 * (*(t_array+i) - (*(t_array+i-1))) / 1e+9;
    }
    Mdot = psi - gmlr;
  }

  return Mdot;
}




double Mdot_Moster(double Mstar,
                   double t,
                   double *sfr,
                   double z,
                   double *t_array,
                   int idx,
                   int len){

  //
  //  * Star formation rate and mass loss rate
  //  * Moster et al. 2018
  //

  int i;
  double Mdot, psi, gmlr;

  double M = pow(10., 10.65 + 0.33*z - 0.08*z*z);
  double N = pow(10., 0.69 + 0.71*z - 0.088*z*z);
  double a = 1. - 0.022*z + 0.009*z*z;
  double b = 1.8 - 1.*z - 0.1*z*z;

  gmlr = 0.;
  psi = 2*N / (pow(pow(10., Mstar) / M, -a) + pow(pow(10., Mstar)/M, b));

  for (i=idx; i<len; i++){
    gmlr += sfr[i] * (f_dot_Moster((*(t_array+i) - t)/1e+9) + f_dot_Moster((*(t_array+i-1) - t)/1e+9) ) * 0.5 * (*(t_array+i) - (*(t_array+i-1))) / 1e+9;
  }

  Mdot = psi - gmlr;

  return Mdot;
}



#endif
