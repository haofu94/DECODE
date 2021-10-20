/** @ file ellipticals.c
  *
  * Written by Hao Fu
  *
  * The main goal of this module is to compute the fraction of elliptical galaxies.
  **/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/dream.h"
#include "../include/cosmological_model.h"

#include "../numerical/num_c.c"
#include "../numerical/allocate.c"
#include "../dark_matter/centrals.c"
#include "../star_formation/star_formation_satellites.c"



int compute_fraction_ellipticals(char *file_name,
                                 char *output_folder,
                                 DM_catalogue *halo_catalogue,
                                 double observed_redshift,
                                 stellar_mass_halo_mass *SMHM,
                                 SMHM_matrix *smhm_data,
                                 int include_sat_SF,
                                 int include_mass_loss,
                                 int include_quenching,
                                 int include_stripping,
                                 cosmological_parameters *cosmo_params,
                                 cosmological_time *cosmo_time) {

  int i, j, k;

  double evolved_satellite_mass;
  double Mstar_at_z_merge;
  int num_this_bin_centrals;
  int num_this_bin_ellipticals;

  FILE *file_pointer;

  char filename[28] = "data/output_ellipticals.txt";
  const size_t len1 = strlen(output_folder);
  const size_t len2 = strlen(filename);
  char *output_filename;
  dream_call(char_malloc(len1 + len2 + 1,
                          &output_filename),
             _alloc_error_message_);
  strcpy(output_filename, output_folder);
  strcpy(output_filename + len1, filename);

  double major_merger_mass_frac = 0.25;
  double frac;

  double z_range_MIN = 0.;
  double z_range_MAX = 4.;
  double z_bin = 0.1;

  int len_z_range = (int) ( (z_range_MAX - z_range_MIN) / z_bin );

  double *z_range;

  dream_call(double_malloc(len_z_range,
                           &z_range),
             _alloc_error_message_);

  dream_call(arange(z_range_MIN, z_range_MAX, z_bin, len_z_range, &z_range), _arange_error_message);

  double Mstar_MIN = 10.;
  double Mstar_MAX = 12.5;
  double Mstar_bin = 0.1;

  int len_Mstar_range = (int) ( (Mstar_MAX - Mstar_MIN) / Mstar_bin );

  double *Mstar_range;

  dream_call(double_malloc(len_Mstar_range,
                           &Mstar_range),
             _alloc_error_message_);

  dream_call(arange(Mstar_MIN, Mstar_MAX, Mstar_bin, len_Mstar_range, &Mstar_range),
             _arange_error_message);


  printf("\n\nComputing fraction of elliptical galaxies...\n");

  double *mass_centrals;

  dream_call(double_malloc(halo_catalogue->len_parents,
                           &mass_centrals),
             _alloc_error_message_);

  for (i=0; i<halo_catalogue->len_parents; i++){

    dream_call(SMHM_numerical_interp(*(halo_catalogue->mass_parents+i), SMHM, smhm_data, 0., &(*(mass_centrals+i))),
               _SMHM_error_message_);

  }

  /*file_pointer = fopen("galaxies.txt", "w");
  for (i=0; i<halo_catalogue->len_parents; i++){
    fprintf(file_pointer, "%lf\n", mass_centrals[i]);
  }
  fclose(file_pointer);*/

  file_pointer = fopen(output_filename, "w");

  for (i=0; i<len_Mstar_range; i++){

    num_this_bin_centrals = 0;
    num_this_bin_ellipticals = 0;

    for (j=0; j<halo_catalogue->len_parents; j++){

      //printf("%lf %lf\n", halo_catalogue->mass_parents[j], mass_centrals[j]);

      if ( (Mstar_range[i]-Mstar_bin/2. <= mass_centrals[j]) && (mass_centrals[j] < Mstar_range[i]+Mstar_bin/2.) ){

        num_this_bin_centrals += 1;

        // MHALO MEAN GROWTH

        double Mhalo;

        dream_call(SMHM_numerical_interp_inverse(Mstar_range[i], SMHM, smhm_data, 0., &Mhalo),
                   _SMHM_error_message_);

        double *Mhalo_mean_growth;

        dream_call(double_malloc(len_z_range,
                                 &Mhalo_mean_growth),
                   _alloc_error_message_);

        dream_call(Mass_acc_history_VDB(Mhalo, z_range, len_z_range, cosmo_params, &Mhalo_mean_growth),
                   _mass_acc_error_message_);

        // MSTAR MEAN GROWTH

        double *Mstar_mean_growth;

        dream_call(double_malloc(len_z_range,
                                 &Mstar_mean_growth),
                   _alloc_error_message_);

        dream_call(Mhalo_track_to_Mstar_track(Mhalo_mean_growth, SMHM, smhm_data, z_range, len_z_range, &Mstar_mean_growth),
                   _Mhalo_track_to_Mstar_track_error_message_);

        /*printf("\n%lf\n", Mstar_mean_growth[0]);
        for (k=0; k<len_z_range; k++){
          printf("%lf %lf\n", z_range[k], Mstar_mean_growth[k]);
        }*/

        for (k=0; k<halo_catalogue->len_mergers; k++){

         if ( (halo_catalogue->id_mergers[k] == halo_catalogue->id_parents[j]) && (halo_catalogue->z_at_merge[k] >= observed_redshift) && (halo_catalogue->z_at_merge[k] <= z_range_MAX) && (halo_catalogue->order_mergers[k] == 1) ) {

           //printf("%d %d\n", halo_catalogue->id_mergers[k], halo_catalogue->id_parents[j]);

           Mstar_at_z_merge = linear_interp(halo_catalogue->z_at_merge[k], z_range, Mstar_mean_growth, len_z_range);

           /*dream_call(get_evolved_satellite_mass(SMHM, smhm_data,
                                                 halo_catalogue->mass_parents[j],
                                                 *(halo_catalogue->mass_mergers+k),
                                                 halo_catalogue->z_at_merge[k],
                                                 *(halo_catalogue->z_infall+k),
                                                 include_sat_SF,
                                                 include_mass_loss,
                                                 include_stripping,
                                                 include_quenching,
                                                 cosmo_params,
                                                 cosmo_time,
                                                 &evolved_satellite_mass),
                      _mstar_at_z_error_message_);*/

            //printf("%lf %lf\n", halo_catalogue->mass_mergers[k], halo_catalogue->mass_parents[j]);

            dream_call(SMHM_numerical_interp(halo_catalogue->mass_mergers[k], SMHM, smhm_data, halo_catalogue->z_infall[k], &evolved_satellite_mass),
                       _SMHM_error_message_);

            //printf("%lf %lf %lf %lf\n", evolved_satellite_mass, Mstar_at_z_merge, pow(10., evolved_satellite_mass) / pow(10., Mstar_at_z_merge), Mstar_range[i]);

            if ( (pow(10., evolved_satellite_mass) / pow(10., Mstar_at_z_merge) >= major_merger_mass_frac ) &&
                (pow(10., evolved_satellite_mass) / pow(10., Mstar_at_z_merge) < 1. ) ) {

              //printf("%lf %lf %lf %lf\n", evolved_satellite_mass, Mstar_at_z_merge, pow(10., evolved_satellite_mass) / pow(10., Mstar_at_z_merge), Mstar_range[i]);

              num_this_bin_ellipticals += 1;

              break;

            }
          }

        }

        dream_call(dealloc(Mhalo_mean_growth),
                   _dealloc_error_message_);

        dream_call(dealloc(Mstar_mean_growth),
                   _dealloc_error_message_);

      }

    }

    frac = (double)num_this_bin_ellipticals / (double)num_this_bin_centrals;

    fprintf(file_pointer, "%lf %lf %d %d\n", Mstar_range[i], frac, num_this_bin_ellipticals, num_this_bin_centrals);

    printf("Progress: %.2f %%\r", (double)(i+1) / (double)len_Mstar_range * 100.);
    fflush(stdout);

  }

  fclose(file_pointer);


  dream_call(dealloc(z_range),
             _dealloc_error_message_);

  dream_call(dealloc(Mstar_range),
             _dealloc_error_message_);

  dream_call(dealloc(mass_centrals),
             _dealloc_error_message_);

  dream_call(dealloc(output_filename),
             _dealloc_error_message_);

  return _success_;

}
