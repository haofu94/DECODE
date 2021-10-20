/** @ file satellite_galaxies.c
  *
  * Written by Hao Fu
  *
  * The main goal of this module is to compute the population of
  * satellite galaxies at the given redshift of observation.
  **/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/cosmological_model.h"

#include "../numerical/num_c.c"
#include "../numerical/allocate.c"
#include "../dark_matter/mergers.c"
#include "../dark_matter/Mhalo_to_Mstar.c"
#include "ellipticals.c"
#include "../star_formation/star_formation_satellites.c"



void print_included_evolution_physics(int include_sat_SF,
                                      int include_mass_loss,
                                      int include_quenching,
                                      int include_stripping) {

  if ( (include_sat_SF == _False_)
    && (include_mass_loss == _False_)
    && (include_quenching == _False_)
    && (include_stripping == _False_) ) {
      printf("Frozen model\n");
  }

  if (include_sat_SF == _True_){
    printf("Including star formation: yes\n");
  } else if (include_sat_SF == _False_){
    printf("Including star formation: no\n");
  }

  if (include_mass_loss == _True_){
    printf("Including mass loss: yes\n");
  } else if (include_mass_loss == _False_){
    printf("Including mass loss: no\n");
  }

  if (include_quenching == _True_){
    printf("Including quenching: yes\n");
  } else if (include_quenching == _False_){
    printf("Including quenching: no\n");
  }

  if (include_stripping == _True_){
    printf("Including stripping: yes\n");
  } else if (include_stripping == _False_){
    printf("Including stripping: no\n");
  }

  return;
}



void get_satellite_galaxies_at_z(char *logfile_name,
                                 char *file_name,
                                 char *output_folder,
                                 int len_mergers,
                                 int len_parents,
                                 stellar_mass_halo_mass *SMHM,
                                 double observed_redshift,
                                 int ignore_high_orders,
                                 int include_sat_SF,
                                 int include_mass_loss,
                                 int include_quenching,
                                 int include_stripping,
                                 int compute_ellipticals,
                                 mergers_parameters *mergers_params,
                                 cosmological_time *cosmo_time,
                                 cosmological_parameters *cosmo_params){

  int i, j, k;
  int *id_parents;
  double *mass_parents;
  int *parID_satellites;
  int *upID_satellites;
  int *id_satellites;
  int *order;
  double *satellites;
  double *redshift_infall;
  double *merging_timescale;
  double *z_at_merge;
  double age_today, age_at_merge;
  double satellite_mass;

  double new_tau_merge_2nd, new_age_at_merge_2nd;

  SMHM_matrix *smhm_data;
  smhm_data = (SMHM_matrix *)malloc(sizeof(SMHM_matrix));

  dream_call(SMHM_read_matrix(SMHM, &smhm_data), _SMHM_read_matrix_error_message_);

  FILE *file_pointer;

  char filename[27] = "data/output_satellites.txt";
  const size_t len1 = strlen(output_folder);
  const size_t len2 = strlen(filename);
  char *output_filename;
  dream_call(char_malloc(len1 + len2 + 1,
                          &output_filename),
              _alloc_error_message_);
  strcpy(output_filename, output_folder);
  strcpy(output_filename + len1, filename);

  char filename_parents[24] = "data/output_parents.txt";
  const size_t len3 = strlen(filename_parents);
  char *path_parents;
  dream_call(char_malloc(len1 + len3 + 1,
                          &path_parents),
              _alloc_error_message_);
  strcpy(path_parents, output_folder);
  strcpy(path_parents + len1, filename_parents);
  //double *id_parents = read_data(0, logfile_name, path_parents, 0);
  //double *mass_parents = read_data(0, logfile_name, path_parents, 1);
  dream_call(int_malloc(len_parents, &id_parents), _alloc_error_message_);
  dream_call(double_malloc(len_parents, &mass_parents), _alloc_error_message_);
  dream_call(load_data_int(0, logfile_name, path_parents, 0, &id_parents), _load_data_error_message_);
  dream_call(load_data(0, logfile_name, path_parents, 1, &mass_parents), _load_data_error_message_);
  dealloc(path_parents);

  dream_call(int_malloc(len_mergers, &parID_satellites), _alloc_error_message_);
  dream_call(int_malloc(len_mergers, &upID_satellites), _alloc_error_message_);
  dream_call(int_malloc(len_mergers, &id_satellites), _alloc_error_message_);
  dream_call(int_malloc(len_mergers, &order), _alloc_error_message_);
  dream_call(double_malloc(len_mergers, &satellites), _alloc_error_message_);
  dream_call(double_malloc(len_mergers, &redshift_infall), _alloc_error_message_);
  dream_call(double_malloc(len_mergers, &merging_timescale), _alloc_error_message_);
  dream_call(load_data_int(1, logfile_name, file_name, 0, &parID_satellites), _load_data_error_message_);
  dream_call(load_data_int(1, logfile_name, file_name, 1, &upID_satellites), _load_data_error_message_);
  dream_call(load_data_int(1, logfile_name, file_name, 2, &id_satellites), _load_data_error_message_);
  dream_call(load_data_int(1, logfile_name, file_name, 3, &order), _load_data_error_message_);
  dream_call(load_data(1, logfile_name, file_name, 4, &satellites), _load_data_error_message_);
  dream_call(load_data(1, logfile_name, file_name, 5, &redshift_infall), _load_data_error_message_);
  dream_call(load_data(1, logfile_name, file_name, 6, &merging_timescale), _load_data_error_message_);


  dream_call(double_malloc(len_mergers,
                            &z_at_merge),
              _alloc_error_message_);

  dream_call(compute_redshift_at_merging(redshift_infall,
                                         merging_timescale,
                                         len_mergers,
                                         cosmo_time,
                                         &z_at_merge),
             _redshift_merging_error_message_);

  file_pointer = fopen(output_filename, "a");

  printf("\nComputing satellites' evolved mass...\n");

  print_included_evolution_physics(include_sat_SF, include_mass_loss, include_quenching, include_stripping);

  int *idx;
  int parID_counter = 1;
  dream_call(int_malloc(parID_counter, &idx), _alloc_error_message_);
  idx[0] = 0;
  for (i=1; i<len_mergers; i++){
    if (parID_satellites[i] != parID_satellites[i-1]){
      parID_counter += 1;
      idx = realloc(idx, parID_counter*sizeof(int));
      idx[parID_counter-1] = i;
    }
  }

  idx = realloc(idx, (parID_counter+1)*sizeof(int));
  idx[parID_counter] = len_mergers;

  double sat_2nd_mass;
  int idx_par;

  for (i=0; i<parID_counter; i++){

    for (j=0; j<len_parents; j++){
      if ( parID_satellites[idx[i]] == id_parents[j] ) {
        idx_par = j;
        break;
      }
    }


    for (j=idx[i]; j<idx[i+1]; j++){

      if (order[j] == 1){

        if (z_at_merge[j] > observed_redshift) {

          /*
          CASE 1
          The first order subhalo merges
          */

          for (k=idx[i]; k<idx[i+1]; k++){

            if ( (parID_satellites[k] == parID_satellites[j]) && (upID_satellites[k] == id_satellites[j]) ){

              if ( (z_at_merge[k] <= z_at_merge[j]) && (z_at_merge[k] > observed_redshift) ) {

                /*
                The second order subhalo does not merge before the first order does
                It is released to the parent halo
                */

                dream_call(compute_merging_timescale(mass_parents[idx_par],
                                                     get_mass_at_t(satellites[k],
                                                                   mass_parents[idx_par],
                                                                   redshift_infall[j],
                                                                   redshift_infall[k],
                                                                   cosmo_time,
                                                                   cosmo_params),
                                                     redshift_infall[j],
                                                     mergers_params,
                                                     cosmo_params,
                                                     &new_tau_merge_2nd),
                           _merging_timescale_error_message_);

                 new_age_at_merge_2nd =  linear_interp(redshift_infall[j], cosmo_time->redshift, cosmo_time->age, cosmo_time->length) + new_tau_merge_2nd;

                 if (new_tau_merge_2nd > age_today) {

                   dream_call(get_evolved_satellite_mass(SMHM, smhm_data,
                                                         mass_parents[idx_par],
                                                         *(satellites+k),
                                                         observed_redshift,
                                                         *(redshift_infall+j),
                                                         include_sat_SF,
                                                         include_mass_loss,
                                                         include_stripping,
                                                         include_quenching,
                                                         cosmo_params,
                                                         cosmo_time,
                                                         &sat_2nd_mass),
                              _mstar_at_z_error_message_);

                   fprintf(file_pointer, "%d %d %lf\n", *(id_satellites+k), *(order+k), sat_2nd_mass);
                 }
               }
             }
           }

        } else if (z_at_merge[j] <= observed_redshift) {

          /*
          CASE 2
          The first order subhalo does not merge
          */

          dream_call(get_evolved_satellite_mass(SMHM, smhm_data,
                                                *(mass_parents+i),
                                                *(satellites+j),
                                                observed_redshift,
                                                *(redshift_infall+j),
                                                include_sat_SF,
                                                include_mass_loss,
                                                include_stripping,
                                                include_quenching,
                                                cosmo_params,
                                                cosmo_time,
                                                &satellite_mass),
                     _mstar_at_z_error_message_);

          for (k=idx[i]; k<idx[i+1]; k++){

            if ( (parID_satellites[k] == parID_satellites[j]) && (upID_satellites[k] == id_satellites[j]) ){

              dream_call(get_evolved_satellite_mass(SMHM, smhm_data,
                                                    *(mass_parents+i),
                                                    *(satellites+k),
                                                    z_at_merge[j],
                                                    *(redshift_infall+k),
                                                    include_sat_SF,
                                                    include_mass_loss,
                                                    include_stripping,
                                                    include_quenching,
                                                    cosmo_params,
                                                    cosmo_time,
                                                    &sat_2nd_mass),
                         _mstar_at_z_error_message_);

              if (z_at_merge[k] > z_at_merge[j]) {

                /*
                The second order merges
                */

                satellite_mass = log10( pow(10., satellite_mass) + pow(10., sat_2nd_mass) );

              } else if ( (z_at_merge[k] <= z_at_merge[j]) && (z_at_merge[k] <= observed_redshift) ) {

                /*
                The second order does not merge
                */

                fprintf(file_pointer, "%d %d %lf\n", *(id_satellites+k), *(order+k), sat_2nd_mass);

              }

            }

          }

          fprintf(file_pointer, "%d %d %lf\n", *(id_satellites+j), *(order+j), satellite_mass);

        }

      }

    }

    printf("Progress: %.2f %%\r", (double)(i+1) / (double)parID_counter * 100.);
    fflush(stdout);

  }


  /*
  COMPUTE ELLIPTICAL (CENTRAL) GALAXIES
  */

  if (compute_ellipticals == _True_) {

    DM_catalogue *halo_catalogue = malloc(sizeof(DM_catalogue));
    halo_catalogue->id_parents = id_parents;
    halo_catalogue->mass_parents = mass_parents;
    halo_catalogue->len_parents = len_parents;
    halo_catalogue->len_mergers = len_mergers;
    halo_catalogue->id_mergers = parID_satellites;
    halo_catalogue->order_mergers = order;
    halo_catalogue->mass_mergers = satellites;
    halo_catalogue->z_infall = redshift_infall;
    halo_catalogue->z_at_merge = z_at_merge;

    dream_call(compute_fraction_ellipticals(file_name,
                                            output_folder,
                                            halo_catalogue,
                                            observed_redshift,
                                            SMHM, smhm_data,
                                            include_sat_SF,
                                            include_mass_loss,
                                            include_quenching,
                                            include_stripping,
                                            cosmo_params,
                                            cosmo_time),
               _compute_ellipticals_error_message_);

    dream_call(dealloc(halo_catalogue), _dealloc_error_message_);

  }

  printf("\n");


  fclose(file_pointer);

  dream_call(dealloc(output_filename),
             _dealloc_error_message_);

  dream_call(dealloc(id_parents),
             _dealloc_error_message_);

  dream_call(dealloc(mass_parents),
             _dealloc_error_message_);

  dream_call(dealloc(parID_satellites),
             _dealloc_error_message_);

  dream_call(dealloc(upID_satellites),
             _dealloc_error_message_);

  dream_call(dealloc(id_satellites),
             _dealloc_error_message_);

  dream_call(dealloc(order),
             _dealloc_error_message_);

  dream_call(dealloc(satellites),
             _dealloc_error_message_);

  dream_call(dealloc(redshift_infall),
             _dealloc_error_message_);

  dream_call(dealloc(merging_timescale),
             _dealloc_error_message_);

  dream_call(dealloc(z_at_merge),
             _dealloc_error_message_);

  dream_call(dealloc(idx),
             _dealloc_error_message_);

  dream_call(dealloc(smhm_data->redshift),
             _dealloc_error_message_);
  dream_call(dealloc(smhm_data->matrix),
             _dealloc_error_message_);
  dream_call(dealloc(smhm_data->Mstar),
             _dealloc_error_message_);
  dream_call(dealloc(smhm_data),
             _dealloc_error_message_);

  return;
}
