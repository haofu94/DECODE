/** @ file satellite_galaxies.c
  *
  * Written by Hao Fu
  *
  *
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
#include "../star_formation/star_formation_satellites.c"



void get_satellite_galaxies_at_z(char *logfile_name,
                                 char *file_name,
                                 char *output_folder,
                                 int len_mergers,
                                 int len_parents,
                                 stellar_mass_halo_mass *SMHM,
                                 double satellites_redshift,
                                 int ignore_high_orders,
                                 int include_quenching,
                                 int include_stripping,
                                 mergers_parameters *mergers_params,
                                 cosmological_time *cosmo_time,
                                 cosmological_parameters *cosmo_params){

  int i, j;
  double *id_parents;
  double *mass_parents;
  double *id_satellites;
  double *order;
  double *satellites;
  double *redshift_infall;
  double *merging_timescale;
  double *z_at_merge;
  //double SMHM_params[8];
  double age_today, age_at_merge;
  double satellite_mass;

  double random_mass_1st, z_inf_1st, new_tau_merge_2nd, new_age_at_merge_2nd;
  double *shmf_1st_order, *mass_range, *cumulative_mass_function;
  double resolution; int len;

  int smhm_data_rows, smhm_data_cols;
  double *smhm_data_matrix, *smhm_data_redshift, *smhm_data_Mstar;

  if (SMHM->is_analytical == _False_){
    dream_call(SMHM_read_matrix(SMHM, &smhm_data_rows, &smhm_data_cols, &smhm_data_matrix, &smhm_data_redshift, &smhm_data_Mstar), _SMHM_read_matrix_error_message_);
  }

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
  dream_call(double_malloc(len_parents, &id_parents), _alloc_error_message_);
  dream_call(double_malloc(len_parents, &mass_parents), _alloc_error_message_);
  dream_call(load_data(0, logfile_name, path_parents, 0, &id_parents), _load_data_error_message_);
  dream_call(load_data(0, logfile_name, path_parents, 1, &mass_parents), _load_data_error_message_);
  dealloc(path_parents);

  /*id_satellites = read_data(1, logfile_name, file_name, 0);
  order = read_data(1, logfile_name, file_name, 1);
  satellites = read_data(1, logfile_name, file_name, 2);
  redshift_infall = read_data(1, logfile_name, file_name, 3);
  merging_timescale = read_data(1, logfile_name, file_name, 4);*/
  dream_call(double_malloc(len_mergers, &id_satellites), _alloc_error_message_);
  dream_call(double_malloc(len_mergers, &order), _alloc_error_message_);
  dream_call(double_malloc(len_mergers, &satellites), _alloc_error_message_);
  dream_call(double_malloc(len_mergers, &redshift_infall), _alloc_error_message_);
  dream_call(double_malloc(len_mergers, &merging_timescale), _alloc_error_message_);
  dream_call(load_data(1, logfile_name, file_name, 0, &id_satellites), _load_data_error_message_);
  dream_call(load_data(1, logfile_name, file_name, 1, &order), _load_data_error_message_);
  dream_call(load_data(1, logfile_name, file_name, 2, &satellites), _load_data_error_message_);
  dream_call(load_data(1, logfile_name, file_name, 3, &redshift_infall), _load_data_error_message_);
  dream_call(load_data(1, logfile_name, file_name, 4, &merging_timescale), _load_data_error_message_);

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

  /*if (SMHM->is_analytical == _True_){

    dream_call(get_SMHM_params(SMHM_model,
                                &SMHM_params[0]),
                _get_SMHM_params_error_message_);

  }*/


  if (ignore_high_orders == _False_){

    resolution = 1.e-3; len = 50;

    age_today = linear_interp(0., cosmo_time->redshift, cosmo_time->age, cosmo_time->length);

    dream_call(double_malloc(len,
                              &shmf_1st_order),
                _alloc_error_message_);

    dream_call(double_malloc(len,
                              &mass_range),
                _alloc_error_message_);

    dream_call(double_malloc(len,
                              &cumulative_mass_function),
                _alloc_error_message_);

  }


  for (i=0; i<len_parents; i++){


    if (ignore_high_orders == _False_){

      dream_call(linear_space(resolution*mass_parents[i],
                               mass_parents[i],
                               len,
                               &mass_range),
                  _linear_space_error_message);

      dream_call(vdB_USHMF_1st_order(params_jiang_vdb_1st_order,
                                      mass_parents[i],
                                      mass_range,
                                      len,
                                      &shmf_1st_order),
                  _mass_function_error_message_);

      dream_call(cumsum(shmf_1st_order,
                         len,
                         &cumulative_mass_function),
                  _cumsum_error_message_);

    }


    for (j=0; j<len_mergers; j++){

      if (order[j]>2){
        continue;
      }

      if (*(id_parents+i) == *(id_satellites+j)){

        if (*(order+j) == 1){

          if (z_at_merge[j] < satellites_redshift){

            dream_call(get_evolved_satellite_mass(SMHM, smhm_data_rows, smhm_data_cols, smhm_data_matrix, smhm_data_redshift, smhm_data_Mstar,
                                                   *(mass_parents+i),
                                                   *(satellites+j),
                                                   satellites_redshift,
                                                   *(redshift_infall+j),
                                                   include_stripping,
                                                   include_quenching,
                                                   cosmo_time,
                                                   &satellite_mass),
                        _mstar_at_z_error_message_);

            fprintf(file_pointer, "%d %d %lf\n", (int)(*(id_satellites+j)), (int)(*(order+j)), satellite_mass);

          }

        } else if ((*(order+j) == 2) && (ignore_high_orders == _False_)) {

          dream_call(generate_halo_from_mass_function(cumulative_mass_function,
                                                       mass_range,
                                                       len,
                                                       &random_mass_1st),
                      _mass_from_pdf_error_message_);

          dream_call(get_infall_redshift_from_pdf(random_mass_1st,
                                                   1, 0.,
                                                   *(redshift_infall+j),
                                                   &z_inf_1st),
                      _infall_redshift_from_pdf_error_message_);

          if (*(z_at_merge+j) < z_inf_1st){

            dream_call(compute_merging_timescale(random_mass_1st,
                                                  get_mass_at_t(satellites[j],
                                                                random_mass_1st,
                                                                z_inf_1st,
                                                                redshift_infall[j],
                                                                cosmo_time,
                                                                cosmo_params),
                                                  z_inf_1st,
                                                  mergers_params,
                                                  cosmo_params,
                                                  &new_tau_merge_2nd),
                        _merging_timescale_error_message_);

            new_age_at_merge_2nd =  linear_interp(z_inf_1st, cosmo_time->redshift, cosmo_time->age, cosmo_time->length) + new_tau_merge_2nd;

            if (new_tau_merge_2nd > age_today) {

              dream_call(get_evolved_satellite_mass(SMHM, smhm_data_rows, smhm_data_cols, smhm_data_matrix, smhm_data_redshift, smhm_data_Mstar,
                                                     *(mass_parents+i),
                                                     *(satellites+j),
                                                     satellites_redshift,
                                                     *(redshift_infall+j),
                                                     include_stripping,
                                                     include_quenching,
                                                     cosmo_time,
                                                     &satellite_mass),
                          _mstar_at_z_error_message_);

              fprintf(file_pointer, "%d %d %lf\n", (int)(*(id_satellites+j)), (int)(*(order+j)), satellite_mass);

            }
          }

        }

      }

    }
  }

  fclose(file_pointer);

  dream_call(dealloc(output_filename),
              _dealloc_error_message_);

  dream_call(dealloc(id_parents),
              _dealloc_error_message_);

  dream_call(dealloc(mass_parents),
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

  if (SMHM->is_analytical == _False_){
    dream_call(dealloc_SMHM_matrix(&smhm_data_matrix, &smhm_data_redshift, &smhm_data_Mstar),
                _dealloc_error_message_);
  }

  return;
}
