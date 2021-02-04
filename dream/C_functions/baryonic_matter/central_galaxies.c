/** @ file central_galaxies.c
  *
  * Written by Hao Fu
  *
  *
  **/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/dream.h"

#include "../numerical/num_c.c"
#include "../numerical/allocate.c"
#include "../dark_matter/Mhalo_to_Mstar.c"


void get_central_galaxies(char *logfile_name,
                          char *file_name,
                          char *output_folder,
                          int len_centrals,
                          stellar_mass_halo_mass *SMHM){

  int i, j;
  double *centrals;
  double *id_centrals;
  FILE *file_pointer;

  SMHM_matrix *smhm_data;
  smhm_data = (SMHM_matrix *)malloc(sizeof(SMHM_matrix));

  dream_call(SMHM_read_matrix(SMHM, &smhm_data), _SMHM_read_matrix_error_message_);

  char filename[25] = "data/output_centrals.txt";
  const size_t len1 = strlen(output_folder);
  const size_t len2 = strlen(filename);
  char *output_filename;
  dream_call(char_malloc(len1 + len2 + 1,
                          &output_filename),
              _alloc_error_message_);
  strcpy(output_filename, output_folder);
  strcpy(output_filename + len1, filename);

  //centrals = read_data(0, logfile_name, file_name, 1);
  //id_centrals = read_data(0, logfile_name, file_name, 0);
  dream_call(double_malloc(len_centrals, &centrals), _alloc_error_message_);
  dream_call(double_malloc(len_centrals, &id_centrals), _alloc_error_message_);
  dream_call(load_data(0, logfile_name, file_name, 1, &centrals), _load_data_error_message_);
  dream_call(load_data(0, logfile_name, file_name, 0, &id_centrals), _load_data_error_message_);

  file_pointer = fopen(output_filename, "a");

  for (i=0; i<len_centrals; i++){

    dream_call(SMHM_numerical_interp(*(centrals+i), SMHM, smhm_data, 0.1, &(*(centrals+i))),
               _SMHM_error_message_);

    fprintf(file_pointer, "%d %lf\n", (int)(*(id_centrals+i)), *(centrals+i));

  }


  fclose(file_pointer);

  dream_call(dealloc(centrals),
             _dealloc_error_message_);

  dream_call(dealloc(id_centrals),
             _dealloc_error_message_);

  dream_call(dealloc(output_filename),
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
