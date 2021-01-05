/** @ file Mhalo_to_Mstar.c
  *
  * Same as file DM_to_SM.py
  **/

#ifndef Mhalo_to_Mstar
# define Mhalo_to_Mstar

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/dream.h"

#include "../numerical/num_c.c"



int SMHM_read_matrix(stellar_mass_halo_mass *SMHM,
                     int *smhm_data_rows,
                     int *smhm_data_cols,
                     double **smhm_data_matrix,
                     double **smhm_data_redshift,
                     double **smhm_data_Mstar){

  int i, j;

  FILE *file_pointer;
  file_pointer = fopen(SMHM->file, "r");
  fscanf(file_pointer, "%d %d", &(*smhm_data_rows), &(*smhm_data_cols));

  dream_call(double_malloc((*smhm_data_rows-1)*(*smhm_data_cols-1),
                            &(*smhm_data_matrix)),
              _alloc_error_message_);

  dream_call(double_malloc(*smhm_data_cols-1,
                            &(*smhm_data_redshift)),
              _alloc_error_message_);

  dream_call(double_malloc(*smhm_data_rows-1,
                            &(*smhm_data_Mstar)),
              _alloc_error_message_);

  for (j=0; j<*smhm_data_cols-1; j++){

    fscanf(file_pointer, "%lf", &(*smhm_data_redshift)[j]);
  }

  fscanf(file_pointer, "%*lf");

  for (i=0; i<*smhm_data_rows-1; i++){

    for (j=0; j<*smhm_data_cols-1; j++){

      fscanf(file_pointer, "%lf", &(*smhm_data_matrix)[i*(*smhm_data_cols-1)+j]);
    }

    fscanf(file_pointer, "%lf", &(*smhm_data_Mstar)[i]);
  }

  fclose(file_pointer);

  return _success_;
}



int SMHM_numerical_interp(double DM,
                          stellar_mass_halo_mass *SMHM,
                          int smhm_data_rows,
                          int smhm_data_cols,
                          double *smhm_data_matrix,
                          double *smhm_data_redshift,
                          double *smhm_data_Mstar,
                          double redshift,
                          double *Mstar_out){

  int i, j;

  double *Mhalo, *Mhalo_of_z;

  dream_call(double_calloc(smhm_data_rows-1,
                            &Mhalo),
              _alloc_error_message_);

  dream_call(double_calloc(smhm_data_cols-1,
                            &Mhalo_of_z),
              _alloc_error_message_);

  for (i=0; i<smhm_data_rows-1; i++){

    if (redshift <= *(smhm_data_redshift+smhm_data_cols-2)) {
    //if (redshift >= 0.) {

      for (j=0; j<smhm_data_cols-1; j++){

        Mhalo_of_z[j] = smhm_data_matrix[i*(smhm_data_cols-1)+j];
      }

      Mhalo[i] = linear_interp(redshift, smhm_data_redshift, Mhalo_of_z, smhm_data_cols-2); //2?
      //if (redshift >= 2.6){printf("%lf %lf\n", Mhalo[i], smhm_data_Mstar[i]);}

    } else {

      //Mhalo[i] = smhm_data_matrix[i*(smhm_data_cols-1)+smhm_data_cols-2];

      for (j=0; j<smhm_data_cols-1; j++){
        Mhalo_of_z[j] = smhm_data_matrix[i*(smhm_data_cols-1)+j];
      }
      Mhalo[i] = linear_interp(redshift, smhm_data_redshift, Mhalo_of_z, smhm_data_cols-2);

    }

  }


  /*if (redshift > 3.5) {
    printf("redshift = %lf\n", redshift);
    FILE *fp;
    fp = fopen("smhm.txt", "w");
    for (i=0; i<smhm_data_rows-1; i++){
      fprintf(fp, "%lf %lf\n", Mhalo[i], smhm_data_Mstar[i]);
    }
    fclose(fp);
    exit(0);
  }*/


  *Mstar_out = linear_interp(DM, Mhalo, smhm_data_Mstar, smhm_data_rows-1);

  if (SMHM->scatter > 0.) {

    *Mstar_out += get_random_gaussian(0., SMHM->scatter, -2., 2.);
  }

  dream_call(dealloc(Mhalo),
              _dealloc_error_message_);

  dream_call(dealloc(Mhalo_of_z),
              _dealloc_error_message_);

  return _success_;
}



int get_SMHM_params(char *SMHM_model,
                    double *SMHM_params) {

  int i, j;

  for (i=0; i<SMHM_models_analytical_number; i++){

    if (strcmp(SMHM_model, *(SMHM_models_analytical+i)) == 0) {

      for (j=0; j<8; j++){

        *(SMHM_params+j) = *(*(SMHM_params_analytical+i)+j);
      }

      break;
    }
  }

  if (i >= SMHM_models_analytical_number) {

    return _failure_;
  }

  return _success_;
}



int SMHM_Grylls_param(double DM,
                      double *Params,
                      double scatter,
                      double z,
                      double *Mstar) {

  double M = *Params;
  double Mz = *(Params+1);
  double N = *(Params+2);
  double Nz = *(Params+3);
  double b = *(Params+4);
  double bz = *(Params+5);
  double g = *(Params+6);
  double gz = *(Params+7);

  M = M + Mz * ((z-0.1)/(z+1.));
  N = N + Nz * ((z-0.1)/(z+1.));
  b = b + bz * ((z-0.1)/(z+1.));
  g = g + gz * ((z-0.1)/(z+1.));

  *Mstar = log10( pow(10., DM) * (2.*N*pow( (pow(pow(10.,DM-M), -b) + pow(pow(10.,DM-M), g)), -1.)) );

  if (scatter > 0.) {

    *Mstar += get_random_gaussian(0., scatter, -1., 1.);
  }

  return _success_;
}



int Mhalo_track_to_Mstar_track(double *Mhalo_track,
                               stellar_mass_halo_mass *SMHM,
                               int smhm_data_rows,
                               int smhm_data_cols,
                               double *smhm_data_matrix,
                               double *smhm_data_redshift,
                               double *smhm_data_Mstar,
                               double *z,
                               int len,
                               double **Mstar_track) {

  int i;

  double SMHM_params[8];

  if (SMHM->is_analytical == _True_){

    dream_call(get_SMHM_params(SMHM->model,
                                &SMHM_params[0]),
                _get_SMHM_params_error_message_);

    for (i=0; i<len; i++) {

      dream_call(SMHM_Grylls_param(*(Mhalo_track+i),
                                    SMHM_params,
                                    SMHM->scatter,
                                    *(z+i),
                                    &(*(*Mstar_track+i))),
                  _SMHM_error_message_);
    }

  } else if (SMHM->is_analytical == _False_){

    for (i=0; i<len; i++){

      dream_call(SMHM_numerical_interp(*(Mhalo_track+i),
                                        SMHM, smhm_data_rows, smhm_data_cols, smhm_data_matrix, smhm_data_redshift, smhm_data_Mstar,
                                        *(z+i), &(*(*Mstar_track+i))),
                  _SMHM_error_message_);

                  //printf("%d/%d, %lf %lf %lf\n", i+1, len, z[i], *(Mhalo_track+i), *(*Mstar_track+i));

    }

  }

  /*if ((*Mstar_track[0]>10.99) && (*Mstar_track[0]<11.1) ){
    printf("printing\n");
    FILE *fp;
    fp = fopen("van_den_Bosch_tracks/track_12.txt", "w");
    for (i=0; i<len; i++){
      fprintf(fp, "%lf %lf %lf\n", z[i], Mhalo_track[i], (*Mstar_track)[i]);
    }
    fclose(fp);
    exit(0);
  }*/

  return _success_;

}





double SMHM_Moster(double DM,
                   double *Params,
                   double z) {

  /**
    * Moster SMHM relation
  **/

  double Mstar;

  double M = *Params;
  double Mz = *(Params+1);
  double N = *(Params+2);
  double Nz = *(Params+3);
  double b = *(Params+4);
  double bz = *(Params+5);
  double g = *(Params+6);
  double gz = *(Params+7);
  double scatter = *(Params+8);

  M = M + Mz * ((z-0.1)/(z+1.));
  N = N + Nz * ((z-0.1)/(z+1.));
  b = b + bz * ((z-0.1)/(z+1.));
  g = g + gz * ((z-0.1)/(z+1.));

  Mstar = log10( pow(10., DM) * (2.*N*pow( (pow(pow(10.,DM-M), -b) + pow(pow(10.,DM-M), g)), -1.)) );

  if (scatter > 0.) {
    Mstar += get_random_gaussian(0., scatter, -10., 10.);
  }

  return Mstar;
}



#endif
