/** @ file centrals.c
  *
  * Written by:
  * Francesco Shankar
  * Chris Marsden
  * Hao Fu
  *
  * Same as file halo_growth.py
  **/

#ifndef centrals
# define centrals


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../include/dream.h"

#include "../numerical/num_c.c"



double fit_var_sem(double Mlog,
                   double h,
                   double omega_m,
                   double omega_b,
                   double sigma_8,
                   double spectral_index){

  /**
    * Direct fitting formula by van den Bosch [note that the masses are in units of Msun/h]
  **/

  int i;
  int N = 50;

  double *Mh, *sig;
  Mh = (double *)malloc(N * sizeof(double));
  sig = (double *)malloc(N * sizeof(double));
  double Mh_bin = (17. - 5.) / (double)(N-1);
  double gammam, x, g1, g2, f, s, res;
  gammam = omega_m * h * exp(-omega_b * (1. + sqrt(2.*h) / omega_m));

  double c = 3.804e-4;

  *Mh = 17.;
  for (i=0; i<N; i++){
    *(Mh+i) = *Mh - (double)i * Mh_bin;
    x = (c * gammam * pow(10., *(Mh+i)/3.)) / pow(omega_m, 1./3.);
    g1 = 64.087 * ( pow(1. + 1.074*pow(x, 0.3) - 1.581*pow(x, 0.4) + 0.954*pow(x, 0.5) - 0.185*pow(x, 0.6), -10.) );
    x = 32. * gammam;
    g2 = 64.087 * ( pow(1. + 1.074*pow(x, 0.3) - 1.581*pow(x, 0.4) + 0.954*pow(x, 0.5) - 0.185*pow(x,0.6), -10.) );
    f = (g1*g1) / (g2*g2);

    s = sqrt(f * sigma_8*sigma_8);

    *(sig+i) = s * pow(10, (*(Mh+i)-14.09) * (1.-spectral_index) / 9.2) /(1.+(1.-spectral_index)/9.2);
    //printf("%lf %lf\n", *(Mh+i), *(sig+i));
  }

  res = linear_interp(Mlog, Mh, sig, N);

  free(Mh); Mh = NULL;
  free(sig); sig = NULL;
  //exit(0);

  return res;
}



double D_z_White(double z,
                 double omega_m){

  double omega_l = 1. - omega_m;
  double omega_k = 0.;
  double Ez, omega_m_z, omega_l_z, gz, gz0, dz;

  Ez = pow(omega_l + omega_k * (1.+z)*(1.+z) + omega_m * (1.+z)*(1.+z)*(1.+z), 0.5);
  omega_m_z = omega_m * (1.+z)*(1.+z)*(1.+z) / (Ez*Ez);
  omega_l_z = omega_l / (Ez*Ez);
  gz=2.5 * omega_m_z * pow(pow(omega_m_z, 4./7.) - omega_l_z + (1. + omega_m_z / 2.) * (1. + omega_l_z / 70.), -1.);

  gz0 = 2.5 * omega_m * pow(pow(omega_m, 4./7.) - omega_l + (1. + omega_m / 2.) * (1. + omega_l / 70.), -1.);

  dz = (gz / gz0) / (1. + z);

  return dz;
}



int Mass_acc_history_VDB(double M0,
                         double *zc,
                         int len,
                         cosmological_parameters *cosmo_params,
                         double **logMz){

  int i, j;
  //M0 = 12.5;
  //printf("%lf %lf %lf %lf %lf\n", cosmo_params->h, cosmo_params->Om0, cosmo_params->Ob0, cosmo_params->sigma8, cosmo_params->ns);
  // Parameters for average
  double apar1 = 3.2954;
  double apar2 = 0.1975;
  double apar3 = 0.7542;
  double apar4 = 0.0898;
  double apar5 = 0.4415;

  double z0 = 0.;

  double logM0h = M0+log10(cosmo_params->h); //[Msun/h]

  double var = fit_var_sem(logM0h, cosmo_params->h, cosmo_params->Om0, cosmo_params->Ob0, cosmo_params->sigma8, cosmo_params->ns);
  double s0 = pow(var, 2.);

  //printf("%lf\n\n", var);

  double dc0 = 1.686/D_z_White(0., cosmo_params->Om0);

  int mass_resolution = 1000;
  double *logMzc;
  logMzc = (double *)malloc(mass_resolution* sizeof(double));

  dream_call(linear_space(M0-0.0001,
                           M0-7.,
                           mass_resolution,
                           &logMzc),
              _linear_space_error_message); //Possible range of masses

  double *logpsiz;
  double *psiz;
  double *Fpsi;
  double *logMzch;
  double *varz;
  double *sz;
  double *Gz;
  logpsiz = (double *)malloc(mass_resolution* sizeof(double));
  psiz = (double *)malloc(mass_resolution* sizeof(double));
  Fpsi = (double *)malloc(mass_resolution* sizeof(double));
  logMzch = (double *)malloc(mass_resolution* sizeof(double));
  varz = (double *)malloc(mass_resolution* sizeof(double));
  sz = (double *)malloc(mass_resolution* sizeof(double));
  Gz = (double *)malloc(mass_resolution* sizeof(double));


  for (i=0; i<mass_resolution; i++){
    *(logpsiz+i) = *(logMzc+i) - M0; //In log, the rato of the possible range of masses to M0
    *(psiz+i) = pow(10., *(logpsiz+i)); //The value un-logged
    *(Fpsi+i) = apar1 * pow(1. - apar2 * (*(logpsiz+i)), apar3) * pow(1. - pow(*(psiz+i), apar4), apar5);
    *(logMzch+i) = *(logMzc+i) + log10(cosmo_params->h);
    //printf("%d", i);
    *(varz+i) = fit_var_sem(*(logMzch+i), cosmo_params->h, cosmo_params->Om0, cosmo_params->Ob0, cosmo_params->sigma8, cosmo_params->ns);
    *(sz+i) = pow(*(varz+i), 2.);
    *(Gz+i) = 0.57 * pow((*(sz+i))/s0, 0.19) * pow(dc0 / (*(varz+i)), -0.01);
  }

  int nz = len;
  *(*logMz) = M0; //The first element is the last element
  double dcz, logpsi;
  double gam = 0.4;
  double *wz, *wtz, *x_array;
  wz = (double *)malloc(mass_resolution* sizeof(double));
  wtz = (double *)malloc(mass_resolution* sizeof(double));
  x_array = (double *)malloc(mass_resolution* sizeof(double));

  //Loop over every redshift interval beyond the first
  /*FILE *fp;
  fp = fopen("track_DM.txt", "w");
  fprintf(fp, "%lf %lf\n", zc[0], *logMz[0]);*/
  //printf("%lf %lf\n", zc[0], **logMz);

  for (i=0; i<nz-1; i++){
    dcz = 1.686/D_z_White(*(zc+i+1), cosmo_params->Om0);
    for (j=0; j<mass_resolution; j++){
      *(wz+j) = (dcz-dc0)/sqrt((*(sz+j))-s0);
      *(wtz+j) = *(wz+j) * pow(*(Gz+j), gam);
      *(x_array+j) = *(Fpsi+j) - (*(wtz+j));
    }
    logpsi = linear_interp(0., x_array, logpsiz, mass_resolution);
    /*if(i==0){
      //for(j=0;j<mass_resolution;j++){printf("%lf %lf\n", x_array[j], logpsiz[j]);}
      printf("%lf\n", logM0h);
    }*/
    *(*logMz+i+1) = logpsi + M0;
    //printf("%lf %lf\n", zc[i+1], *(*logMz+i+1));
    //fprintf(fp, "%lf %lf\n", zc[i+1], *(*logMz+i+1));
  }
  //fclose(fp); exit(0);

  free(logMzc); logMzc = NULL;
  free(logpsiz); logpsiz = NULL;
  free(psiz); psiz = NULL;
  free(Fpsi); Fpsi = NULL;
  free(logMzch); logMzch = NULL;
  free(varz); varz = NULL;
  free(sz); sz = NULL;
  free(Gz); Gz = NULL;
  free(wz); wz = NULL;
  free(wtz); wtz = NULL;
  free(x_array); x_array = NULL;

  return _success_;
}



double * Mass_acc_history_VDB_(double M0,
                              double *zc,
                              double h,
                              double omega_m,
                              double omega_b,
                              double sigma_8,
                              double spectral_index,
                              int len){

  int i, j;
  // Parameters for average
  double apar1 = 3.2954;
  double apar2 = 0.1975;
  double apar3 = 0.7542;
  double apar4 = 0.0898;
  double apar5 = 0.4415;

  double z0 = 0.;

  double logM0h = M0+log10(h); //[Msun/h]

  double var = fit_var_sem(logM0h, h, omega_m, omega_b, sigma_8, spectral_index);
  double s0 = var*var;

  double dc0 = 1.686/D_z_White(0., omega_m);

  int mass_resolution = 1000;
  double *logMzc;
  logMzc = (double *)malloc(mass_resolution* sizeof(double));

  dream_call(linear_space(M0-0.0001,
                           M0-7.,
                           mass_resolution,
                           &logMzc),
              _linear_space_error_message); //Possible range of masses

  double *logpsiz;
  double *psiz;
  double *Fpsi;
  double *logMzch;
  double *varz;
  double *sz;
  double *Gz;
  logpsiz = (double *)malloc(mass_resolution* sizeof(double));
  psiz = (double *)malloc(mass_resolution* sizeof(double));
  Fpsi = (double *)malloc(mass_resolution* sizeof(double));
  logMzch = (double *)malloc(mass_resolution* sizeof(double));
  varz = (double *)malloc(mass_resolution* sizeof(double));
  sz = (double *)malloc(mass_resolution* sizeof(double));
  Gz = (double *)malloc(mass_resolution* sizeof(double));


  for (i=0; i<mass_resolution; i++){
    *(logpsiz+i) = *(logMzc+i) - M0; //In log, the rato of the possible range of masses to M0
    *(psiz+i) = pow(10., *(logpsiz+i)); //The value un-logged
    *(Fpsi+i) = apar1 * pow(1. - apar2 * (*(logpsiz+i)), apar3)* pow(1. - pow(*(psiz+i), apar4), apar5);
    *(logMzch+i) = *(logMzc+i) + log10(h);
    //printf("%d", i);
    *(varz+i) = fit_var_sem(*(logMzch+i), h, omega_m, omega_b, sigma_8, spectral_index);
    *(sz+i) = (*(varz+i)) * (*(varz+i));
    *(Gz+i) = 0.57 * pow(*(sz+i) / s0, 0.19) * pow(dc0 / (*(varz+i)), -0.01);
  }

  int nz = len;
  double *logMz;
  logMz = (double *)malloc(nz* sizeof(double));
  logMz[0] = M0; //The first element is the last element
  double dcz, logpsi;
  double gam = 0.4;
  double *wz, *wtz, *x_array;
  wz = (double *)malloc(mass_resolution* sizeof(double));
  wtz = (double *)malloc(mass_resolution* sizeof(double));
  x_array = (double *)malloc(mass_resolution* sizeof(double));

  //Loop over every redshift interval beyond the first
  for (i=0; i<nz-1; i++){
    dcz = 1.686/D_z_White(*(zc+i+1), omega_m);
    for (j=0; j<mass_resolution; j++){
      *(wz+j) = (dcz-dc0)/sqrt(*(sz+j)-s0);
      *(wtz+j) = *(wz+j) * pow(*(Gz+j), gam);
      *(x_array+j) = *(Fpsi+j) - (*(wtz+j));
    }
    logpsi = linear_interp(0., x_array, logpsiz, mass_resolution);
    *(logMz+i) = logpsi + M0;
  }

  free(logMzc); logMzc = NULL;
  free(logpsiz); logpsiz = NULL;
  free(psiz); psiz = NULL;
  free(Fpsi); Fpsi = NULL;
  free(logMzch); logMzch = NULL;
  free(varz); varz = NULL;
  free(sz); sz = NULL;
  free(Gz); Gz = NULL;
  free(wz); wz = NULL;
  free(wtz); wtz = NULL;
  free(x_array); x_array = NULL;

  return logMz;
}



int *get_halo_index(double *halo_catalog,
                    int catalog_len,
                    double *mass_range,
                    int mass_range_len){

  /**
    * Functions for identifing the index of the
      corresponding mean halo track of each halo.
  **/

  int i,j;
  int *index;
  index = (int *)malloc(catalog_len* sizeof(int));

  for (i=0; i<catalog_len; i++){
    for (j=0; j<mass_range_len; j++){
      if ((*(mass_range+j)<=(*(halo_catalog+i))) && (*(halo_catalog+i)<(*(mass_range+j+1)))){
        *(index+i) = j;
        break;
      }
    }
  }

  return index;

}



int *get_halo_IDs(int Num_haloes,
                  int min,
                  int max){

  /**
    * Functions for getting the halo IDs
  **/

  int i;//,j;
  int id_trial;
  //int exist = 0; // 0 if does not exist, 1 if exists

  int *id_list;
  id_list = (int *)malloc(Num_haloes* sizeof(int));

  *(id_list) = min;
  for (i=1; i<Num_haloes; i++){
    *(id_list+i) = *(id_list+i-1)+1;
  }

  /**(id_list) = get_int_random_uniform(min, max);
  for (i=1; i<Num_haloes; i++){

    //printf("%d out of %d, ID = %d\n", i, Num_haloes, *(id_list+i-1));

    do {
      id_trial = get_int_random_uniform(min, max);

      for (j=0; j<i; j++){

        if (*(id_list+j) == (*(id_list+i))){
          exist = 1;
          break;
        }

      }

    } while(exist==1);

    *(id_list+i) = id_trial;

  }*/

  return id_list;

}



#endif
