/** @ file analyze.c
  *
  * Written by Hao Fu
  *
  * UNDER DEVELOPMENT
  **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "include/dream.h"
#include "include/cosmological_model.h"

//#include "analytical_SHMFs.c"
#include "dark_matter/mergers.c"



int get_mask_in_range(double *data,
                      double lower_limit,
                      double upper_limit,
                      int len,
                      int **mask) {

  /**
    * Function for identifying elements in array data which
    * lies in the internal [lower_limit, upper_limit].
    *
    * Returns a mask array of booleans.
  **/

  int i;

  for (i=0; i<len; i++){
    if ( (lower_limit <= (*(data+i))) && (*(data+i) < upper_limit) ) {
      *(*mask+i) = _True_;
    } else {
      *(*mask+i) = _False_;
    }
  }

  return _success_;
}
bool *get_mask_in_range_(double *data,
                        double lower_limit,
                        double upper_limit,
                        int len) {

  /**
    * Function for identifying elements in array data which
    * lies in the internal [lower_limit, upper_limit].
    *
    * Returns a mask array of booleans.
  **/

  int i;
  bool *mask;
  mask = (bool *)calloc(len, sizeof(bool));

  for (i=0; i<len; i++){
    if ( (lower_limit <= (*(data+i))) && (*(data+i) < upper_limit) ) {
      *(mask+i) = true;
    } else {
      *(mask+i) = false;
    }
  }

  return mask;
}



int array1_elements_in_array2(int *array1,
                              int *array2,
                              int len1,
                              int len2,
                              int **mask) {

  /**
    * Function for checking whether elements in array1 are present in array2.
    *
    * Returns an array of booleans of the length of array1.
  **/

  int i, j;

  for (i=0; i<len1; i++) {

    *(*mask+i) = _False_;

    for (j=0; j<len2; j++) {

      if (*(array1+i) == (*(array2+j))) {

        *(*mask+i) = _True_;
        break;

      }
    }
  }

  return _success_;
}
bool *array1_elements_in_array2_(int *array1,
                                int *array2,
                                int len1,
                                int len2) {

  /**
    * Function for checking whether elements in array1 are present in array2.
    *
    * Returns an array of booleans of the length of array1.
  **/

  int i,j;
  bool *mask;
  mask = (bool *)calloc(len1, sizeof(bool));

  for (i=0; i<len1; i++) {

    *(mask+i) = false;

    for (j=0; j<len2; j++) {

      if (*(array1+i) == (*(array2+j))) {

        *(mask+i) = true;
        break;

      }
    }
  }

  return mask;
}



int get_evolved_mergers_number(double z_,
                               double *mergers,
                               double *z_at_merge,
                               double *z_infall,
                               int len){

  /**
    * Function for counting the evolved mergers at given z_
    *
    * Inputs:
    *   1) z_ --> redshift for mergers rate
    *   2) mergers --> mergers array
    *   3) z_at_merge --> redshift at full merging
    *   4) z_infall --> redshift at first accretion
    *   5) len --> mergers array length
  **/

  int i, count;

  count = 0;

  for (i=0; i<len; i++) {

    if ( (*(z_at_merge+i) < z_) && (z_ < (*(z_infall+i))) ) {
    //if ( (z_ < (*(z_at_merge+i))) && (*(z_at_merge+i) < (*(z_infall+i))) ) {

      count += 1;

    }
  }

  return count;
}



double *get_evolved_mergers(double Mh,
                            double z_,
                            double *mergers,
                            double *z_at_merge,
                            double *z_infall,
                            double age_at_z,
                            double *age_at_z_inf,
                            double h,
                            double *Om,
                            double *Hz,
                            double H0,
                            int len,
                            int evolved_mergers_number) {

  /**
    * Function for calculating the the evolved mergers mass at z_
    *
    * Inputs:
    *   1) Mh --> parent halo mass
    *   2) z_ --> redshift for mergers rate
    *   3) mergers --> mergers mass
    *   4) z_at_merge --> redshift at full merging
    *   5) z_infall --> redshift at first accretion
    *   6) age_at_z --> age of the Universe at z_
    *   7) age_at_z_inf --> age of the Universe at first accretion
    *   8) h, 0m, Hz, H0 --> cosmological parameters
    *   9) len --> length of the mergers array
    *   10) evolved_mergers_number --> number of evolved mergers
  **/

  cosmological_parameters *cosmo_params = malloc (sizeof (cosmological_parameters));
  cosmo_params->h = h;
  cosmo_params->H0 = H0;

  int i, count;

  double *evolved_mergers;
  evolved_mergers = (double *)calloc(evolved_mergers_number, sizeof(double));

  count = 0;

  cosmo_params->Om = Om;
  cosmo_params->Hz = Hz;

  for (i=0; i<len; i++) {

    /*cosmo_params->Om = *(Om+i);
    cosmo_params->Hz = *(Hz+i);*/

    if ( (*(z_at_merge+i) < z_) && (z_ < (*(z_infall+i))) ) {

      *(evolved_mergers+count) = mass_at_t(*(mergers+i), Mh, z_, *(z_infall+i), age_at_z, *(age_at_z_inf+i), cosmo_params);
      count += 1;

    }

    if (count == evolved_mergers_number) {

      break;

    }
  }

  return evolved_mergers;

}
