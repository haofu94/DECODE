/** @ file num_c.c
  *
  * Written by Hao Fu
  *
  * The main goal of this module is to provide some useful basic numerical tools.
  **/


#ifndef num_c
# define num_c


#include <time.h>

#include "../include/dream.h"


int find_index(double *array,
               double element,
               int len) {

  /**
    * Ordered array
  **/

  int i, idx;

  idx = 0;

  if ((element < (*array)) || (element >= (*(array+len-1)))) {

    idx = -1;

  } else {

    for (i=0; i<len-1; i++) {
      if ((*(array+i) <= element) && (element < (*(array+i+1)))) {
        idx = i;
        break;
      }
    }
  }

  return idx;

}



int arange(double min,
           double max,
           double bin,
           int len,
           double **res){

  int i;
  //int len = (int)((max-min)/bin);
  //double *res;
  //res = (double *)calloc(len, sizeof(double));
  *res[0] = min;
  for (i=1; i<len; i++){
    *(*res+i) = *(*res+i-1) + bin;
  }

  return _success_;
}


int linear_space(double min,
                 double max,
                 int len,
                 double **res){

  /**
    * Returns a linear space.
  **/

  double bin = (max - min) / (double)(len-1);
  int i;
  *res[0] = min;
  for (i=1; i<len; i++){
    *(*res+i) = *(*res+i-1) + bin;
  }

  return _success_;
}


double * flip(double *array,
              int len){

  /**
    * Returns the inverted array.
  **/

  double *res;
  res = (double *)malloc(len* sizeof(double));
  int i;

  for (i=0; i<len; i++){
    *(res+i) = *(array+len-1-i);
  }

  return res;
}


int cumsum(double *array, int len, double **res){

  /**
    * Return the cumulative sum of an array's element.
  **/

  int i, j;

  for (i=0; i<len; i++){
    for (j=0; j<=i; j++){
      *(*res+i) += *(array+j);
    }
  }

  return _success_;
}


double sum(double *array,
           int len) {

  /**
    * Return the sum of an array's element.
  **/

  int i;
  double res;
  res = 0.;
  for (i=0; i<len; i++){
    res += *(array+i);
  }
  return res;
}


int int_sum(int *array,
            int len) {

  /**
    * Return the sum of an integer array's element.
  **/

  int i;
  int res;
  res = 0.;
  for (i=0; i<len; i++){
    res += *(array+i);
  }
  return res;
}


double mean(double *array, int len){

  int i;

  double res = 0.;
  for (i=0; i<len; i++){

    res += *(array+i);
  }
  res /= (double)len;

  return res  ;
}


double min(double *array,
           int len){

  /**
    * Return the minimum of an array.
  **/

  double Min;
  int i;
  Min = array[0];
  for (i=1; i<len; i++){
    if (*(array+i)<Min) {
      Min = *(array+i);
    }
  }

  return Min;
}


double max(double *array,
           int len){

  /**
    * Return the maximum of an array.
  **/

  double Max;
  int i;
  Max = array[0];
  for (i=1; i<len; i++){
    if (*(array+i)>Max) {
      Max = *(array+i);
    }
  }

  return Max;
}


int double_append(int length, double **pointer, double value) {
  *pointer = (double *)realloc(*pointer, (length+1)* sizeof(double));
  *(*pointer+length) = value;
  return _success_;
}


int int_append(int length, int **pointer, int value) {
  *pointer = (int *)realloc(*pointer, (length+1)* sizeof(int));
  *(*pointer+length) = value;
  return _success_;
}


double integrate_trapz(double *y, double *x, int len){

  /**
    * Function for integration with trapezoidal rule
  **/

  int i;
  double res;

  res = 0.;
  for (i=0; i<len-1; i++){
    //res += 0.5 * fabs(*(x+i+1) - (*(x+i))) * fabs(*(y+i+1) + (*(y+i)));
    res += 0.5 * (*(x+i+1) - (*(x+i))) * (*(y+i+1) + (*(y+i)));
  }

  return res;
}


double line_through_two_points(double x, double x1, double y1, double x2, double y2){
  return (x-x1) * (y1-y2) / (x1-x2) + y1;
}


double linear_interp(double x,
                     double *x_array,
                     double *y_array,
                     int len){

  int i;
  double y;

  /*int j;
  if ((x >= 1.3) && (x<=1.5)){
    for(j=0;j<len;j++){printf("%lf %lf %lf\n", x, x_array[j], y_array[j]);}
  }*/

  if (x_array[0] < x_array[len-1]){ // x crescente

    if (x<x_array[0]){
      y = line_through_two_points(x, x_array[0], y_array[0], x_array[1], y_array[1]);
    } else if (x_array[len-1]<=x){
      y = line_through_two_points(x, x_array[len-2], y_array[len-2], x_array[len-1], y_array[len-1]);
    } else {
      for (i=0; i<len-1; i++){
        if ( (x_array[i]<=x) && (x<x_array[i+1]) ){
          //if ((x >1.3) && (x<1.5)){printf("%d %lf <= x < %lf\n", i, x_array[i], x_array[i+1]);}
          y = line_through_two_points(x, x_array[i], y_array[i], x_array[i+1], y_array[i+1]);
          break;
        }
      }
    }

  } else if (x_array[0] > x_array[len-1]){

    if (x_array[0]<x){
      //printf("ciao 1\n");
      y = line_through_two_points(x, x_array[0], y_array[0], x_array[1], y_array[1]);
    } else if (x<=x_array[len-1]){
      //printf("ciao 2\n");
      y = line_through_two_points(x, x_array[len-2], y_array[len-2], x_array[len-1], y_array[len-1]);
    } else {
      //printf("ciao 3\n");
      for (i=0; i<len-1; i++){
        //printf("%lf %lf\n", x_array[i], x_array[i+1]);
        if ( (x_array[i]>=x) && (x>x_array[i+1]) ){
          y = line_through_two_points(x, x_array[i], y_array[i], x_array[i+1], y_array[i+1]);
          break;
        }
      }
    }
  }

  //if ((x >= 1.3) && (x<=1.5)){printf("y= %lf\n", y);}

  return y;
}
double linear_interp_(double x,
                     double *x_array,
                     double *y_array,
                     int len){

  /**
    * Function for linear interpolation
  **/

  double y;
  int i, idx;

  if (*x_array < (*(x_array+len-1))){
    if (x == *x_array){
      y = *y_array;
    } else if (x < *x_array) {
      y = *y_array + (x - (*x_array)) * (*(y_array+1) - (*y_array)) / (*(x_array+1) - (*x_array));
    } else if (x > *(x_array+len-1)) {
      y = *(y_array+len-2) + (x - (*(x_array+len-2))) * (*(y_array+len-1) - (*(y_array+len-2))) / (*(x_array+len-1) - (*(x_array+len-2)));
    } else if ((x > *x_array) && (x <= *(x_array+len-1))) {
      for (i=0; i<len-1; i++){
        if ((*(x_array+i) > x) && (x <= *(x_array+i+1))) {
          y = *(y_array+i-1) + (x - (*(x_array+i-1))) * (*(y_array+i) - (*(y_array+i-1))) / (*(x_array+i) - (*(x_array+i-1)));
          break;
        }
      }
    }

  } else if (*x_array > (*(x_array+len-1))){

    if (x == *x_array){
      y = *y_array;
    } else if (x > *x_array) {
      y = *y_array + (x - (*x_array)) * (*(y_array+1) - (*y_array)) / (*(x_array+1) - (*x_array));
    } else if (x < *(x_array+len-1)){
      y = *(y_array+len-2) + (x - (*(x_array+len-2))) * (*(y_array+len-1) - (*(y_array+len-2))) / (*(x_array+len-1) - (*(x_array+len-2)));
    } else if ((x < *x_array) && (x >= *(x_array+len-1))) {
      for (i=0; i<len-1; i++){
        if ((*(x_array+i) < x) && (x >= *(x_array+i+1))) {
          y = *(y_array+i-1) + (x - (*(x_array+i-1))) * (*(y_array+i) - (*(y_array+i-1))) / (*(x_array+i) - (*(x_array+i-1)));
          break;
        }
      }
    }

  }

  return y;
}


double get_random_uniform(double min,
                          double max) {

  /**
    * Returns a random number with uniform probability distribution
  **/

  double num;
  //srand(time(NULL));
  //srand(0);
  num = ((double)rand() / (double)RAND_MAX);
  num = num * (max-min) + min;
  return num;
}


int get_int_random_uniform(int min,
                           int max){

  /**
    * Returns a random integer number with uniform probability distribution
  **/

  double num;
  num = ((double)rand() / (double)RAND_MAX);
  num = num * ((double)max-(double)min) + (double)min;
  return (int)(round(num));
}


double gaussian_distribution(double x,
                             double mean,
                             double sigma,
                             double norm) {

  /**
    * Returns Gaussian function
  **/

  return norm * exp(-0.5 * ((x - mean) / sigma) * ((x - mean) / sigma)) / (sigma * sqrt(2.*pi));

}


double get_random_gaussian(double mean,
                           double sigma,
                           double lower_bound,
                           double upper_bound) {

  /**
    * Returns a random number with Gaussian probability distribution
  **/

  int i;
  double num;

  int len = 100;

  double *x;
  x = (double *)malloc(len* sizeof(double));

  dream_call(linear_space(lower_bound,
                           upper_bound,
                           len,
                           &x),
              _linear_space_error_message);

  double *gauss;
  gauss = (double *)malloc(len* sizeof(double));
  for (i=0; i<len; i++) {
    *(gauss+i) = gaussian_distribution(*(x+i), mean, sigma, 1.);
  }

  double *cumu_gauss;
  cumu_gauss = (double *)calloc(len, sizeof(double));

  dream_call(cumsum(gauss,
                     len,
                     &cumu_gauss),
              _cumsum_error_message_);

  double low = min(cumu_gauss, len);
  double high = max(cumu_gauss, len);
  num = linear_interp(get_random_uniform(low, high), cumu_gauss, x, len);

  free(x); x = NULL;
  free(gauss); gauss = NULL;

  return num;
}



double get_random_from_distribution(double *x, double *PDFx, int len){

  // get_random_from_distribution(z_array, PDFz, halo_accretion->length-1);
  //int i; for(i=0; i<len; i++) { printf("%d %lf %lf\n", i, x[i], PDFx[i]); }
  /**
    * Returns a random number with the given probability distribution
  **/

  double num;

  double *cumu_pdf;
  cumu_pdf = (double *)calloc(len, sizeof(double));

  dream_call(cumsum(PDFx,
                    len,
                    &cumu_pdf),
              _cumsum_error_message_);

  double low = min(cumu_pdf, len);
  double high = max(cumu_pdf, len);
  num = linear_interp(get_random_uniform(low, high), cumu_pdf, x, len);

  free(cumu_pdf); cumu_pdf=NULL;

  return num;

}



int compute_histogram(double *data, int len_data, double *bins, int len_bins, double **histogram){

  /**
    * Inputs
    * histogram: array to write histogram on, must be of lenght len_bins-1
    *
    * Returns the histogram of data
    **/

  int i, j;

  for (i=0; i<len_bins-1; i++){

    *(*histogram+i) = 0.;

    if (i == len_bins-2) {

      for (j=0; j<len_data; j++){

        if ( (bins[i] <= data[j]) && (data[j] <= bins[i+1]) ){
          *(*histogram+i) += 1.;
        }

      }

    } else {

      for (j=0; j<len_data; j++){

        if ( (bins[i] <= data[j]) && (data[j] < bins[i+1]) ){
          *(*histogram+i) += 1.;
        }

      }

    }

  }

  return _success_;

}



int load_data(int halo_type,
              char *logfile_name,
              char *file_name,
              int column,
              double **data){

  int i, j, stop, len, Num_columns;
  FILE *file_pointer;

  int skiprows = 10;
  int max_line_length = 1000;
  char line[max_line_length];
  int headlines = 0;
  char * left;

  file_pointer = fopen(logfile_name, "r");

  stop=0;
  while (fgets(line, max_line_length, file_pointer) != NULL) {
    headlines++;

    /* Eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
        left++;
    }

    if ((left[0] != '#') && (left[0] != '\n')) {
      sscanf(line, "%d %d %d", &stop, &len, &Num_columns);
      if (stop == halo_type) {
        break;}
    }
  }

  fclose(file_pointer);

  double *data_trial;
  data_trial = (double *)calloc(Num_columns, sizeof(double));

  file_pointer = fopen(file_name, "r");

  if (Num_columns == 2){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf", &data_trial[0], &data_trial[1]);
        *(*data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  } else if (Num_columns == 3){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf %lf", &data_trial[0], &data_trial[1], &data_trial[2]);
        *(*data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  } else if (Num_columns == 4){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf %lf %lf", &data_trial[0], &data_trial[1], &data_trial[2], &data_trial[3]);
        *(*data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  } else if (Num_columns == 5){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf %lf %lf %lf", &data_trial[0], &data_trial[1], &data_trial[2], &data_trial[3], &data_trial[4]);
        *(*data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  } else if (Num_columns == 6){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf %lf %lf %lf %lf", &data_trial[0], &data_trial[1], &data_trial[2], &data_trial[3], &data_trial[4], &data_trial[5]);
        *(*data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  } else if (Num_columns == 7){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &data_trial[0], &data_trial[1], &data_trial[2], &data_trial[3], &data_trial[4], &data_trial[5], &data_trial[6]);
        *(*data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  }

  fclose(file_pointer);

  dream_call(dealloc(data_trial), _dealloc_error_message_);

  return _success_;
}
double *read_data(int halo_type,
                  char *logfile_name,
                  char *file_name,
                  int column) {

  /**
    * Functions for reading data from file.
    * Returns an array with the desired column read from the file.
  **/

  int i, j, stop, len, Num_columns;
  FILE *file_pointer;

  int skiprows = 10;
  int max_line_length = 1000;
  char line[max_line_length];
  int headlines = 0;
  char * left;

  file_pointer = fopen(logfile_name, "r");

  stop=0;
  //printf("ciao\n");
  //printf("%s\n", logfile_name);
  //printf("%p\n", file_pointer);printf("ciao\n");
  while (fgets(line, max_line_length, file_pointer) != NULL) {
    headlines++;

    /* Eliminate blank spaces at beginning of line */
    left=line;
    while (left[0]==' ') {
        left++;
    }

    if ((left[0] != '#') && (left[0] != '\n')) {
      sscanf(line, "%d %d %d", &stop, &len, &Num_columns);
      if (stop == halo_type) {
        //printf("%d %d %d\n", stop, len, Num_columns);
        break;}
    }
  }

  fclose(file_pointer);

  double *data; // NOTE TO MYSELF: memory leak to be solved
  data = (double *)calloc(len, sizeof(double));

  double *data_trial;
  data_trial = (double *)calloc(Num_columns, sizeof(double));

  file_pointer = fopen(file_name, "r");

  if (Num_columns == 2){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf", &data_trial[0], &data_trial[1]);
        *(data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  } else if (Num_columns == 3){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf %lf", &data_trial[0], &data_trial[1], &data_trial[2]);
        *(data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  } else if (Num_columns == 4){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf %lf %lf", &data_trial[0], &data_trial[1], &data_trial[2], &data_trial[3]);
        *(data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  } else if (Num_columns == 5){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf %lf %lf %lf", &data_trial[0], &data_trial[1], &data_trial[2], &data_trial[3], &data_trial[4]);
        *(data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  } else if (Num_columns == 6){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf %lf %lf %lf %lf", &data_trial[0], &data_trial[1], &data_trial[2], &data_trial[3], &data_trial[4], &data_trial[5]);
        *(data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  } else if (Num_columns == 7){
    j = 0;
    for (i=0; i<len+skiprows; i++){
      fgets(line,max_line_length, file_pointer);
      if ((line[0] != '#') && (line[0] != '\n')){
        sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &data_trial[0], &data_trial[1], &data_trial[2], &data_trial[3], &data_trial[4], &data_trial[5], &data_trial[6]);
        *(data+j) = (*(data_trial+column));
        j += 1;
      }
      if (j==len){
        break;
      }
    }
  }

  fclose(file_pointer);

  free(data_trial); data_trial = NULL;

  return data;

}


#endif
