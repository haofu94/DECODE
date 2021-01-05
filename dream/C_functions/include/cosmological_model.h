/**
  * Includes for cosmological models set up.
  *
  * If the user is not familiar with the definition of "cosmological model",
  * a good reference could be the following:
  * Dedelson Scott - Modern Cosmology
  * https://ui.adsabs.harvard.edu/abs/2003moco.book.....D/abstract
  **/

#ifndef cosmological_model
# define cosmological_model


typedef struct cosmological_parameters{

  double Om0;
  double Ob0;
  double sigma8;
  double ns;
  double h;
  double H0;
  double *Om;
  double *Ob;
  double *Hz;
  double *z;
  int length;

} cosmological_parameters;


typedef struct cosmological_time{

  double *redshift;
  double *lookback_time;
  double *age;
  int length;

} cosmological_time;


#endif
