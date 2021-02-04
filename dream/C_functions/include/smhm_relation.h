#ifndef smhm_relation
# define smhm_relation

#define MAX_SMHM_matrix_redshift_length 32
#define MAX_SMHM_matrix_Mstar_length 256


typedef struct stellar_mass_halo_mass{

  char *file;
  double *Mhalo;
  double *Mstar;
  int length;
  int constant_scatter;
  double scatter;

} stellar_mass_halo_mass;


typedef struct SMHM_matrix{

  int rows;
  int cols;
  double *matrix; // length = (cols-1)*(rows-1)
  double *redshift; // length = (cols-1)
  double *Mstar; // length = (rows-1)

} SMHM_matrix;



#endif
