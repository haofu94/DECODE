#ifndef smhm_relation
# define smhm_relation

#define MAX_SMHM_matrix_redshift_length 32
#define MAX_SMHM_matrix_Mstar_length 256


typedef struct stellar_mass_halo_mass{

  int is_analytical;
  char *model;
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
  double *matrix;
  double *redshift;
  double *Mstar;
  /*double matrix[MAX_SMHM_matrix_redshift_length*MAX_SMHM_matrix_Mstar_length];
  double redshift[MAX_SMHM_matrix_redshift_length];
  double Mstar[MAX_SMHM_matrix_Mstar_length];*/

} SMHM_matrix;



#endif
