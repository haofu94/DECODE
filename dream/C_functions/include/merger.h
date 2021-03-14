/**
  * Includes for (sub)haloes catalogue.
  **/

#ifndef merger
# define merger


typedef struct DM_halo_accretion{

  double *track;
  double *redshift;
  int length;

} DM_halo_accretion;



/**
  * @ id -> ID of the parent halo
  * @ halo_mass_at_z0 -> mass of the parent halo at z=0 [log(M/Msun)]
  * @ halo_mass_at_z -> mass of the parent halo at z [log(M/Msun)]
  * @ halo_mass_at_z_plus_dz -> mass of the parent halo at z+dz [log(M/Msun)]
  * @ redshift -> redshift at which the parent halo mass is given
  * @ redshift_bin -> redshift bin resolution
  * @ redshift_max -> maximum infall redshift to assign to subhaloes
  * @ subhalo_mass_range -> mass range for subhaloes [log(M/Msun)]
  * @ subhalo_mass_bin -> subhalo mass bin resolution [log(M/Msun)]
  * @ length -> length of subhalo_mass_range array
  * @ fudge -> fudge factor
  * @ max_order -> maximum order required for subhaloes
  **/
typedef struct mergers_parameters{

  int id;
  double halo_mass_at_z0;
  double halo_mass_at_z;
  double halo_mass_at_z_plus_dz;
  double redshift;
  double redshift_bin;
  double redshift_max;
  double *subhalo_mass_range;
  double subhalo_mass_bin;
  int length;
  int type_orbital_circularity; // 0 for constant, 1 for gaussian
  double orbital_circularity;
  double fudge;
  int max_order;

} mergers_parameters;



/**
  * @ total -> subhalo mass function with all orders subhaloes [dex^-1]
  * @ first -> first order subhalo mass function [dex^-1]
  * @ second -> second order subhalo mass function [dex^-1]
  * @ third -> third order subhalo mass function [dex^-1]
  * @ fourth -> forth order subhalo mass function [dex^-1]
  * @ fifth -> fifth order subhalo mass function [dex^-1]
  * @ psi -> ratio of the subhalo mass and parent halo mass
  * @ length_psi -> length of psi array
  **/
typedef struct subhalo_mass_functions{

  double *total;
  double *first;
  double *second;
  double *third;
  double *fourth;
  double *fifth;
  double *psi;
  int length_psi;

} subhalo_mass_functions;



/**
  * @ len_parents -> number of parent haloes
  * @ id_parents -> ID of the parent haloes
  * @ mass_parents -> mass of the parent haloes [log(M/Msun)]
  * @ len_mergers -> number of subhaloes
  * @ id_mergers -> ID of the subhaloes
  * @ mass_mergers -> mass of the subhaloes [log(M/Msun)]
  * @ z_infall -> infall redshift of the subhaloes
  * @ tau_merge -> merging timescale [Gyr]
  * @ z_at_merge -> redshift at merging
  **/
typedef struct DM_catalogue{

  int len_parents;
  int *id_parents;
  double *mass_parents;
  int len_mergers;
  int *id_mergers;
  int *order_mergers;
  double *mass_mergers;
  double *z_infall;
  double *tau_merge;
  double *z_at_merge;

} DM_catalogue;




#endif
