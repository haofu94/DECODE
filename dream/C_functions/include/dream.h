#ifndef dream
# define dream


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cosmological_model.h"
#include "error_messages.h"
#include "merger.h"
#include "smhm_relation.h"


#define _True_ 1
#define _False_ 0

#define _success_ 1
#define _failure_ 0


#define pi M_PI


static double params_jiang_vdb_total_unevolved[5] = {0.22, -0.91, 6., 3., 1.}; //Jiang & van den Bosch 2016 Table A1, Unevolved, total
static double params_jiang_vdb_1st_order[6] = {0.13, -0.83, 1.33, -0.02, 5.67, 1.19}; //Jiang & van den Bosch 2016 Eq. (14), Unevolved, 1st order


static int SMHM_models_analytical_number = 1;
static char *SMHM_models_analytical[1] = {"Grylls_2019"};
static double SMHM_params_analytical[1][8] = {{11.92, 0.58, 0.032, -0.014, 1.64, -0.69, 0.53, 0.03}};

static int SMHM_models_numerical_number = 1;
static char *SMHM_models_numerical[1] = {"Fu_Dickson_2020"};


#include "../numerical/test.c"
#include "../numerical/allocate.c"
#include "../numerical/num_c.c"


#endif
