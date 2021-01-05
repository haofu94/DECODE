#ifndef allocate
# define allocate


#include "../include/dream.h"


int test_pointer(void *pointer){
  if (pointer == NULL) {
    return _failure_;
  } else {
    return _success_;
  }
}


//////////////////////////////
void *check_alloc(void *pointer){
  dream_call(test_pointer(pointer), _alloc_error_message_);
  return pointer;
}
//////////////////////////////


int double_malloc(int length, double **pointer){
  *pointer = (double *)malloc(length* sizeof(double));
  return test_pointer(pointer);
}


int double_calloc(int length, double **pointer){
  *pointer = (double *)calloc(length, sizeof(double));
  return test_pointer(pointer);
}


int int_malloc(int length, int **pointer){
  *pointer = (int *)malloc(length* sizeof(int));
  return test_pointer(pointer);
}


int int_calloc(int length, int **pointer){
  *pointer = (int *)calloc(length, sizeof(int));
  return test_pointer(pointer);
}


int char_malloc(int length, char **pointer){
  *pointer = malloc(length);
  return test_pointer(pointer);
}


int dealloc(void *pointer){
  free(pointer);
  pointer = NULL;
  if (pointer != NULL) {
    return _failure_;
  } else {
    return _success_;
  }
}



int dealloc_SMHM_matrix(double **smhm_data_matrix,
                        double **smhm_data_redshift,
                        double **smhm_data_Mstar){

  dream_call(dealloc(*smhm_data_matrix), _dealloc_error_message_);
  dream_call(dealloc(*smhm_data_redshift), _dealloc_error_message_);
  dream_call(dealloc(*smhm_data_Mstar), _dealloc_error_message_);

  return _success_;
}
int dealloc_SMHM_matrix_(SMHM_matrix **smhm_data){

  //int i;

  dream_call(dealloc((*smhm_data)->redshift), _dealloc_error_message_);
  printf("redshift freed\n");
  dream_call(dealloc((*smhm_data)->Mstar), _dealloc_error_message_);
  printf("Mstar freed\n");

  /*for(i=0; i<smhm_data->rows-1; i++){
    dream_call(dealloc(*(smhm_data->matrix+i)), _dealloc_error_message_);
  }*/
  dream_call(dealloc((*smhm_data)->matrix), _dealloc_error_message_);
  printf("matrix freed\n");
  //dream_call(dealloc(smhm_data), _dealloc_error_message_);

  return _success_;
}



#endif
