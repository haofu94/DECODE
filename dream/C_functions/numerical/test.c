#ifndef test
# define test


void dream_call(int outcome, char *error_message){
  if (outcome == _failure_){
    printf("%s\n", error_message);
    exit(0);
  } else if (outcome == _success_){
    return;
  }
}


#endif
