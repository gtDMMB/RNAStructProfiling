#include <stdlib.h>
#include <stdio.h>
#include "Options.h"

Options* make_options() {
  Options *opt = malloc(sizeof(Options));
  opt->OUTPUT = DEF_OUTPUT;
  opt->INPUT = NULL;
  opt->NATIVE = NULL;
  opt->MIN_HEL_LEN = DEF_MIN_HEL_LEN;
  opt->HC_FREQ = 0;
  opt->VERBOSE = 0;
  opt->PROF_FREQ = 0;
  opt->NUMSTRUCTS = 0;
  opt->GRAPH = 1;
  opt->REP_STRUCT = 0;
  return opt;
}

void free_options(Options* opt) {
  free(opt);
}
