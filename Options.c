#include <stdlib.h>
#include <stdio.h>
#include "Options.h"

Options* make_options() {
  Options *opt = malloc(sizeof(Options));
  opt->OUTPUT = DEF_OUTPUT;
  opt->INPUT = NULL;
  opt->NATIVE = NULL;
  opt->MIN_HEL_LEN = DEF_MIN_HEL_LEN;
  opt->NUM_FHC = 0;
  opt->NUM_SPROF = 0;
  opt->HC_FREQ = 0;
  opt->VERBOSE = 0;
  opt->PROF_FREQ = 0;
  opt->COVERAGE = 0.5;
  opt->NUMSTRUCTS = 0;
  opt->PNOISE = 5;
  opt->CYCLES = 10;
  opt->GRAPH = 1;
  opt->REP_STRUCT = 0;
  opt->TOPDOWN = 0;
  opt->ALTTHRESH = 1;
  return opt;
}

void free_options(Options* opt) {
  free(opt);
}
