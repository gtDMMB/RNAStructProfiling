#include <stdlib.h>
#include <stdio.h>
#include "Options.h"

Options* make_options() {
  Options *opt = (Options*) malloc(sizeof(Options));
  opt->OUTPUT = (char*) DEF_OUTPUT;
  opt->INPUT = NULL;
  opt->NATIVE = NULL;
  opt->MIN_HEL_LEN = DEF_MIN_HEL_LEN;
  opt->NUM_FHC = 0;
  opt->NUM_SPROF = 0;
  opt->HC_FREQ = -1;
  opt->VERBOSE = 0;
  opt->SFOLD = 0;
  opt->PROF_FREQ = -1;
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

void print_options() {
  puts("OPTIONS");
  puts("-e [FILE]\tExternal structure as input, following gtboltzmann format (dot bracket, energy, set of triplets on same line)");
  puts("-sfold [FILE]\tExternal structure as input, following Sfold format ('Structure xx' followed by a triplet per line");
  puts("-h DBL\t\tSet frequency threshold for features(in percentage, e.g. 10.5 for 10.5% threshold)");
  puts("-p DBL\t\tSet frequency threshold for selected profiles (in percentage, e.g. 10.5 for 10.5% threshold)");
  puts("-c DBL\t\tSet minimum coverage requirement for featured helix classes and selected profiles (in percentage)");
  puts("-f INT\t\tSet number of featured helix classes");
  puts("-s INT\t\tSet number of selected profiles");
  puts("-l INT\t\tSet minimum helix length");
  puts("-u INT\t\tSet number of structures to profile");
  /*puts("-m INT\t\tSet PNOISE");*/
  puts("-o NAME\t\tPrefix of output files");
  puts("-i FILE\t\tFile containing external structure to be inserted into summary profile graph");
  puts("-n FILE\t\tFile containing native structure to be inserted into summary profile graph");
  puts("-v \t\tRun in verbose mode");
  puts("-g \t\tRun without generating summary profile graph");
  //puts("-t \t\tRun with top-down alternate algorithm");
  //puts("-a \t\tRun with alternate threshold");
  puts("gtboltzmann options (passed to gtboltzmann):");
  //puts("-d, --dangle INT Restricts treatment of dangling energies (INT=0,2)");
  puts("--limitcd  INT Set a maximum base pair contact distance to INT. If no limit is given, base pairs can be over any distance");
  puts("--paramdir DIR   Path to directory from which parameters are to be read");
  puts("--useSHAPE FILE  Use SHAPE constraints from FILE");
  puts("-w, --workdir DIR    Path of directory where output files will be written");
  puts("--sample INT Number of structures to sample");
}

void free_options(Options* opt) {
  free(opt);
}
