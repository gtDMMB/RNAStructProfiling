#ifndef OPTIONS
#define OPTIONS

/*
//file name to write output
extern char *OUTPUT;
//file name contains input structure(s)
extern char *INPUT;
//file name contains native structure
extern char *NATIVE;
//if filter is 1, keep number of freq helices under 10
extern int FILTER; 
//defines h and p st they cut off no more than NOISE percentage of area  
extern int NOISE;
//if 1, prints additional info
extern int VERBOSE;
//minimum size a helix must be to be considered significant
extern int MIN_HEL_LEN;
//if freq > thresh_freq, is significant (in percentage)
extern double THRESH_FREQ;
//if freq > thresh_common, assume is in everything
extern double THRESH_COMMON;
//limit to number of frequent helices
extern int NUMFREQ;
//limit to number of profiles made
extern int NUMPROF;
//profiles must have freq at least PROF_FREQ
extern double PROF_FREQ;
//number of structures processed; sfold generates 1000
extern int NUMSTRUCTS;
//lower bound for percent of structures represented by final profiles
extern double THRESH_STRUCT;
//sets minimum length a freq helix must be
extern int LENGTH;
//makes triplet stats in process_structs
extern int STATS;
//turns on or off making representative structures by consensus
extern int REP_STRUCT;
*/

/*default minimum significant helix size*/
#define DEF_MIN_HEL_LEN 1
/*default output file name*/
#define DEF_OUTPUT "profile.dot"

typedef struct options {
   /*all the options*/
  char *OUTPUT;
  char *INPUT;
  char *NATIVE;
  int VERBOSE;
  int SFOLD;
  int MIN_HEL_LEN;
  int NUM_FHC;
  int NUM_SPROF;
  double HC_FREQ;
  double PROF_FREQ;
  double COVERAGE;
  int NUMSTRUCTS;
  int PNOISE;
  int CYCLES;
  int GRAPH;
  int REP_STRUCT;
  int TOPDOWN;
  int ALTTHRESH;
} Options;

Options* make_options();
void print_options();
void free_options(Options* opt);

#endif
