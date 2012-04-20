#ifndef SET
#define SET

#include "helix_class.h"
#include "Profile.h"
#include "hashtbl.h"
#include "graph.h"
#include "Options.h"

//default value that looking for h dropoff starts at (in percent)
#define H_START 100
//default value that looking for p dropoff starts at (in percent)
#define P_START 10
#define HASHSIZE 31
#define ARRAYSIZE 20
#define INIT_SIZE 2

typedef struct {
  char *seq;
  char *structfile;
  int hc_size;
  int hc_num;
  int num_fhc;
  HC **helices;
  //HC **freqhelices;
  int prof_size;
  int prof_num;
  int num_sprof;
  Profile **profiles;
  //Profile **freqprof;  //?
  double h_cutoff;     //? in options already  set->inputprof = NULL;

  double p_cutoff;     //? in options already
  Options *opt;
  node *inputnode;
  node *graph;
} Set;

Set* make_Set(char *name);
void input_seq(Set *set,char *seqfile);
void process_structs(Set *set);
char* longest_possible(int id,int i,int j,int k,char *seq);
int match(int i,int j,char *seq);
void addHC(Set *set, HC *hc, int idcount);
void reorder_helices(Set *set);
int freqcompare(const void *v1, const void *v2);
double set_threshold(Set *set, int start);
int compare(const void *v1, const void *v2);
void print_all_helices(Set *set);
double set_num_fhc(Set *set);
void find_freq(Set *set);
void make_profiles(Set *set);
char* process_profile(HASHTBL *halfbrac,int *profile,int numhelix,Set *set);
void make_bracket_rep(HASHTBL *brac,Profile *prof);
void make_brackets(HASHTBL *brac, int i, int j, int id);
void make_rep_struct(HASHTBL *consensus,char *profile, char* trips);
void print_profiles(Set *set);
int profsort(const void *v1, const void *v2);
double set_num_sprof(Set *set);
double set_p_threshold(Set *set, int start);
void select_profiles(Set *set);
void process_one_input(Set *set);
int* process_native(Set *set,int i, int j, int k);
void find_consensus(Set *set);
int print_consensus(Set *set);
void free_Set(Set *set);

#endif
