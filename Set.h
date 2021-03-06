/* RNAStructProfiling -- Profiles RNA structures and produces a summary graph in graphviz format.
 * Copyright 2013, 2014, 2018 Emily Rogers
 *
 * This file is part of RNAStructProfiling.
 *
 * RNAStructProfiling is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RNAStructProfiling is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RNAStructProfiling.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef SET
#define SET

#include "helix_class.h"
#include "Profile.h"
#include "hashtbl.h"
//#include "graph.h"
#include "Options.h"
#include "Profnode.h"

//default value that looking for h dropoff starts at (in percent)
#define H_START 10
//default value that looking for p dropoff starts at (in percent)
#define P_START 5
#define HASHSIZE 31
#define ARRAYSIZE 20
#define INIT_SIZE 2

typedef struct node
{
  char *label;
  struct node **neighbors;
  int numNeighbors;
  int nsize;
  int sfreq;
  int gfreq;
  char *bracket;
  char **diff;
  int DFS; 
  unsigned long sum;
} node;

typedef struct {
  char *seq;
  char *structfile;
  int hc_size;
  int hc_num;
  int helsum;
  int num_fhc;
  HC **helices;
  int **joint;
  int prof_size;
  int prof_num;
  int num_sprof;
  Profile **profiles;
  Profnode ***proftree;
  int *treeindex;
  int treesize;
  //Profile **freqprof;  //?
  double h_cutoff;     //? in options already  set->inputprof = NULL;

  double p_cutoff;     //? in options already
  Options *opt;
  node *inputnode;
  node *graph;
} Set;

node* createNode(char *name);
Set* make_Set(char *name);
void input_seq(Set *set,char *seqfile);
void process_structs(Set *set);
void process_structs_sfold(Set *set);
char* longest_possible(int id,int i,int j,int k,char *seq);
int match(int i,int j,char *seq);
void addHC(Set *set, HC *hc, int idcount);
void reorder_helices(Set *set);
int freqcompare(const void *v1, const void *v2);
double set_threshold_entropy(Set *set);
void init_joint(Set *set);
double set_threshold(Set *set, int start);
int compare(const void *v1, const void *v2);
int print_all_helices(Set *set);
double set_num_fhc(Set *set);
void find_freq(Set *set);
int top_down_h(Set *set,int minh);
int find_kink(double *opt, int j);
void translate(Profile *prof);
double split(Set *set, int index);
int nodecompare(const void *v1, const void *v2);
int split_one(Set *set,Profnode *node,int index);
int check_hc(char *prof, int index);
char* convert(int *array,int length);
int top_down_p(Set *set,int h);
int find_kink_p(Profnode **profs,int start, int stop);
void print_topdown_prof(Set *set, int h, int p);

void make_profiles(Set *set);
void make_profiles_sfold(Set *set);
void calc_joint(Set *set, int *prof, int num);
char* process_profile(HASHTBL *halfbrac,int *profile,int numhelix,Set *set);
void make_bracket_rep(HASHTBL *brac,Profile *prof);
void make_brackets(HASHTBL *brac, int i, int j, int id);
void make_rep_struct(HASHTBL *consensus,char *profile, char* trips);
void print_meta(Set *set);
void print_profiles(Set *set);
int profsort(const void *v1, const void *v2);
double set_num_sprof(Set *set);
double set_p_threshold_entropy(Set *set);
void find_general_freq(Set *set);
int subset(Set *set,char *one, char *two);
unsigned long binary_rep(Set *set,char *profile);
double set_p_threshold(Set *set, int start);
void select_profiles(Set *set);
void process_one_input(Set *set);
int* process_native(Set *set,int i, int j, int k);
void find_consensus(Set *set);
int print_consensus(Set *set);
void free_Set(Set *set);

#endif
