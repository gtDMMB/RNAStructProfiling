#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hashtbl.h"
#include "Set.h"
#include "Options.h"
#include "memoryDFS.h"

//input first the fasta file, then the sample_1000.out file run on the fasta, then options
int main(int argc, char *argv[]) {
  int i, h, minh,p;
  HASHTBL *deleteHash;
  FILE *fp;
  Set *set;
  Options *opt;

  if (argc < 3) {
    fprintf(stderr,"Not enough arguments\n");
    exit(EXIT_FAILURE);
  }

  set = make_Set(argv[2]);
  opt = set->opt;
  for (i = 3; i < argc; i++) {
    //printf("argv[%d] is %s\n",i,argv[i]);
    if (!strcmp(argv[i],"-h")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%f",&(opt->HC_FREQ))) {
	opt->HC_FREQ = atof(argv[i+1]);
	if (opt->HC_FREQ < 0 || opt->HC_FREQ > 100) {
	  fprintf(stderr,"Error: invalid input %f for frequency threshold\n",opt->HC_FREQ);
	  opt->HC_FREQ = 0;
	}
	i++;
      }
    }
    else if (!strcmp(argv[i],"-p")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%f",&(opt->PROF_FREQ))) {
	opt->PROF_FREQ = atof(argv[i+1]);
	if (opt->PROF_FREQ < 0 || opt->PROF_FREQ > 100) {
	  fprintf(stderr,"Error: invalid input %f for frequency threshold\n",opt->PROF_FREQ);
	  opt->PROF_FREQ = 0;
	}
	i++;
      }
    }
    else if (!strcmp(argv[i],"-c")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%f",&(opt->COVERAGE))) {
	opt->COVERAGE = atof(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-f")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&(opt->NUM_FHC))) {
	opt->NUM_FHC = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-s")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&(opt->NUM_SPROF))) {
	opt->NUM_SPROF = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-l")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&(opt->MIN_HEL_LEN))) {
	opt->MIN_HEL_LEN = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-u")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&(opt->NUMSTRUCTS))) {
	opt->NUMSTRUCTS = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-m")) {
      if (i + 1 <= argc - 1) {
	opt->PNOISE = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-o")) {
      if (i + 1 <= argc - 1) {
	opt->OUTPUT = argv[i+1];
	i++;
      }
    }
    else if (!strcmp(argv[i],"-i")) {
      if (i + 1 <= argc - 1) {
	opt->INPUT = argv[i+1];
	i++;
      }
    }
    else if (!strcmp(argv[i],"-n")) {
      if (i + 1 <= argc - 1) {
	opt->NATIVE = argv[i+1];
	i++;
      }
    }
    else if (!strcmp(argv[i],"-k")) {
      if (i + 1 <= argc - 1) {
	opt->CYCLES = argv[i+1];
	i++;
      }
    }
    else if (!strcmp(argv[i],"-v"))
      opt->VERBOSE = 1;
    else if (!strcmp(argv[i],"-g"))
      opt->GRAPH = 0;
    else if (!strcmp(argv[i],"-r"))
      opt->REP_STRUCT = 1;
    else if (!strcmp(argv[i],"-t"))
      opt->TOPDOWN = 1;
    else if (!strcmp(argv[i],"-a"))
      opt->ALTTHRESH = 0;
  }

  input_seq(set,argv[1]);
  process_structs(set);
  reorder_helices(set);
  minh = print_all_helices(set);
  printf("Total number of helix classes: %d\n",set->hc_num);
  
  if (set->opt->TOPDOWN) {
    printf("Total number of extended profiles: %d\n",set->prof_num);
    h = top_down_h(set,minh);
    //if (set->opt->VERBOSE)
    printf("Number of featured helix classes: %d\n",h+1);
    find_freq(set);
    p = top_down_p(set,h);
    //if (set->opt->VERBOSE)
    printf("Number of selected profiles: %d\n",p+1);
    print_topdown_prof(set,h,p);
  } else {
    if (set->opt->NUM_FHC)
      set->opt->HC_FREQ = set_num_fhc(set);
    else if (set->opt->HC_FREQ==0) 
      set->opt->HC_FREQ = set_threshold(set,H_START);
    
    if (set->opt->VERBOSE) {
      printf("Threshold to find frequent helices: %.1f\%\n",set->opt->HC_FREQ);
      printf("Number of structures processed: %d\n",set->opt->NUMSTRUCTS);
    }
    
    find_freq(set);
    
    printf("Total number of featured helix classes: %d\n",set->num_fhc);
    make_profiles(set);
    
    printf("Total number of profiles: %d\n",set->prof_num);
    print_profiles(set);
    
    if (set->opt->NUM_SPROF)
      set->opt->PROF_FREQ = set_num_sprof(set);
    else if (set->opt->PROF_FREQ == 0) {
      set->opt->PROF_FREQ = set_p_threshold(set,P_START);
    }
    if (set->opt->VERBOSE)
      printf("setting p to %.1f\n",set->opt->PROF_FREQ);
    select_profiles(set);
    printf("Total number of selected profiles: %d\n",set->num_sprof);
  }
  
  if (set->opt->INPUT)
    process_one_input(set);
  if (set->opt->REP_STRUCT) {
    find_consensus(set);
    print_consensus(set);
  }
  
  if (set->opt->GRAPH) {
    fp = fopen(set->opt->OUTPUT,"w");
    init_graph(fp,set);
    initialize(set);
    if (set->opt->INPUT)
      print_input(fp,set);
    find_LCAs(fp,set);
    calc_gfreq(fp,set);
    //printGraph();
    deleteHash = MemoryDFS(set->graph);
    removeEdges(deleteHash);
    //start_trans_reductn(set->graph);
    //printGraph();
    print_edges(fp,set);
    fputs("}",fp);
    fclose(fp);
    hashtbl_destroy(deleteHash);
  }
  return 0;
}
