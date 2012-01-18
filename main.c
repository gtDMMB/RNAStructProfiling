#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hashtbl.h"
#include "profile.h"

HASHTBL *bp;
HASHTBL *marginals;
HASHTBL *idhash;
HASHTBL *freq;
HASHTBL *cluster;
HASHTBL *binary;
HASHTBL *bracket;
HASHTBL *avetrip;
HASHTBL *freqprofs;
char *OUTPUT;
char *INPUT;
char *NATIVE;
int FILTER;
int VERBOSE;
int NOISE;
int MIN_HEL_LEN;
double THRESH_FREQ;
double THRESH_COMMON;
int NUMFREQ;
int NUMPROF;
double PROF_FREQ;
int NUMSTRUCTS;
double THRESH_STRUCT;
int LENGTH;
int STATS;
int REP_STRUCT;

//turns on or off graphing section of code
static int GRAPH;

//input first the fasta file, then the sample_1000.out file run on the fasta, then options
int main(int argc, char *argv[]) {
  int i,total,notcommon,most;
  char **mostfreq;
  FILE *fp;

  OUTPUT = DEF_OUTPUT;
  INPUT = NULL;
  NATIVE = NULL;
  FILTER = 0;
  NOISE = DEF_NOISE;
  MIN_HEL_LEN = DEF_MIN_HEL_LEN;
  //  THRESH_FREQ = DEF_THRESH_FREQ;
  THRESH_FREQ = 0;
  THRESH_COMMON = DEF_THRESH_COMMON;
  VERBOSE = 0;
  NUMFREQ = 0;
  NUMPROF = 0;
  PROF_FREQ = 0;
  //  PROF_FREQ = DEF_PROF_FREQ;
  NUMSTRUCTS = 0;
  THRESH_STRUCT = DEF_THRESH_STRUCT;
  LENGTH = 0;
  STATS = 0;
  GRAPH = 1;
  REP_STRUCT = 0;

  if (argc < 3) {
printf("./main [sequence FASTA file] [sfold sample structures file] [-options]\n");
      printf("\t-o\toutput file name: default is profile.dot\n");
      printf("\t-h\tset helix class frequency threshold; in percentage\n");
      printf("\t-p\tset profile frequency threshold; in percentage\n");
      printf("\t-v\tverbose output\n");
      printf("\t-g\tturn off graphing option\n");

    exit(EXIT_FAILURE);
  }

  //OUTPUT = malloc(strlen(argv[1])+5);
  //sprintf(OUTPUT,"%s.dot",argv[1]);

  for (i = 3; i < argc; i++) {
    //printf("argv[%d] is %s\n",i,argv[i]);
    if (!strcmp(argv[i],"--help")) {
      printf("./main [sequence FASTA file] [sfold sample structures file] [-options]\n");
      printf("\t-o\toutput file name: default is profile.dot\n");
      printf("\t-h\tset helix class frequency threshold; in percentage\n");
      printf("\t-p\tset profile frequency threshold; in percentage\n");
      printf("\t-v\tverbose output\n");
      printf("\t-g\tturn off graphing option\n");
    }
    else if (!strcmp(argv[i],"-f")) {
      FILTER = 1;
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&NUMFREQ)) {
	//	sscanf(argv[i+1],"%s",val);
	//printf("val is %s and argv %s\n",val,argv[i+1]);
	NUMFREQ = atoi(argv[i+1]);
	i++;
      }
      else
	NUMFREQ = DEF_NUMFREQ;
    }
    else if (!strcmp(argv[i],"-z")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&NOISE)) {
	NOISE = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-h")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%f",&THRESH_FREQ)) {
	THRESH_FREQ = atof(argv[i+1]);
	if (THRESH_FREQ < 0 || THRESH_FREQ > 100) {
	  fprintf(stderr,"Error: invalid input %f for frequency threshold\n",THRESH_FREQ);
	  THRESH_FREQ = DEF_THRESH_FREQ;
	}
	i++;
      }
    }
    else if (!strcmp(argv[i],"-c")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%f",&THRESH_COMMON)) {
	THRESH_COMMON = atof(argv[i+1]);
	if (THRESH_COMMON < 0.0 || THRESH_COMMON > 100.0) {
	  fprintf(stderr,"Error: invalid input %f for common threshold\n",THRESH_COMMON);
	  THRESH_COMMON = DEF_THRESH_COMMON;
	}
	i++;
      }
    }
    else if (!strcmp(argv[i],"-q")) {
      //PRUNE = 1;
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&NUMPROF)) {
	NUMPROF = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-p")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%f",&PROF_FREQ)) {
	PROF_FREQ = atof(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-l")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&MIN_HEL_LEN)) {
	MIN_HEL_LEN = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-s")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&NUMSTRUCTS)) {
	NUMSTRUCTS = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-t")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%f",&THRESH_STRUCT)) {
	THRESH_STRUCT = atof(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-o")) {
      if (i + 1 <= argc - 1) {
	OUTPUT = argv[i+1];
	i++;
      }
    }
    else if (!strcmp(argv[i],"-i")) {
      if (i + 1 <= argc - 1) {
	INPUT = argv[i+1];
	i++;
      }
    }
    else if (!strcmp(argv[i],"-n")) {
      if (i + 1 <= argc - 1) {
	NATIVE = argv[i+1];
	i++;
      }
    }
    else if (!strcmp(argv[i],"-v"))
      VERBOSE = 1;
    else if (!strcmp(argv[i],"-a"))
      STATS = 1;
    else if (!strcmp(argv[i],"-g"))
      GRAPH = 0;
    else if (!strcmp(argv[i],"-r"))
      REP_STRUCT = 1;
  }

  if (!(bp = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for bp failed");
    exit(EXIT_FAILURE);
  }
  
  if (!(marginals = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for marginals failed");
    exit(EXIT_FAILURE);
  }

  if (!(idhash = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for idhash failed");
    exit(EXIT_FAILURE);
  }

  if (!(binary = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for binary failed");
    exit(EXIT_FAILURE);
  }

  total = process_structs(argv[1],argv[2]);
  reorder_helices(total);
  if (THRESH_FREQ==0) 
    THRESH_FREQ = set_h_dropoff(marginals, H_START);
  if (VERBOSE) {
    printf("Threshold to find frequent helices: %.1f\%\n",THRESH_FREQ);
    printf("Maximum number of frequent helices: ");
    if (NUMFREQ == 0)
      puts("no limit");
    else
      printf("%d\n",NUMFREQ);
    printf("Number of structures processed: %d\n",NUMSTRUCTS);
  }
  printf("Total number of equivalence helix classes: %d\n",total-1);
  if (VERBOSE)
    print_all_helices(total);
  
  //  make_graph(marginals,max,id,total,argv[1],fp);
  
  mostfreq = find_freq(total);
  printf("Total number of selected helices: %d\n",hashtbl_numkeys(freq));
  notcommon = make_profiles(argv[2]);
  printf("Total number of profiles: %d\n",hashtbl_numkeys(cluster));
  print_profiles();

  if (PROF_FREQ == 0) {
    PROF_FREQ = set_h_dropoff(cluster,P_START);
    if (VERBOSE)
      printf("setting p to %.1f\n",PROF_FREQ);
  }

  most = select_profiles(mostfreq,notcommon)-1;
  printf("Total number of selected profiles: %d\n",hashtbl_numkeys(freqprofs));
  if (hashtbl_numkeys(cluster) == 0)
    GRAPH = 0;
  if (REP_STRUCT) {
    //fp = fopen("structures.out","w");
    //fprintf(fp,"Processing %s\n",argv[2]);
    find_consensus();
    print_consensus(argv[1]);
    //print_cluster(argv[1]);
    //fclose(fp);
  }
  if (GRAPH) {
    fp = fopen(OUTPUT,"w");
    insert_graph(fp,argv[1],most);  
    fputs("}",fp);
    fclose(fp);
  }
  hashtbl_destroy(bp);
  hashtbl_destroy(marginals);
  hashtbl_destroy(idhash);
  hashtbl_destroy(binary);
  hashtbl_destroy(freq);
  hashtbl_destroy(cluster);
  hashtbl_destroy(bracket);
  hashtbl_destroy(freqprofs);
  return 0;
}
