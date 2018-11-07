#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Set.h"
#include "hashtbl.h"
#include "Options.h"
#include "helix_class.h"
#include "Profile.h"
#include "Profnode.h"
#include <math.h>
#include <float.h>

static HASHTBL *bp;
static HASHTBL *translate_hc;
static HASHTBL *consensus;

node* createNode(char *name)
{
  node* newNode;
  newNode = (node*)malloc(sizeof(node));
  newNode->label = name;
  newNode->neighbors = NULL;
  //newNode->neighbors = (node**)malloc(sizeof(node*)*GRAPHSIZE);
  newNode->nsize = 1;
  newNode->numNeighbors = 0;
  newNode->sfreq = 0;
  newNode->gfreq = 0;
  newNode->bracket = NULL;
  newNode->diff = NULL;
  newNode->DFS = 0;
  newNode->sum = 0;
  return newNode;
}

Set* make_Set(char *name) {
  Set *set = (Set*) malloc(sizeof(Set));
  set->seq = NULL;
  set->structfile = (char *) malloc((strlen(name) + 1) * sizeof(char));
  strcpy(set->structfile, name);
  set->hc_size = 5;
  set->hc_num = 0;
  set->helsum = 0;
  set->num_fhc = 0;
  set->helices = (HC**) malloc(sizeof(HC*)*ARRAYSIZE*5);
/*
  set->joint = (HASHTBL**) malloc(sizeof(HASHTBL*)*ARRAYSIZE*5);
  for (i=0; i < ARRAYSIZE*5; i++)
    set->joint[i] = NULL;
  */
  set->prof_size = 5;
  set->prof_num = 0;
  set->num_sprof = 0;
  set->profiles = (Profile**) malloc(sizeof(Profile*)*ARRAYSIZE*5);
  set->proftree = (Profnode***) malloc(sizeof(Profnode**)*ARRAYSIZE*2);
  set->treeindex = (int*) malloc(sizeof(int)*ARRAYSIZE*2);
  set->treesize = 2;
  set->h_cutoff = 0;
  set->p_cutoff = 0;
  set->opt = make_options();
  set->inputnode = NULL;
  set->graph = NULL;
  return set;
}

void input_seq(Set *set,char *seqfile) {
  FILE * fp;
  int size = 5,fasta = 0;
  char temp[100],*blank = (char*)" \n",*part,*final;

  fp = fopen(seqfile,(char*)"r");
  if (fp == NULL) {
    fprintf(stderr, "can't open %s\n",seqfile);
  }
  final = (char*) malloc(sizeof(char)*ARRAYSIZE*size);
  final[0] = '\0';
  while (fgets(temp,100,fp)) {    
    //put error handling in case first line has more than 100 chars
    //skipping first line if fasta header and not seq
    if (temp[0] == '>' || fasta) {
      //printf("found fasta header %s",temp);
      if (strlen(temp) < 99 || (strlen(temp) == 99 && temp[98] == '\n')) {
	fasta = 0;
	continue;
      }
      else 
	fasta = 1;
      continue;
    } 
    for (part = strtok(temp,blank); part; part = strtok(NULL,blank)) {
      if (strlen(final)+strlen(part) > (unsigned int)ARRAYSIZE*size-1) {
	while ((unsigned int)(++size*ARRAYSIZE - 1) < strlen(final)+strlen(part)) ;
	final = (char*) realloc(final,sizeof(char)*ARRAYSIZE*size);
      }
      final = strcat(final,part);
    }
  }
  if (set->opt->VERBOSE) 
    printf("seq in %s is %s with length %d\n",seqfile,final,(signed int)strlen(final));
  //printf("final char is %c\n",final[strlen(final)-1]);
  set->seq = final;
}

void process_structs(Set *set) {
  FILE *fp;
  int i,j,k,*helixid,idcount=1,*lg,last = 0, toosmall = 0,numhelix=0,*profile=NULL,size=1,seqsize;
  double *trip;
  char *tmp,*key,dbl[ARRAYSIZE],*max, *delim = (char*)" ,\t\n", *val;
  HASHTBL *halfbrac,*extra,*avetrip;
  HC *hc;

  fp = fopen(set->structfile,(char*)"r");
  if (fp == NULL) {
    fprintf(stderr, "can't open %s\n",set->structfile);
  }
  if (!(bp = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for bp failed");
    exit(EXIT_FAILURE);
  }
  if (!(avetrip = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for avetrip failed");
    exit(EXIT_FAILURE);
  }
  if (!(extra = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for extra failed");
    exit(EXIT_FAILURE);
  }
  if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
    exit(EXIT_FAILURE);
  }
  profile = (int*) malloc(sizeof(int)*ARRAYSIZE);
  key = (char*) malloc(sizeof(char)*ARRAYSIZE);
  seqsize = strlen(set->seq)*10;
  tmp = (char*) malloc(sizeof(char)*seqsize);
  while (fgets(tmp,seqsize,fp)) {
    numhelix=0;
    val = strtok(tmp,delim);
    val = strtok(NULL,delim);
    while ((val = strtok(NULL,delim))) {
      i = atoi(val);
      if ((val = strtok(NULL,delim)))
        j = atoi(val);
      else
        fprintf(stderr, "Error in file input format: struct %d\n",set->opt->NUMSTRUCTS);
      if ((val = strtok(NULL,delim)))
        k = atoi(val);
      else
        fprintf(stderr, "Error in file input format: struct %d\n",set->opt->NUMSTRUCTS);
      if (k < set->opt->MIN_HEL_LEN) {
        toosmall++;
      //printf("too small %d %d %d\n",i,j,k);
        continue;
      }
      sprintf(dbl,"%d %d",i,j);
      helixid = (int*)hashtbl_get(bp,dbl);
      numhelix++;
      //printf("%d %d %d\n",i,j,k);
      if (!helixid) {
	max = longest_possible(idcount,i,j,k,set->seq);
	hc = create_HC(idcount,max);
	addHC(set,hc,idcount);
	//triplet stats
	trip = (double*) malloc(sizeof(double)*3);
	trip[0] = i;
	trip[1] = j;
	trip[2] = k;
	sprintf(key,"%d",idcount);
	hashtbl_insert(avetrip,key,trip);
	last = idcount++;
      }
      else {
	//printf("Found %d %d with id %d\n",i,j,*helixid);
	sprintf(key,"%d",*helixid);
	if (last != *helixid) {
	  hc = set->helices[*helixid-1];
	  hc->freq++;
	} else {
	  if ((lg = (int*) hashtbl_get(extra,key)))
	    ++*lg;
	  else {
	    lg = (int*) malloc(sizeof(int));
	    *lg = 1;
	    hashtbl_insert(extra,key,lg);
	  }
	  //if (VERBOSE) 
	  //printf("Found repeat id %d:%s\n",last,hashtbl_get(idhash,key));
	}
	//average stats
	trip = (double*) hashtbl_get(avetrip,key);
	trip[0] += i;
	trip[1] += j;
	trip[2] += k;
	
	last = *helixid;
      }
      if (numhelix >= ARRAYSIZE*size) 
	profile = (int*) realloc(profile,sizeof(int)*ARRAYSIZE*++size);
      profile[numhelix-1] = last;
      //calc_joint(set,prof,numhelix);
    }
    set->opt->NUMSTRUCTS++;
  }
  for (i = 1; i < idcount; i++) {
    j = set->helices[i-1]->freq;
    sprintf(key,"%d",i);
    if ((lg = (int*) hashtbl_get(extra,key))) 
      k = j + *lg;
    else
      k = j;
    trip = (double*) hashtbl_get(avetrip,key);
    sprintf(key,"%.1f %.1f %.1f", trip[0]/k,trip[1]/k,trip[2]/k);
    hc = set->helices[i-1];
    hc->avetrip = mystrdup(key);
  }
  if (fclose(fp))
    fprintf(stderr, "File %s not closed successfully\n",set->structfile);
  hashtbl_destroy(extra);
  hashtbl_destroy(avetrip);
}


/*calculates joint for last value in prof vs everything else

void calc_joint(Set *set, int *prof, int num) {
  int k, *freq;  
  char val[ARRAYSIZE];
  HASHTBL *hash;
  
  if (!set->joint[prof[num-1]-1])
    set->joint[prof[num-1]-1] = hashtbl_create(HASHSIZE,NULL);
  if (num < 2)
    return;
  for (k = 0; k < num-1; k++) {
    if (prof[k] < prof[num-1]) {
      sprintf(val,"%d",prof[num-1]);
      hash = set->joint[prof[k]-1];
      if (!hash) fprintf(stderr,"hash doesn't exist in calc_joint\n");
      freq = (int*) hashtbl_get(hash,val);
      if (freq)
	*freq++;
      else {
	freq = (int*) malloc(sizeof(int));
	*freq = 1;
	hashtbl_insert(hash,val,freq);
      }
    }
    else {
      sprintf(val,"%d",prof[k]);
      hash = set->joint[prof[num-1]-1];
      if (!hash) fprintf(stderr,"hash doesn't exist in calc_joint\n");
      freq = (int*) hashtbl_get(hash,val);
      if (freq)
	*freq++;
      else {
	freq = (int*) malloc(sizeof(int));
	*freq = 1;
	hashtbl_insert(hash,val,freq);
      }
    }
  }
  return;
}
*/

void process_structs_sfold(Set *set) {
  FILE *fp;
  int i,j,k,*helixid,idcount=1,*lg,last = 0, toosmall = 0,numhelix=0,*profile=NULL,size=1;
  double *trip;
  char tmp[100],*key,dbl[ARRAYSIZE],*max;
  HASHTBL *halfbrac,*extra,*avetrip;
  HC *hc;

  fp = fopen(set->structfile,"r");
  if (fp == NULL) {
    fprintf(stderr, "can't open %s\n",set->structfile);
  }
  if (!(bp = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for bp failed");
    exit(EXIT_FAILURE);
  }
  if (!(avetrip = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for avetrip failed");
    exit(EXIT_FAILURE);
  }
  if (!(extra = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for extra failed");
    exit(EXIT_FAILURE);
  }
  if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
    exit(EXIT_FAILURE);
  }
  key = (char*) malloc(sizeof(char)*ARRAYSIZE);
  while (fgets(tmp,100,fp) != NULL) {
    if (sscanf(tmp,"%d %d %d",&i,&j,&k) == 3) {
      if (i == 0) continue;
      if (k < set->opt->MIN_HEL_LEN) {
	toosmall++;
	//printf("too small %d %d %d\n",i,j,k);
	continue;
      }
      sprintf(dbl,"%d %d",i,j);
      helixid = (int*) hashtbl_get(bp,dbl);
      //printf("%d %d %d\n",i,j,k);
      if (!helixid) {
	max = longest_possible(idcount,i,j,k,set->seq);
	hc = create_HC(idcount,max);
	addHC(set,hc,idcount);
	//triplet stats
	trip = (double*)malloc(sizeof(double)*3);
	trip[0] = i;
	trip[1] = j;
	trip[2] = k;
	sprintf(key,"%d",idcount);
	hashtbl_insert(avetrip,key,trip);
	if (set->opt->TOPDOWN) {
	  if (numhelix >= ARRAYSIZE*size) 
	    profile = (int*) realloc(profile,sizeof(int)*ARRAYSIZE*++size);
	  profile[numhelix++]=idcount;
	  make_brackets(halfbrac,i,j,idcount);
	}
	last = idcount++;
      }
      else {
	//printf("Found %d %d with id %d\n",i,j,*helixid);
	sprintf(key,"%d",*helixid);
	if (last != *helixid) {
	  hc = set->helices[*helixid-1];
	  hc->freq++;
	  if (set->opt->TOPDOWN) {
	    if (numhelix >= ARRAYSIZE*size) 
	      profile = (int*)realloc(profile,sizeof(int)*ARRAYSIZE*++size);
	    profile[numhelix++]=*helixid;
	    make_brackets(halfbrac,i,j,*helixid);
	  }
	} else {
	  if ((lg = (int*) hashtbl_get(extra,key)))
	    ++*lg;
	  else {
	    lg = (int*) malloc(sizeof(int));
	    *lg = 1;
	    hashtbl_insert(extra,key,lg);
	  }
	  //if (VERBOSE) 
	  //printf("Found repeat id %d:%s\n",last,hashtbl_get(idhash,key));
	}
	//average stats
	trip = (double*) hashtbl_get(avetrip,key);
	trip[0] += i;
	trip[1] += j;
	trip[2] += k;
	
	last = *helixid;
      }
    }
    else if (sscanf(tmp,"Structure %d",&i) == 1) {
      if (set->opt->TOPDOWN) {
	if (profile) {
	  process_profile(halfbrac,profile,numhelix,set);
	  if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
	    fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
	    exit(EXIT_FAILURE);
	  }
	} else
	  profile = (int*) malloc(sizeof(int)*ARRAYSIZE);
	numhelix = 0;
      }
      set->opt->NUMSTRUCTS++;
    }
  }
  if (set->opt->TOPDOWN) 
    process_profile(halfbrac,profile,numhelix,set);
  for (i = 1; i < idcount; i++) {
    j = set->helices[i-1]->freq;
    sprintf(key,"%d",i);
    if ((lg = (int*) hashtbl_get(extra,key))) 
      k = j + *lg;
    else
      k = j;
    trip = (double*) hashtbl_get(avetrip,key);
    sprintf(key,"%.1f %.1f %.1f", trip[0]/k,trip[1]/k,trip[2]/k);
    hc = set->helices[i-1];
    hc->avetrip = mystrdup(key);
  }
  if (fclose(fp))
    fprintf(stderr, "File %s not closed successfully\n",set->structfile);
  hashtbl_destroy(extra);
  hashtbl_destroy(avetrip);
}

char* longest_possible(int id,int i,int j,int k,char *seq) {
  int m = 1,*check,*num,diff;
  char *val;

  val = (char*) malloc(sizeof(char)*ARRAYSIZE);
  for (diff = j-i-2*(k+1); diff >= 2 && match(i+k,j-k,seq); diff = j-i-2*(k+1)) k++;
  //if (diff < 2 && match(i+k,j-k)) printf("found overlap for %d %d %d\n",i,j,k+1);
  while (match(i-m,j+m,seq)) m++;
  m--;
  i -= m;
  j += m;
  k+= m;
  sprintf(val,"%d %d %d",i,j,k);

  num = (int*) malloc(sizeof(int));
  *num = id;
  for (m = 0; m < k; m++) {
    sprintf(val,"%d %d",i+m,j-m);
    if ((check = (int*) hashtbl_get(bp,val)))
      printf("%s (id %d) already has id %d\n",val,id,*check);
    hashtbl_insert(bp,val,num);
  }
  sprintf(val,"%d %d %d",i,j,k);
  return val;
}

//finds longest possible helix based on i,j
//populates all bp in longest possible with id in hash
int match(int i,int j,char *seq) {
  char l,r;
  if (i >= j) return 0;
  if (i < 1) return 0;
  if ((unsigned int) j > strlen(seq)) return 0;
  l = seq[i-1];
  r = seq[j-1];
  //printf("l(%d) is %c and r(%d) is %c\n",i,l,j,r);
  if ((l == 'a' || l == 'A') && (r == 'u' || r == 'U' || r == 't' || r == 'T'))
    return 1;
  else if ((l == 'u' || l == 'U' || l == 't' || l == 'T') && (r == 'a' || r == 'A' || r == 'g' || r == 'G'))
    return 1;
  else if ((l == 'c' || l == 'C') && (r == 'g' || r == 'G'))
    return 1;
  else if ((l == 'g' || l == 'G') && (r == 'c' || r == 'C' || r == 'u' || r == 'U' || r == 't' || r == 'T' ))
    return 1;
  else
    return 0;
}

void addHC(Set *set, HC *hc, int idcount) {
  int i,k;

  if (idcount > ARRAYSIZE*set->hc_size) {
    set->hc_size++;
    k = set->hc_size;
    set->helices = (HC**) realloc(set->helices, sizeof(HC*)*ARRAYSIZE*k);
    /*set->joint = (int**) realloc(set->joint, sizeof(int*)*ARRAYSIZE*k);
    for (i=ARRAYSIZE*(k - 1); i < ARRAYSIZE*k; i++) {
      set->joint[i] = (int*) malloc(sizeof(int)*(k-1-i));
    }
    */
  } 
  
  set->helices[idcount-1] = hc;
  set->hc_num = idcount;
}

void reorder_helices(Set *set) {
  int i,*nw, total, sum=0;
  char *old;
  HC **helices = set->helices;

  total = set->hc_num;
  if (!(translate_hc = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for translate_hc failed");
    exit(EXIT_FAILURE);
  }

  qsort(helices,total,sizeof(HC*),freqcompare);
  for (i = 0; i < total; i++) {
    old = helices[i]->id;
    nw = (int*) malloc(sizeof(int)*ARRAYSIZE);
    *nw = i+1;
    hashtbl_insert(translate_hc,old,nw);
    sprintf(old,"%d",i+1);
    sum += helices[i]->freq;
  }
  set->helsum = sum;
}

//will sort to have descending freq
int freqcompare(const void *v1, const void *v2) {
  HC *h1,*h2;
  h1 = *((HC**)v1);
  int i1 = h1->freq;
  h2 = *((HC**)v2);
  int i2 = h2->freq;
  return (i2 - i1);
}

double set_threshold_entropy(Set *set) {
  int i=0;
  double ent =0,frac,last=0,ave,norm;
  HC **list = set->helices;

  norm = (double)list[i]->freq;
  for (i=0; i < set->hc_num; i++) {
    frac = (double)list[i]->freq/norm;
    ent -= frac*log(frac);
    if (frac != 1)
      ent -= (1-frac)*log(1-frac);
    ave = ent/(double)(i+1);
    if ((ave > last) || (fabs(ave-last) < FLT_EPSILON*2)) {
      last = ave;      
    } 
    else {
      //printf("%f is lower than %f\n",ent/(i+1), last);
      set->num_fhc = i;
      //init_joint(set);
      return (100*(double) list[i-1]->freq/(double) set->opt->NUMSTRUCTS);
    }
  }
  return 0;
}

/*LU array
void init_joint(Set *set) {
  int i,j,k=set->num_fhc;
  
  set->joint = (int**) malloc(sizeof(int*)*(k-1));
  for (i = 0; i < k-1; i++) {
    set->joint[i] = (int*) malloc(sizeof(int)*(k-1-i));
    for (j = 0; j < k-1-i; j++) {
      set->joint[i][j] = 0;
    }
  }
}
*/ 

int compare(const void *v1, const void *v2) {
  return (*(int*)v1 - *(int*)v2);
}

//looks up and prints actual helices for all id's
int print_all_helices(Set *set) {
  HC **helices = set->helices;
  char *val,*trip;
  int i,k,m,total = set->hc_num,target=0,cov=0;

  target = set->opt->COVERAGE*set->helsum;
  for (i = 0; i < total; i++) {
    val = helices[i]->maxtrip;
    m = helices[i]->freq;
    trip = helices[i]->avetrip;
    if (set->opt->VERBOSE) {
      if (val != NULL)
	printf("Helix %d is %s (%s) with freq %d\n",i+1,val,trip,m);
      else
	printf("No entry for %d\n",i);
    }
    if (cov < target) {
      cov += helices[i]->freq;
      k = i;
    }
   }
  printf("%.2f coverage of helices after first %d HC\n",set->opt->COVERAGE,k+1);
  return k+1;
 }

double set_num_fhc(Set *set) {
  int marg;
  double percent;
  
  if (set->opt->NUM_FHC > set->hc_num) {
    printf("Number of requested fhc %d greater than total number of hc %d\n",set->opt->NUM_FHC, set->hc_num);
    set->opt->NUM_FHC = set->hc_num;
  }
  if (set->opt->NUM_FHC > 63) 
    printf("Warning: requested num fhc %d exceeds limit of 63\n",set->opt->NUM_FHC);

  marg = set->helices[set->opt->NUM_FHC-1]->freq;
  percent = ((double) marg)*100.0/((double)set->opt->NUMSTRUCTS);
  return percent;
}
/*
void find_bools(Set *set) {
  int i,j,k,joint;
  double cond1,cond2,freqthresh=0.1,thresh=0.1;
  
  //find 10% thresh
  for (k = set->num_fhc; k < set->hc_num && set->helices[k]->freq > NUMSTRUCTS*freqthresh; k++) ;
  //i only ranges in entropy features, j to 10%
  for (i=0; i < set->num_fhc; i++) {    
    for (j=i+1; j < k; j++) {
      joint = set->joint[i];
      cond1 = (double)set->joint[i][j]/(double)set->helices[i+j+1]->freq;
      cond2 = (double)set->joint[i][j]/(double)set->helices[i]->freq;
      if (set->opt->VERBOSE) {
	printf("p(%d|%d) = %.3f\n",i+1,i+j+2,cond1);
	printf("p(%d|%d) = %.3f\n",i+j+2,i+1,cond2);
      }
      if (cond1 > 1-thresh && cond2 > 1-thresh) {
	printf("%d AND %d\n", i+1,i+j+2);
      } else if (cond1 < thresh && cond2 < thresh)
	printf("%d XOR %d\n",i+1,i+j+2);
    }
  }
}
*/
void print_meta(Set *set) {
  int i,j,k=set->num_fhc-1;
  double cond1,cond2,thresh=0.1;

  for (i=0; i < k; i++) {
    for (j=0; j < k-i; j++) {
      cond1 = (double)set->joint[i][j]/(double)set->helices[i+j+1]->freq;
      cond2 = (double)set->joint[i][j]/(double)set->helices[i]->freq;
      if (set->opt->VERBOSE) {
	printf("p(%d|%d) = %.3f\n",i+1,i+j+2,cond1);
	printf("p(%d|%d) = %.3f\n",i+j+2,i+1,cond2);
      }
      if (cond1 > 1-thresh && cond2 > 1-thresh) {
	printf("%d AND %d\n", i+1,i+j+2);
      } else if (cond1 < thresh && cond2 < thresh)
	printf("%d XOR %d\n",i+1,i+j+2);
    }
  }
}

/*assumes threshold found, now setting actual helices
  calculates conditionals, augments features with appropriate
*/
void find_freq(Set *set) {
  int marg,i,total;
  double percent,cov=0;
  HC *hc;

  total = set->hc_num;
  set->num_fhc = set->hc_num;
  for (i = 0; i < total; i++) {
    hc = set->helices[i];
    marg = hc->freq;
    percent = ((double) marg*100.0)/((double)set->opt->NUMSTRUCTS);
    if (percent >= set->opt->HC_FREQ) {
      if (set->opt->VERBOSE)
	printf("Featured helix %d: %s with freq %d\n",i+1,hc->maxtrip,marg);
      hc->isfreq = 1;      
      hc->binary = 1<<i;
      cov += (double)marg/(double)set->helsum;;
    }
    else {
      set->num_fhc = i;
      i = total;
    }
  } 
  if (set->num_fhc > 63) {
    //fprintf(stderr,"number of helices greater than allowed in find_freq()\n");
    set->num_fhc = 63;
    marg = set->helices[62]->freq;
    printf("Capping at 63 fhc with freq %d\n",marg);
    set->opt->HC_FREQ = ((double) marg)*100.0/((double)set->opt->NUMSTRUCTS);
  }
  printf("Coverage by featured helix classes: %.3f\n",cov);
}

void make_profiles(Set *set) {
  FILE *fp,*file;
  int num=1,*id = 0,i,j,k,last = -1, lastfreq = -1,*profile,totalhc=0,allhelix=0;
  int numhelix = 0,size=1,tripsize = INIT_SIZE,allbp=0,fbp=0,seqsize;
  double coverage=0,bpcov=0;
  char *temp,val[ARRAYSIZE],*trips,*name,*prof,*sub,*delim = (char*)" ,\t\n";
  HASHTBL *halfbrac;

  name = set->structfile;
  profile = (int*) malloc(sizeof(int)*ARRAYSIZE);

  if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
    exit(EXIT_FAILURE);
  }
  if (set->opt->REP_STRUCT) {
    if (!(consensus = hashtbl_create(HASHSIZE,NULL))) {
      fprintf(stderr, "ERROR: hashtbl_create() for consensus failed");
      exit(EXIT_FAILURE);
    }
    trips = (char*) malloc(sizeof(char)*tripsize*ARRAYSIZE);
    trips[0] = '\0';
  }
  fp = fopen(name,"r");
  if (fp == NULL) {
    fprintf(stderr, "can't open %s\n",name);
    return;
  }
  file = fopen("structure.out","w");
  if (file == NULL)
    fprintf(stderr,"Error: can't open structure.out\n");
  fprintf(file,"Processing %s\n",name);
  seqsize = strlen(set->seq)*10;
  temp = (char*) malloc(sizeof(char)*seqsize);
  while (fgets(temp,seqsize,fp)) {
    fprintf(file,"Structure %d: ",num++);
    sub = strtok(temp,delim);
    sub = strtok(NULL,delim);
    while ((sub = strtok(NULL,delim))) {
      i = atoi(sub);
      if ((sub = strtok(NULL,delim)))
        j = atoi(sub);
      else
        fprintf(stderr, "Error in file input format, processing %s\n",sub);
      if ((sub = strtok(NULL,delim)))
        k = atoi(sub);
      else
        fprintf(stderr, "Error in file input format, processing %s\n",sub);
      if (k < set->opt->MIN_HEL_LEN) 
	continue;
      sprintf(val,"%d %d",i,j);
      id = (int*) hashtbl_get(bp,val);
      sprintf(val,"%d",*id);
      totalhc++;
      allhelix++;
      allbp+=k;
      id = (int*) hashtbl_get(translate_hc,val);
      if (!id) fprintf(stderr,"no valid translation hc exists for %s\n",val);
      if (*id != -1 && *id != last) {
	fprintf(file,"%d ",*id);
	if (*id <= set->num_fhc && *id != lastfreq) {
	  numhelix++;
	  fbp+=k;
	  if (numhelix >= ARRAYSIZE*size) 
	    profile = (int*) realloc(profile,sizeof(int)*ARRAYSIZE*++size);
	  profile[numhelix-1] = *id;
	  //calc_joint(set,profile,numhelix);
	  make_brackets(halfbrac,i,j,*id);
	  lastfreq = *id;
	}
	last = *id;
	//printf("assigning %d to %d %d\n",*id,i,j);
      }
      if (set->opt->REP_STRUCT) {
	sprintf(val,"%d %d %d ",i,j,k);
	while (strlen(trips)+strlen(val) > (unsigned int)(ARRAYSIZE*tripsize-1))
	  trips = (char*) realloc(trips,sizeof(char)*++tripsize*ARRAYSIZE);
	strcat(trips,val);
      }
    }
    prof = process_profile(halfbrac,profile,numhelix,set);
    //printf("processing %d with profile %s\n",num,prof);
    fprintf(file,"\n\t-> %s\n",prof);
     
    if (set->opt->REP_STRUCT) {
      make_rep_struct(consensus,prof,trips);
    }
     
    if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
      fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
      exit(EXIT_FAILURE);
    }
    last = 0;
    lastfreq = 0;
    coverage += ((double)numhelix/(double)allhelix);
    bpcov += ((double)fbp/(double)allbp);
    numhelix = 0;
    allhelix = 0;
    if (set->opt->REP_STRUCT)
      trips[0] = '\0';
  }
  if (set->opt->VERBOSE) {
  printf("Ave number of HC per structure: %.1f\n",(double)totalhc/(double) set->opt->NUMSTRUCTS);
  printf("Ave structure coverage by fhc: %.2f\n",coverage/(double)set->opt->NUMSTRUCTS);
  printf("Ave structure coverage by fbp: %.2f\n",bpcov/(double)set->opt->NUMSTRUCTS);
	}
  free(profile);
  fclose(fp);
  fclose(file);
}

/*calculates joint for last value in prof vs everything else
Moved up to process_structs

void calc_joint(Set *set, int *prof, int num) {
  int k;  

  if (num < 2)
    return;
  for (k = 0; k < num-1; k++) {
    if (prof[k] < prof[num-1])
      set->joint[prof[k]-1][prof[num-1]-prof[k]-1]++;
    else
      set->joint[prof[num-1]-1][prof[k]-prof[num-1]-1]++;
    //if (j)
  }
  return;
}
*/

void make_profiles_sfold(Set *set) {
  FILE *fp,*file;
  int num=0,*id = 0,i,j,k,last = -1, lastfreq = -1,*profile,totalhc=0,allhelix=0;
  int numhelix = 0,size=1,tripsize = INIT_SIZE,allbp=0,fbp=0;
  double coverage=0,bpcov=0;
  char temp[100],val[ARRAYSIZE],*trips,*name,*prof;
  HASHTBL *halfbrac;

  name = set->structfile;
  profile = (int*) malloc(sizeof(int)*ARRAYSIZE);

  if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
    exit(EXIT_FAILURE);
  }
  if (set->opt->REP_STRUCT) {
    if (!(consensus = hashtbl_create(HASHSIZE,NULL))) {
      fprintf(stderr, "ERROR: hashtbl_create() for consensus failed");
      exit(EXIT_FAILURE);
    }
    trips = (char*) malloc(sizeof(char)*tripsize*ARRAYSIZE);
    trips[0] = '\0';
  }
  fp = fopen(name,"r");
  if (fp == NULL) {
    fprintf(stderr, "can't open %s\n",name);
    return;
  }
  file = fopen("structure.out","w");
  if (file == NULL)
    fprintf(stderr,"Error: can't open structure.out\n");
  fprintf(file,"Processing %s\n",name);
  while (fgets(temp,100,fp) != NULL) {
    if (sscanf(temp,"Structure %d",&num) == 1) {
     if (last == -1) {
       fprintf(file,"Structure %d: ",num);	
       continue;
     }
     prof = process_profile(halfbrac,profile,numhelix,set);
     //printf("processing %d with profile %s\n",num,prof);
     fprintf(file,"\n\t-> %s\nStructure %d: ",prof,num);
     
     if (set->opt->REP_STRUCT) {
       make_rep_struct(consensus,prof,trips);
     }
     
     if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
       fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
       exit(EXIT_FAILURE);
     }
     last = 0;
     lastfreq = 0;
     coverage += ((double)numhelix/(double)allhelix);
     bpcov += ((double)fbp/(double)allbp);
     numhelix = 0;
     allhelix = 0;
     if (set->opt->REP_STRUCT)
       trips[0] = '\0';
    } 
    else if (sscanf(temp,"%d %d %d",&i,&j,&k) == 3) {
      if (k < set->opt->MIN_HEL_LEN)
	continue;
      sprintf(val,"%d %d",i,j);
      id = (int*) hashtbl_get(bp,val);
      sprintf(val,"%d",*id);
      totalhc++;
      allhelix++;
      allbp+=k;
      id = (int*) hashtbl_get(translate_hc,val);
      if (!id) fprintf(stderr,"no valid translation hc exists for %s\n",val);
      if (*id != -1 && *id != last) {
	fprintf(file,"%d ",*id);
	if (*id <= set->num_fhc && *id != lastfreq) {
	  numhelix++;
	  fbp+=k;
	  if (numhelix >= ARRAYSIZE*size) 
	    profile = (int*) realloc(profile,sizeof(int)*ARRAYSIZE*++size);
	  profile[numhelix-1] = *id;
	  make_brackets(halfbrac,i,j,*id);
	  lastfreq = *id;
	}
	last = *id;
	//printf("assigning %d to %d %d\n",*id,i,j);
      }
      if (set->opt->REP_STRUCT) {
	sprintf(val,"%d %d %d ",i,j,k);
	while (strlen(trips)+strlen(val) > (unsigned int)(ARRAYSIZE*tripsize-1))
	  trips = (char*) realloc(trips,sizeof(char)*++tripsize*ARRAYSIZE);
	strcat(trips,val);
      }
    }
  }
  fprintf(file,"\n\t-> %s ",prof);
  process_profile(halfbrac,profile,numhelix,set);
 //fprintf(file,"Structure %d: %s\n",num,profile);
  printf("Ave number of HC per structure: %.1f\n",(double)totalhc/(double) set->opt->NUMSTRUCTS);
  printf("Ave structure coverage by fhc: %.2f\n",coverage/(double)set->opt->NUMSTRUCTS);
  printf("Ave structure coverage by fbp: %.2f\n",bpcov/(double)set->opt->NUMSTRUCTS);
  free(profile);
  fclose(fp);
  fclose(file);
}

char* process_profile(HASHTBL *halfbrac,int *profile,int numhelix,Set *set) {
  int i,size=INIT_SIZE;
  char val[ARRAYSIZE],*dup;
  Profile **profiles;  

  profiles = set->profiles;
  qsort(profile,numhelix,sizeof(int),compare);

  dup = (char*) malloc(sizeof(char)*ARRAYSIZE*size);
  dup[0] = '\0';
  for (i = 0; i < numhelix; i++) {
    sprintf(val,"%d ",profile[i]);
    if (strlen(dup) + strlen(val) >= (unsigned int) ARRAYSIZE*size-1) {
      dup = (char*) realloc(dup,sizeof(char)*ARRAYSIZE*++size);
    }
    dup = strcat(dup,val);
  }
  for (i = 0; i < set->prof_num; i++) {
    if (!strcmp(profiles[i]->profile,dup)) {
      profiles[i]->freq++;
      break;
    }
  }
  if (i == set->prof_num) {
    if (i >= ARRAYSIZE*set->prof_size) {
      set->prof_size++;
      profiles = (Profile**) realloc(profiles,sizeof(Profile*)*ARRAYSIZE*set->prof_size);
      set->profiles = profiles;
    }
    profiles[i] = create_profile(dup);
    //printf("adding profile[%d] %s\n",i,profiles[i]->profile);
    set->prof_num++;
    make_bracket_rep(halfbrac,profiles[i]);
  }
  hashtbl_destroy(halfbrac);
  return dup;
}

//makes the bracket representation of dup, using values in hashtbl brac
//dup is a (mod) profile in graph
//called by process_profile()
void make_bracket_rep(HASHTBL *brac,Profile *prof) {
  int num,*array,k=0,size = INIT_SIZE,total;
  char *profile,*val;
  KEY *node = NULL;

  num = hashtbl_numkeys(brac);
  array = (int*) malloc(sizeof(int)*num);
  for (node = hashtbl_getkeys(brac); node; node=node->next) 
    array[k++] = atoi(node->data);
  //sort by i,j position  
  qsort(array,num,sizeof(int),compare);
  profile = (char*) malloc(sizeof(char)*ARRAYSIZE*size);
  profile[0] = '\0';
  val = (char*) malloc(sizeof(char)*ARRAYSIZE);
  for (k = 0; k < num; k++) {
    sprintf(val,"%d",array[k]);
    val = (char*) hashtbl_get(brac,val);
    if ((total = strlen(profile)+strlen(val)) > ARRAYSIZE*size-1)
      profile = (char*) realloc(profile,sizeof(char)*ARRAYSIZE*++size);
    strcat(profile,val);
  }
  prof->bracket = profile;
  //free(val);
  free(array);
}

//inserts bracket representation for i,j into a hash
void make_brackets(HASHTBL *brac, int i, int j, int id) {
  char key[ARRAYSIZE],*val;

  sprintf(key,"%d",i);
  val = (char*) malloc(sizeof(char)*ARRAYSIZE);
  sprintf(val,"[%d",id);
  //  printf("making bracket %s for %d\n",val,i);
  hashtbl_insert(brac,key,val);
  sprintf(key,"%d",j);
  val = (char*) malloc(sizeof(char)*2);
  val[0] = ']';
  val[1] = '\0';
  hashtbl_insert(brac,key,val);
}

void make_rep_struct(HASHTBL *consensus,char *profile, char* trips) {
  int *bpfreq,i,j,k;
  char *val,*blank = (char*)" ",bpair[ARRAYSIZE];
  HASHTBL *ij;
  
  ij = (HASHTBL*) hashtbl_get(consensus,profile);
  if (!ij) {
    if (!(ij = hashtbl_create(HASHSIZE,NULL))) {
      fprintf(stderr, "ERROR: hashtbl_create() for ij failed");
      exit(EXIT_FAILURE);
    }
    hashtbl_insert(consensus,profile,ij);
  }
  for (val = strtok(trips,blank); val; val = strtok(NULL,blank)) {
    i = atoi(val);
    j = atoi(strtok(NULL,blank));
    k = atoi(strtok(NULL,blank));
    for (k--; k >= 0; k--) {
      sprintf(bpair,"%d %d",i+k,j-k);
      bpfreq = (int*)hashtbl_get(ij,bpair);
      if (bpfreq)
	(*bpfreq)++;
      else {
	bpfreq = (int*) malloc(sizeof(int));
	*bpfreq = 1;
	hashtbl_insert(ij,bpair,bpfreq);
      }
      //printf("in rep struct for %s, inserting %d %d\n",profile,i+k,j-k);      
    }
  }
}

void print_profiles(Set *set) {
  int i;
  Profile **profiles = set->profiles;
  
  find_general_freq(set);
  qsort(profiles,set->prof_num,sizeof(Profile*),profsort); 
  for (i = 0; i < set->prof_num; i++) 
    if (set->opt->VERBOSE)
      printf("Profile %s with freq %d (%d)\n",profiles[i]->profile,profiles[i]->freq,profiles[i]->genfreq);
}

int profsort(const void *v1, const void *v2) {
  Profile *p1,*p2;
  p1 = *((Profile**)v1);
  p2 = *((Profile**)v2);
  //return (p2->genfreq + p2->freq - p1->genfreq - p1->freq);
  //return (p2->genfreq - p1->genfreq);
  return (p2->freq - p1->freq);
}

double set_num_sprof(Set *set) {
  int marg;

  if (set->opt->NUM_SPROF > set->prof_num) {
    printf("Number of requested sprof %d greater than total number of hc %d\n",set->opt->NUM_SPROF, set->prof_num);
    set->opt->NUM_SPROF = set->prof_num;
  }
  /*
  for (k = set->opt->NUM_SPROF; set->profiles[set->opt->NUM_SPROF-1]->freq == set->profiles[k]->freq; k++) 
    if (k == set->prof_num -1)
      break;
  if (k != set->opt->NUM_SPROF) {
    if (set->opt->VERBOSE)
      printf("Increasing number of sprof to %d\n",k);
    set->opt->NUM_SPROF = k;
  }
  */
  set->num_sprof = set->opt->NUM_SPROF;
  marg = set->profiles[set->opt->NUM_SPROF-1]->freq;
  return (100*(double) marg/(double) set->opt->NUMSTRUCTS);
}

/*If all profiles have frequency of 1, we default to displaying first 10
 */
double set_p_threshold_entropy(Set *set) {
  int i=0;
  double ent =0,frac,last=0,norm,ave;
  Profile **list = set->profiles;
  
  norm = (double)list[0]->freq;
  //  find_general_freq(set);
  if (norm == 1) {
    set->num_sprof = 10;
    return -2;
  }
  for (i=0; i < set->prof_num; i++) {
    if (list[i]->freq == 1) {
      set->num_sprof = i;
      return (100*(double) list[i-1]->freq/(double) set->opt->NUMSTRUCTS);
    }
    frac = (double)list[i]->freq/norm;
    ent -= frac*log(frac);
    if (frac != 1)
      ent -= (1-frac)*log(1-frac);
    ave = ent/(double)(i+1);
    if ((ave > last) || (fabs(ave-last) < FLT_EPSILON*2)) {
      last = ave;      
    } 
    else {
      set->num_sprof = i;
      return (100*(double) list[i-1]->freq/(double) set->opt->NUMSTRUCTS);
    }
  }
  set->num_sprof = i;
  return (100*(double) list[i-1]->freq/(double) set->opt->NUMSTRUCTS);
}

void find_general_freq(Set *set) {
  int i,j,val;

  //  genfreq = (int*)malloc(sizeof(int)*set->prof_num);
  for (i = 0; i < set->prof_num-1; i++) {
    set->profiles[i]->genfreq += set->profiles[i]->freq;
    for (j = i+1; j < set->prof_num; j++) {
      val = subset(set,set->profiles[i]->profile,set->profiles[j]->profile); 
      if (val == 1)
	set->profiles[i]->genfreq += set->profiles[j]->freq;
      else if (val == 2)
	set->profiles[j]->genfreq += set->profiles[i]->freq;
    }
  }
}

/*tests whether one profile is a subset of profile one
returns 1 if first profile is a subset, 2 if second is ,0 if none are*/
int subset(Set *set,char *one, char *two) {
  unsigned long rep1,rep2;

  rep1 = binary_rep(set,one);
  rep2 = binary_rep(set,two);
  if ((rep1 & rep2) == rep1)
    return 1;
  else if ((rep1 & rep2) == rep2)
    return 2;

  return 0;
}

unsigned long binary_rep(Set *set,char *profile) {
  int i;
  unsigned long sum = 0;
  char *copy = mystrdup(profile),*helix;

  for (helix = strtok(copy," "); helix; helix = strtok(NULL," ")) {
    for (i = 0; i < set->num_fhc; i++) 
      if (!strcmp(set->helices[i]->id,helix))
	break;
    sum += set->helices[i]->binary;
  }
  free(copy);
  return sum;
}

double set_p_threshold(Set *set, int start) {
  int i=0,freq_target,ave = 0,diff,index,dropoff=0,partial=0,total;
  double frac=0;
  //where you start in drop off calc has to cover at least 50% of area
  //double coverage = 50.0;
  Profile **list = set->profiles;
  
  total = set->prof_num;
  //translate start percentage to actual num of elements
  freq_target = start*set->opt->NUMSTRUCTS/100;
  /*
  profiles = malloc(sizeof(int)*total);
  for (i = 0; i < total; i++) {
    profiles[i] = list[i]->freq;
    //printf("profiles[%d] = %d\n",i,profiles[i]);
    sum += profiles[i];
  }
  printf("sum is %d, freq target %d\n",sum,freq_target);
  qsort(profiles,total,sizeof(int),compare);
  */

  //altthres brings you to freq target first
  if (set->opt->ALTTHRESH) {
    for (i=0; i<total && list[i]->freq > freq_target; i++) {
      partial += list[i]->freq;
    }
    //if everything above start, return the start
    if (i == total) {
      set->num_sprof = i;
      return (double) start;
    }
  }
  //if stopping at freq_target results in less than 50% coverage, go to at least 50%
  frac = ((double) partial)/((double)set->opt->NUMSTRUCTS);
  if (frac < set->opt->COVERAGE) {
    //printf("not up to coverage with %d/%d\n",partial,sum);
    while (frac < set->opt->COVERAGE && i < total) {
      partial += list[i++]->freq;
      frac = ((double) partial)/((double)set->opt->NUMSTRUCTS);
    }
    //printf("now at %d/%d = %.1f at helices[%d] = %d\n",partial,sum,frac,i+1,helices[i+1]);
  }
  if (i == 0) 
    i++;
  index = i;
  for ( ; i < total-1; i++) {
    ave = (list[i+1]->freq+list[i-1]->freq)/2;
    diff = ave-list[i]->freq;
    //printf("ave is %d and diff is %d for %d\n",ave,diff,i);
    if (diff > dropoff) {
      //printf("bumping off %d for %d\n",dropoff,diff);
      dropoff = diff;
      index = i;
    }
  }
  set->num_sprof = index;
  return (100*(double) list[index-1]->freq/(double) set->opt->NUMSTRUCTS);
}

void select_profiles(Set *set) {
  int i,coverage=0,cov=0,target;
  //double percent;
  Profile *prof;

  if (set->num_sprof == 0) 
	set->num_sprof = set->prof_num;
  target = set->opt->COVERAGE*set->opt->NUMSTRUCTS;
  for (i = 0; i < set->num_sprof; i++) {
    prof = set->profiles[i];
    //percent = ((double) prof->freq)*100.0 / ((double)set->opt->NUMSTRUCTS);
    //printf("%s has percent %.1f\n",node->data,percent);
    if (((double)prof->freq*100.0/((double)set->opt->NUMSTRUCTS)) < set->opt->PROF_FREQ) {
	set->num_sprof = i;
	break;
    }
    prof->selected = 1;
    coverage += prof->freq;
    if (coverage < target)
      cov = i+1;
    if (set->opt->VERBOSE)
      printf("Selected profile %swith freq %d (%d)\n",prof->profile,prof->freq,prof->genfreq);
  }
  printf("Coverage by selected profiles: %.3f\n",(double)coverage/(double)set->opt->NUMSTRUCTS);
  if (coverage < target) {
    while (coverage < target) {
      coverage += set->profiles[i++]->freq;
      if (i >= set->prof_num)
	fprintf(stderr,"in select profiles() exceeding number of profiles\n");
    }
    cov = i;
  }
  printf("Number of profiles needed to get %.2f coverage: %d\n",set->opt->COVERAGE, cov+1);
  //printf("Number of structures represented by top %d profiles: %d/%d\n",stop,cov15,set->opt->NUMSTRUCTS);
}

void process_one_input(Set *set) {
  HASHTBL *halfbrac;
  FILE *file;
  char temp[100],tmp[ARRAYSIZE],*profile=NULL,*fullprofile=NULL,*diff=NULL,*native,*diffn=NULL;
  int i,j,k,*id,last=0,*prof,hc_num;
  int numhelix = 0,fullnum = 0,natnum = 0;
  int size = INIT_SIZE;
  Profile *natprof;
  node *inode;
  
  if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
    exit(EXIT_FAILURE);
  }
  prof = (int*) malloc(sizeof(int)*ARRAYSIZE*size);
  if (!(file = fopen(set->opt->INPUT,"r")))
    fprintf(stderr,"Cannot open %s\n",set->opt->INPUT);
  hc_num = set->hc_num;
  while (fgets(temp,100,file)) {
    //    if (sscanf(temp,"Structure %d (%d)",&i,&prob) == 2) 
    if (sscanf(temp,"%d %d %d",&i,&j,&k) == 3) {
      sprintf(tmp,"%d %d",i,j);
      id = (int*) hashtbl_get(bp,tmp);
      if (!id) {
	id = process_native(set,i,j,k);
      } else {
	sprintf(tmp,"%d",*id);
	id = (int*) hashtbl_get(translate_hc,tmp);
      }
      printf("found %d for %d %d %d\n",*id,i,j,k);
      if (*id != last) {
	if (natnum >= ARRAYSIZE*size) 
	  prof = (int*) realloc(prof,sizeof(int)*ARRAYSIZE*++size);
	prof[natnum++] = *id;
	last = *id;
	make_brackets(halfbrac,i,j,*id);
      }
      else if (set->opt->VERBOSE)
	printf("helix %d %d %d is duplicate with id %d: process_input()\n",i,j,k,*id);
    }
  }
  set->hc_num = hc_num;
  qsort(prof,natnum,sizeof(int),compare);

  native = (char*) malloc(sizeof(char)*ARRAYSIZE*++size);
  native[0] = '\0';
  numhelix = natnum;
  fullnum = natnum;
  for (i = 0; i<natnum; i++) {
    if (prof[i] > set->num_fhc && !profile) {
      profile = mystrdup(native);
      numhelix = i;
      native[0] = '\0';
    }
    if (prof[i] > set->hc_num && !diffn) {
      diffn = mystrdup(native);
      fullnum = i;
      native[0]='\0';
    }
    sprintf(tmp,"%d ",prof[i]);
    if (strlen(tmp)+strlen(native)+1 > (unsigned int) ARRAYSIZE*size)
      native = (char*) realloc(native,sizeof(char)*ARRAYSIZE*++size);
    native = strcat(native,tmp);
  }
  diff = mystrdup(native);
  if (!profile) {
    profile = (char*) malloc(sizeof(char));
    profile = (char*)"";
  }
  if (!diffn) {
    diffn = diff;
    diff = (char*)"";
  }

  //printf("native %s%s%s(%d), fullprofile %s%s(%d), profile %s(%d)\n",profile,diffn,diff,natnum,profile,diffn,fullnum,profile,numhelix);
  
  /*making data structure*/
  while (strlen(profile)+strlen(diffn)+strlen(diff)+1 > (unsigned int) ARRAYSIZE*size)
    native = (char*) realloc(native,sizeof(char)*ARRAYSIZE*++size);
  sprintf(native,"%s%s%s",profile,diffn,diff);
  fullprofile = (char*) malloc(strlen(native));
  sprintf(fullprofile,"%s%s",profile,diffn);

  //printf("now native %s, full profile %s\n",native,fullprofile);
  natprof = create_profile(native);
  make_bracket_rep(halfbrac,natprof);
  hashtbl_destroy(halfbrac);

  set->inputnode = createNode(native);
  set->inputnode->bracket = natprof->bracket;

  if (fullnum != natnum) {
    //printf("creating input node %s with diff %s\n",fullprofile,diff);
    inode = createNode(fullprofile);
    inode->neighbors = (node**) malloc(sizeof(node*)*ARRAYSIZE);
    inode->neighbors[0] = set->inputnode;
    inode->diff = (char**) malloc(sizeof(char*)*ARRAYSIZE);
    inode->diff[0] = diff;
    inode->numNeighbors++;
    set->inputnode = inode;
  }
  if (numhelix != fullnum) {
    //printf("creating input node %s with diffn %s\n",profile,diffn);
    inode = createNode(profile);
    inode->neighbors = (node**) malloc(sizeof(node*)*ARRAYSIZE);
    inode->neighbors[0] = set->inputnode;
    inode->diff = (char**) malloc(sizeof(char*)*ARRAYSIZE);
    inode->diff[0] = diffn;
    inode->numNeighbors++;
    set->inputnode = inode;
  }
}

//processes native helices if necessary; called by process_input
//returns id for i,j,k helix
int* process_native(Set *set,int i, int j, int k) {
  int *id = NULL,l;
  char *tmp,key[ARRAYSIZE];

  tmp = (char*) malloc(sizeof(char)*ARRAYSIZE);
  for (l=1; l < k; l++) {
    sprintf(tmp,"%d %d",i+l,j-l);
    id = (int*) hashtbl_get(bp,tmp);
    if (id) {
      sprintf(tmp,"%d",*id);
      id = (int*) hashtbl_get(translate_hc,tmp);
      for (l-- ; l >= 0; l--) {
	sprintf(tmp,"%d %d",i+l,j+l);
	hashtbl_insert(bp,tmp,id);
      }
      return id;
    }
  }

  id = (int*) malloc(sizeof(int));
  *id = ++(set->hc_num);
  sprintf(key,"%d %d",i,j);
  hashtbl_insert(bp,key,id);
  return id;
}

void find_consensus(Set *set) {
  int freq,*bpfreq,i;
  KEY *node,*bpnode;
  HASHTBL *ij,*final;

  if (!(final = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for final failed");
    exit(EXIT_FAILURE);
  }

  //cycle thru all profiles, either calc consensus for selected or destroying
  for (node = hashtbl_getkeys(consensus); node; node = node->next) {
    freq = 0;
    for (i = 0; i < set->num_sprof; i++) 
      if (!strcmp(set->profiles[i]->profile,node->data))
	freq = set->profiles[i]->freq;
    ij = (HASHTBL*) hashtbl_get(consensus,node->data);
    if (freq) {
      if (!ij)
	fprintf(stderr, "ij not found in find_consensus()\n");
      for (bpnode = hashtbl_getkeys(ij); bpnode; bpnode = bpnode->next) {
	bpfreq = (int*) hashtbl_get(ij,bpnode->data);
	if (*bpfreq*100/freq <= 50)
	  hashtbl_remove(ij,bpnode->data);
      }
    } else {
      hashtbl_destroy(ij);
      //insert dummy pointer so remove won't seg fault
      hashtbl_insert(consensus,node->data, malloc(sizeof(char)));
      hashtbl_remove(consensus,node->data);
    }
  }
}

//in ct format
int print_consensus(Set *set) {
  int l,k=0,m,seqlen;
  char outfile[ARRAYSIZE],key[ARRAYSIZE],*pair,*i,*j,*blank = (char*)" ";
  KEY *bpnode;
  HASHTBL *bpairs,*temp;
  FILE *fp;

  seqlen = strlen(set->seq);

  //foreach profile
  //for (node = hashtbl_getkeys(consensus); node; node = node->next) {
  for (l = 0; l < set->num_sprof; l++) {
    if (!(temp = hashtbl_create(HASHSIZE,NULL))) {
      fprintf(stderr, "ERROR: hashtbl_create() for temp failed");
      exit(EXIT_FAILURE);
    }
    sprintf(outfile,"Structure_%d.ct",++k);
    fp = fopen(outfile,"w");

    fprintf(fp,"Profile: %s\n",set->profiles[l]->profile);
    fprintf(fp,"Freq: %d\n",set->profiles[l]->freq);
    fprintf(fp,"%d dG = n/a\n",seqlen);
    bpairs = (HASHTBL*) hashtbl_get(consensus,set->profiles[l]->profile);
    if (!bpairs)
      fprintf(stderr,"no bpairs found\n");
    for (bpnode = hashtbl_getkeys(bpairs); bpnode; bpnode = bpnode->next) {
      pair = mystrdup(bpnode->data);
      //printf("processing pair %s\n",pair);
      i = strtok(pair,blank);
      j = strtok(NULL,blank);
      hashtbl_insert(temp,i,j);
      hashtbl_insert(temp,j,i);

    }
    for (m = 0; m < seqlen; m++) {
      sprintf(key,"%d",m+1);
      if ((j = (char*) hashtbl_get(temp,key)) )
	fprintf(fp,"\t%d %c\t%d   %d   %s   %d\n",m+1,set->seq[m],m,m+2,j,m+1);
      else
	fprintf(fp,"\t%d %c\t%d   %d   0   %d\n",m+1,set->seq[m],m,m+2,m+1);
    }
    hashtbl_destroy(temp);
    fclose(fp);
  }
  hashtbl_destroy(consensus);
  return 0;
}

void free_Set(Set* set) {
  free(set);
}
