/* Part of cluster graph
Program assigns all helices characterized by a maximal helix an ID and counts them
Finds all freq helices and clusters 1000 structs based on them
*/

#include "hashtbl.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "profile.h"
#include <math.h>

static HASHTBL *triplet;
static HASHTBL *common;
static HASHTBL *consensus;
static HASHTBL *translate_hc;
static char *seq;

/*
Input: the fasta file sfold ran on. First line must be description, rest of lines is sequence
Does: Stores nucleotide at position i in array
Output: string containing sequence
Used to find longest possible
*/
char* input_seq(char *seqfile) {
FILE * fp;
  int size = 7,fasta = 0;
  char temp[100],*blank = " \n",*part,*final;

  fp = fopen(seqfile,"r");
  if (fp == NULL) {
    fprintf(stderr, "can't open %s\n",seqfile);
    return 0;
  }
  final = malloc(sizeof(char)*ARRAYSIZE*size);
  final[0] = '\0';
  while (fgets(temp,100,fp)) {    
    //put error handling in case first line has more than 100 chars
    //skipping first line if fasta header and not seq
    if (temp[0] == '>' || fasta) {
      //printf("found fasta header %s",temp);
      if (strlen(temp) < 99 || (strlen(temp) == 99 && temp[98] == '\n'))
	continue;
      else 
	fasta = 1;
      continue;
    } 
    for (part = strtok(temp,blank); part; part = strtok(NULL,blank)) {
      if (strlen(final)+strlen(part) > ARRAYSIZE*size-1) {
	while (++size*ARRAYSIZE - 1 < strlen(final)+strlen(part)) ;
	final = realloc(final,ARRAYSIZE*size);
      }
      final = strcat(final,part);
      //printf("adding %s of length %d for total length of %d\n", part,strlen(part),strlen(final));
    }
  }
  if (VERBOSE) 
    printf("seq in %s is %s with length %d\n",seqfile,final,strlen(final));
  //printf("final char is %c\n",final[strlen(final)-1]);
  return final;
}

/*
name is the file name containing structures to process
marginals is the hash: id->frequency of id
return idcount = num of helix equivalence classes + 1
*/
int process_structs(char *seqfile,char *name) {
  FILE *fp;
  int i,j,k,*helixid,idcount=1,*lg,*id=NULL,last = 0, toosmall = 0;
  double *trip;
  char tmp[100],key[ARRAYSIZE],dbl[ARRAYSIZE];
  HASHTBL *temp,*extra;

  fp = fopen(name,"r");
  if (fp == NULL) {
    fprintf(stderr, "can't open %s\n",name);
    return 0;
  }
  seq = input_seq(seqfile);
  if (STATS) 
    triplet = hashtbl_create(HASHSIZE,NULL);

  if (!(avetrip = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for avetrip failed");
    exit(EXIT_FAILURE);
  }
  if (!(extra = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for extra failed");
    exit(EXIT_FAILURE);
  }

  while (fgets(tmp,100,fp) != NULL) {
    if (sscanf(tmp,"%d %d %d",&i,&j,&k) == 3) {
      if (i == 0) continue;
      //if (longest < k) 
      //longest = k;
      if (k < MIN_HEL_LEN) {
	toosmall++;
	//printf("too small %d %d %d\n",i,j,k);
	continue;
      }
      sprintf(dbl,"%d %d",i,j);
      helixid = hashtbl_get(bp,dbl);
      //printf("%d %d %d\n",i,j,k);
      if (!helixid) {
	//printf("%s not found in bp\n",dbl);
	longest_possible(i,j,k,idcount);
	//increment marginal
	sprintf(key,"%d",idcount++);
	id = malloc(sizeof(int));
	*id = 1;
	hashtbl_insert(marginals,key,id);
	//printf("inserting %s for %d %d %d\n",key,i,j,k);
	//triplet stats
	trip = malloc(sizeof(double)*3);
	trip[0] = i;
	trip[1] = j;
	trip[2] = k;
	hashtbl_insert(avetrip,key,trip);
	if (STATS) {
	  temp = hashtbl_create(HASHSIZE,NULL);
	  sprintf(dbl,"%d %d",i,k);
	  lg = malloc(sizeof(int));
	  *lg = 1;
	  hashtbl_insert(temp,dbl,lg);
	  hashtbl_insert(triplet,key,temp);
	}
	last = idcount - 1;
      }
      else {
	//printf("Found %d %d with id %d\n",i,j,helixid);
	sprintf(key,"%d",*helixid);
	//increment marginal
	if (last != *helixid) {
	  id = hashtbl_get(marginals,key);
	  ++*id;
	} else {
	  if ((lg = hashtbl_get(extra,key)))
	    ++*lg;
	  else {
	    lg = malloc(sizeof(int));
	    *lg = 1;
	    hashtbl_insert(extra,key,lg);
	  }
	  //if (VERBOSE) 
	  //printf("Found repeat id %d:%s\n",last,hashtbl_get(idhash,key));
	}
	//triplet stats
	trip = hashtbl_get(avetrip,key);
	trip[0] += i;
	trip[1] += j;
	trip[2] += k;
	//	if (*helixid == 2)
	//printf("%.1f ",trip[1]);
	if (STATS) {
	  temp = hashtbl_get(triplet,key);
	  sprintf(dbl,"%d %d",i,k);
	  lg = hashtbl_get(temp,dbl);
	  if (lg)
	    (*lg)++;
	  else {
	    lg = malloc(sizeof(int));
	    *lg = 1;
	    hashtbl_insert(temp,dbl,lg);
	  }
	}
	last = *helixid;
      }
    }
    else if (sscanf(tmp,"Structure %d",&i) == 1)
      NUMSTRUCTS++;
  }

  for (i = 1; i < idcount; i++) {
    sprintf(key,"%d",i);
    id = hashtbl_get(marginals,key);
    if ((lg = hashtbl_get(extra,key))) 
      k = *id + *lg;
    else
      k = *id;
    trip = hashtbl_get(avetrip,key);
    for (j = 0; j < 3; j++)
      //printf(" %.1f", trip[j]);
      trip[j] /= k;
  }
  if (fclose(fp))
    fprintf(stderr, "File %s not closed successfully\n",name);
  hashtbl_destroy(extra);
  return idcount;   //idcount augmented, so can use < in loops
}


//finds longest possible helix based on i,j
//populates all bp in longest possible with id in hash
void longest_possible(int i,int j,int k,int id) {
  int m = 1,*check,*num,diff;
  char val[ARRAYSIZE],key[ARRAYSIZE];
  //printf("for %d %d %d\n",i,j,k);
  for (diff = j-i-2*(k+1); diff >= 2 && match(i+k,j-k); diff = j-i-2*(k+1)) k++;
  //if (diff < 2 && match(i+k,j-k)) printf("found overlap for %d %d %d\n",i,j,k+1);
  while (match(i-m,j+m)) m++;
  m--;
  i -= m;
  j += m;
  k+= m;
  sprintf(val,"%d %d %d",i,j,k);
  sprintf(key,"%d",id);
  hashtbl_insert(idhash,key,mystrdup(val));
  //printf("inserting %s -> %s into idhash\n",key,val);

  num = malloc(sizeof(int));
  *num = id;

  for (m = 0; m < k; m++) {
    sprintf(val,"%d %d",i+m,j-m);
    if ((check = hashtbl_get(bp,val)))
      printf("%s (id %d) already has id %d\n",val,id,*check);

    hashtbl_insert(bp,val,num);
    //printf("inserting %s\n",val);
  }
}

int match(int i,int j) {
  char l,r;
  if (i >= j) return 0;
  if (i < 1) return 0;
  if (j > strlen(seq)) return 0;
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

void reorder_helices(int total) {
  int i,*new,*val;
  char *old,**helices,*key;
  HASHTBL *new_marg, *new_id, *new_ave,*temp;

  if (!(translate_hc = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for translate_hc failed");
    exit(EXIT_FAILURE);
  }
  if (!(new_marg = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for new_marg failed");
    exit(EXIT_FAILURE);
  }
  if (!(new_id = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for new_id failed");
    exit(EXIT_FAILURE);
  }
  if (!(new_ave = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for new_ave failed");
    exit(EXIT_FAILURE);
  }
  
  if (total > 99998)
    fprintf(stderr,"Need to allocate more space in reorder_helices(): error\n");
  helices = malloc(sizeof(char*)*total);
  key = malloc(sizeof(char)*ARRAYSIZE);
  for (i = 1; i< total; i++) {
    old = malloc(sizeof(char)*5);
    sprintf(old,"%d",i);
    helices[i-1] = old;
  }
  qsort(helices,total-1,sizeof(char*),freqcompare);
  for (i = 1; i < total; i++) {
    old = helices[i-1];
    //printf("helix %s becomes %d with freq %d\n",old,i,*((int*)hashtbl_get(marginals,old)));
    new = malloc(sizeof(int));
    *new = i;
    hashtbl_insert(translate_hc,old,new);
    sprintf(key,"%s",old);
    sprintf(old,"%d",i);
    hashtbl_insert(new_marg,old,hashtbl_get(marginals,key));
    hashtbl_insert(new_id,old,hashtbl_get(idhash,key));
    hashtbl_insert(new_ave,old,hashtbl_get(avetrip,key));
  }
  temp = marginals;
  marginals = new_marg;
  hashtbl_destroy(temp);
  temp = idhash;
  idhash = new_id;
  hashtbl_destroy(temp);
  temp = avetrip;
  avetrip = new_ave;
  hashtbl_destroy(temp);
  free(helices);
}

//current algorithm
double set_h_dropoff(HASHTBL *hash, int start) {
  int *helices,i,sum=0,freq_target,ave = 0,*dropoffs,diff,*val,total,partial=0;
  char key[ARRAYSIZE];
  double h, frac;
  //where you start in drop off calc has to cover at least 50% of area
  double coverage = 50.0;
  KEY *node;
  HASHTBL *diff_to_key;

  total = hashtbl_numkeys(hash);
  //if (total < 4)
  //return 1.0;

  if (!(diff_to_key = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for diff_to_key failed");
    exit(EXIT_FAILURE);
  }

  //translate start percentage to actual num of elements
  freq_target = start*NUMSTRUCTS/100;

  helices = malloc(sizeof(int)*total);
  //initialize top three dropoff locations
  dropoffs = malloc(sizeof(int)*3);
  for (i = 0; i<3; dropoffs[i++] = 0) ;
  i = 0;
  for (node = hashtbl_getkeys(hash); node; node = node->next) {
    helices[i] = *((int*)hashtbl_get(hash,node->data));
    //printf("helices[%d] = %d\n",i,helices[i]);
    sum += helices[i++];
  }
  sum = NUMSTRUCTS;
  qsort(helices,total,sizeof(int),compare);
  //printf("freq target is %d with sum %d last num %d at %d\n",freq_target,sum,helices[i-1],i);
  //for (i = 0; i < total; i++)
  //printf("helices[%d] is %d\n",i,helices[i]);
  for (--i; helices[i] > freq_target && i >= 0; i--)
    partial += helices[i];
  //if everything above start, return the start
  if (i < 0)
    return (double) start;

  //if stopping at freq_target results in less than 50% coverage, go to at least 50%
  frac = ((double) partial*100)/((double)sum);
  if (frac < coverage) {
    //printf("not up to coverage with %d/%d\n",partial,sum);
    while (frac < coverage) {
      partial += helices[i--];
      frac = ((double) partial*100)/((double)sum);
    }
    //printf("now at %d/%d = %.1f at helices[%d] = %d\n",partial,sum,frac,i+1,helices[i+1]);
  }

  //create entry for 0, in case no drop off exists
  sprintf(key,"%d",0);
  val = malloc(sizeof(int));
  *val = helices[i];
  hashtbl_insert(diff_to_key,key,val);

  for ( ; i > 0; i--) {
    ave = (helices[i+1]+helices[i-1])/2;
    diff = ave - helices[i];
    //printf("ave is %d and diff is %d for %d\n",ave,diff,i);
    if (diff > dropoffs[0]) {
      //printf("bumping off %d for %d\n",dropoffs[0],diff);
      dropoffs[0] = diff;
      qsort(dropoffs,3,sizeof(int),compare);
      sprintf(key,"%d",diff);
      val = malloc(sizeof(int));
      *val = helices[i];
      hashtbl_insert(diff_to_key,key,val);
      //printf("inserting %s with %d\n",key,*val);
    }
  }
  printf("Possible cutoffs: ");
  for (i = 0; i<3; i++) {
    sprintf(key,"%d",dropoffs[i]);
    val = hashtbl_get(diff_to_key,key);
    if (val) {
      h = ((double)(*val+1)*100)/(double)NUMSTRUCTS;
      printf("%.1f ",h);
    }
  }
  printf("\n");
  hashtbl_destroy(diff_to_key);
  return h;
}

//looks up and prints actual helices for all id's
int print_all_helices(int total) {
  FILE *fp;
  HASHTBL *temp;
  KEY *node;
  char key[ARRAYSIZE],*val;
  double *trip;
  int i,*m;

  if (STATS)
    fp = fopen("triplet.out","w");

  for (i = 1; i < total; i++) {
    sprintf(key,"%d",i);
    val = hashtbl_get(idhash,key);
    m = hashtbl_get(marginals,key);
    trip = hashtbl_get(avetrip,key);
    if (val != NULL)
      printf("Helix %d is %s (%.1f %.1f %.1f) with freq %d\n",i,val,trip[0],trip[1],trip[2],*m);
    else
      printf("No entry for %d\n",i);
    if (STATS) {
      //printing triplet info
      val = hashtbl_get(idhash,key);
      fprintf(fp,"For id %d with frequency %d, represented by %s:\n",i,*m,val);
      temp = hashtbl_get(triplet,key);
      if (!temp)
	fprintf(stderr,"error: no entry for %d: print_all_helices\n",i);
      else {
	for (node = hashtbl_getkeys(temp); node; node = node->next) {
	  m = hashtbl_get(temp,node->data);
	  fprintf(fp, "\t%s %d\n",node->data,*m);
	}   
	hashtbl_destroy(temp);   
      }      
    }
  }
  if (STATS) {
    fclose(fp);
    hashtbl_destroy(triplet);
  }
  return 0;
}

//finds all frequent helices > threshold and prints them out
//inserts into hash with key = id, val = 1
//returns linked list of frequent helices ordered in ascending frequency
char** find_freq(int total) {
  int *marg = NULL,i,count = 0,h,j,k;
  double percent;
  char key[ARRAYSIZE],**mostfreq,*val;
  KEY *node = NULL,*begin = NULL;

  if (!(freq = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() of freq failed");
    exit(EXIT_FAILURE);
  }
  if (!(common = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() of common failed");
    exit(EXIT_FAILURE);
  }
  if (FILTER)
    mostfreq = malloc(sizeof(char*)*NUMFREQ);
  for (i = 1; i < total; i++) {
    sprintf(key,"%d",i);
    //LENGTH depreciated for now, filtering on size happening in process_structs
    if (LENGTH) {
      val = hashtbl_get(idhash,key);
      sscanf(val,"%d %d %d",&h,&j,&k);
      if (k < LENGTH) continue;
    }
    marg = hashtbl_get(marginals,key);
    if (!marg)
      fprintf(stderr,"error: marginal in find_freq is null");
    percent = ((double)*marg)*100.0/((double)NUMSTRUCTS);
    if (percent >= THRESH_COMMON) {
      if (VERBOSE)
	printf("Common helix %s with freq %d/%d\n",key,*marg,NUMSTRUCTS);
      hashtbl_insert(common,key,"1");
    }
    else if (percent >= THRESH_FREQ) {
      //printf("freq helix %s: with freq %d\n",key,*marg);
      if (FILTER) {
	filter(key,mostfreq,count);
      } else {
	//freq_insert(key,*marg,length++);
	//length++;
	node = malloc(sizeof(KEY*));
	node->data = mystrdup(key);
	node->next = begin;
	begin = node;
	//printf("found freq helix %s\n",key);
      }
      count++;
    }
  }

  if (FILTER) {  
    if (count > NUMFREQ)
      count = NUMFREQ;    
  }
  else {  //no filter
    mostfreq = malloc(sizeof(char*)*(count+hashtbl_numkeys(common)));
    for (i = 0; begin; begin = node) {
      if (i >= count) fprintf(stderr,"mismatch in freq helices: find_freq\n");
      mostfreq[i++] = begin->data;
      node = begin->next;
      free(begin);
    }
  }

  //add in common helices to array of most common; 

  //printf("found %d non-common helices and %d common ones; numfreq now %d\n",count, hashtbl_numkeys(common),NUMFREQ);
  if ((i = hashtbl_numkeys(common)) > 0) {
    if (FILTER && (count+i > NUMFREQ)) {
      mostfreq = realloc(mostfreq,(count+i)*sizeof(char*));
    }
    //printf("new size %d\n",count+i);
    NUMFREQ = count + hashtbl_numkeys(common);    
    i = count;
    count += hashtbl_numkeys(common);
    for (node = hashtbl_getkeys(common); node; node = node->next)
      mostfreq[i++] = node->data;
  }

  //insert most freq helices in ascending order
  qsort(mostfreq,count,sizeof(char*),charcompare);
  if (count > 63) fprintf(stderr,"number of helices greater than allowed in find_freq()\n");
  for (i = 0; i < count; i++)
    freq_insert(mostfreq[i],*((int*)hashtbl_get(marginals,mostfreq[i])),i);

    //sort it back for pruning purposes
  //if (count == NUMFREQ)
  qsort(mostfreq,NUMFREQ,sizeof(char*),freqcompare);
    //else
    //for (i = count; i < NUMFREQ; i++)
    //mostfreq[i] = "0";

  return mostfreq;
}

int charcompare(const void *v1, const void *v2) {
  return (atoi(*((char**)v1))-atoi(*((char**)v2)));
}
//maintains most frequent helices, at most NUMFREQ
//opted for this rather than keep a long list of all freq helices, then sort and take top
//because long sequences may have very long list of frequent helices
void filter(char *key,char **mostfreq, int count) {
  int i,k,*least,*keyfreq;
  char *temp,*last;

  if (count < NUMFREQ) {
    mostfreq[count] = strdup(key);
    //printf("in mostfreq: %s as %d\n",mostfreq[count],count);
    if (count == NUMFREQ -1) {
      qsort(mostfreq,NUMFREQ,sizeof(char*),freqcompare);
    }
  } else {
    //if not among top 10 most freq, return
    least = hashtbl_get(marginals,mostfreq[NUMFREQ-1]);
    keyfreq = hashtbl_get(marginals,key);
    if (*keyfreq <= *least) return;
    //else, find where to insert...
    //    printf("testing key %s with freq %d\n",key,*((int*)hashtbl_get(marginals,key)));
    k = binsearch(mostfreq,key);
    if (*((int*)hashtbl_get(marginals,mostfreq[k])) > *((int*)hashtbl_get(marginals,key))) 
      k++;
    //printf("in mostfreq replacing %s with %s at %d\n",mostfreq[k],key,k);
    last = strdup(key);
    //..and bump everything down one
    for (i = k; i < NUMFREQ; i++ ) {
      temp = mostfreq[i];
      mostfreq[i] = last;
      last = temp;
    }
    free(last);
  }
}

//will sort to have descending freq
int freqcompare(const void *v1, const void *v2) {
  int *i1 =  hashtbl_get(marginals,*(char**)v1);
  if (!i1) fprintf(stderr,"%s not found in marginals\n",*(char**)v1);
  int *i2 =  hashtbl_get(marginals,*(char**)v2);
  if (!i2) fprintf(stderr,"%s not found in marginals\n",*(char**)v2);
  //  printf("i1 is %d and i2 is %d\n",*i1,*i2);
  return (*i2 - *i1);
}
//returns the index i next to where key should be (either side possible)
int binsearch(char **mostfreq, char *key) {
  int left = 0;
  int right = NUMFREQ - 1;
  int mid,val;
  char **l,**r;

  while (left < right) {
    mid = (left+right)/2;
    //printf("left is %d, right is %d, and mid %d\n",left,right,mid);
    l = &(mostfreq[mid]);
    r = &key;
    val = freqcompare(l,r);
    //    printf("val is %d\n",val);
    if ( val > 0)
      right = mid - 1;
    else if (val < 0)
      left = mid + 1;
    else
      return mid;
  }
  return left;
}

//inserts frequent helix and makes binary rep
void freq_insert(char *key,int marg,int length) {
  unsigned long *bin;
  char *val;

  //  printf("inserting into freq, key %s\n",key);
  hashtbl_insert(freq,key,"1");
  bin = malloc(sizeof(unsigned long));
  *bin = (1<<length);
  hashtbl_insert(binary,key,bin);
  //printf("inserting %s with value %d\n",key,(1<<length));
  val = hashtbl_get(idhash,key);
  if (val) {
    if (VERBOSE)
      printf("Freq helix %s: %s with freq %d\n",key,val,marg);
  }
  else 
    printf("error: no triplet found for helix ID %s\n",key);
}

//like cluster.pl
//hash freq: key = helix ID. val = 1 (indicator)
//hash cluster: key = (num helices) list of helices. val = freq
int make_profiles(char *name) {
  FILE *fp,*file;
  int num=0,*id = 0,i,j,k,last = -1, lastfreq = -1,iscommon = 0;
  int numhelix = 0,size=INIT_SIZE,notcommon = 0,tripsize = INIT_SIZE;
  char temp[100],val[ARRAYSIZE],*l=NULL,*profile=NULL,*trips;
  HASHTBL *halfbrac;

  profile = malloc(sizeof(char)*ARRAYSIZE*size);
  //printf("length of profile is %d\n",(int)strlen(profile));
  profile[0] = '\0';

  if (!(cluster = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for cluster failed");
    exit(EXIT_FAILURE);
  }
  if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
    exit(EXIT_FAILURE);
  }
  if (!(bracket = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for bracket failed");
    exit(EXIT_FAILURE);
  }
  if (REP_STRUCT) {
    if (!(consensus = hashtbl_create(HASHSIZE,NULL))) {
      fprintf(stderr, "ERROR: hashtbl_create() for consensus failed");
      exit(EXIT_FAILURE);
    }
    trips = malloc(sizeof(char)*tripsize*ARRAYSIZE);
    trips[0] = '\0';
  }
  fp = fopen(name,"r");
  file = fopen("structure.out","w");
  fprintf(file,"Processing %s\n",name);
  if (fp == NULL) {
    fprintf(stderr, "can't open %s\n",name);
    return 0;
  }
  while (fgets(temp,100,fp) != NULL) {
    //printf("temp is %s\n",temp);
    if (sscanf(temp,"Structure %d",&num) == 1) {
      if (last == -1) {
	fprintf(file,"Structure %d: ",num);	
	continue;
      }
      fprintf(file,"\n\t-> %s\nStructure %d: ",profile,num);
      if (iscommon < hashtbl_numkeys(common)) {
	if (VERBOSE) 
	  printf("Found profile %snot having common helices\n",profile);
	notcommon++;
      }
      //else
      profile = process_profile(halfbrac,profile,numhelix,&size);
      if (REP_STRUCT) {
	make_rep_struct(profile,trips);
      }
      //if (VERBOSE && (count = hashtbl_get(cluster,profile)) && (*count == 1))
      //printf("First struct %d with profile %s\n",num,profile);
      if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
	fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
	exit(EXIT_FAILURE);
      }
      last = 0;
      numhelix = 0;
      profile[0] = '\0';
      if (REP_STRUCT)
	trips[0] = '\0';
      iscommon = 0;
    }
    else if (sscanf(temp,"%d %d %d",&i,&j,&k) == 3) {
      sprintf(val,"%d %d",i,j);
      id = hashtbl_get(bp,val);
      if (!id) continue;
      sprintf(val,"%d",*id);
      id = hashtbl_get(translate_hc,val);
      if (!id) continue;
      if (*id != -1 && *id != last) {
	fprintf(file,"%d ",*id);
	sprintf(val,"%d",*id);
	if (hashtbl_get(common,val)) {
	  //printf("found helix %s out of %d common\n",val,hashtbl_numkeys(common));
	  iscommon++;
	} 
	l = hashtbl_get(freq,val);
	if (l != NULL && *id != lastfreq) {   //is a freq helix, so save
	  numhelix++;
	  if (strlen(profile)+strlen(val) > (ARRAYSIZE*size-2)) 
	    profile = resize(&size,strlen(profile)+strlen(val)+2,profile);
	  //printf("adding %d to profile\n",id);
	  strcat(profile,val);
	  strcat(profile," ");
	  make_brackets(halfbrac,i,j,*id);
	  lastfreq = *id;
	} 
	last = *id;
      } 

      if (REP_STRUCT) {
	sprintf(val,"%d %d %d ",i,j,k);
	while (strlen(trips)+strlen(val) > (ARRAYSIZE*tripsize-1))
	  trips = realloc(trips,++tripsize*ARRAYSIZE);
	strcat(trips,val);
      }
    }
  }
  fprintf(file,"\n\t-> %s ",profile);
  if (iscommon != hashtbl_numkeys(common)) {
    if (VERBOSE)
      printf("Found profile %snot having common helices\n",profile);
    notcommon++;
  }
  profile = process_profile(halfbrac,profile,numhelix,&size);
  //fprintf(file,"Structure %d: %s\n",num,profile);
  
  free(profile);
  fclose(fp);
  fclose(file);
  return notcommon;
}

void make_rep_struct(char *profile, char* trips) {
  int *bpfreq,i,j,k;
  char *val,*blank = " ",bpair[ARRAYSIZE];
  HASHTBL *ij;
  
  ij = hashtbl_get(consensus,profile);
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
      bpfreq = hashtbl_get(ij,bpair);
      if (bpfreq)
	(*bpfreq)++;
      else {
	bpfreq = malloc(sizeof(int));
	*bpfreq = 1;
	hashtbl_insert(ij,bpair,bpfreq);
      }
      //printf("in rep struct for %s, inserting %d %d\n",profile,i+k,j-k);      
    }
  }
}

int print_profiles() {
  int *val;
  KEY *node;

  for (node = hashtbl_getkeys(cluster); node; node = node->next) {
    val = hashtbl_get(cluster,node->data);
    if (VERBOSE)
      printf("Profile %s with freq %d\n",node->data,*val);
  }
  return 0;
}

int select_profiles(char **mostfreq,int notcommon) {
  KEY *node = NULL;
  int *count = NULL,toosmall = 0,i,j,k,num = 0,*val,most=0;
  double percent,rep;
  char **prof,*copy,*blank = " ";

  if (!(freqprofs = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for freqprofs failed");
    exit(EXIT_FAILURE);
  }

  if (PROF_FREQ) {
    //if (VERBOSE)
    //printf("Original number of profiles before filtering: %d\n",hashtbl_numkeys(cluster));
    for (node = hashtbl_getkeys(cluster); node; node = node->next) {
      count = hashtbl_get(cluster,node->data);
      percent = ((double) *count)*100.0 / ((double)NUMSTRUCTS);
      //printf("%s has percent %.1f\n",node->data,percent);
      if (percent < PROF_FREQ) {
	toosmall += *count;
      }
      else {
	if (VERBOSE)
	  printf("Freq profile %swith freq %d\n",node->data,*count);
	val = malloc(sizeof(int));
	*val = *count;
	hashtbl_insert(freqprofs,node->data,val);
	
	copy = mystrdup(node->data);
	i = atoi(strtok(copy,blank));
	if (most < i)
	  most = i;
	free(copy);	
      }
    }
  }
  if (VERBOSE) {
    printf("Number of structs with infrequent profile: %d\n",toosmall);
    //printf("Number of profiles before pruning: %d\n",hashtbl_numkeys(cluster));
  }

  if (NUMPROF) {
    j = hashtbl_numkeys(cluster);
    if (j > NUMPROF) {
      prof = malloc(sizeof(char*)*j);
      for (i = 0,node = hashtbl_getkeys(cluster); node; node = node->next)
	prof[i++] = mystrdup(node->data);
      qsort(prof,j,sizeof(char*),profcompare);
      //for (k = 0; k < j; k++) 
      //printf("%s with freq %d\n",prof[k],*((int*)hashtbl_get(cluster,prof[k])));      
      for (k = NUMPROF; k < j; k++) {
	i = *((int*)hashtbl_get(cluster,prof[k]));
	//if (VERBOSE)
	//printf("removing %s with freq %d\n",prof[k],i);
	num += i;
	if (hashtbl_remove(cluster,prof[k]) == -1)
	  fprintf(stderr,"Failed remove of %s in cluster\n",prof[k]);
      }
      if (VERBOSE)
	printf("Total number of structures without profile in top %d: %d\n",NUMPROF,num);
      for (i = 0; i < j; i++)
	free(prof[i]);
      free(prof);
    }
  }
  rep = ((double)(NUMSTRUCTS - (toosmall+num)))/((double) NUMSTRUCTS);
  /*
  if (rep < THRESH_STRUCT) {
    for (; hashtbl_numkeys(cluster) > NUMPROF; prune_profiles(mostfreq))
      if (VERBOSE)
	printf("pruning profiles, now at %d profiles\n",hashtbl_numkeys(cluster));    
  }
  */
  if (VERBOSE)
    printf("Number of structures without common helices: %d\n",notcommon);
  printf("Number of structures with direct representation in graph: %d/%d\n",NUMSTRUCTS - (toosmall+num),NUMSTRUCTS);

  if (VERBOSE)
    for (node = hashtbl_getkeys(freqprofs); node; node = node->next) 
      printf("Node %s with freq %d\n",node->data,*((int*)hashtbl_get(freqprofs,node->data)));
  //count = malloc(sizeof(int));
  //*count = most;
  //hashtbl_insert(cluster,"most",count);
  return most;
}

//will sort to have descending freq
int profcompare(const void *v1, const void *v2) {
  int *i1 =  hashtbl_get(cluster,*(char**)v1);
  if (!i1) fprintf(stderr,"%s not found in cluster\n",*(char**)v1);
  int *i2 =  hashtbl_get(cluster,*(char**)v2);
  if (!i2) fprintf(stderr,"%s not found in cluster\n",*(char**)v2);
  //  printf("i1 is %d and i2 is %d\n",*i1,*i2);
  return (*i2 - *i1);
}

//if profile is unique, insert into cluster, and make bracket representation
char* process_profile(HASHTBL *halfbrac,char *profile,int numhelix,int *size) {
  int *count,size2;
  char val[ARRAYSIZE],*dup;

  //printf("profile %s\n",profile);
  for (size2 = INIT_SIZE; strlen(profile) > (ARRAYSIZE *size2-1); size2++);  
  dup = malloc(sizeof(char)*ARRAYSIZE*size2);
  //numhelix += hashtbl_numkeys(common);
  sprintf(val,"%d ",numhelix);

  if (strlen(profile)+strlen(val)> ARRAYSIZE*(*size)-1)
    profile = resize(size,strlen(profile)+strlen(val)+1,profile);
  //  strcat(profile,comm);
  profile = strcat_front(profile,val);
  //printf("profile is %s, val is %s\n",profile,val);
  profile = quicksort(profile,dup);
  if ((count = hashtbl_get(cluster,profile)) == NULL) {
    count = malloc(sizeof(int));
    *count = 1;
    hashtbl_insert(cluster,profile,count);
    make_bracket_rep(halfbrac,dup);
    //printf("inserting %s for level %d\n",profile,numhelix-1);
  }
  else {
    ++*count;
    //printf("augmenting count of %s to %d\n",profile,*count-1);
  }
  //if (*most < numhelix)
  //*most = numhelix;
  free(dup);
  hashtbl_destroy(halfbrac);
  return profile;
}

//inserts bracket representation for i,j into a hash
void make_brackets(HASHTBL *brac, int i, int j, int id) {
  char key[ARRAYSIZE],*val;

  sprintf(key,"%d",i);
  val = malloc(sizeof(char)*ARRAYSIZE);
  sprintf(val,"[%d",id);
  //  printf("making bracket %s for %d\n",val,i);
  hashtbl_insert(brac,key,val);
  sprintf(key,"%d",j);
  val = malloc(sizeof(char)*2);
  val[0] = ']';
  val[1] = '\0';
  hashtbl_insert(brac,key,val);
}

//makes the bracket representation of dup, using values in hashtbl brac
//dup is a (mod) profile in graph
//called by process_profile()
void make_bracket_rep(HASHTBL *brac,char *dup) {
  int num,*array,k=0,size = INIT_SIZE,total;
  char *profile,*val;
  KEY *node = NULL;

  num = hashtbl_numkeys(brac);
  array = malloc(sizeof(int)*num);
  for (node = hashtbl_getkeys(brac); node; node=node->next) 
    array[k++] = atoi(node->data);
  //sort by i,j position  
  qsort(array,num,sizeof(int),compare);
  profile = malloc(sizeof(char)*ARRAYSIZE*size);
  profile[0] = '\0';
  val = malloc(sizeof(char)*ARRAYSIZE);
  for (k = 0; k < num; k++) {
    sprintf(val,"%d",array[k]);
    val = hashtbl_get(brac,val);
    if ((total = strlen(profile)+strlen(val)) > ARRAYSIZE*size-1)
      profile = resize(&size,++total,profile);
    strcat(profile,val);
  }
  //if (VERBOSE)
  //printf("%s for %s\n",profile,dup);
  hashtbl_insert(bracket,dup,profile);
  free(val);
  free(array);
}

/*
void process_inputs(FILE *fp) {
  char temp[100];
  FILE *file;

  if (!(file = fopen(INPUT,"r")))
    fprintf(stderr,"Cannot open %s\n",INPUT);
  while (fgets(temp,100,file)) {
    if (sscanf(temp,"%d %d %d",&i,&j,&k) == 3) {
      sprintf(tmp,"%d %d",i,j);
      id = hashtbl_get(bp,tmp);
      //printf("id is %d for %d %d %d\n",id,i,j,k);
      if (!id)
	id = process_native(i,j,k);
      if (*id != last) {
	sprintf(tmp,"%d",*id);
	if (hashtbl_get(freq,tmp)) {   //is a freq helix, save to profile
	}
	else {}
      }
    }
  }
}
*/

//resizes a dynamically allocated string s to appropriate size
char* resize(int *size,int total,char *s) {
  int old = *size;
  char *temp = NULL;
  for (; total > ARRAYSIZE * (*size);(*size)++) ;
  //if (old != *size)
  //printf("resizing to %d from %d for size %d\n",*size,old,total);
  temp = realloc(s,sizeof(char)*ARRAYSIZE*(*size));
  if (!temp)
    fprintf(stderr, "unable to resize %s\n",s);
  /*
  temp = malloc(sizeof(char)*ARRAYSIZE*(*size));
  strcpy(temp,s);
  free(s);
  s = temp;
  */
  return temp;
}

//implements quicksort to sort profile ID's into ascending order
char* quicksort(char *profile,char *dup) {
  int *array,i,length;
  char *blank = " ",*k=NULL,val[ARRAYSIZE],*copy = strdup(profile);

  k = strtok(copy,blank);
  length = atoi(k);
  array = malloc(sizeof(int)*length);
  for (i = 0; i < length; i++) {
    array[i] = atoi(strtok(NULL,blank));
    //printf("setting array[%d] to %d\n",i,array[i]);
  }
  qsort(array,length,sizeof(int),compare);
  dup[0] = '\0';
  for (i = 0; i < length; i++) {
    sprintf(val,"%d",array[i]);
    //printf("now adding %s\n",val);
    strcat(dup,val);
    strcat(dup," ");
  }
  //  strcat(k," ");
  //  profile = strcat_front(profile,k);
  sprintf(profile,"%s %s",k,dup);
  free(array);
  return profile;
}

//sorts numerically into ascending order
int compare(const void *v1, const void *v2) {
  return (*(int*)v1 - *(int*)v2);
}

//concatenate ct to the front of s; return s
//assumes s has enough space to add ct
char* strcat_front(char *s, char *ct) {
  char *temp = mystrdup(s);
  strcpy(s,ct);
  strcat(s,temp);
  free(temp);
  return s;
  /*
  char *hold,*temp = malloc(sizeof(char)*(strlen(s)+strlen(ct)+1));
  //  printf("concatenating %s to %s\n",ct,s);
  strcpy(temp,ct);
  strcat(temp,s);
  hold = s;
  s = temp;
  printf("destroying %d and making new %d\n",hold,s);
  free(hold);
  return s;*/
}

//gets number of profiles to below NUMPROF
//eliminates profiles by removing least freq helix from consideration
void prune_profiles(char **mostfreq) {
  int k,*frq,*frq2,m;
  char *least,*profile;
  char *modprofile;
  KEY *node = NULL;
  HASHTBL *hash;

  for (k = NUMFREQ-1; !strcmp(mostfreq[k],"0"); k--) ;
  least = mostfreq[k];
  mostfreq[k] = "0";
  for (node = hashtbl_getkeys(cluster); node; node = node->next) {
    modprofile = malloc(strlen(node->data));
    profile = delete_helix(node->data,least,modprofile,&m);
    if (!strcmp(profile,node->data)) continue;
    frq = hashtbl_get(cluster,node->data);
    //modified profile exists or not
    if (!(frq2 = hashtbl_get(cluster,profile))) {
      //printf("adding %s t cluster\n",profile);
      hashtbl_insert(cluster,profile,frq);
      if (!(hash = hashtbl_create(HASHSIZE,NULL))) {
	fprintf(stderr, "ERROR: hashtbl_create() failed");
	exit(EXIT_FAILURE);
      }
      hashtbl_insert(hash,least,"1");
      edge_label(hash,modprofile,node->data,m);
    }
    else
      *frq2 += *frq;
    //printf("removing %s from cluster; num is %d\n",node->data,hashtbl_numkeys(cluster));
    hashtbl_remove(cluster,node->data);

    hashtbl_remove(bracket,modprofile);
  }
  
}

//takes cluster entry origprof and removes least if present
char *delete_helix(char *origprof, char *least,char *modprofile,int *m) {
  int length,found = 0;
  char *k,*blank = " ";
  char *copy = mystrdup(origprof);
  char *profile = malloc(strlen(origprof));

  modprofile[0] = '\0';
  profile[0] = '\0';
  length = atoi(strtok(copy,blank));
  *m = length;
  for (k = strtok(NULL,blank); k; k = strtok(NULL,blank)) {
    sprintf(modprofile,"%s%s ",modprofile,k);
    if (strcmp(k,least))
      sprintf(profile,"%s%s ",profile,k);
    else
      found = 1;
  }
  if (found) {
    length--;
    //assuming num of helices in profile will never go past 4 digits
    k = malloc(sizeof(char)*5);
    sprintf(k,"%d ",length);
    profile = strcat_front(profile,k);
    if (VERBOSE)
      printf("profile %safter deleting %s is %s\n",origprof,least,profile);
    free(k);
    free(copy);
    return profile;
  }
  free(copy);
  return origprof;
}

void find_consensus() {
  int *freq,*bpfreq;
  KEY *node,*bpnode;
  HASHTBL *ij,*final;

  if (!(final = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for final failed");
    exit(EXIT_FAILURE);
  }

  node = hashtbl_getkeys(consensus);
  //cycle thru all profiles, either calc consensus or destroying
  for (node = node->next; node; node = node->next) {
    freq = hashtbl_get(cluster,node->data);
    ij = hashtbl_get(consensus,node->data);
    if (freq) {
      if (!ij)
	fprintf(stderr, "ij not found in find_consensus()\n");
      for (bpnode = hashtbl_getkeys(ij); bpnode; bpnode = bpnode->next) {
	bpfreq = hashtbl_get(ij,bpnode->data);
	if (*bpfreq*100/(*freq) < 50)
	  hashtbl_remove(ij,bpnode->data);
	else
	  ;
	  //printf("for node %s, found bp %s with %d/%d\n",node->data,bpnode->data,*bpfreq,*freq);
      }
    } else {
      hashtbl_destroy(ij);
      //insert dummy pointer so remove won't seg fault
      hashtbl_insert(consensus,node->data, malloc(sizeof(char)));
      hashtbl_remove(consensus,node->data);
    }
  }
}

int print_consensus(char *seqfile) {
  int k=0,m,seqlen;
  char outfile[ARRAYSIZE],key[ARRAYSIZE],*pair,*i,*j,*blank = " ";
  KEY *node,*bpnode;
  HASHTBL *bpairs,*temp;
  FILE *fp;

  seqlen = strlen(seq);

  //foreach profile
  for (node = hashtbl_getkeys(consensus); node; node = node->next) {
    if (!(temp = hashtbl_create(HASHSIZE,NULL))) {
      fprintf(stderr, "ERROR: hashtbl_create() for temp failed");
      exit(EXIT_FAILURE);
    }
    sprintf(outfile,"Structure_%d.ct",++k);
    fp = fopen(outfile,"w");

    fprintf(fp,"Filename: %s\n",seqfile);
    fprintf(fp,"%d dG = n/a\n",seqlen);
    //    printf("processing %s\n",node->data);
    bpairs = hashtbl_get(consensus,node->data);
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
      if ((j = hashtbl_get(temp,key)) )
	fprintf(fp,"\t%d %c\t%d   %d   %s   %d\n",m+1,seq[m],m,m+2,j,m+1);
      else
	fprintf(fp,"\t%d %c\t%d   %d   0   %d\n",m+1,seq[m],m,m+2,m+1);
    }
    hashtbl_destroy(temp);
    fclose(fp);
  }
  hashtbl_destroy(consensus);
  return 0;
}

int print_cluster(char *seqfile) {
  int *count, length,i,j,seqlen,num=0;
  char key[ARRAYSIZE], *blank = " ",outfile[ARRAYSIZE];
  KEY *node = hashtbl_getkeys(cluster);
  FILE *fp;

  seqlen = strlen(seq);
  *count = hashtbl_numkeys(cluster)-1;
  
  for (i = 0; i < *count; i++) {
    strcpy(outfile,"Structure_");
    sprintf(key,"%d",i+1);
    strcat(outfile,key);
    strcat(outfile,".ct");
    fp = fopen(outfile,"w");

    fprintf(fp,"Filename: %s\n",seqfile);
    fprintf(fp,"%d dG = n/a\n",seqlen);
    for (j = 0; j < seqlen; j++) {
      fprintf(fp,"\t%d %c\t%d   %d   %d   %d",j+1,seq[j],j,j+2,num,j+1);
    }
    fclose(fp);
  }
  //fprintf(fp,"Using a threshold of %.1f:\n",THRESH_FREQ);
  for (node = node->next; node; node = node->next) {
    //key = mystrdup(node->data);
    count = hashtbl_get(cluster,key);
    if (count)
      printf("Key %s has %d elements\n",key,*count);
    else
      printf("No entry for %s",key);

    length = atoi(strtok(key,blank));
    fprintf(fp,"Representative structure:\n");
    for (i = 0; i < length; i++) {
      fprintf(fp,"\t%s\n",(char*)hashtbl_get(idhash,strtok(NULL,blank)));
    }
    
  }
  return 0;
}

