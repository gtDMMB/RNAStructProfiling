/*Used with main.c,profile.c in profile graph
Finds freq helices in 1000 structs; finds their LCA
Creates Hasse diagram with these vertices
*/

#include "hashtbl.h"
#include "summaryGraph.h"
#include "profile.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "memoryDFS.h"

static HASHTBL *edges;
static HASHTBL **graph1;
static HASHTBL *input;
static char **table;
static char **profileID;
static char **modprofileID;
static unsigned long *sums;
static int most;
static int graphsize;

//runs through nodes in cluster, finding their LCA's
int insert_graph(FILE *fp,char *file,int gsize) {
  int i,*freq,k = 0,numkeys,total,size=1;
  char *profile,*copy,*diff,*parent;
  KEY *node;
  struct hashnode_s *begin;

  if (!(edges = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() failed");
    exit(EXIT_FAILURE);
  }

  fputs("digraph G {\n",fp);
  fprintf(fp,"\tlabel = \"%s\";\n",file);
  fprintf(fp,"\tpad = 0.5;\n");
  fprintf(fp,"\tnodesep = 0.5;\n");
  //fprintf(fp,"\"\" [label=\"(0/%d)\"]\n",NUMSTRUCTS);
  //fputs("{ node [shape = plaintext]; ",fp);
  graphsize = gsize;
  if (INPUT) {
    i = process_one_input(fp);
    if (i > gsize)
      graphsize = i;
  }
  //printf("graphsize is %d for %s\n",*count,node->data);
  graph1 = malloc(sizeof(HASHTBL*)*(graphsize+1));
  for (i = 0; i <= graphsize; i++) graph1[i] = NULL;

  numkeys = hashtbl_numkeys(freqprofs);
  while (numkeys >= ARRAYSIZE*size) size++;
  sums = malloc(sizeof(unsigned long)*ARRAYSIZE*size);
  k = numkeys-1;
  profileID = malloc(sizeof(char*)*ARRAYSIZE*size);
  modprofileID = malloc(sizeof(char*)*ARRAYSIZE*size);
  most = 0;
  //for each profile, insert with freq
  
  for (node = hashtbl_getkeys(freqprofs); node; node = node->next) {
    //printf("node data is %s with k = %d\n",node->data,k);
    freq = hashtbl_get(freqprofs,node->data);
    //need to insert into graph
    if (freq) {
      if (most < *freq)
	most = *freq;
      profile = malloc(strlen(node->data)+1);
      sums[k] = insert_and_binary(node->data,profile,*freq);      
      profileID[k] = node->data;
      modprofileID[k--] = profile;
      if (*freq > 0)
	fprintf(fp,"\"%s\" [shape = box];\n",profile);
    }
    else
      fprintf(stderr,"No entry for %s\n",node->data);
    //printf("for profile %s binary rep is %d with ID %d\n",node->data,sums[k+1],k+1);
  }
  //make_key();
  sums[numkeys] = 0;
  modprofileID[numkeys] = " ";

  total = find_LCAs(fp,numkeys+1,&size);
  //printf("found %d vertices\n",total);
  calc_gfreq(fp,total);
  
  //total = print_vertices(fp);
  printf("Total number of vertices: %d\n",total);
  if (hashtbl_numkeys(edges) > 1) {
    start_trans_reductn(edges,total);
    print_edges(fp);
    //printGraph(); 
  } else {
    node = hashtbl_getkeys(edges);
    copy = mystrdup(node->data);
    printf("copy is %s\n",copy);
    diff = hashtbl_get(edges,copy);
    parent = strtok(copy,"->");
    copy = strtok(NULL,"->");
    fprintf(fp,"\"%s\" -> \"%s\" [label = \"%s\", arrowhead = vee];\n",parent,copy,diff);
  }
  if (INPUT)
    hashtbl_destroy(input);
  for (i = 0; i < graphsize; i++) {
    if (graph1[i])
      hashtbl_destroy(graph1[i]);
  }
  freeGraph();
  free(sums);
  free(table);
  free(profileID);
  free(modprofileID);

  hashtbl_destroy(edges);
  return 0;
}

//inserts key into graph and converts from rep struct to binary notation
unsigned long insert_and_binary(char *key,char *profile,int freq) {
  char *blank = " ",*helix,*copy = strdup(key);
  unsigned long sum = 0;
  int i,length = atoi(strtok(copy,blank)),*val;
  HASHTBL *hash = graph1[length-1];

  profile[0] = '\0';
  for (i = 0; i < length; i++) {
    helix = strtok(NULL,blank);
    strcat(profile,helix);
    strcat(profile," ");
    sum += *((unsigned long*) hashtbl_get(binary,helix));
    //printf("sum is now %d\n",sum);
  }
  //printf("profile is %s\n",profile);
  if (!hash) {
    if (!(hash = hashtbl_create(HASHSIZE,NULL))) {
      fprintf(stderr, "ERROR: hashtbl_create() failed");
      exit(EXIT_FAILURE);
    }
    graph1[length-1] = hash;
  }
  val = malloc(sizeof(int));
  *val = 0;
  hashtbl_insert(hash,profile,val);
  //printf("profile %sinserted with length %d and freq %d\n",profile,length,*val);
  return sum;
}

//table[freq id] = helix id
//produces table,where ith element is val
//where (1<<i) = hash_get(binary,val)
void make_key() {
  int count;
  KEY *node =  hashtbl_getkeys(binary);

  count = hashtbl_numkeys(binary);
  //for (count = 0; (*num & (1<<count)) == 0; count++) ;
  table = malloc(sizeof(char*)*count);
  for (count--; count >= 0; count--) {
    printf("table[%d] id is %s\n",count,node->data);
    table[count] = node->data;
    node = node->next;
  }
}

/* finished profiles = [0 start-1]
   present LCA to be intersected = [start oldk]
   newly generated LCA = [oldk k]
   returns k, the number of vertices in graph
 */
int find_LCAs(FILE *fp,int k,int *size) {
  int new,oldk,go,start,count,firsttime = 0;
  unsigned long num;
  char *profile,*prof;

  start = 0;
  for (oldk = k; start != k; oldk = k) {
    for (new = start; new != oldk; new++) {
      for (go = advance(new,oldk); go != start; go = advance(go,oldk)) {
	num = sums[new] & sums[go];
	//printf("num is %u of s[%d] = %u and s[%d] = %u\n",num,new,sums[new],go,sums[go]);
	if (not_in_sums(num,k)) {
	  //printf("found new profile for %s and %s\n",modprofileID[new],modprofileID[go]);

	  profile = convert_binary(num,&count);
	  fprintf(fp,"\"%s\" [style = dashed];\n",profile);
	  insert_prof(count,profile);
	  if (k >= ARRAYSIZE*(*size)) {
	    sums = realloc(sums,sizeof(unsigned long)*ARRAYSIZE*++*size);
	    modprofileID = realloc(modprofileID,sizeof(char*)*ARRAYSIZE*(*size));
	    profileID = realloc(profileID,sizeof(char*)*ARRAYSIZE*(*size));
	  }
	  //printf("k is %d\n",k);
	  modprofileID[k] = profile;
	  if (count > 999) fprintf(stderr,"count > 999 in find_LCAs\n");
	  prof = malloc(sizeof(char)*strlen(profile)+5);
	  sprintf(prof,"%d %s",count,profile);
	  profileID[k] = prof;
	  sums[k] = num;     
	  firsttime = found_edge(fp,new,k);	  
	  firsttime = found_edge(fp,go,k);	  
	  k++;
	} else if (num == sums[new]) {
	  firsttime = found_edge(fp,go,new);	  
	}
	else if (num == sums[go])
	  firsttime = found_edge(fp,new,go);	  
      } 
    }
    start = oldk;
  }
  return k;
}

//wraps around like mod function
int advance(int new, int oldk) {
  if (new+1 != oldk)
    return new+1;
  else
    return 0;
}

//returns 1 if num doesn't match anything in sums up to oldk
//returns 0 otherwise
int not_in_sums(unsigned long num, int k) {
  int i;
  for (i = 0; i < k; i++) {
    if (sums[i] == num)
      return 0;
  }
  return 1;
}


//converts binary rep to string of helices (profile)
char* convert_binary(unsigned long binary,int *count) {
  int k,size = 1,num=0;
  char val[ARRAYSIZE],*profile;

  profile = malloc(sizeof(char)*ARRAYSIZE);
  profile[0] = '\0';
  for (k = 0; binary > 0; binary >>= 1, k++) {
    //printf("binary is %u\n",binary);
    if ((binary & 1) == 1) {
      sprintf(val,"%d",k+1);
      if (strlen(profile)+strlen(val) > ARRAYSIZE*size-2) 
	profile = realloc(profile,sizeof(char)*ARRAYSIZE*++size);
      //printf("adding %s, with k %d, binary is %u, shifted is %u\n",table[k],k,binary,binary>>1);
      strcat(profile,val);
      strcat(profile," ");
      num++;
    }
  }
  *count = num;
  return profile;
}

void insert_prof(int k, char *profile) {
  int *val;
  HASHTBL *hash = graph1[k-1];

  if (!hash) {
    if (!(hash = hashtbl_create(HASHSIZE,NULL))) {
      fprintf(stderr, "ERROR: hashtbl_create() failed");
      exit(EXIT_FAILURE);
    }
    graph1[k-1] = hash;
  }
  val = malloc(sizeof(int));
  *val = 0;
  hashtbl_insert(hash,profile,val);  
  //return val;
}

  /*make edge
  update parents' general freq in graph[] with child's spec. freq if exists
  helix set difference = edge label
  child and parent are indices
  returns 1 if inserted into edge hash, 0 otherwise
  */
int found_edge(FILE *fp,int child,int parent) {
  int count;
  unsigned long xor;
  char *edge,*childprof,*parentprof,*diff;

  childprof = modprofileID[child];
  parentprof = modprofileID[parent];
  //printf("child is %s and parent is %s\n",childprof,parentprof);
  edge = malloc(sizeof(char)*(strlen(childprof)+strlen(parentprof)+3));
  sprintf(edge,"%s->%s",parentprof,childprof);
  if (!hashtbl_get(edges,edge)) {
    xor = sums[child] ^ sums[parent];
    diff = convert_binary(xor,&count);
    hashtbl_insert(edges,edge,diff);
    //fprintf(fp,"\"%s\" -> \"%s\" [label = \"%s\", arrowhead = vee];\n",parentprof,childprof,diff);
    return 1;
  }
  //printf("Found %d is child of %d\n",child,parent);
  return 0;
}

void calc_gfreq(FILE *fp,int total) {
  int i,j=0,length,*val,k,*freq;
  unsigned long *sum;
  char *fprof,*orig,*blank = " ",*helix,*copy;
  KEY *node,*vert;
  HASHTBL *hash;

  sum = malloc(sizeof(unsigned long)*hashtbl_numkeys(cluster));
  for (node = hashtbl_getkeys(cluster); node; node = node->next) {
    orig = mystrdup(node->data);
    length = atoi(strtok(orig,blank));
    sum[j] = 0;
    //find its binary rep
    for (k = 0; k < length; k++) {
      helix = strtok(NULL,blank);
      //all helices are freq helices in binary hash
      sum[j] += *((unsigned long*) hashtbl_get(binary,helix));
    }
    j++;
  }
  for (i = 0; i < total; i++) {
    if (!strcmp(modprofileID[i]," ")) {
      fprintf(fp,"\" \" [label =\"0\"];\n");
      continue;
    }
    j=0;
    for (node = hashtbl_getkeys(cluster); node; node = node->next) {
      //printf("investigating %s and prof[%d] = %s\n",node->data,i,profileID[i]);
      //printf("%u & %u = %u\n",sum[j],sums[i],sum[j]&sums[i]);
      if ((sum[j++] & sums[i]) == sums[i]) {
	fprof = profileID[i];
	//printf("found parent[%d] %s of %s,adding %d\n",i,fprof,node->data,*((int*)hashtbl_get(cluster,node->data)));
	copy = mystrdup(fprof);
	length = atoi(strtok(copy,blank));
	hash = graph1[length-1];
	fprof = modprofileID[i];
	val = hashtbl_get(hash,fprof);
	if (val)
	  *val += *((int*)hashtbl_get(cluster,node->data));
	else
	  fprintf(stderr,"error:val not found in hash()\n");
	free(copy);
      }
    }
    freq = hashtbl_get(cluster,profileID[i]);
    helix = hashtbl_get(bracket,modprofileID[i]);
    if (!freq) 
      fprintf(fp,"\"%s\" [label = \"0/%d\"];\n",modprofileID[i],*val);
    else
      fprintf(fp,"\"%s\" [label = \"%s\\n%d/%d\"];\n",modprofileID[i],helix,*freq,*val);
  }
}

void print_edges(FILE *fp) {
  int i, j,size = 2,oldsize;
  char *edge,*diff=NULL;
  edge = malloc(sizeof(char)*ARRAYSIZE*size);
  for(i = 0; i < GRAPHSIZE; i++) {
    if(graph[i].label != NULL) {
      if (VERBOSE)
	printf("node: '%s'\n", graph[i].label);
      for(j = 0; j < graph[i].numNeighbors; j++) {
	oldsize = size;
	while (strlen(graph[i].label) + strlen(graph[i].neighbors[j]->label) + 3 > ARRAYSIZE*size)
	  size++;
	if (oldsize != size)
	  edge = realloc(edge,sizeof(char)*ARRAYSIZE*size);
	sprintf(edge,"%s->%s",graph[i].label,graph[i].neighbors[j]->label);
	diff = hashtbl_get(edges,edge);
	if (diff) {
	  fprintf(fp,"\"%s\" -> \"%s\" [label = \"%s\", arrowhead = vee];\n",graph[i].label,graph[i].neighbors[j]->label,diff);
	} else {
	  fprintf(stderr, "edge %s found in graph not found in hash: print_edges()\n",edge);
	}
	//printf("'%s', ", graph[i].neighbors[j]->label);
      }
    }
  }	
}

//processes native helices if necessary; called by process_input
//returns id for i,j,k helix
int* process_native(int i, int j, int k) {
  int *id = NULL,l;
  char *key,tmp[ARRAYSIZE];

  for (l=1; l < k; l++) {
    sprintf(tmp,"%d %d",i+l,j-l);
    id = hashtbl_get(bp,tmp);
    if (id) {
      for (l-- ; l >= 0; l--) {
	sprintf(tmp,"%d %d",i+l,j+l);
	hashtbl_insert(bp,tmp,id);
      }
      return id;
    }
  }
  //printf("helix %d %d %d doesn't exist in sfold sample\n",i,j,k);
  id = malloc(sizeof(int));
  *id = hashtbl_numkeys(idhash)+1;
  hashtbl_insert(bp,tmp,id);
  sprintf(tmp,"%d",*id);
  key = malloc(sizeof(char)*ARRAYSIZE);
  sprintf(key,"%d %d %d",i,j,k);
  hashtbl_insert(idhash,tmp,key);
  return id;
}

//process input into three profiles Hs,Hu,Hn
//s = selected profiles, u = unselected in sample,n=in native only
//connects them, puts Hs in cluster with freq = 0 if not already there
int process_one_input(FILE *fp) {
  HASHTBL *halfbrac;
  FILE *file;
  char temp[100],tmp[ARRAYSIZE],*profile,*fullprofile,*diff,*native,*diffn,*dup;
  int i,j,k,*id,last=0,insample;
  int numhelix = 0,fullnum = 0,natnum = 0;
  int size = INIT_SIZE,size2 = INIT_SIZE,size3 = INIT_SIZE,size4 = INIT_SIZE,size5 = INIT_SIZE;
  
  if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
    exit(EXIT_FAILURE);
  }
  if (!(input = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for input failed");
    exit(EXIT_FAILURE);
  }
  //longest = hashtbl_get(max,"longest");
  profile = malloc(sizeof(char)*ARRAYSIZE*size);
  fullprofile = malloc(sizeof(char)*ARRAYSIZE*size2);
  diff = malloc(sizeof(char)*ARRAYSIZE*size3);
  native = malloc(sizeof(char)*ARRAYSIZE*size4);
  diffn = malloc(sizeof(char)*ARRAYSIZE*size5);
  profile[0] = '\0';
  fullprofile[0] = '\0';
  diff[0] = '\0';
  native[0] = '\0';
  diffn[0] = '\0';
  if (!(file = fopen(INPUT,"r")))
    fprintf(stderr,"Cannot open %s\n",INPUT);
  while (fgets(temp,100,file)) {
    //    if (sscanf(temp,"Structure %d (%d)",&i,&prob) == 2) 
    if (sscanf(temp,"%d %d %d",&i,&j,&k) == 3) {
      insample = 1;
      sprintf(tmp,"%d %d",i,j);
      id = hashtbl_get(bp,tmp);
      if (!id) {
	id = process_native(i,j,k);
	//printf("number in marginals %d\n",hashtbl_numkeys(marginals));
	if (*id > hashtbl_numkeys(marginals))
	  insample = 0;
      }
      if (*id != last) {
	sprintf(tmp,"%d",*id);
	if (strlen(native)+strlen(tmp) > (ARRAYSIZE*size4-2))
	  native = realloc(native,ARRAYSIZE*(++size4));
	sprintf(native,"%s%s ",native,tmp);
	natnum++;
	if (insample) {
	  fullnum++;
	  if (strlen(fullprofile)+strlen(tmp) > (ARRAYSIZE*size2-2)) 
	    fullprofile = realloc(fullprofile,ARRAYSIZE*(++size2));
	  sprintf(fullprofile,"%s%s ",fullprofile,tmp);
	  if (hashtbl_get(freq,tmp)) {   //is a freq helix, save to profile
	    numhelix++;
	    if (strlen(profile)+strlen(tmp) > (ARRAYSIZE*size-2)) 
	      profile = realloc(profile,ARRAYSIZE*(++size));
	    //printf("adding %d to profile\n",*id);
	    sprintf(profile,"%s%s ",profile,tmp);
	  } 
	  else { //if not freq record diff
	    if (strlen(diff)+strlen(tmp)+2 > ARRAYSIZE*size3)
	      diff = realloc(diff,ARRAYSIZE*(++size3));
	    sprintf(diff,"%s%s ",diff,tmp);	  
	    //printf("printing diff %s\n",diff);
	  }
	} 
	else {
	  if (strlen(diffn)+strlen(tmp)+2 > ARRAYSIZE*size5)
	    diffn = realloc(diffn,ARRAYSIZE*(++size5));
	  sprintf(diffn,"%s%s ",diffn,tmp);
	}	
	last = *id;
	make_brackets(halfbrac,i,j,*id);
      }
      else if (VERBOSE)
	printf("helix %d %d %d is duplicate with id %d: process_input()\n",i,j,k,*id);
    }
  }
  
  native = sort_input(native,natnum);
  //printf("native is now %s\n",native);
  make_bracket_rep(halfbrac,native);
  hashtbl_destroy(halfbrac);

  fullprofile = sort_input(fullprofile,fullnum);
  profile = sort_input(profile,numhelix);
  //printf("native %s, fullprofile %s, profile %s, diff %s, diffn %s\n",native,fullprofile,profile,diff,diffn);
  if (fullnum != natnum)
    make_edge_and_node(fp,fullprofile,native,diffn,natnum);
  if (numhelix != fullnum)
    make_edge_and_node(fp,profile,fullprofile,diff,fullnum);

  fprintf(fp,"\"%s\" [style = filled, fillcolor = gray60];\n",profile);  
  sprintf(tmp,"%d ",numhelix);
  if (strlen(tmp)+strlen(profile) > ARRAYSIZE*size+1)
    profile = realloc(profile,ARRAYSIZE*(++size));
  dup = mystrdup(profile);
  sprintf(profile,"%s%s",tmp,dup);
  id = hashtbl_get(cluster,profile);
  free(dup);
  
  if (!id) {			    
    id = malloc(sizeof(int));
    *id = 0;
    hashtbl_insert(cluster,profile,id);
    //printf("inserting input %s into cluster\n",profile);
  }
  
  free(profile);
  return numhelix;
}

//writes the destination node and edge to output file
void make_edge_and_node(FILE *fp,char *from, char *to,char *diff,int fullnum) {
  char *brac;
  HASHTBL *temp;

  if (!(temp = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for temp in make_edge_and_node() failed");
    exit(EXIT_FAILURE);
  }

  fprintf(fp,"\"%s\" [style = filled,fillcolor = gray70];\n",to);
  //how to handle getting general frequency?
  
  diff = insert_diff(temp,diff);
  //printf("diff is %s\n",diff);
  brac = edge_label(temp,from,to,fullnum);      
  fprintf(fp,"\"%s\"-> \"%s\" [label =\" %s\\n%s \",fontsize=8];\n",from,to,diff,brac); 
  hashtbl_destroy(temp);
}

//processes input file containing structures of interest, in triplet form
//similar code to make_cluster: converts triplet to profile,
void process_input(FILE *fp) {
  HASHTBL *halfbrac;
  FILE *file;
  char temp[100],tmp[ARRAYSIZE],*profile,*fullprofile,*diff;
  int i,j,k,*id,last=0,lastprob;
  int numhelix = 0,fullnum = 0,size = INIT_SIZE,size2 = INIT_SIZE,size3 = INIT_SIZE;

  if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
    exit(EXIT_FAILURE);
  }
  if (!(input = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for input failed");
    exit(EXIT_FAILURE);
  }
  //longest = hashtbl_get(max,"longest");
  profile = malloc(sizeof(char)*ARRAYSIZE*size);
  fullprofile = malloc(sizeof(char)*ARRAYSIZE*size2);
  diff = malloc(sizeof(char)*ARRAYSIZE*size3);
  profile[0] = '\0';
  fullprofile[0] = '\0';
  diff[0] = '\0';
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
	  numhelix++;
	  if (strlen(profile)+strlen(tmp) > (ARRAYSIZE*size-2)) 
	    profile = resize(&size,strlen(profile)+strlen(tmp)+1,profile);
	  //printf("adding %d to profile\n",*id);
	  sprintf(profile,"%s%s ",profile,tmp);
	}
	else {
	  if (strlen(diff)+strlen(tmp)+2 > ARRAYSIZE*size3)
	    diff = resize(&size3,strlen(diff)+strlen(tmp)+2,diff);	  
	  sprintf(diff,"%s%s ",diff,tmp);	  
	  //printf("printing diff %s\n",diff);
	}
	fullnum++;
	if (strlen(fullprofile)+strlen(tmp) > (ARRAYSIZE*size-2)) 
	  fullprofile = resize(&size2,strlen(fullprofile)+strlen(tmp)+2,fullprofile);
	sprintf(fullprofile,"%s%s ",fullprofile,tmp);
	//printf("helix %d added is %s\n",fullnum,tmp);
	last = *id;
	make_brackets(halfbrac,i,j,*id);
      }
      else //if id == last
	printf("helix %d %d %d is duplicate with id %d: process_input()\n",i,j,k,*id);
    } 
    else if (sscanf(temp,"Structure %d (%d)",&i,&j) == 2) {
      if (strlen(fullprofile) == 0) {
	lastprob = j;
	continue;
      } 
      //printf("profile is %s, fullprofile %s with diff %s\n\n",profile,fullprofile,diff);
      halfbrac = process_input_profile(fp,halfbrac,fullprofile,fullnum,profile,numhelix,diff,lastprob);
      numhelix = 0;
      fullnum = 0;
      profile[0] = '\0';
      fullprofile[0] = '\0';
      diff[0] = '\0';
      lastprob = j;
    }
  }
  //printf("input profile is %s with fullprofile %s and diff %s\n",profile,fullprofile,diff);
  halfbrac = process_input_profile(fp,halfbrac,fullprofile,fullnum,profile,numhelix,diff,lastprob);
  free(profile);
  hashtbl_destroy(halfbrac);
  
  //finds edges between centroids
  find_centroid_edges(fp);
}

//takes input profile and checks if in graph
//if so, changes node shape to hexagon
//if not, insert new vertex
HASHTBL* process_input_profile(FILE *fp,HASHTBL *brac,char *fullprofile, int fullnum,char *profile,int numhelix, char *diff, int prob) {
  HASHTBL *hash, *temp=NULL;
  char *diff1,*bracket,*difftrip;
  int *val,k1=0,k2=0;
  KEY *parent,*next;

  if (!(temp = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() failed");
    exit(EXIT_FAILURE);
  }

  make_bracket_rep(brac,fullprofile);
  hashtbl_destroy(brac);
  
  profile = sort_input(profile,numhelix);
  //printf("(sorted) profile is %s with fullprofile %s and diff %s\n",profile,fullprofile,diff);

  if ((hash = graph1[numhelix-1]) && (hashtbl_get(hash,profile))) {
    if (numhelix == fullnum) {
      //printf("case 1: full profile found in graph\n");
      fullprofile = sort_input(fullprofile,numhelix);
    } else {
    //cannot use find_diff because fullprofile has helices not in table[]
      //puts("case 2: profile found in graph");
      diff = insert_diff(temp,diff);
      bracket = edge_label(temp,profile,fullprofile,fullnum);      
      fprintf(fp,"\"%s\"-> \"%s\" [label =\" %s\\n%s \",fontsize=8,style = dotted];\n",profile,fullprofile,diff,bracket); 
    }
  }
  else {
    /*          
    if (numhelix == fullnum)
      puts("case 3: full profile not found");
    else
      puts("case 4: profile not found");
    */
    for (parent = find_parents(profile); parent; parent = next) {
      diff1 = find_diff(temp,parent->data,profile,&k1,&k2);      
      if (numhelix != fullnum) {
	difftrip = insert_diff(temp,diff);
	diff1 = realloc(diff1,strlen(diff1)+strlen(difftrip)+4);
	//printf("for parent %s, diff1 is now %s, diff is %s and difftrip is %s\n",parent->data,diff1,diff,difftrip);
	sprintf(diff1,"%s\\n%s",diff1,difftrip);
	//printf("Diff is %s for parent %s of profile %s; diff %s for %s\n",diff1,parent->data,profile,difftrip,fullprofile);	
      }
      bracket = edge_label(temp,parent->data,fullprofile,fullnum);
      fprintf(fp,"\"%s\"-> \"%s\" [label =\" %s\\n%s \",fontsize=8,style = dotted];\n",parent->data,fullprofile,bracket,diff1); 
      next = parent->next;
      free(parent);
    }
  }
  //printf("%s has size %d and prob %d\n",fullprofile,fullnum,prob);
  fprintf(fp,"\"%s\" [shape = hexagon];\n",fullprofile);
  val = malloc(sizeof(int)*2);
  val[0] = fullnum;
  val[1] = prob;
  hashtbl_insert(input,fullprofile,val);
  hashtbl_destroy(temp);  

  if (!(brac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() failed");
    exit(EXIT_FAILURE);
  }
  return brac;
}

//find all parents of profile in graph; profile has only freq helices
//transitive reduction to the rescue again
KEY* find_parents(char *profile) {
  int i,b1,b2,N,num = 0,k = 0;
  HASHTBL *hash;
  KEY *node,*parent,*begin = NULL;

  begin = malloc(sizeof(KEY));
  begin->data = "";
  begin->next = NULL;
  b2 = make_binary(profile,&num);
  //printf("finding parents for %s of length %d and binary %d\n",profile,num,b2);
  for (i = 0; i < num-1; i++) {
    if (!(hash = graph1[i])) continue;
    for (node = hashtbl_getkeys(hash); node; node = node->next) {
      b1 = make_binary(node->data,&k);
      //printf("investigating %sof size %d with binary %d\n",node->data,i+1,b1);
      N = b1 & b2;
      if (N == b1) {
	//printf("found parent %sof size %d with binary %d\n",node->data,i+1,b1);
	parent = malloc(sizeof(KEY));
	parent->data = node->data;
	parent->next = begin;
	begin = parent;
      }
    }
  }
  return begin;
}

//diff contains helix nums that are inserted into a hash
//to be used before edge_label
//add triplet info to diff
char* insert_diff(HASHTBL *temp,char *diff) {
  int size = INIT_SIZE;
  char *blank = " ", *copy,*k,*val;

  if (strlen(diff) == 0) return "";
  copy = malloc(sizeof(char)*ARRAYSIZE*size);
  copy[0] = '\0';
  for (k = strtok(diff,blank); k; k = strtok(NULL,blank)) {
    hashtbl_insert(temp,k,"1");
    //printf("insert_diff: inserting %s into hash\n",k);
    val = hashtbl_get(idhash,k);
    if (strlen(copy)+strlen(k)+strlen(val)+5 > ARRAYSIZE*size)
      copy = resize(&size,strlen(copy)+strlen(k)+strlen(val)+5,copy);	  
    if (strlen(copy)> 0) 
      strcat(copy,"\\n");
    sprintf(copy,"%s%s: %s",copy,k,val);
  }
  //free(diff);
  return copy;
}

//sorts input profile; code similar to quicksort()
char* sort_input(char *profile,int length) {
  int *array,i=0;
  char *blank = " ", *k,*copy = mystrdup(profile);

  array = malloc(sizeof(int)*length);
  for (k = strtok(copy,blank); k ; k = strtok(NULL,blank)) 
    array[i++] = atoi(k); 
  qsort(array,length,sizeof(int),compare);
  profile[0] = '\0';
  for (i = 0; i < length; i++) {
    sprintf(profile,"%s%d ",profile,array[i]);
  }
  //  printf("input profile sorted now %s\n",profile);
  free(copy);
  free(array);
  return profile;
}

//inserts centroids into graph
//finds and prints edges between them
void find_centroid_edges(FILE *fp) {
  int *i,*zero;
  KEY *node;
  HASHTBL *hash = NULL;

  zero = malloc(sizeof(int));
  *zero = 0;
  for (node = hashtbl_getkeys(input); node; node = node->next) {
    i = hashtbl_get(input,node->data);
    if (*i-1 > graphsize)
      *i = graphsize+1;
    //printf("inserting %s into graph[%d]\n",node->data,*i-1);
    if (!(hash = graph1[*i-1])) {
      if (!(hash = hashtbl_create(HASHSIZE,NULL))) {
	fprintf(stderr, "ERROR: hashtbl_create() for input failed");
	exit(EXIT_FAILURE);
      }
      graph1[*i-1] = hash;
    }
    if (!hashtbl_get(hash,node->data))
      hashtbl_insert(hash,node->data,zero);    
  }
}

int print_vertices(FILE *fp) {
  int i,*val,size = 5,total = 0,size2=INIT_SIZE,*frq = NULL,zero = 0,start,end;
  char *rank,*v;;
  HASHTBL *hash;
  KEY *node = NULL;

  rank = malloc(sizeof(char)*ARRAYSIZE*size);
  v = malloc(sizeof(char)*ARRAYSIZE*size2);
  for (i = 0; !graph1[i]; i++) ;
  start = i;
  //for (node = hashtbl_getkeys(graph[i]); node; node = node->next)
  //check_insert_edge(fp,"",node->data);
  //print ranks
  fputs("{ node [shape = plaintext]; ",fp);
  if (graph1[graphsize])
    end = graphsize;
  else
    end = graphsize-1;
  for ( ; i <= end; i++) {
    if (!graph1[i]) continue;
    fprintf(fp,"%d",i+1);
    if (i == end)
      fprintf(fp,"; }\n");
    else
      fprintf(fp,"->");
  }

  for (i = start; i <= end; i++) {
    if (!(hash = graph1[i])) continue;
    //printf("printing level %d\n",i+1);
    node = hashtbl_getkeys(hash);
    sprintf(rank,"{ rank = same; %d;",i+1);
    for (; node; node = node->next) {
      if (strlen(node->data)+ 4 > ARRAYSIZE*size2-1) {
	v = resize(&size2,strlen(node->data)+5,v);
	//printf("resizing v of size %d to %d\n",strlen(v),size2);
      }
      sprintf(v,"%d %s",i+1,node->data);
      frq = hashtbl_get(cluster,v);
      if (!frq) 
	frq = &zero;
      sprintf(v," \"%s\";",node->data);
      val = hashtbl_get(hash,node->data);
      //      printf("found %s for %s\n",hashtbl_get(bracket,node->data),node->data);
      if (VERBOSE)
	printf("Vertex %swith frequency %d, originally %d\n",node->data,*val,*frq);
      if (strlen(rank)+strlen(v) > ARRAYSIZE*size-1) {
	//	printf("resizing rank %s and v %s of size %d to %d\n",rank,v,strlen(rank)+strlen(v),size);
	rank = resize(&size,strlen(rank)+strlen(v)+1,rank);
      }
      strcat(rank,v);
      //fprintf(fp,"\"%s\" [label = \"%s ",node->data,hashtbl_get(bracket,node->data));
      fprintf(fp,"\"%s\" [label = \"(%d/%d)\"];",node->data,*frq,*val);
      /*
      if (*frq == most)
	fprintf(fp,"**");
      fprintf(fp,"%s",hashtbl_get(bracket,node->data));
      fprintf(fp,"(%d/%d)",*frq,*val);
      if (*frq == most)
	fprintf(fp,"**");
      if (INPUT && (val = hashtbl_get(input,node->data)))
      fprintf(fp,"\\n(%d)",val[1]);

      fprintf(fp,"\",shape = box,style=filled,color=black,fillcolor=grey%d];",(1000-*frq)/20+49);
      fprintf(fp,"\"%s\" [shape = box, label = \"%s (%d)\",style=filled,color=black,fillcolor=grey%d];\n",node->data,*val,(1000-*frq)/20+49);
      */
    }
    fprintf(fp,"%s }\n",rank);
    total += hashtbl_numkeys(hash);
    //v]0] = '\0';
  }
  return total;
}

//edges[LCA ID] = profile ID
//profile is LCA found, found[] is array of contributing originals, count their number
int find_edges(FILE *fp,char *profile, int *found, int count) {
  int i;
  char *origprof;

  //for each contributing original profile...
  for (i = 0; i < count; i++) {
    origprof = modprofileID[found[i]];
    //printf("original profile is %s, and subprofile is %swith total length %d\n",origprof,profile,strlen(profile)+strlen(origprof));
    if (!strcmp(profile,origprof)) continue;
    check_insert_edge(fp,profile,origprof);
  }

  return 0;
}

//checks whether edge exists and inserts if not
//ie insert if profile -> origprof doesn't exist yet
void check_insert_edge(FILE *fp,char *profile,char *origprof) {
  int k1=0,k2=0; //*f1,*f2;
  //double ratio;
  char *diff,*copy,*brac;
  HASHTBL *hash = NULL;

  /*
  hash = graph[k1-1];
  f1 = hashtbl_get(hash,profile);
  hash = graph[k2-1];
  f2 = hashtbl_get(hash,origprof);
  ratio = ((double)*f2)/((double)*f1);
  //  printf("ratio for %s and %s is %d/%d = %.2f\n",origprof,profile,*f2,*f1,ratio);
  */

  if (strlen(profile) == 0)
    profile = " ";
  copy = malloc(sizeof(char)*(strlen(profile)+strlen(origprof)+4));
  //  if (strlen(profile)+strlen(origprof) > (ARRAYSIZE*size-4)) 
  //copy = resize(&size,strlen(profile)+strlen(origprof)+4,copy);

  sprintf(copy,"%s->%s",profile,origprof);
  //printf("edge %s->%s\n",profile,origprof);

  if (hashtbl_get(edges,copy)) return;

  hashtbl_insert(edges,copy,"1");
  if (!(hash = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() failed");
    exit(EXIT_FAILURE);
  }
  diff = find_diff(hash,profile,origprof,&k1,&k2);
  brac = edge_label(hash,profile,origprof,k2);
  //color=gray%.0f,100-(ratio*100);
  if (strlen(brac) > 1)
    fprintf(fp,"\"%s\"-> \"%s\" [label=\"  %s\\n%s  \",fontsize=8];\n",profile,origprof,brac,diff);
  //if (VERBOSE)
  //printf("inserting edge %s\n",copy);

  free(copy);
  free(diff);
  free(brac);
  hashtbl_destroy(hash);
}  

//finds helical difference between profile -> origprof
//k1 and k2 are num of helices in profile and origprof
//return difference with triplets
char* find_diff(HASHTBL *hash,char *profile, char *origprof, int *k1, int *k2) {
  int size = INIT_SIZE,b1,b2,xor,i;
  char *diff,*id;

  diff = malloc(sizeof(char)*ARRAYSIZE*size);
  diff[0] = '\0';
  b1 = make_binary(profile,k1);
  b2 = make_binary(origprof,k2);
  xor = b1 ^ b2;
  //printf("b1 is %d, b2 is %d, and xor is %d\n",b1,b2,xor);
  for (i = 0; xor > 0; xor>>=1,i++) {
    if ((xor & 1)==1) {
      hashtbl_insert(hash,table[i],"1");
      id = hashtbl_get(idhash,table[i]);
      
      if (strlen(diff)+strlen(table[i])+strlen(id)+6 > ARRAYSIZE*size)
	diff = resize(&size,strlen(diff)+strlen(table[i])+strlen(id)+6,diff);
      if (strlen(diff)>1)
	strcat(diff,"\\n");
      sprintf(diff,"%s%s: %s",diff,table[i],id);
    }
  }
  //printf("diff is %s between %s and %s\n",diff,origprof,profile);

  return diff;
}

//finds the bracket notation for profile->origprof based on originating profile of size k
//difference between profile and origprof is in hash
//make sure origprof has a bracket rep (should be done in make_bracket_rep)
char* edge_label(HASHTBL *hash,char *profile, char *origprof,int k) {
  int i=0,count = 0,j=0,num=0,*diff,*save,ind=0,m=0;
  char *origbrac,*brac,*blank = "[]",*copy,*val;
  char **array;

  if (!(origbrac = hashtbl_get(bracket,origprof))) {
    fprintf(stderr,"Error: origprof %s has no bracket representation in bracket: edge_label()\n",origprof);
    return "";
  }
  //  if (strlen(profile) == 0) 
  copy = mystrdup(origbrac);
  brac = malloc(strlen(origbrac));
  //printf("finding edge label between %s and %s\n",origbrac,profile);
  array = malloc(sizeof(char*)*k);
  diff = malloc(sizeof(int)*hashtbl_numkeys(hash));
  save = malloc(sizeof(int)*hashtbl_numkeys(hash));
  //put helices of origprof into array; array in chron order
  for (val = strtok(copy,blank); val; val = strtok(NULL,blank)) {
    //printf("val is %s\n",val);
    if (i >= k) fprintf(stderr,"mismatch between k=%d and number of helices %d for %s: edge_label\n",k,i,origprof);
    array[i++] = val;
    //printf("val %s is at %d\n",val,i-1);
  }
  //save index in origprof of all diff helices; in ascending order
  for (count = 0; count < i; count++) {
    if (hashtbl_get(hash,array[count])) {
      //printf("saving %s at %d to index %d\n",array[count],count,j);
      diff[j++] = count;
    }
  }
  copy[0] = '\0';
  brac[0] = '\0';
  count = 0;
  j = -1;
  //i is index for origbrac, j is index for what level we need to match ']'
  //ind is index for level that increases for '[' and decreases for ']'
  //num is number of '[' encountered
  //count is index of different helix being proecessed
  //m is index for array of origprof helices, to print out ones not in diff
  for (i = 0; origbrac[i] != '\0'; i++) {
    //keep track of how many '['s
    if (origbrac[i] == '[') {
      ind++;
      val = mystrdup(array[m]);
      if (diff[count] == num++) {
	count++;
	save[++j] = ind;
	//printf("\nsaving %d to j=%d\n",ind,j);
	strcat(copy,"{");

      }
      else {
	//strcat(brac,"[");
	//strcat(brac,array[m]);
	sprintf(brac,"%s[%s",brac,array[m]);
	strcat(copy,"[");
      }
      strcat(copy,val);
      free(val);
      m++;
    }
    //keep track of which level brackets are at
    else if (origbrac[i] == ']') {
      //printf("\nchecking ind %d against %d at %d\n",ind,save[j]  puts("after check insert edge");,j);
      if (j >= 0 && save[j] == ind--) {
	j--;
	strcat(copy,"}");
      }
      else {
	strcat(brac,"]");
	strcat(copy,"]");
      }
    }
  }
  if (!(hashtbl_get(bracket,profile))) {
    hashtbl_insert(bracket,profile,brac);
    //printf("new brac is %s for (%s ->) %s with copy %s\n",brac,origprof,profile,copy);
  }
  free(array);
  free(diff);
  free(save);
  return copy;
}

int make_binary(char *profile,int *k) {
  int sum=0,*bin = NULL;
  char *blank = " ";
  char *copy = strdup(profile);
  char *helix;
  for (helix = strtok(copy,blank); helix; helix = strtok(NULL,blank)) {
    bin = hashtbl_get(binary,helix);
    (*k)++;
    if (!bin) continue;
    sum += *bin;

  }
  free(copy);
  return sum;
}

char* print_edge(KEY *node,char **table,int v,int* sum) {
  int count = 0,val;
  char *v1 = malloc(sizeof(char)*30);
  
  v1[0] = '\0';
  if (v == 0) {
    printf("%s-- %s\n",node->data,"E ");
    return v1;
  }

  for (val = (v & 1); v > 0; v >>= 1, count++, val = (v & 1))
    if (val == 1) {
      *sum += *((int*) hashtbl_get(binary,table[count]));
      strcat(v1,table[count]);
      strcat(v1," ");
    }
  
  printf("%s-- %s\n",node->data,v1);
  return v1;
}

