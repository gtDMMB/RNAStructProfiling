#include "graph.h"
#include "hashtbl.h"
#include "Set.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int GRAPHSIZE;
node** graph;


void init_graph(FILE *fp, Set *set) {
  int i;

  fputs("digraph G {\n",fp);
  //fprintf(fp,"\tlabel = \"%s\";\n",set->opt->OUTPUT);
  fprintf(fp,"\tpad = 0.5;\n");
  fprintf(fp,"\tnodesep = 0.5;\n");
  fprintf(fp,"\"legend\" [label = < <table border=\"0\" cellborder=\"1\" cellspacing=\"0\"><tr><td>Helix</td><td>Triplet</td><td>Frequency</td></tr>\n");
  for (i = 0; i < set->num_fhc; i++) {
    fprintf(fp,"<tr><td>%d</td><td>%s</td><td>%d</td></tr>\n",i+1,set->helices[i]->maxtrip,set->helices[i]->freq);
  }
  fprintf(fp,"</table>>, shape = plaintext, fontsize=11];\n");

  for (i = 0; i < set->num_sprof; i++) 
    fprintf(fp,"\"%s\" [shape = box];\n",set->profiles[i]->profile);
}

int initialize(Set *set) {
  int i,j,k,found=0,size=1;
  node *root;
  char **diff,*rt,*num;
  
  for (j = 0; j < set->hc_num; j++) {
    if (set->helices[j]->freq < set->opt->NUMSTRUCTS)
      break;
  }
  if (j>0) {
    //assuming <100 always helices
    rt = (char*)malloc(sizeof(char)*j*3);
    sprintf(rt,"%d ",1);
    for (k=1; k<j; k++) {
      num = (char*)malloc(sizeof(char)*3);
      sprintf(num,"%d ",k+1);
      strcat(rt,num);
    }
    root = createNode(rt);
    root->sum = binary_rep(set,root->label);
  } else
    root = createNode((char*)"");
  for (i = 0; i < set->prof_num; i++) 
    if (!strcmp(root->label,set->profiles[i]->profile)) 
      root->sfreq = set->profiles[i]->freq;
  
  while (set->num_sprof +1 > ARRAYSIZE*size) size++;
  node **neighbors = (node**)malloc(sizeof(node*)*ARRAYSIZE*size);
  diff = (char**)malloc(sizeof(char*)*ARRAYSIZE*size);
  for (i = 0; i < set->num_sprof; i++) {
    /*
    if (!strcmp(root->label,set->profiles[i]->profile)) {
      root->sfreq = set->profiles[i]->freq;
      continue;
    }
    */
    if (set->opt->INPUT && !strcmp(set->profiles[i]->profile,set->inputnode->label)) {
      neighbors[i] = set->inputnode;
      found = 1;
    }
    else
      neighbors[i] = createNode(set->profiles[i]->profile);
    neighbors[i]->sum = binary_rep(set,neighbors[i]->label);
    if (j>0) { //create labels from root, offset past root
      diff[i] = rm_root(j,neighbors[i]->label);
    } else
      diff[i] = neighbors[i]->label;
  }
  //  root = createNode((char*)" ");
  if (set->opt->INPUT && !found) {
    neighbors[i] = set->inputnode;
    neighbors[i]->sum = binary_rep(set,neighbors[i]->label);
    diff[i] = set->inputnode->label;
    root->numNeighbors = set->num_sprof+1;
  }
   else
    root->numNeighbors = set->num_sprof;

  root->neighbors = neighbors;
  root->nsize = size;
  //root->numNeighbors = set->num_sprof;
  root->diff = diff;
  set->graph = root;
  return j;
}

//removes first j helices from label
char* rm_root(int j, char* label) {
  int k,size;
  char *num, *val;
  k = 1;
  size = 0;
  num = mystrdup(label);
  for (val = strtok(num," "); val; val = strtok(NULL," ")) {
    size += strlen(val)+1;
    //printf("%s: val is %s and num %s\n",label,val,num+size);
    if (++k > j) break;
  }
  return num+size;
}

void print_input(FILE *fp,Set *set) {
  node *vert;
  
  fprintf(fp,"\"%s\" [style=filled, fillcolor=gray60];\n",set->inputnode->label);
  if (set->inputnode->numNeighbors == 0) return;
  vert = set->inputnode->neighbors[0];
  //printf("making bracket for %s\n",vert->label);
  if (!(vert->bracket)) {
    make_oval_bracket(vert);
  }
  fprintf(fp,"\"%s\" [label = \"%s\\n0/0\",style=filled, fillcolor=gray60];\n",vert->label,vert->bracket);
  if (vert->numNeighbors == 0) return;
  fprintf(fp,"\"%s\" -> \"%s\" [label = \"%s\", arrowhead = vee];\n",vert->label,vert->neighbors[0]->label,vert->diff[0]);
  vert = vert->neighbors[0];
  fprintf(fp,"\"%s\" [label = \"%s\\n0/0\",style=filled, fillcolor=gray60];\n",vert->label,vert->bracket);
}
/*
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
  return sum;
}
*/
/* finished profiles = [0 start-1]
   present LCA to be intersected = [start oldk]
   newly generated LCA = [oldk k]
   returns k, the number of vertices in graph
 */
void find_LCAs(FILE *fp,Set *set, int i) {
  int newk, oldk, go, start, size, k, cycles=0;
  unsigned long num;
  char *profile,**diff;
  node **vertices = set->graph->neighbors;

  k = set->graph->numNeighbors;
  size = set->graph->nsize;
  diff = set->graph->diff;
  start = 0;
  for (oldk = k; start != k; oldk = k) {
    for (newk = start; newk != oldk; newk++) {
      for (go = advance(newk,oldk); go != start; go = advance(go,oldk)) {
	//printf("start is %d oldk is %d and go is %d\n",start,oldk,go);
	num = vertices[newk]->sum & vertices[go]->sum;
	//printf("num is %u of s[%d] = %u and s[%d] = %u\n",num,new,sums[new],go,sums[go]);
	if (not_in_sums(num,k,vertices)) {
	  //printf("found new profile for %s and %s\n",modprofileID[new],modprofileID[go]);

	  profile = convert_binary(num);
	  fprintf(fp,"\"%s\" [style = dashed];\n",profile);
	  if (k >= ARRAYSIZE*(size)) {
	    vertices = (node**) realloc(vertices,sizeof(node*)*ARRAYSIZE*++size);
	    diff = (char**) realloc(diff,sizeof(char*)*ARRAYSIZE*size);
	  }
	  vertices[k] = createNode(profile);
	  vertices[k]->sum = num;
	  if (i>0) 
	    diff[k] = rm_root(i,profile);
	  else 
	    diff[k] = profile;
	  //printf("k is %d\n",k);

	  found_edge(vertices[newk],vertices[k]);	  
	  found_edge(vertices[go],vertices[k]);	  
	  k++;
	} else if (num == vertices[newk]->sum) {
	  found_edge(vertices[go],vertices[newk]);	  
	}
	else if (num == vertices[go]->sum)
	  found_edge(vertices[newk],vertices[go]);	  
      }
      //printf("comparings against %d with end %d, vertices %d\n",new,oldk,k);
    }
    start = oldk;
    if (++cycles == set->opt->CYCLES)
      break;
    //printf("Cycle %d with %d vertices\n",cycles,k);
  }
  /*  if (set->opt->VERBOSE) 
    printf("Total cycles %d\n",cycles);
  */
  set->graph->neighbors = vertices;
  set->graph->diff = diff;
  set->graph->numNeighbors = k;
  if (set->graph->sfreq == 0) k++;
  printf("Total number of vertices: %d\n",k);
  set->graph->nsize = size;
}

//wraps around like mod function
int advance(int newk, int oldk) {
  int next = newk+1;
  if (next == oldk) {
    //qputs("returning 0\n");
    return 0;
  }
  else
    return next;
}

//returns 1 if num doesn't match anything in sums up to oldk
//returns 0 otherwise
int not_in_sums(unsigned long num, int k,node **vertices) {
  int i;
  for (i = 0; i < k; i++) {
    if (vertices[i]->sum == num)
      return 0;
  }
  return 1;
}

//converts binary rep to string of helices (profile)
char* convert_binary(unsigned long binary) {
  int k,size = 1,num=0;
  char val[ARRAYSIZE],*profile;

  profile = (char*) malloc(sizeof(char)*ARRAYSIZE);
  profile[0] = '\0';
  for (k = 0; binary > 0; binary >>= 1, k++) {
    //printf("binary is %u\n",binary);
    if ((binary & 1) == 1) {
      sprintf(val,"%d",k+1);
      if (strlen(profile)+strlen(val) > (unsigned int)ARRAYSIZE*size-2) {
	profile = (char*) realloc(profile,sizeof(char)*ARRAYSIZE*++size);
      }
      //printf("adding %s, with k %d, binary is %u, shifted is %u\n",table[k],k,binary,binary>>1);
      strcat(profile,val);
      strcat(profile," ");
      num++;
    }
  }
  return profile;
}

  /*make edge
  update parents' general freq in graph[] with child's spec. freq if exists
  helix set difference = edge label
  child and parent are indices
  returns 1 if inserted into edge hash, 0 otherwise
  */
void found_edge(node *child,node *parent) {
  int i;
  unsigned long xr;
  char *diff;

  //printf("child is %s and parent is %s\n",childprof,parentprof);
  for (i = 0; i < parent->numNeighbors; i++) 
    if (!strcmp(parent->neighbors[i]->label,child->label)) 
      return;
  if (parent->numNeighbors >= ARRAYSIZE*parent->nsize) {
    parent->nsize++;
    parent->neighbors = (node**) realloc(parent->neighbors,sizeof(node*)*ARRAYSIZE*parent->nsize);
    parent->diff = (char**) realloc(parent->diff,sizeof(char*)*ARRAYSIZE*parent->nsize);
  }
  else if (parent->numNeighbors == 0) {
    parent->neighbors = (node**) malloc(sizeof(node*)*ARRAYSIZE*parent->nsize);
    parent->diff = (char**) malloc(sizeof(char*)*ARRAYSIZE*parent->nsize);
  }
  xr = child->sum ^ parent->sum;
  diff = convert_binary(xr);
  parent->neighbors[parent->numNeighbors] = child;
  parent->diff[parent->numNeighbors] = diff;
  parent->numNeighbors++;

  //fprintf(fp,"\"%s\" -> \"%s\" [label = \"%s\", arrowhead = vee];\n",parentprof,childprof,diff);
  //printf("Found %d is child of %d\n",child,parent);
}

 void calc_gfreq(FILE *fp,Set *set) {
  int i,j;
  unsigned long *sum;
  node *vert;

  GRAPHSIZE = set->graph->numNeighbors + 1;
  graph = (node**)malloc(sizeof(node*)*GRAPHSIZE);
  sum = (long unsigned int*) malloc(sizeof(unsigned long)*set->prof_num);
  for (i = 0; i < set->prof_num; i++) {
    sum[i] = binary_rep(set,set->profiles[i]->profile);
    if (!strcmp(set->profiles[i]->profile," "))
      set->graph->sfreq = set->profiles[i]->freq;
  }

  for (i = 0; i < set->graph->numNeighbors; i++) {
    vert = set->graph->neighbors[i];
    graph[i] = vert;
    for (j = 0; j < set->prof_num; j++) {
      if (sum[j] == vert->sum) {
	vert->sfreq = set->profiles[j]->freq;
	vert->bracket = set->profiles[j]->bracket;
      }	
      if ((sum[j] & vert->sum) == vert->sum) 
	vert->gfreq += set->profiles[j]->freq;
    }
    if (vert->sfreq == 0)
      make_oval_bracket(vert);
    fprintf(fp,"\"%s\" [label = \"%s\\n%d/%d\"];\n",vert->label,vert->bracket,vert->sfreq,vert->gfreq);
  }
  //if (set->graph->sfreq == 0) {
    if (!strcmp(set->graph->label,""))
      set->graph->bracket = (char*)"[]";
    else
      make_oval_bracket(set->graph);
    fprintf(fp,"\"%s\" [label = \"%s\\n%d/%d\"];\n",set->graph->label,set->graph->bracket,set->graph->sfreq,set->opt->NUMSTRUCTS);
    graph[i] = set->graph;
    //}
}

void make_oval_bracket(node *vert) {
  int i = 0,j=0,k,h=0,m=0,count = 0,*df,*skip,*helices;
  char *pbrac, *cbrac,*diff,*val;
  node *child;

  if (vert->bracket) return;
  child = (node*) malloc(sizeof(node));
  diff = find_child_bracket(vert,child);

  df = (int*) malloc(sizeof(int)*(strlen(diff)/2 + 1));
  for (val = strtok(mystrdup(diff)," "); val; val = strtok(NULL," ")) {
    df[i++] = atoi(val);
  }
  cbrac = mystrdup(child->bracket);
  pbrac = (char*) malloc(sizeof(char)*strlen(cbrac));
  skip = (int*) malloc(sizeof(int)*i);
  helices = (int*) malloc(sizeof(int)*(strlen(cbrac)/3 + 1));
  for (val = strtok(cbrac,"[]"); val; val = strtok(NULL,"[]")) {
    helices[j++] = atoi(val);
  }
  val = (char*) malloc(sizeof(char)*ARRAYSIZE);
  cbrac = child->bracket;
  pbrac[0] = '\0';
  for (j = 0; (unsigned int) j < strlen(cbrac); j++) {
    if (cbrac[j] == '[') {
      for (k=0; k < i; k++)
	if (df[k] == helices[h])
	  //printf("skipping %d\n",df[k]);
	  break;
      if (k == i) {
	sprintf(val,"[%d",helices[h]);
	pbrac = strcat(pbrac,val);
      } else
	skip[m++] = count;
      h++;
      count++;
    }
    
    else if (cbrac[j] == ']') {
      count--;
      if (m == 0 || count != skip[m-1])
	pbrac = strcat(pbrac,"]");
      else 
	m--;
    }
  }
  //printf("bracket based on %s for %s is %s\n",child->bracket,vert->label,pbrac);
  vert->bracket = pbrac;
  free(df);
  free(helices);
  free(skip);
  free(val);
  free(child);
}

char* find_child_bracket(node *vert, node *child) {
  int i = 0;
  char *nowdiff,*concat,*diff;

  if (vert->numNeighbors == 0)
    fprintf(stderr,"Oval should have children in find_child_bracket()\n");

  //printf("finding bracket child for %s\n",vert->label);
  for (i = 0; i < vert->numNeighbors; i++) {
    if (vert->neighbors[i]->bracket) {
      //printf("found bracket %s for %s\n",vert->neighbors[i]->bracket,vert->neighbors[i]->label);
      *child = *(vert->neighbors[i]);
      return vert->diff[i];
    }
  }
  //printf("going to child %s\n",vert->neighbors[0]->label);
  diff = find_child_bracket(vert->neighbors[0],child);
  nowdiff = vert->diff[0];
  concat = (char*) malloc(sizeof(char)*(strlen(nowdiff)+strlen(diff)+1));
  sprintf(concat,"%s%s",diff,nowdiff);
  return concat;
}

void print_edges(FILE *fp,Set *set) {
  int i, j;
  char *diff=NULL;

  for (i = 0; i < GRAPHSIZE; i++) {
    if (graph[i]->label != NULL) {
      if (set->opt->VERBOSE)
	printf("node: '%s'\n", graph[i]->label);
      for (j = 0; j < graph[i]->numNeighbors; j++) {
	diff = graph[i]->diff[j];
	if (diff) {
	  fprintf(fp,"\"%s\" -> \"%s\" [label = \"%s\", arrowhead = vee];\n",graph[i]->label,graph[i]->neighbors[j]->label,diff);
	} else {
	  fprintf(stderr, "no diff for %s and %s in print_edges()\n",graph[i]->label,graph[i]->neighbors[j]->label);
	}
	//printf("'%s', ", graph[i].neighbors[j]->label);
      }
    }
  }	
}

void removeEdges(HASHTBL *deleteHash) {
  int i,j;
  KEY *parent,*child;

  for (parent = hashtbl_getkeys(deleteHash); parent; parent = parent->next) {
    for (i = 0; i < GRAPHSIZE; i++) 
      if (!strcmp(graph[i]->label,parent->data)) break;
    if (i == GRAPHSIZE) fprintf(stderr,"didn't find %s in removeEdges\n",parent->data);
    for (child = (KEY*) hashtbl_get(deleteHash,parent->data); child; child = child->next) {
      for (j = 0; j < graph[i]->numNeighbors; j++)
	if (!strcmp(graph[i]->neighbors[j]->label,child->data)) break;
      if (j == graph[i]->numNeighbors) fprintf(stderr,"didn't find %s in neighbors of %s in removeEdges()\n",child->data,parent->data);
      removeEdge(i,j);
    }
  }
}

void removeEdge(int i, int j) {
  node *root = graph[i];

  //  printf("removing %s -> %s\n",root->label,root->neighbors[j]->label);
  if (j < root->numNeighbors-1) {
    int probe;
    probe = j;
    while(probe < root->numNeighbors - 1) {
      root->neighbors[probe] = root->neighbors[probe+1];
      root->diff[probe] = root->diff[probe+1];
      probe++;
    }
  }
  root->numNeighbors--;
}

void printGraph()
{
	int i, j;
	for(i = 0; i < GRAPHSIZE; i++)
	{
        if(graph[i]->label != NULL)
        {
            printf("node: '%s'\nneighbors: ", graph[i]->label);
            for(j = 0; j < graph[i]->numNeighbors; j++)
            {
                printf("'%s', ", graph[i]->neighbors[j]->label);
            }
            printf("\n");
            /*printf("total neighbors: %d\n", graph[i]->numNeighbors);*/
        }
	}	
}

void freeGraph()
{
    int i;
    for(i = 0; i < GRAPHSIZE; i++) {
      free(graph[i]->neighbors);
      free(graph[i]->bracket);
      free(graph[i]->diff);
    }
    free(graph);
}
