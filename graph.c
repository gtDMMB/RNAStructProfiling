#include "graph.h"
#include "hashtbl.h"
#include "Set.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int GRAPHSIZE;
node** graph;

node* createNode(char *name)
{
  node* newNode;
  newNode = malloc(sizeof(node));
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

void initialize(Set *set) {
  int i,size=1;
  node* root;
  char **diff;

  root = createNode(" ");
  while (set->num_sprof > ARRAYSIZE*size) size++;
  node **neighbors = malloc(sizeof(node*)*ARRAYSIZE*size);
  diff = malloc(sizeof(char*)*ARRAYSIZE*size);
  for (i = 0; i < set->num_sprof; i++) {
    neighbors[i] = createNode(set->profiles[i]->profile);
    neighbors[i]->sum = binary_rep(set,neighbors[i]->label);
    diff[i] = neighbors[i]->label;
    //printf("making %s node\n", neighbors[i]->label);
  }
  root->neighbors = neighbors;
  root->nsize = size;
  root->numNeighbors = set->num_sprof;
  root->diff = diff;
  set->graph = root;
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
  return sum;
}

/* finished profiles = [0 start-1]
   present LCA to be intersected = [start oldk]
   newly generated LCA = [oldk k]
   returns k, the number of vertices in graph
 */
void find_LCAs(FILE *fp,Set *set) {
  int new,oldk,go,start,size,k;
  unsigned long num;
  char *profile,**diff;
  node **vertices = set->graph->neighbors;

  k = set->num_sprof;
  size = set->graph->nsize;
  diff = set->graph->diff;
  start = 0;
  for (oldk = k; start != k; oldk = k) {
    for (new = start; new != oldk; new++) {
      for (go = advance(new,oldk); go != start; go = advance(go,oldk)) {
	num = vertices[new]->sum & vertices[go]->sum;
	//printf("num is %u of s[%d] = %u and s[%d] = %u\n",num,new,sums[new],go,sums[go]);
	if (not_in_sums(num,k,vertices)) {
	  //printf("found new profile for %s and %s\n",modprofileID[new],modprofileID[go]);

	  profile = convert_binary(num);
	  fprintf(fp,"\"%s\" [style = dashed];\n",profile);
	  if (k >= ARRAYSIZE*(size)) {
	    vertices = realloc(vertices,sizeof(node*)*ARRAYSIZE*++size);
	    diff = realloc(diff,sizeof(char*)*ARRAYSIZE*size);
	  }
	  vertices[k] = createNode(profile);
	  vertices[k]->sum = num;
	  diff[k] = profile;
	  //printf("k is %d\n",k);

	  found_edge(vertices[new],vertices[k]);	  
	  found_edge(vertices[go],vertices[k]);	  
	  k++;
	} else if (num == vertices[new]->sum) {
	  found_edge(vertices[go],vertices[new]);	  
	}
	else if (num == vertices[go]->sum)
	  found_edge(vertices[new],vertices[go]);	  
      } 
    }
    start = oldk;
  }
  set->graph->neighbors = vertices;
  set->graph->diff = diff;
  set->graph->numNeighbors = k;
  printf("graph size is now %d\n",k);
  set->graph->nsize = size;
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
  unsigned long xor;
  char *diff;

  //printf("child is %s and parent is %s\n",childprof,parentprof);
  for (i = 0; i < parent->numNeighbors; i++) 
    if (!strcmp(parent->neighbors[i]->label,child->label)) 
      return;
  if (parent->numNeighbors >= ARRAYSIZE*parent->nsize) {
    parent->neighbors = realloc(parent->neighbors,sizeof(node*)*ARRAYSIZE*(parent->nsize++));
    parent->diff = realloc(parent->diff,sizeof(char*)*ARRAYSIZE*parent->nsize);
  }
  else if (parent->numNeighbors == 0) {
    parent->neighbors = malloc(sizeof(node*)*ARRAYSIZE*parent->nsize);
    parent->diff = malloc(sizeof(char*)*ARRAYSIZE*parent->nsize);
  }
  xor = child->sum ^ parent->sum;
  diff = convert_binary(xor);
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

  sum = malloc(sizeof(unsigned long)*set->prof_num);
  for (i = 0; i < set->prof_num; i++) {
    sum[i] = binary_rep(set,set->profiles[i]->profile);
    if (!strcmp(set->profiles[i]->profile," "))
      set->graph->sfreq = set->profiles[i]->freq;
  }

  for (i = 0; i < set->graph->numNeighbors; i++) {
    vert = set->graph->neighbors[i];
    graph[i] = vert;
    for (j = 0; j < set->prof_num; j++) {
      //printf("investigating %s and prof[%d] = %s\n",node->data,i,profileID[i]);
      //printf("%u & %u = %u\n",sum[j],sums[i],sum[j]&sums[i]);
      if (sum[j] == vert->sum) {
	vert->sfreq = set->profiles[j]->freq;
	vert->bracket = set->profiles[j]->bracket;
      }	
      if ((sum[j] & vert->sum) == vert->sum) 
	vert->gfreq += set->profiles[j]->freq;
    }
    if (vert->sfreq == 0)
      fprintf(fp,"\"%s\" [label = \"0/%d\"];\n",vert->label,vert->gfreq);
    else
      fprintf(fp,"\"%s\" [label = \"%s\\n%d/%d\"];\n",vert->label,vert->bracket,vert->sfreq,vert->gfreq);
  }
  fprintf(fp,"\" \" [label = \"%d/%d\"];\n",set->graph->sfreq,set->opt->NUMSTRUCTS);
  graph[i] = set->graph;
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
    for (child = hashtbl_get(deleteHash,parent->data); child; child = child->next) {
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
      root->diff[probe] = root->diff[++probe];
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
    for(i = 0; i < GRAPHSIZE; i++)
    {
        free(graph[i]->neighbors);
    }
    free(graph);
}
