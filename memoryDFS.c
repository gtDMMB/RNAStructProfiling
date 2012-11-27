#include "memoryDFS.h"
#include "graph.h"
#include "hashtbl.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

node** currentPath;
int position;
HASHTBL *deleteHash;

HASHTBL* MemoryDFS(node* root) {
  int idx, num, counter, counter2;
  node* myroot;
  currentPath = (node**)malloc(sizeof(node*)*GRAPHSIZE);
  position = 0;
  if (!(deleteHash = hashtbl_create(31,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for deleteHash failed");
    exit(EXIT_FAILURE);
  }
  /*myroot = root->neighbors[0];*/
  myroot = root;
  addToPath(myroot);
  num = myroot->numNeighbors;
  /*printNeighbors(myroot);*/
  counter = counter2 = 0;
  while(counter < num) {
    node* v;
    counter2 = 0;
    /*printf("while counter\n");*/
    for(idx = 0; idx < myroot->numNeighbors; idx++) {	
      /*printf("for numneighbors\n");*/
      /*printf("numNeighbors: %d, counter: %d, counter2: %d\n", myroot->numNeighbors, counter, counter2);*/
      if(!myroot->neighbors[idx]->DFS) {
	/*printf("if myroot -> neighbors -> DFS\n");*/
	myroot->neighbors[idx]->DFS = 1;
	counter++;
	v = &*(myroot->neighbors[idx]);
	if (findEdge(myroot, v, 0)) {
	  /*printf("if find edge\n");*/
	  /*printf("v: %s\n", v->label);*/
	  addToPath(v);
	  MDFSHelper(myroot, v);
	  removeFromPath(v);
	}
	break;
      }
      if (myroot->neighbors[idx]->DFS) {
	counter2++;
      }
    }
    if (counter2 == myroot->numNeighbors) {
      break;
    }
  }
  removeFromPath(myroot);
  free(currentPath);
  return deleteHash;
}

void MDFSHelper(node* root, node* v)
{
  int i, j, index;
  /*printf("MDFSHepler root: %s, v: %s\n", root->label, v->label);*/
  /*printCurrentPath();*/
  /*printNeighbors(v);*/
  index = v->numNeighbors;
  //printf("examining %s with %d neighbors\n",v->label,v->numNeighbors);
  for(i = 0; i < index; i++)
    {
      node* u;
      u = &*(v->neighbors[i]);
      //printf("looking at %s and n[%d]=%s\n",v->label,i,u->label);
      if(findEdge(v, u, 0))
        {
 	  //printf("found edge %s and n[%d]=%s\n",v->label,i,u->label);
	  addToPath(u);
	  /*printCurrentPath();*/
	  for(j = 0; j < position-1; j++)
            {
	      findEdge(currentPath[j], u, 1);
	      MDFSHelper(root, u);
	      removeFromPath(u);
            }
        }
    }
}

int findEdge(node* root, node* v, int remove) {
  int i;
  KEY *child;

  /*printf("root: %s, v: %s, remove: %d\n", root->label, v->label, remove);*/
  for(i = 0; i < root->numNeighbors; i++) {
    if(!strcmp(root->neighbors[i]->label, v->label)) {
      /*printf("edge exists\n");*/
      if(remove) {
	for (child = (KEY*)hashtbl_get(deleteHash,root->label); child; child = child->next)
	  if (!strcmp(v->label,child->data)) 
	    return 1;
	//printf("removing edge (%s, %s)\n", root->label, v->label);
	child = (KEY*)malloc(sizeof(KEY));
	child->data = v->label;
	child->next = (keynode*)hashtbl_get(deleteHash,root->label);
	hashtbl_insert(deleteHash,root->label,child);
      }
      return 1;
    }
  }
  return 0;
}

void removeFromPath(node* n)
{
  int i;
    for(i = 0; i < position; i++)
    {
        if(!strcmp(currentPath[i]->label,n->label))
        {
            if(i < position-1)
            {
                int probe;
                probe = i + 1;
                while(probe < position)
                {
                    currentPath[i] = currentPath[probe];
                    probe++;
                }
            }
            position--;
        }
    }
    /*printCurrentPath();*/
}

void addToPath(node* n)
{
    currentPath[position] = n;
    /* printf("currentPath[%d]: %s\n", 
           position,currentPath[position]->label); */
    position++;
}

void printCurrentPath()
{
    int i;
    printf("current path: ");
    for(i = 0; i < position; i++)
    {
        printf("%s ", currentPath[i]->label);
    }
    printf("\n");
}

void printNeighbors(node* n)
{
    int i;
    printf("neighbors of %s: ", n->label);
    for(i = 0; i < n->numNeighbors; i++)
    {
        printf("%s ", n->neighbors[i]->label);
    }
    printf("\n");
    
}
