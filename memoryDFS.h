#ifndef MEMORYDFS_H
#define MEMORYDFS_H

#include "graph.h"
#include "hashtbl.h"

HASHTBL* MemoryDFS(node* root);
void MDFSHelper(node* root, node* v);

/* Returns 1 if the node was found, 0 otherwise */
int inCurrentPath(node* n);

/* Looks for edge between root and v
   if remove is 1, removes edge.
   Returns 1 if edge was removed, 0 otherwise */
int findEdge(node* root, node* v, int remove);

/* Adds node n to current path */
void addToPath(node* n);

/* Removes node n from current path */
void removeFromPath(node* n);

/* Prints the current path */
void printCurrentPath();

/* Prints neighbors of node n */
void printNeighbors(node* n);

extern node** currentPath;
extern int position;


#endif
