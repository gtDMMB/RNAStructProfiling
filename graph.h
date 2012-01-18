#ifndef GRAPH_H
#define GRAPH_H

#include "hashtbl.h"

#define SIZE 4
#define HASHSIZE 31

/* Node Struct */
typedef struct node
{
    /*int label;*/
    char *label;
    struct node **neighbors;
    int numNeighbors;
    int DFS; /* boolean value */
} node;

void generateGraph(HASHTBL *hash);
node* createNode();
void addNode(node* n);
void addNeighbors(HASHTBL *hash);
void printGraph();
void freeGraph();

extern node* graph;
extern HASHTBL* hashgraph;
extern int GRAPHSIZE;

#endif
