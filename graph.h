#ifndef GRAPH_H
#define GRAPH_H

#define SIZE 4

/* Node Struct */
typedef struct node
{
    /*int label;*/
  char *label;
  struct node **neighbors;
  int numNeighbors;
  int nsize;
  int sfreq;
  int gfreq;
  char *bracket;
  char **diff;
  int DFS; /* boolean value */
  unsigned long sum;
} node;

node* createNode();
void init_graph();
void initialize();
unsigned long binary_rep();
void find_LCAs();
int advance(int new, int oldk);
int not_in_sums(unsigned long num, int k,node **vertices);
char* convert_binary(unsigned long binary);
void found_edge(node *child,node *parent) ;
void calc_gfreq();
void removeEdges();
void removeEdge();
void make_oval_bracket(node *vert);
char* find_child_bracket(node *vert, node *child);
void print_edges();
void printGraph();
void freeGraph();

extern node** graph;
extern int GRAPHSIZE;

#endif
