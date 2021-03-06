/*RNAStructProfiling -- Profiles RNA structures and produces a summary graph in graphviz format.
 * Copyright 2013, 2014, 2018 Emily Rogers
 *
 * This file is part of RNAStructProfiling.
 *
 * RNAStructProfiling is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RNAStructProfiling is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RNAStructProfiling.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef GRAPH_H
#define GRAPH_H

#define SIZE 4

#include "Set.h"

/* Node Struct 
typedef struct node
{
  char *label;
  struct node **neighbors;
  int numNeighbors;
  int nsize;
  int sfreq;
  int gfreq;
  char *bracket;
  char **diff;
  int DFS; 
  unsigned long sum;
} node;
*/ 

void init_graph(FILE *fp, Set *set); 
int initialize(Set *set);
char* rm_root(int j, char* label);
void print_input(FILE *fp,Set *set);
//unsigned long binary_rep(Set *set,char *profile);
void find_LCAs(FILE *fp,Set *set,int i);
int advance(int newk, int oldk);
int not_in_sums(unsigned long num, int k,node **vertices);
char* convert_binary(unsigned long binary);
void found_edge(node *child,node *parent) ;
void calc_gfreq(FILE *fp,Set *set);
void removeEdges(HASHTBL *deleteHash);
void removeEdge(int i, int j);
void make_oval_bracket(node *vert);
char* find_child_bracket(node *vert, node *child);
void print_edges(FILE *fp,Set *set);
void printGraph();
void freeGraph();

extern node** graph;
extern int GRAPHSIZE;

#endif
