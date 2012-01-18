#include <stdio.h>

int insert_graph(FILE *fp,char *file,int gsize);
unsigned long insert_and_binary(char *key,char *profile,int freq);
void make_key();
int find_LCAs(FILE *fp,int k, int *size);
int advance(int new, int oldk);
int not_in_sums(unsigned long num, int k);
void insert_prof(int k, char *profile);
int found_edge(FILE *fp,int child,int parent);
void calc_gfreq(FILE *fp,int total);
char* convert_binary(unsigned long binary,int *count);
void print_edges(FILE *fp);
int* process_native(int i, int j, int k);
int process_one_input(FILE *fp);
void make_edge_and_node(FILE *fp,char *from, char *to,char *diff,int fullnum);
void process_input(FILE *fp);
HASHTBL* process_input_profile(FILE *fp,HASHTBL *brac,char *fullprofile, int fullnum,char *profile,int numhelix,char *diff,int prob);
KEY* find_parents(char *profile);
char* insert_diff(HASHTBL *temp,char *diff);
char* sort_input(char *profile,int length);
void find_centroid_edges(FILE *fp);
int print_vertices(FILE *fp);
int find_edges(FILE *fp,char *profile, int *found, int count);
void check_insert_edge(FILE *fp,char *profile,char *origprof);
char* find_diff(HASHTBL *hash,char *profile, char *origprof, int *k1, int *k2);
char* edge_label(HASHTBL *hash,char *profile, char *origprof,int k);
int make_binary(char *profile,int *k);
char* print_edge(KEY *node,char **table,int v,int* sum);
