#ifndef Helix_class
#define Helix_class

typedef struct helix_class {
  char *id;
  char *maxtrip;
  char *avetrip;
  int freq;
  int isfreq;
  unsigned long binary;
} HC;

HC* create_HC(int id,char *max);
void free_hc(HC* hc);

#endif
