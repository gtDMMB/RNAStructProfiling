#ifndef PROFNODE
#define PROFNODE

#include "Profile.h"

typedef struct profnode{
  //array of hc that make the prof
  int *prof;
  //number of hc in prof
  int profnum;
  //array of extended prof that make the prof
  int *extended;
  //length of extended
  int extnum;
  //sum of freq of extended profiles
  int coverage;
  struct profnode *withNext;
  struct profnode *woNext;
  struct profnode *parent;
} Profnode;

struct profnode* makeProfnode(int *prof);

#endif
