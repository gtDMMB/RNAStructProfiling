#ifndef PROFILE_H
#define PROFILE_H

#include <stdio.h>

typedef struct profile {
  int freq;
  int genfreq;
  int selected;
  char *profile;
  char *bracket;
  //char *repstruct;
  //struct profile **children;
} Profile;

Profile* create_profile(char *profile);

#endif
