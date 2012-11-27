#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "helix_class.h"
#include "Set.h"

HC* create_HC(int id, char* max) {
  HC *hc = (HC*) malloc(sizeof(HC));
  char *key = (char*) malloc(sizeof(char)*ARRAYSIZE);
  sprintf(key,"%d",id);
  hc->id = key;
  hc->maxtrip = max;
  hc->avetrip = NULL;
  hc->freq = 1;
  hc->isfreq = 0;
  hc->binary = 0;
  return hc;
}


