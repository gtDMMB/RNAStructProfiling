#include <stdlib.h>
#include "Profile.h"

Profile* create_profile(char *profile) {
  Profile *prof = (Profile*) malloc(sizeof(Profile));
  prof->freq = 1;
  prof->genfreq = 0;
  prof->selected = 0;
  prof->profile = profile;
  prof->bracket = NULL;
  //prof->repstruct = NULL;
  //prof->children = NULL;
  return prof;
}

