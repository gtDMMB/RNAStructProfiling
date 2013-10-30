#include <stdlib.h>
#include "Profnode.h"

struct profnode* makeProfnode(int *prof) {
  struct profnode *profile;

  profile = (struct profnode*) malloc(sizeof(Profnode));
  profile->prof = prof;
  profile->extended = NULL;
  profile->extnum = 0;
  profile->coverage = 0;
  profile->withNext = NULL;
  profile->woNext = NULL;
  profile->parent = NULL;
  return profile;
}
