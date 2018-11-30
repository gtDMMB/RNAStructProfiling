/* RNAStructProfiling -- Profiles RNA structures and produces a summary graph in graphviz format.
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
