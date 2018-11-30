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
 * */

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


