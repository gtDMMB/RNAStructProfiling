/* Copyright (c) 2011 the authors listed at the following URL, and/or
the authors of referenced articles or incorporated external code:
http://en.literateprograms.org/Hash_table_(C)?action=history&offset=20100620072342

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Retrieved from: http://en.literateprograms.org/Hash_table_(C)?oldid=16749
*/

#ifndef HASHTBL_H_INCLUDE_GUARD
#define HASHTBL_H_INCLUDE_GUARD

#include<stdlib.h>

typedef size_t hash_size;

struct hashnode_s {
  char *key;
  void *data;
  struct hashnode_s *next;
};

typedef struct keynode {
  char *data;
  struct keynode *next;
} KEY;

typedef struct hashtbl {
  hash_size size;
  struct hashnode_s **nodes;
  hash_size (*hashfunc)(char *);
  KEY *begin;
  int numkeys;
} HASHTBL;

char *mystrdup(char *s);
HASHTBL *hashtbl_create(hash_size size, hash_size (*hashfunc)(char *));
void hashtbl_destroy(HASHTBL *hashtbl);
int hashtbl_insert(HASHTBL *hashtbl, char *key, void *data);
int hashtbl_remove(HASHTBL *hashtbl, char *key);
void *hashtbl_get(HASHTBL *hashtbl, char *key);
KEY* hashtbl_getkeys(HASHTBL *hashtbl);
int hashtbl_numkeys(HASHTBL *hashtbl);
int hashtbl_resize(HASHTBL *hashtbl, hash_size size);


#endif

