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

#include "hashtbl.h"

#include <string.h>
#include <stdio.h>

char *mystrdup(char *s)
{
	char *b;
	if(!(b=(char*)malloc(strlen(s)+1))) return NULL;
	strcpy(b, s);
	return b;
}


static hash_size def_hashfunc(char *key)
{
	hash_size hash=0;
	
	while(*key) hash+=(unsigned char)*key++;

	return hash;
}



HASHTBL *hashtbl_create(hash_size size, hash_size (*hashfunc)(char *))
{
	HASHTBL *hashtbl;

	if(!(hashtbl=(HASHTBL*) malloc(sizeof(HASHTBL)))) return NULL;

	if(!(hashtbl->nodes=(hashnode_s**) calloc(size, sizeof(struct hashnode_s*)))) {
		free(hashtbl);
		return NULL;
	}

	hashtbl->size=size;

	if(hashfunc) hashtbl->hashfunc=hashfunc;
	else hashtbl->hashfunc=def_hashfunc;

	hashtbl->begin = NULL;
        hashtbl->numkeys = 0;
	return hashtbl;
}


void hashtbl_destroy(HASHTBL *hashtbl)
{
	hash_size n;
	struct hashnode_s *node, *oldnode;
	KEY *next,*last;

	for (last = hashtbl->begin; last != NULL; last = next) {
	  next = last->next;
	  free(last);
	}
	for(n=0; n<hashtbl->size; ++n) {
	  node=hashtbl->nodes[n];
	  while(node) {
	    free(node->key);
	    oldnode=node;
	    node=node->next;
	    free(oldnode);
	  }
	}
	free(hashtbl->nodes);
	free(hashtbl);
}


int hashtbl_insert(HASHTBL *hashtbl, char *key, void *data)
{
	struct hashnode_s *node;
	hash_size hash=hashtbl->hashfunc(key)%hashtbl->size;
	KEY *temp = NULL;

/*	fprintf(stderr, "hashtbl_insert() key=%s, hash=%d, data=%s\n", key, hash, (char*)data);*/

	node=hashtbl->nodes[hash];
	while(node) {
		if(!strcmp(node->key, key)) {
			node->data=data;
			return 0;
		}
		node=node->next;
	}

	if(!(node=(struct hashnode_s*)malloc(sizeof(struct hashnode_s)))) return -1;
	if(!(node->key=mystrdup(key))) {
		free(node);
		return -1;
	}
	node->data=data;
	node->next=hashtbl->nodes[hash];
	hashtbl->nodes[hash]=node;
	
	temp = (KEY*) malloc(sizeof(KEY));
	temp->data = mystrdup(key);
	temp->next = hashtbl->begin;
	hashtbl->begin = temp;
	
        hashtbl->numkeys++;
	return 0;
}


int hashtbl_remove(HASHTBL *hashtbl, char *key)
{
	struct hashnode_s *node, *prevnode=NULL;
	hash_size hash=hashtbl->hashfunc(key)%hashtbl->size;
	KEY *current = NULL, *prev = NULL;

	node=hashtbl->nodes[hash];
	while(node) {
	  if(!strcmp(node->key, key)) {
	    free(node->key);
	    if(prevnode) prevnode->next=node->next;
	    else hashtbl->nodes[hash]=node->next;
	    free(node);

	    for (current = hashtbl->begin; current != NULL; current = current->next) {
	      if (!strcmp(current->data, key)) {
		free(current->data);
		if (prev)
		  prev->next = current->next;
		else
		  hashtbl->begin = current->next;
		free(current);
	      }
	      prev = current;
	    }
	    hashtbl->numkeys--;
	    return 0;
	  }
	  prevnode=node;
	  node=node->next;
	}
	
	return -1;
}


void *hashtbl_get(HASHTBL *hashtbl, char *key)
{
	struct hashnode_s *node;
	if (!hashtbl) fprintf(stderr, "hashtbl doesn't exist\n");
	hash_size hash=hashtbl->hashfunc(key)%hashtbl->size;

/*	fprintf(stderr, "hashtbl_get() key=%s, hash=%d\n", key, hash);*/

	node=hashtbl->nodes[hash];
	while(node) {
		if(!strcmp(node->key, key)) return node->data;
		node=node->next;
	}

	return NULL;
}

KEY* hashtbl_getkeys(HASHTBL *hashtbl) {
  return hashtbl->begin;
}

int hashtbl_numkeys(HASHTBL *hashtbl) {
    return hashtbl->numkeys;
}

int hashtbl_resize(HASHTBL *hashtbl, hash_size size)
{
	HASHTBL newtbl;
	hash_size n;
	struct hashnode_s *node,*next;

	newtbl.size=size;
	newtbl.hashfunc=hashtbl->hashfunc;

	if(!(newtbl.nodes=(struct hashnode_s**) calloc(size, sizeof(struct hashnode_s*)))) return -1;

	for(n=0; n<hashtbl->size; ++n) {
		for(node=hashtbl->nodes[n]; node; node=next) {
			next = node->next;
			hashtbl_insert(&newtbl, node->key, node->data);
			hashtbl_remove(hashtbl, node->key);
			
		}
	}

	free(hashtbl->nodes);
	hashtbl->size=newtbl.size;
	hashtbl->nodes=newtbl.nodes;

	return 0;
}

