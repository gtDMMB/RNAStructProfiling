#include "graph.h"
#include "hashtbl.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int GRAPHSIZE;
node* graph;
node* createNode(char *name)
{
	node* newNode;
	newNode = malloc(sizeof(node));
	newNode->label = name;
	newNode->neighbors = (node**)malloc(sizeof(node*)*GRAPHSIZE);
	newNode->numNeighbors = 0;
	newNode->DFS = 0;
	return newNode;
}

void generateGraph(HASHTBL *hash)
{	int i, j, boolean;
    KEY *key;

    //GRAPHSIZE = hashtbl_numkeys(hash);
    graph = (node*)malloc(sizeof(node)*GRAPHSIZE);
    key = hashtbl_getkeys(hash);
    i = 0;
    for(key=key; key; key=key->next)
    {       
      node *temp;
        char *name, *endpt;
	//printf("processing %s of length %d\n",key->data,strlen(key->data));
        name = (char*)malloc(sizeof(char)*strlen(key->data)+1);
        strcpy(name, key->data);
        name = strtok(name, "->");
        
        j = 0; boolean = 0;
	//printf("looking for %s\n",name);
        while(j < i && graph[j].label != NULL)
        {
	  //printf("checking for graph[%d].label = %s\n",j,graph[j].label);
            if(!strcmp(graph[j].label, name))
            {
                boolean = 1;
                break;
            }
            j++;
        }
        if(!boolean)
        {
	  temp = createNode(name);
            graph[i] = *temp;
	    //printf("inserting %s at position %d\n",name,i);
            i++;
	    
        }        
	
	while(name != NULL) {
	  endpt = name;
	  name = strtok( NULL, "->");
	}
	//printf("looking for %s\n",endpt);
        j = 0; boolean = 0;
        while(j < i && graph[j].label != NULL)
	  {
            if(!strcmp(graph[j].label, endpt))
            {
	      //printf("found %s at graph[%d]\n",endpt,j);
                boolean = 1;
                break;
            }
            j++;
        }
        if(!boolean)
        {
	  temp = createNode(endpt);
	  graph[i] = *temp;
	  //printf("inserting %s at position %d\n",endpt,i);
	  i++;
        }
    }
    addNeighbors(hash);    
}

void addNeighbors(HASHTBL *hash)
{
    char *name, *endpt;
    int i, j, position;
    KEY *key;
    
    key = hashtbl_getkeys(hash);
    for(key=key; key != NULL; key = key->next)
    {
        /*printf("key->data: %s\n", key->data);*/
        name = (char*)malloc(sizeof(char)*strlen(key->data)+1);
        strcpy(name, key->data);
        name = strtok(name, "->");
        /*printf("name: %s\n", name);*/
        for(i = 0; i < GRAPHSIZE; i++)
        {   if(graph[i].label != NULL)
        {
            if(!strcmp(graph[i].label, name))
            {
                while(name != NULL) {
                    endpt = name;
                    name = strtok( NULL, "->");
                }
                /*printf("endpt: %s\n", endpt);*/
                for(j = 0; j < GRAPHSIZE; j++)
                {
					
                    if(graph[j].label != NULL)
                    {
                        if(!strcmp(graph[j].label, endpt))
                        {
                            /*printf("graph[%d].label: %s\n", j, graph[j].label);*/
                            position = graph[i].numNeighbors;
                            graph[i].neighbors[position] = &(graph[j]);
                            /*printf("graph[%d].neighbors[%d].label: %s\n", i, j, graph[i].neighbors[position]->label);*/
                            graph[i].numNeighbors = graph[i].numNeighbors + 1;
                            /*printf("graph[%d].numNeighbors: %d\n",i , graph[i].numNeighbors);*/
                            break;
                        } 
                    }
                }
                break;
            }
        }
        }
    }
}

void printGraph()
{
	int i, j;
	for(i = 0; i < GRAPHSIZE; i++)
	{
        if(graph[i].label != NULL)
        {
            printf("node: '%s'\nneighbors: ", graph[i].label);
            for(j = 0; j < graph[i].numNeighbors; j++)
            {
                printf("'%s', ", graph[i].neighbors[j]->label);
            }
            printf("\n");
            /*printf("total neighbors: %d\n", graph[i].numNeighbors);*/
        }
	}	
}

void freeGraph()
{
    int i;
    for(i = 0; i < GRAPHSIZE; i++)
    {
        free(graph[i].neighbors);
    }
    free(graph);
}
