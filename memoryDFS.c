#include "memoryDFS.h"
#include "graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

node** currentPath;
int position;

void start_trans_reductn(HASHTBL *hash,int numNodes)
{   
    int i;
    GRAPHSIZE = numNodes;
    generateGraph(hash);
    
    for(i = 0 ; i < GRAPHSIZE; i++)
    {
        if(!strcmp(graph[i].label," "))
            break;
    }
    /*printGraph();*/
    MemoryDFS(&graph[i]);
}


void MemoryDFS(node* root)
{
    int idx, num, counter, counter2;
    node* myroot;
    currentPath = (node**)malloc(sizeof(node*)*GRAPHSIZE);
    position = 0;
    /*myroot = root->neighbors[0];*/
    myroot = root;
    /*printf("begin\n");*/
	addToPath(myroot);
    num = myroot->numNeighbors;
    /*printNeighbors(myroot);*/
    counter = counter2 = 0;
    while(counter < num)
    {
        node* v;
        counter2 = 0;
		/*printf("while counter\n");*/
        for(idx = 0; idx < myroot->numNeighbors; idx++)
        {	
			/*printf("for numneighbors\n");*/
            /*printf("numNeighbors: %d, counter: %d, counter2: %d\n", myroot->numNeighbors, counter, counter2);*/
            if(!myroot->neighbors[idx]->DFS)
            {
				/*printf("if myroot -> neighbors -> DFS\n");*/
                myroot->neighbors[idx]->DFS = 1;
                counter++;
                v = &*(myroot->neighbors[idx]);
                if(findEdge(myroot, v, 0))
                {
					/*printf("if find edge\n");*/
                    /*printf("v: %s\n", v->label);*/
                    addToPath(v);
                    MDFSHelper(myroot, v);
                    removeFromPath(v);
                }
                break;
            }
            if(myroot->neighbors[idx]->DFS)
            {
                counter2++;
            }
        }
        if(counter2 == myroot->numNeighbors)
        {
            break;
        }
    }
    removeFromPath(myroot);
    free(currentPath);
}

void MDFSHelper(node* root, node* v)
{
    int i, j, index;
    /*printf("MDFSHepler root: %s, v: %s\n", root->label, v->label);*/
    /*printCurrentPath();*/
    /*printNeighbors(v);*/
    index = v->numNeighbors;
    for(i = 0; i < index; i++)
    {
        node* u;
        u = &*(v->neighbors[i]);
        if(findEdge(v, u, 0))
        {
            addToPath(u);
            /*printCurrentPath();*/
            for(j = 0; j < position-1; j++)
            {
                findEdge(currentPath[j], u, 1);
                MDFSHelper(root, u);
                removeFromPath(u);
            }
        }
    }
}

int findEdge(node* root, node* v, int remove)
{
    int i;
    /*printf("root: %s, v: %s, remove: %d\n", root->label, v->label, remove);*/
    for(i = 0; i < root->numNeighbors; i++)
    {
        if(!strcmp(root->neighbors[i]->label, v->label))
        {
            /*printf("edge exists\n");*/
            if(remove)
            {
                /*printNeighbors(root);*/
                /*printf("removing edge (%s, %s)\n", root->label, v->label);*/
                if(i < root->numNeighbors-1)
                {
                    int probe;
                    probe = i;
                    while(probe < root->numNeighbors - 1)
                    {
                        /*printf("BEFORE probe: %d, root->neighors[probe]: %s\n", probe, root->neighbors[probe]->label);*/
                        root->neighbors[probe] = root->neighbors[probe+1];
                        /*printf("AFTER probe: %d, root->neighors[probe]: %s\n", probe, root->neighbors[probe]->label);*/
                        probe++;
                    }
                }
                root->numNeighbors = root->numNeighbors - 1;
                /*printNeighbors(root);*/
            }
            return 1;
        }
    }
    return 0;
    
}

void removeFromPath(node* n)
{
    int i;
    for(i = 0; i < position; i++)
    {
        if(!strcmp(currentPath[i]->label,n->label))
        {
            if(i < position-1)
            {
                int probe;
                probe = i + 1;
                while(probe < position)
                {
                    currentPath[i] = currentPath[probe];
                    probe++;
                }
            }
            position--;
        }
    }
    /*printCurrentPath();*/
}

void addToPath(node* n)
{
    currentPath[position] = n;
    /* printf("currentPath[%d]: %s\n", 
           position,currentPath[position]->label); */
    position++;
}

void printCurrentPath()
{
    int i;
    printf("current path: ");
    for(i = 0; i < position; i++)
    {
        printf("%s ", currentPath[i]->label);
    }
    printf("\n");
}
void printNeighbors(node* n)
{
    int i;
    printf("neighbors of %s: ", n->label);
    for(i = 0; i < n->numNeighbors; i++)
    {
        printf("%s ", n->neighbors[i]->label);
    }
    printf("\n");
    
}
