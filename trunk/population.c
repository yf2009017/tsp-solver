/*
 * population.c
 *
 *  Created on: May 21, 2011
 *      Author: helder
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "common.h"
#include "population.h"

int generateRandomPopulation(tspsPopulation_t *pop, tspsMap_t *map, tspsConfig_t *config){
	int i;

	pop->numIndividuals = config->populationSize;

	if((pop->individual = (tspsIndividual_t*)malloc(sizeof(tspsIndividual_t) * pop->numIndividuals)) == NULL){
		return TSPS_RC_FAILURE;
	}

	for(i=0; i<pop->numIndividuals; i++){
		generateRandomChromossome(&pop->individual[i]);

		calculateFitness(&pop->individual[i], map);
	}


	return TSPS_RC_SUCCESS;
}

/* Generate a random chromossome*/
int generateRandomChromossome(tspsIndividual_t *ind){
	int isAlreadyUsed[NUM_NODES];
	short unsigned int chosenNode = 0;
	int i;

	for(i=0; i<NUM_NODES; i++){
		isAlreadyUsed[i] = 0;
	}

	// initialize random seed:
	//srand ( time(NULL) );

	for(i=0; i<NUM_NODES; i++){

		// choose random node of the map
		//chosenNode = lrand48() % NUM_NODES;
		chosenNode = rand() % NUM_NODES;

		//if this node was chosen already
		if(isAlreadyUsed[chosenNode] != 0){

			//choose the next one of the list
			while(isAlreadyUsed[chosenNode] != 0){
				if(chosenNode + 1 < NUM_NODES){
					chosenNode++;
				}else{
					chosenNode = 0;
				}
			}
		}

		isAlreadyUsed[chosenNode] = 1;
		ind->chromossome[i] = chosenNode;
	}

	//srand(1);
	return TSPS_RC_SUCCESS;
}

/* calculate the fitness for the entire chromossome*/
int calculateFitness(tspsIndividual_t *ind, tspsMap_t *map){
	int fitness = 0;
	int i;
	int currentNode, nextNode;

	for(i=0; i<NUM_NODES-1; i++){
		currentNode = ind->chromossome[i];
		nextNode = ind->chromossome[i+1];

		fitness += map->weights[currentNode][nextNode];
	}

	currentNode = ind->chromossome[NUM_NODES-1];
	nextNode = ind->chromossome[0];
	fitness += map->weights[currentNode][nextNode];

	ind->fitness = fitness;

	return TSPS_RC_SUCCESS;
}

/* sort population by fitness (greatest to lowest)*/
int sortPopulation(tspsPopulation_t *pop){

	qsort(pop->individual, pop->numIndividuals, sizeof(tspsIndividual_t), comparePopulation);

	return TSPS_RC_SUCCESS;
}

/* qsort function*/
int comparePopulation(const void * a, const void * b){
	tspsIndividual_t *a1 = (tspsIndividual_t*)a;
	tspsIndividual_t *b1 = (tspsIndividual_t*)b;

	return (a1->fitness - b1->fitness);
}

/* apply edge crossover over the entire population, generating a new one*/
int crossoverPopulation(tspsPopulation_t *pop, tspsMap_t *map, tspsConfig_t *config){
	short int *alreadyBred = NULL;
	int parent1, parent2;
	int i;

	//initialize flags of individuals who already bred
	if((alreadyBred = (short int*)malloc(sizeof(short int) * pop->numIndividuals) ) == NULL){
		return TSPS_RC_FAILURE;
	}
	for(i=0; i<pop->numIndividuals; i++){
		alreadyBred[i] = 0;
	}

	//the first 'numElitism' individuals will be skipped, will treat them as if they have already bred
	for(i=0; i<config->numElitism; i++){
		alreadyBred[i] = 1;
	}

	// initialize random seed:
	//srand ( time(NULL) );

	for(i=0; i<pop->numIndividuals-1; i++){
		if(alreadyBred[i] != 0)
			continue;

		//randomly select two individuals to breed
		parent1 = i;
		alreadyBred[parent1] = 1;

		parent2 = parent1 + (rand() % (pop->numIndividuals-i));

		//choose the next one of the list
		if(alreadyBred[parent2] != 0){
			while(alreadyBred[parent2] != 0){
				if(parent2 + 1 < pop->numIndividuals){
					parent2++;
				}else{
					parent2 = 0;
				}
			}
		}

		alreadyBred[parent2] = 1;

		if(edgeCrossover(&pop->individual[parent1], &pop->individual[parent2], &pop->individual[parent1], &pop->individual[parent2]) != TSPS_RC_SUCCESS){
			free(alreadyBred);
			return TSPS_RC_FAILURE;
		}

		if(calculateFitness(&pop->individual[parent1], map) != TSPS_RC_SUCCESS ||
			calculateFitness(&pop->individual[parent2], map) != TSPS_RC_SUCCESS){

			free(alreadyBred);
			return TSPS_RC_FAILURE;
		}


	}

	free(alreadyBred);
	return TSPS_RC_SUCCESS;
}

int edgeCrossover(tspsIndividual_t *s1, tspsIndividual_t *s2, tspsIndividual_t *par1, tspsIndividual_t *par2){
	tspsEdgeTable_t edgeTable;
	tspsIndividual_t son1, son2;
	int isAlreadyUsed[NUM_NODES];
	int chosenNode;
	int i, j;

	memset(son1.chromossome, 0, sizeof(son1.chromossome));
	memset(son2.chromossome, 0, sizeof(son2.chromossome));

	for(i=0; i<NUM_NODES; i++){

		//search for the edges in parent 1
		for(j=0; j<NUM_NODES; j++){
			if(par1->chromossome[j] != i)
				continue;

			if(j-1 >= 0){
				edgeTable.node[i].edges[0] = par1->chromossome[j-1];
			}else{
				edgeTable.node[i].edges[0] = par1->chromossome[NUM_NODES-1];
			}

			if(j+1 < NUM_NODES){
				edgeTable.node[i].edges[1] = par1->chromossome[j+1];
			}else{
				edgeTable.node[i].edges[1] = par1->chromossome[0];
			}
		}

		//search for the edges in parent 2
		for(j=0; j<NUM_NODES; j++){
			if(par2->chromossome[j] != i)
				continue;

			if(j-1 >= 0){
				edgeTable.node[i].edges[2] = par2->chromossome[j-1];
			}else{
				edgeTable.node[i].edges[2] = par2->chromossome[NUM_NODES-1];
			}

			if(j+1 < NUM_NODES){
				edgeTable.node[i].edges[3] = par2->chromossome[j+1];
			}else{
				edgeTable.node[i].edges[3] = par2->chromossome[0];
			}
		}

		//uses this same for to initialize this flag list
		isAlreadyUsed[i] = 0;
	}

	//the first node will always be choosed randomly
	chosenNode = rand() % NUM_NODES;

	for(i=0; i<NUM_NODES-1; i++){

		//apply the chosen node to the chromossome
		son1.chromossome[i] = chosenNode;
		isAlreadyUsed[chosenNode] = 1;

		//search through the edges for the next chosen node
		chosenNode = chooseNodeFromEdges(&edgeTable, &edgeTable.node[chosenNode], &isAlreadyUsed[0]);
	}
	son1.chromossome[NUM_NODES-1] = chosenNode;

	memset(isAlreadyUsed, 0 , sizeof(isAlreadyUsed));

	//the first node will always be choosed randomly
	chosenNode = rand() % NUM_NODES;

	for(i=0; i<NUM_NODES-1; i++){

		//apply the chosen node to the chromossome
		son2.chromossome[i] = chosenNode;
		isAlreadyUsed[chosenNode] = 1;

		//search through the edges for the next chosen node
		chosenNode = chooseNodeFromEdges(&edgeTable, &edgeTable.node[chosenNode], &isAlreadyUsed[0]);
	}
	son2.chromossome[NUM_NODES-1] = chosenNode;

	memcpy(s1->chromossome, son1.chromossome, sizeof(s1->chromossome));
	memcpy(s2->chromossome, son2.chromossome, sizeof(s2->chromossome));
	return TSPS_RC_SUCCESS;
}

int chooseNodeFromEdges(tspsEdgeTable_t *edgeTable, tspsNodeEdges_t *nodeEdges, int *isAlreadyUsed){
	int node1, node2;
	int numDiffEdges, shortestDiff;
	int chosenNode = -1;
	tspsNodeEdges_t *aux = NULL;
	int i, j, k;


	//search for common edges
	for(i=0; i<4; i++){
		node1 = nodeEdges->edges[i];

		if(isAlreadyUsed[node1] != 0)
			continue;

		//search for another edge equals to 'node'
		for(j=i+1; j<4; j++){
			node2 = nodeEdges->edges[j];

			if(isAlreadyUsed[node2] != 0)
				continue;

			//if there is a common edge, return it
			if(node1 == node2){
				return node1;
			}
		}
	}

	//search edge with the shortest list
	numDiffEdges = 0;
	shortestDiff = 99;
	//chosenNode = nodeEdges->edges[0];
	for(i=0; i<4; i++){

		node1 = nodeEdges->edges[i];

		if(isAlreadyUsed[node1] != 0)
			continue;

		aux = &edgeTable->node[node1];
		numDiffEdges = 0;

		for(j=0; j<4; j++){
			for(k=j+1; k<4; k++){
				if(aux->edges[j] != aux->edges[k])
					numDiffEdges++;
			}
		}

		if(numDiffEdges < shortestDiff){
			chosenNode = node1;
			shortestDiff = numDiffEdges;
		}
	}

	if(chosenNode != -1)
		return chosenNode;

	//choose the next node randomly
	chosenNode = rand() % NUM_NODES;
	if(isAlreadyUsed[chosenNode] != 0){

		//choose the next one of the list
		while(isAlreadyUsed[chosenNode] != 0){
			if(chosenNode + 1 < NUM_NODES){
				chosenNode++;
			}else{
				chosenNode = 0;
			}
		}
	}

	return chosenNode;
}

int mutatePopulation(tspsPopulation_t *pop, tspsMap_t *map, tspsConfig_t *config){
	tspsIndividual_t *ind = NULL;
	int alreadySwaped[NUM_NODES];
	int mutationRate = config->mutationRate * 100;
	int index1, index2;		//the index of the nodes to be swapped
	int aux;
	int i, j;

	for(i=config->numElitism; i<pop->numIndividuals; i++){
		if(rand()%100 > mutationRate)
			continue;

		memset(alreadySwaped, 0, sizeof(alreadySwaped));
		ind = &pop->individual[i];

		//mutate!
		//swap mutationSize nodes in the chromossome
		for(j=0; j<config->mutationSize; j++){
			index1 = rand() % NUM_NODES;

			//if already swaped, jump to the next of the list
			while(alreadySwaped[index1] !=0){
				if(index1 + 1 < NUM_NODES)
					index1++;
				else
					index1 = 0;
			}
			alreadySwaped[index1] = 1;

			index2 = rand() % NUM_NODES;

			//if already swaped, jump to the next of the list
			while(alreadySwaped[index2] !=0){
				if(index2 + 1 < NUM_NODES)
					index2++;
				else
					index2 = 0;
			}
			alreadySwaped[index2] = 1;

			//swap the nodes
			aux = ind->chromossome[index1];
			ind->chromossome[index1] = ind->chromossome[index2];
			ind->chromossome[index2] = aux;
		}

		//recalculate the fitness
		if(calculateFitness(ind, map) != TSPS_RC_SUCCESS){
			return TSPS_RC_FAILURE;
		}
	}

	return TSPS_RC_SUCCESS;
}

void printIndividual(tspsIndividual_t *ind, unsigned long int numGenerations){
	FILE *file = NULL;
	int i;

	if((file = fopen("tsps.log", "a")) == NULL){
		return;
	}

	fprintf(file, "*** Best Individual\n");
	fprintf(file, "\t --> generat = [%lu]\n", numGenerations);
	fprintf(file, "\t --> fitness = [%d]\n", ind->fitness);
	fprintf(file, "\t --> chromos = ");

	for(i=0; i<NUM_NODES; i++){
		fprintf(file, "[%d]", ind->chromossome[i]);
		if((i+1) % 10 == 0)
			fprintf(file, "\n\t\t\t");
	}
	fprintf(file, "\n");
	fprintf(file, "------------------------------------------------------------------\n");

	fclose(file);
}
