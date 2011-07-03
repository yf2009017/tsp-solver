/*
 * population.c
 *
 *  Created on: May 21, 2011
 *      Author: helder
 */
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "mpi.h"
#include "common.h"
#include "population.h"

tspsImigrant_t imigrant;

int generateRandomPopulation(tspsPopulation_t *pop, tspsMap_t *map, tspsConfig_t *config){
	int i;

	pop->numIndividuals = config->populationSize;

	if((pop->individual = (tspsIndividual_t*)malloc(sizeof(tspsIndividual_t) * pop->numIndividuals)) == NULL){
		return TSPS_RC_FAILURE;
	}

	for(i=0; i<pop->numIndividuals; i++){
		generateRandomChromossome(&pop->individual[i]);

		calculateFitness(&pop->individual[i], map);
		pop->individual[i].breedChance = 0.0;
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
int sortPopulation(tspsPopulation_t *pop, int mpiId){
	int isCompleted = 0;
	MPI_Status status;

	if(imigrant.popIndex != -1){
		MPI_Test(&imigrant.req, &isCompleted, &status);
		if(isCompleted){
			memcpy(&pop->individual[imigrant.popIndex], &imigrant.individual, sizeof(tspsIndividual_t));
			imigrant.popIndex = -1;
		}
	}
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
	//short int *alreadyBred = NULL;
	//int parent1, parent2;
	int i, j;

	tspsBreedersList_t breederList;

	//calculate the breed chance that each individual has, based on the fitness
	calculateBreedChance(pop);

	breederList.numBreeders = 0;
	breederList.breeders = NULL;

	//select those individuals that will breed
	if(selectBreeders(&breederList, pop, config->maxBreeding) != TSPS_RC_SUCCESS){
		return TSPS_RC_FAILURE;
	}

	//crossover the parents and recalculate the fitness
	for(i=0; i<breederList.numBreeders; i++){

		if(edgeCrossover(&breederList.breeders[i].son1, &breederList.breeders[i].son2, &breederList.breeders[i].parent1, &breederList.breeders[i].parent2) != TSPS_RC_SUCCESS){
			free(breederList.breeders);
			return TSPS_RC_FAILURE;
		}

		if(calculateFitness(&breederList.breeders[i].son1, map) != TSPS_RC_SUCCESS ||
			calculateFitness(&breederList.breeders[i].son2, map) != TSPS_RC_SUCCESS){

			free(breederList.breeders);
			return TSPS_RC_FAILURE;
		}
	}

	//replace the worst individuals with the new sons
	j = pop->numIndividuals-1;
	for(i=0; i<breederList.numBreeders; i++){
		if(j < config->numElitism || j < 1)
			break;

		pop->individual[j].fitness = breederList.breeders[i].son1.fitness;
		pop->individual[j-1].fitness = breederList.breeders[i].son2.fitness;

		memcpy(&pop->individual[j].chromossome, &breederList.breeders[i].son1.chromossome, sizeof(pop->individual[j].chromossome));
		memcpy(&pop->individual[j-1].chromossome, &breederList.breeders[i].son2.chromossome, sizeof(pop->individual[j].chromossome));
		j = j-2;
	}

	free(breederList.breeders);
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

void printIndividual(tspsIndividual_t *ind, unsigned long int numGenerations, int mpiId){
	FILE *file = NULL;
	char *filename;
	int i;

	if(asprintf(&filename, "tsps_%d.log", mpiId) < 0){
		return;
	}

	if((file = fopen(filename, "a")) == NULL){
		free(filename);
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
	free(filename);
}

void calculateBreedChance(tspsPopulation_t *pop){
	int i;
	int bestFitness = pop->individual[0].fitness;
	float rate = 0.0;

	for(i=0; i<pop->numIndividuals; i++){
		rate = (float)pop->individual[i].fitness / (float)bestFitness;
		pop->individual[i].breedChance = 1 / rate;
	}
}

int selectBreeders(tspsBreedersList_t *brl, tspsPopulation_t *pop, int maxBreeding){
	tspsIndividual_t *parent1 = NULL, *parent2 = NULL;
	char *alreadyBreed = NULL;
	int numBreeding = maxBreeding;
	int i, j;

	if((alreadyBreed = (char*)malloc(sizeof(char) * pop->numIndividuals)) == NULL){
		return TSPS_RC_FAILURE;
	}
	memset(alreadyBreed, 0, sizeof(char) * pop->numIndividuals);

	for(i=0; i<pop->numIndividuals-1; i++){

		if(numBreeding == 0)
			break;

		if(alreadyBreed[i] == 1)
			continue;

		if(rand()%100 > (pop->individual[i].breedChance * 100))
			continue;

		alreadyBreed[i] = 1;
		parent1 = &pop->individual[i];

		brl->numBreeders += 1;
		if((brl->breeders = (tspsBreeders_t*)realloc(brl->breeders, sizeof(tspsBreeders_t)*brl->numBreeders)) == NULL){
			free(alreadyBreed);
			return TSPS_RC_FAILURE;
		}
		memcpy(&brl->breeders[brl->numBreeders-1].parent1, parent1, sizeof(tspsIndividual_t));

		for(j=i+1; j<pop->numIndividuals; j++){
			if(alreadyBreed[j] == 1)
				continue;

			if(rand()%100 > (pop->individual[j].breedChance * 100))
				continue;

			parent2 = &pop->individual[j];
			alreadyBreed[j] = 1;
			break;
		}

		if(parent2 == NULL){
			j = i+1;
			while(alreadyBreed[j] != 0){
				if(j+1 < pop->numIndividuals)
					j++;
				else
					j=0;

			}
			parent2 = &pop->individual[j];
			alreadyBreed[j] = 1;
		}

		memcpy(&brl->breeders[brl->numBreeders-1].parent2, parent2, sizeof(tspsIndividual_t));

		memset(&brl->breeders[brl->numBreeders-1].son1, 0, sizeof(tspsIndividual_t));
		memset(&brl->breeders[brl->numBreeders-1].son2, 0, sizeof(tspsIndividual_t));

		numBreeding--;
	}

	free(alreadyBreed);

	return TSPS_RC_SUCCESS;
}

/* choose an individual to migrate between other populations through MPI*/
int migrateIndividual(tspsPopulation_t *population, tspsConfig_t *config, int mpiId, int mpiNumProcs){
	tspsIndividual_t emigrant;
	int i;
	MPI_Status status;

	//choose randomly for now the individual to migrate
	i = rand() % population->numIndividuals;
	memcpy(&emigrant, &population->individual[i], sizeof(tspsIndividual_t));

	//migrate the chosen individual
	MPI_Send(&emigrant, sizeof(tspsIndividual_t), MPI_CHAR, (mpiId+1 < mpiNumProcs? mpiId + 1 : 0), MPI_MIGRATION_TAG, MPI_COMM_WORLD);

	//wait to receive the imigrant from other population that was sent in past generations
	if(imigrant.popIndex != -1){
		MPI_Wait(&imigrant.req, &status);
		memcpy(&population->individual[imigrant.popIndex], &imigrant.individual, sizeof(tspsIndividual_t));
	}

	//receive the new imigrant from the other population
	memset(&imigrant, 0, sizeof(tspsImigrant_t));
	imigrant.popIndex = i;
	MPI_Irecv(&imigrant.individual, sizeof(tspsIndividual_t), MPI_CHAR, (mpiId > 0? mpiId-1 : mpiNumProcs-1),
			MPI_MIGRATION_TAG, MPI_COMM_WORLD, &imigrant.req);

	return TSPS_RC_SUCCESS;
}

void printMigrants(tspsIndividual_t *imigrant, tspsIndividual_t *emigrant, int mpiId){
	FILE *file = NULL;
	char *filename;
	int i;

	if(asprintf(&filename, "tsps_%d.log", mpiId) < 0){
		return;
	}

	if((file = fopen(filename, "a")) == NULL){
		free(filename);
		return;
	}

	fprintf(file, "*** Migration Event!!!\n");
	fprintf(file, "*** Imigrant\n");
	fprintf(file, "\t --> fitness = [%d]\n", imigrant->fitness);
	fprintf(file, "\t --> chromos = ");

	for(i=0; i<NUM_NODES; i++){
		fprintf(file, "[%d]", imigrant->chromossome[i]);
		if((i+1) % 10 == 0)
			fprintf(file, "\n\t\t\t");
	}

	fprintf(file, "\n");
	fprintf(file, "*** Emigrant\n");
	fprintf(file, "\t --> fitness = [%d]\n", emigrant->fitness);
	fprintf(file, "\t --> chromos = ");

	for(i=0; i<NUM_NODES; i++){
		fprintf(file, "[%d]", emigrant->chromossome[i]);
		if((i+1) % 10 == 0)
			fprintf(file, "\n\t\t\t");
	}
	fprintf(file, "\n");
	fprintf(file, "------------------------------------------------------------------\n");

	fclose(file);
	free(filename);
}
