/*
 * population.h
 *
 *  Created on: May 21, 2011
 *      Author: helder
 */

#ifndef POPULATION_H_
#define POPULATION_H_

#include "common.h"
#include "map.h"
#include "mpi.h"

#define MPI_MIGRATION_TAG 1

/*--------------------------------------*/

typedef struct{
	int fitness;
	float breedChance;
	short unsigned int chromossome[NUM_NODES];
} tspsIndividual_t;

typedef struct{
	int numIndividuals;
	tspsIndividual_t *individual;
} tspsPopulation_t;

typedef struct{
	int popIndex;
	MPI_Request req;
	tspsIndividual_t individual;
} tspsImigrant_t;
/*--------------------------------------*/

typedef struct{
	tspsIndividual_t parent1;
	tspsIndividual_t parent2;
	tspsIndividual_t son1;
	tspsIndividual_t son2;
} tspsBreeders_t;

typedef struct{
	int numBreeders;
	tspsBreeders_t *breeders;
} tspsBreedersList_t;
/*--------------------------------------*/

typedef struct{
	int edges[4];
} tspsNodeEdges_t;

typedef struct{
	tspsNodeEdges_t node[NUM_NODES];
} tspsEdgeTable_t;
/*--------------------------------------*/

int generateRandomPopulation(tspsPopulation_t *population, tspsMap_t *map, tspsConfig_t *config);
int generateRandomChromossome(tspsIndividual_t *ind);
int calculateFitness(tspsIndividual_t *ind, tspsMap_t *map);
int sortPopulation(tspsPopulation_t *pop, int mpiId);
int comparePopulation(const void * a, const void * b);
int crossoverPopulation(tspsPopulation_t *pop, tspsMap_t *map, tspsConfig_t *config);
int edgeCrossover(tspsIndividual_t *son1, tspsIndividual_t *son2, tspsIndividual_t *par1, tspsIndividual_t *par2);
int chooseNodeFromEdges(tspsEdgeTable_t *edgeTable, tspsNodeEdges_t *nodeEdges, int *isAlreadyUsed);
int mutatePopulation(tspsPopulation_t *pop, tspsMap_t *map, tspsConfig_t *config);
void printIndividual(tspsIndividual_t *ind, unsigned long int numGenerations, int mpiId);
void calculateBreedChance(tspsPopulation_t *pop);
int selectBreeders(tspsBreedersList_t *brl, tspsPopulation_t *pop, int maxBreeding);
int migrateIndividual(tspsPopulation_t *population, tspsConfig_t *config, int mpiId, int mpiNumProcs);
void printMigrants(tspsIndividual_t *imigrant, tspsIndividual_t *emigrant, int mpiId);

#endif /* POPULATION_H_ */
