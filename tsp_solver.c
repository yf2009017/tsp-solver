/*
 * tsp_solver.c
 *
 *  Created on: May 21, 2011
 *      Author: helder
 */

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include "mpi.h"
#include "common.h"
#include "map.h"
#include "population.h"

extern tspsImigrant_t imigrant;

int main(int argc, char **argv){
	tspsMap_t map;
	tspsConfig_t config;
	tspsPopulation_t population;
	tspsIndividual_t mostFit;
	unsigned long int numGenerations = 0;
	int mpiNumProcs = 0, mpiId = 0;

	//starting MPI directives
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD,&mpiNumProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpiId);

	if(readConfig(&config, argc, argv) != TSPS_RC_SUCCESS){
		return TSPS_RC_FAILURE;
	}

	// initialize random seed:
	srand ( time(NULL)*mpiId );

	if(parseMap(&map) != TSPS_RC_SUCCESS){
		printf("Error! Unable to read map 'maps/brazil58.tsp'!\n");
		return TSPS_RC_FAILURE;
	}

	if(generateRandomPopulation(&population, &map, &config) != TSPS_RC_SUCCESS){
		printf("Error! Unable to generate random population!\n");
		return TSPS_RC_FAILURE;
	}

	mostFit.fitness = 999999;

	memset(&imigrant, 0, sizeof(tspsImigrant_t));
	imigrant.req = MPI_REQUEST_NULL;
	imigrant.popIndex = -1;
	while(1){

		numGenerations++;

		if(sortPopulation(&population, mpiId) != TSPS_RC_SUCCESS){
			return TSPS_RC_FAILURE;
		}

		if(population.individual[0].fitness < mostFit.fitness){
			memcpy(&mostFit, &population.individual[0], sizeof(mostFit));
			printIndividual(&mostFit, numGenerations, mpiId);
		}

		if(crossoverPopulation(&population, &map, &config) != TSPS_RC_SUCCESS){
			return TSPS_RC_FAILURE;
		}

		if(mutatePopulation(&population, &map, &config) != TSPS_RC_SUCCESS){
			return TSPS_RC_FAILURE;
		}

		//end execution
		if(config.numGenerations > 0 && numGenerations == config.numGenerations){
			break;
		}

		//migrate
		if(mpiNumProcs > 1 && numGenerations % 100 == 0){
			if(migrateIndividual(&population, &config, mpiId, mpiNumProcs) != TSPS_RC_SUCCESS){
				return TSPS_RC_FAILURE;
			}
		}

		if(numGenerations % 1000 == 0){
			printf("* (Id%d) Generation [%lu]...\n", mpiId, numGenerations);
		}

	}

	printf("* (Id%d) Generation [%lu]...\n", mpiId, numGenerations);
	printf("* (Id%d) Max number of generations reached! Ending... \n", mpiId);

	free(population.individual);
	MPI_Finalize();
	return TSPS_RC_SUCCESS;

}
