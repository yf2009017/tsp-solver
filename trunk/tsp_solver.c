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
#include "common.h"
#include "map.h"
#include "population.h"

int main(int argc, char **argv){
	tspsMap_t map;
	tspsConfig_t config;
	tspsPopulation_t population;
	tspsIndividual_t mostFit;
	unsigned long int numGenerations = 0;

	if(readConfig(&config, argc, argv) != TSPS_RC_SUCCESS){
		return TSPS_RC_FAILURE;
	}

	// initialize random seed:
	srand ( time(NULL) );

	if(parseMap(&map) != TSPS_RC_SUCCESS){
		printf("Error! Unable to read map 'maps/brazil58.tsp\n'!");
		return TSPS_RC_FAILURE;
	}

	if(generateRandomPopulation(&population, &map, &config) != TSPS_RC_SUCCESS){
		printf("Error! Unable to generate random population!\n");
		return TSPS_RC_FAILURE;
	}

	mostFit.fitness = 999999;

	while(1){

		numGenerations++;

		if(sortPopulation(&population) != TSPS_RC_SUCCESS){
			return TSPS_RC_FAILURE;
		}

		if(population.individual[0].fitness < mostFit.fitness){
			memcpy(&mostFit, &population.individual[0], sizeof(mostFit));
			printIndividual(&mostFit, numGenerations);
		}

		if(crossoverPopulation(&population, &map, &config) != TSPS_RC_SUCCESS){
			return TSPS_RC_FAILURE;
		}

		if(mutatePopulation(&population, &map, &config) != TSPS_RC_SUCCESS){
			return TSPS_RC_FAILURE;
		}

		if(config.numGenerations > 0){
			if(numGenerations == config.numGenerations){
				break;
			}
		}

		if(numGenerations % 1000 == 0){
			printf("* Generation [%lu]...\n", numGenerations);
		}

	}

	return TSPS_RC_SUCCESS;

}
