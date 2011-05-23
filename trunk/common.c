/*
 * common.c
 *
 *  Created on: May 21, 2011
 *      Author: helder
 */

#include "common.h"
#include "stdio.h"
#include <stdlib.h>

int readConfig(tspsConfig_t *config, int argc, char **argv){
	int populationSize = 0;
	double mutationRate = 0.0;
	int numGenerations = 0;
	int numElitism = 0;
	int mutationSize = 0;

	if(argc != 6){
		printf("Invalid number of args\n Usage: ./tsp_solver <population size> <mutation rate> <number generations> <elitism number> <mutation size>\n");
		return TSPS_RC_FAILURE;
	}

	populationSize = atoi(argv[1]);
	if(populationSize <= 0){
		printf("Invalid population size (> 0)!\n Usage: ./tsp_solver <population size> <mutation rate><number generations> <elitism number> <mutation size>\n");
		return TSPS_RC_FAILURE;
	}

	mutationRate = atof(argv[2]);
	if(mutationRate < 0 || mutationRate > 1){
		printf("Invalid mutation rate (0 - 1)!\n Usage: ./tsp_solver <population size> <mutation rate> <number generations> <elitism number> <mutation size>\n");
		return TSPS_RC_FAILURE;
	}

	numGenerations = atoi(argv[3]);
	if(numGenerations < 0){
		printf("Invalid number of generations (>= 0)!\n Usage: ./tsp_solver <population size> <mutation rate> <number generations> <elitism number> <mutation size>\n");
		return TSPS_RC_FAILURE;
	}

	numElitism = atoi(argv[4]);
	if(numElitism < 0 || numElitism >= populationSize){
		printf("Invalid elitism number (population size < elitism <= 0)!\n Usage: ./tsp_solver <population size> <mutation rate> <number generations> <elitism number> <mutation size>\n");
		return TSPS_RC_FAILURE;
	}

	mutationSize = atoi(argv[5]);
	if(mutationSize < 1 || mutationSize > NUM_NODES/2){
		printf("Invalid mutation size (%d < mutationSize < 0)!\n Usage: ./tsp_solver <population size> <mutation rate> <number generations> <elitism number> <mutation size>\n", NUM_NODES/2);
		return TSPS_RC_FAILURE;
	}

	config->populationSize = populationSize;
	config->mutationRate = mutationRate;
	config->numGenerations = numGenerations;
	config->numElitism = numElitism;
	config->mutationSize = mutationSize;
	return TSPS_RC_SUCCESS;

}
