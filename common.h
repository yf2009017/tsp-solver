/*
 * common.h
 *
 *  Created on: May 21, 2011
 *      Author: helder
 */

#ifndef COMMON_H_
#define COMMON_H_

#define TSPS_RC_SUCCESS 0
#define TSPS_RC_FAILURE 1

#define NUM_NODES 58

typedef struct{
	int populationSize;
	double mutationRate;
	int numGenerations;
	int numElitism;
	int mutationSize;
	int maxBreeding;
} tspsConfig_t;

int readConfig(tspsConfig_t *config, int argc, char **argv);

#endif /* COMMON_H_ */
