/*
 * map.h
 *
 *  Created on: May 21, 2011
 *      Author: helder
 */

#ifndef MAP_H_
#define MAP_H_

#include"common.h"

typedef struct{
	unsigned short int weights[NUM_NODES][NUM_NODES];

} tspsMap_t;


int parseMap(tspsMap_t *map);
int loadMap(char **file);

#endif /* MAP_H_ */
