/*
 * map.c
 *
 *  Created on: May 21, 2011
 *      Author: helder
 */

#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "map.h"
#include "string.h"

int parseMap(tspsMap_t *map){

	char *file = NULL;
	char *weights = NULL;
	char *aux = NULL;
	int i, j;

	if(loadMap(&file) != TSPS_RC_SUCCESS)
		return TSPS_RC_FAILURE;

	weights = strstr(file, "EDGE_WEIGHT_SECTION");
	weights += 20;

	for(i=0; i<NUM_NODES; i++){
		map->weights[i][i] = 0;

		for(j=i+1; j<NUM_NODES; j++){

			if(i==0 && j==1)
				aux = strtok(weights, " ");
			else
				aux = strtok(NULL, " ");

			if(aux == NULL)
				break;

			if(*aux == '\n'){
				aux++;
			}
			map->weights[i][j] = atoi(aux);
			map->weights[j][i] = atoi(aux);

		}
	}

	free(file);
	return TSPS_RC_SUCCESS;
}

int loadMap(char **file){

	FILE * pFile;
	long lSize;
	char *buffer;
	size_t result;

	if((pFile = fopen("maps/brazil58.tsp", "rb")) == NULL &&
		(pFile = fopen("../maps/brazil58.tsp", "rb")) == NULL &&
		(pFile = fopen("brazil58.tsp", "rb")) == NULL){

		return TSPS_RC_FAILURE;
	}

	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);

	// allocate memory to contain the whole file:
	buffer = (char*) malloc (sizeof(char)*lSize);
	if (buffer == NULL)
		return TSPS_RC_FAILURE;

	// copy the file into the buffer:
	result = fread (buffer,1,lSize,pFile);
	if (result != lSize)
		return TSPS_RC_FAILURE;

	// terminate
	fclose (pFile);
	*file = buffer;
	return TSPS_RC_SUCCESS;

}


