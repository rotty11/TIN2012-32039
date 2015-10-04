/**
 * @file initialization.cpp
 * @author Juan José Escobar Pérez
 * @date 18/06/2015
 * @brief File with the necessary implementation to initialize the population
 *
 */

/********************************* Includes *******************************/

#include "initialization.h"
#include <stdlib.h> // rand...
#include <string.h> // memset...

/********************************* Methods ********************************/

/**
 * @brief Allocates memory for parents and children. Also, they are initialized
 * @param totalIndividuals The total number of individuals (parents + children)
 * @param nObjectives The number of objectives
 * @param nFeatures The number of features (columns) of the database
 * @param maxFeatures The maximum number of features that will be set to "1"
 * @return Pointer to the population
 */
individual* initPopulation(const int totalIndividuals, const unsigned char nObjectives, const int nFeatures, const int maxFeatures) {


	/********** Initialization of the population and the individuals ***********/

	individual *pop = new individual[totalIndividuals];

	// Allocates memory for parents and children
	for (int i = 0; i < totalIndividuals; ++i) {
		memset(pop[i].chromosome, 0, nFeatures * sizeof(unsigned char));
		memset(pop[i].fitness, 0, nObjectives * sizeof(float));
		pop[i].nSelFeatures = 0;
		pop[i].rank = -1;
		pop[i].crowding = 0.0f;
	}

	// Only the parents are initialized
	int sizePopulation = totalIndividuals >> 1;
	for (int i = 0; i < sizePopulation; ++i) {

		// Set the "1" value at most "maxFeatures" decision variables
		for (int j = 0; j < maxFeatures; ++j) {
			int randomFeature = rand() % nFeatures;
			if (!(pop[i].chromosome[randomFeature] & 1)) {
				pop[i].chromosome[randomFeature] = 1;
				pop[i].nSelFeatures++;
			}
		}
	}

	return pop;
}