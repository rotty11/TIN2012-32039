/**
 * @file crossover.cpp
 * @author Juan José Escobar Pérez
 * @date 29/06/2015
 * @brief File with the necessary implementation to perform crossover between two individuals
 *
 */

/********************************* Includes *******************************/

#include "crossover.h"
#include <stdlib.h> // rand...
#include <string.h> // memset...

/********************************* Methods ********************************/

/**
 * @brief Perform binary crossover between two individuals
 * @param pop Current population
 * @param populationSize The size of the population
 * @param pool Position of the selected individuals for the crossover
 * @param poolSize Number of selected individuals for the crossover
 * @param nObjectives The number of objectives
 * @param nFeatures The number of features (columns) of the database
 * @param crossDistribution Parameter of crossing distribution
 * @param mutDistribution Parameter of mutation distribution
 * @return The number of generated children
 */
int crossover(individual *pop, const int populationSize, const int *pool, const int poolSize, const unsigned char nObjectives, const int nFeatures, const int crossDistribution, const int mutDistribution) {

	int childrenSize = 0;
	for (int i = 0; i < poolSize; ++i) {

		float probCross = (rand() / (float) RAND_MAX);

		// 90% probability perform crossover. Two childen are generated
		if (probCross < 0.9f) {
			int parent1 = rand() % poolSize;
			int parent2 = rand() % poolSize;

			// Avoid repeated parents
			while (parent1 == parent2) {
				parent2 = rand() % poolSize;
			}

			// Initialize the two children
			pop[populationSize + childrenSize].nSelFeatures = 0;
			pop[populationSize + childrenSize + 1].nSelFeatures = 0;
			pop[populationSize + childrenSize].rank = -1;
			pop[populationSize + childrenSize + 1].rank = -1;
			pop[populationSize + childrenSize].crowding = 0.0f;
			pop[populationSize + childrenSize + 1].crowding = 0.0f;

			for (unsigned char obj = 0; obj < nObjectives; ++obj) {
				pop[populationSize + childrenSize].fitness[obj] = 0.0f;
				pop[populationSize + childrenSize + 1].fitness[obj] = 0.0f;
			}

			// Perform crossover for each decision variable in the chromosome
			// Crossover between two points
			int point1 = rand() % nFeatures;
			int point2 = rand() % nFeatures;
			if (point1 > point2) {
				int temp = point1;
				point1 = point2;
				point2 = temp;
			}

			// First part
			for (int f = 0; f < point1; ++f) {

				// Generate the fth element of first child
				pop[populationSize + childrenSize].nSelFeatures += (pop[populationSize + childrenSize].chromosome[f] = pop[pool[parent1]].chromosome[f]);
				
				// Generate the fth element of second child
				pop[populationSize + childrenSize + 1].nSelFeatures += (pop[populationSize + childrenSize + 1].chromosome[f] = pop[pool[parent2]].chromosome[f]);
			}

			// Second part
			for (int f = point1; f < point2; ++f) {

				// Generate the fth element of first child
				pop[populationSize + childrenSize].nSelFeatures += (pop[populationSize + childrenSize].chromosome[f] = pop[pool[parent2]].chromosome[f]);
				
				// Generate the fth element of second child
				pop[populationSize + childrenSize + 1].nSelFeatures += (pop[populationSize + childrenSize + 1].chromosome[f] = pop[pool[parent1]].chromosome[f]);
			}

			// Third part
			for (int f = point2; f < nFeatures; ++f) {

				// Generate the fth element of first child
				pop[populationSize + childrenSize].nSelFeatures += (pop[populationSize + childrenSize].chromosome[f] = pop[pool[parent1]].chromosome[f]);

				// Generate the fth element of second child
				pop[populationSize + childrenSize + 1].nSelFeatures += (pop[populationSize + childrenSize + 1].chromosome[f] = pop[pool[parent2]].chromosome[f]);
			}

			// At least one decision variable must be set to "1"
			if (pop[populationSize + childrenSize].nSelFeatures == 0) {
				int randomFeature = rand() % nFeatures;
				pop[populationSize + childrenSize].chromosome[randomFeature] = 1;
				pop[populationSize + childrenSize].nSelFeatures = 1;
			}

			if (pop[populationSize + childrenSize + 1].nSelFeatures == 0) {
				int randomFeature = rand() % nFeatures;
				pop[populationSize + childrenSize + 1].chromosome[randomFeature] = 1;
				pop[populationSize + childrenSize + 1].nSelFeatures = 1;
			}

			childrenSize += 2;
		}

		// 10% probability perform mutation. One child is generated
		// Mutation is based on random mutation
		else {
			int parent = rand() % poolSize;

			// Initialize the child
			pop[populationSize + childrenSize].nSelFeatures = 0;
			pop[populationSize + childrenSize].rank = -1;
			pop[populationSize + childrenSize].crowding = 0.0f;

			for (unsigned char obj = 0; obj < nObjectives; ++obj) {
				pop[populationSize + childrenSize].fitness[obj] = 0.0f;
			}

			// Perform mutation on each element of the selected parent
			for (int f = 0; f < nFeatures; ++f) {
				float mut = (rand() / (float) RAND_MAX);
				if (mut < 0.1f) {
					if (pop[pool[parent]].chromosome[f] & 1) {
						pop[populationSize + childrenSize].chromosome[f] = 0;
					}
					else {
						pop[populationSize + childrenSize].chromosome[f] = 1;
						pop[populationSize + childrenSize].nSelFeatures++;
					}
				}
				else {
					pop[populationSize + childrenSize]. nSelFeatures += (pop[populationSize + childrenSize].chromosome[f] = pop[pool[parent]].chromosome[f]);
				}
			}

			// At least one decision variable must be set to "1"
			if (pop[populationSize + childrenSize].nSelFeatures == 0) {
				int randomFeature = rand() % nFeatures;
				pop[populationSize + childrenSize].chromosome[randomFeature] = 1;
				pop[populationSize + childrenSize].nSelFeatures = 1;
			}

			childrenSize++;
		}
	}

	// The not generated children are reinitialized
	int totalIndividuals = populationSize << 1;
	for (int i = populationSize + childrenSize; i < totalIndividuals; ++i) {
		memset(pop[i].chromosome, 0, nFeatures * sizeof(unsigned char));
		memset(pop[i].fitness, 0, nObjectives * sizeof(float));
		pop[i].nSelFeatures = 0;
		pop[i].rank = -1;
		pop[i].crowding = 0.0f;
	}

	return childrenSize;
}