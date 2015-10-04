/**
 * @file tournament.cpp
 * @author Juan José Escobar Pérez
 * @date 28/06/2015
 * @brief File with the necessary implementation for the realization of the tournament (selection of individuals)
 *
 */

/********************************* Includes *******************************/

#include "tournament.h"
#include <stdlib.h> // rand...

/********************************* Methods ********************************/

/**
 * @brief Competition between randomly selected individuals. The best individuals are stored in the pool
 * @param pool Position of the selected individuals for the crossover
 * @param poolSize Number of selected individuals for the crossover
 * @param tourSize Number of individuals competing in the tournament
 * @param populationSize The size of the population
 */
void fillPool(int *pool, const int poolSize, const int tourSize, const int populationSize) {

	// Fill pool
	for (int i = 0; i < poolSize; ++i) {
		int bestCandidate = rand() % populationSize;
		for (int j = 0; j < tourSize - 1; ++j) {
			bool repeated;

			// Avoid repeated candidates
			do {
				int randomCandidate = rand() % populationSize;
				if (randomCandidate != bestCandidate) {
					repeated = false;

					// At this point, the individuals already are sorted by rank and crowding distance
					// Therefore, lower index is better
					if (randomCandidate < bestCandidate) {
						bestCandidate = randomCandidate;
					}
				}
				else {
					repeated = true;
				}
			} while (repeated);
		}

		pool[i] = bestCandidate;
	}
}