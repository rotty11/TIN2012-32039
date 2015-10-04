/**
 * @file sort.cpp
 * @author Juan José Escobar Pérez
 * @date 26/06/2015
 * @brief File with the necessary implementation to perform "nonDominationSort" according to the "Pareto front" and the crowding distance. The original code in Matlab is owned by Aravind Seshadri
 *
 */

using namespace std;

/********************************* Includes *******************************/

#include "sort.h"
#include <vector> // std::Vector...
#include <algorithm> // sort...
#include <math.h> // INFINITY...

/********************************* Methods ********************************/

/**
 * @brief Perform "nonDominationSort" on the population
 * @param pop Current population
 * @param nIndividuals The number of individuals which will be sorted
 * @param nObjectives The number of objectives
 * @param nInstances The number of instances (rows) of the database
 * @param nFeatures The number of features (columns) of the database
 * @return The number of individuals in the front 0
 */
int nonDominationSort(individual *pop, const int nIndividuals, const unsigned char nObjectives, const int nInstances, const int nFeatures) {

	// Individuals to classify (really they aren't "individual" objects)
	// Each "individual" p of "indivDomination" contains the number of individuals...
	// ...who dominate p and a list of individuals who are dominated by p
	vector< pair< vector<int>, int> > indivDomination(nIndividuals);
	vector< vector<int> > front(1);

	// Search for individuals who belong to the first front
	int nFronts = 0;
	for (int i = 0; i < nIndividuals; ++i) {
		for (int j = i + 1; j < nIndividuals; ++j) {
			unsigned char domLess = 0;
			unsigned char domEqual = 0;
			unsigned char domMore = 0;
			for (unsigned char obj = 0; obj < nObjectives; ++obj) {
				if (pop[i].fitness[obj] < pop[j].fitness[obj]) {
					domLess++;
				}
				else if (pop[i].fitness[obj] == pop[j].fitness[obj]) {
					domEqual++;
				}
				else {
					domMore++;
				}
			}

			if (domLess == 0 && domEqual != nObjectives) {
				indivDomination[i].second++;
				indivDomination[j].first.push_back(i);
			}
			else if (domMore == 0 && domEqual != nObjectives) {
				indivDomination[i].first.push_back(j);
				indivDomination[j].second++;
			}
		}

		if (indivDomination[i].second == 0) {
			pop[i].rank = 0;
			front[nFronts].push_back(i);
		}
	}

	// Find the subsequent fronts
	while (!front[nFronts].empty()) {
		front.push_back(vector<int>());
		int sizeActualFront = front[nFronts].size();
		for (int i = 0; i < sizeActualFront; ++i) {
			int nDomByInd = (int) indivDomination[front[nFronts][i]].first.size();
			for (int j = 0; j < nDomByInd; ++j) {
				int dominateToInd = indivDomination[front[nFronts][i]].first[j];
				int nDomToInd = (--indivDomination[dominateToInd].second);
				if (nDomToInd == 0) {
					pop[dominateToInd].rank = nFronts + 1;
						front[nFronts + 1].push_back(dominateToInd);
				}
			}
		}
		nFronts++;
	}

	front.pop_back();

	// Find the crowding distance for each individual in each front
	sort(pop, pop + nIndividuals, rankCompare());

	for (int f = 0, i = 0; f < nFronts; ++f) {
		int sizeFrontI = (int) front[f].size();
		individual *begin = pop + i;
		individual *end = begin + sizeFrontI;
		for (unsigned char obj = 0; obj < nObjectives; ++obj) {
			sort(begin, end, objectiveCompare(obj));
			float fMin = begin -> fitness[obj];
			float fMax = (end - 1) -> fitness[obj];
			begin -> crowding = INFINITY;
			(end - 1) -> crowding = INFINITY;
			bool fMaxFminZero = (fMax - fMin == 0.0f);

			for (int j = 1; j < sizeFrontI - 1; ++j) {
				individual *current = begin + j;
				if (fMaxFminZero) {
					current -> crowding = INFINITY;
				}
				else if (current -> crowding != INFINITY) {
					float nextObj = (current + 1) -> fitness[obj];
					float previousObj = (current - 1) -> fitness[obj];
					current -> crowding += (nextObj - previousObj) / (fMax - fMin);
				}
			}
		}

		i += sizeFrontI;
	}

	sort(pop, pop + nIndividuals, rankAndCrowdingCompare());
	return front[0].size();
}