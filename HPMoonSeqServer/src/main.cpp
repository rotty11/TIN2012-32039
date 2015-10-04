/**
 * @file main.cpp
 * @author Juan José Escobar Pérez
 * @date 15/06/2015
 * @brief Multiobjective genetic algorithm
 *
 * Multiobjective genetic algorithm running on a general purpose processor
 *
 */

/********************************* Includes ********************************/

#include "tinyxml2.h"
#include "xml.h"
#include "bd.h"
#include "initialization.h"
#include "evaluation.h"
#include "sort.h"
#include "tournament.h"
#include "crossover.h"
#include <stdio.h> // fprintf...
#include <time.h> // clock...

using namespace tinyxml2;


/**
 * @brief Main program
 * @param argc The number of arguments of the program
 * @param argv Arguments of the program
 * @return Returns nothing if successful or a negative number if incorrect
 */
int main(int argc, char** argv) {


	/********** Get the configuration data from the XML file ***********/

	XMLDocument configDoc;
	configDoc.LoadFile(argv[1]);
	const char *dataBaseName = getDataBaseName(&configDoc);
	const int nGenerations = getNGenerations(&configDoc);
	const int maxFeatures = getMaxFeatures(&configDoc);
	const int tourSize = getTourSize(&configDoc);
	const int crossDistribution = getCrossDistribution(&configDoc);
	const int mutDistribution = getMutDistribution(&configDoc);
	const char *dataName = getDataName(&configDoc);
	const char *plotName = getPlotName(&configDoc);
	const char *imageName = getImageName(&configDoc);


	/********** Check program restrictions ***********/

	if (POPULATION_SIZE < 4) {
		fprintf(stderr, "Error: The number of individuals must be 4 or higher\n");
		exit(-1);
	}

	if (N_FEATURES < 2 || N_INSTANCES < 2) {
		fprintf(stderr, "Error: The number of features and number of instances must be 2 or higher\n");
		exit(-1);
	}

	if (N_OBJECTIVES != 2) {
		fprintf(stderr, "Error: The number of objectives must be 2. If you want to increase this number, the module \"evaluation\" should be modified\n");
		exit(-1);
	}

	if (maxFeatures < 1) {
		fprintf(stderr, "Error: The maximum initial number of features must be 1 or higher\n");
		exit(-1);
	}

	if (tourSize < 2) {
		fprintf(stderr, "Error: The number of individuals in the tournament must be 2 or higher\n");
		exit(-1);
	}

	if (crossDistribution < 0 || mutDistribution < 0) {
		fprintf(stderr, "Error: The cross distribution and the mutation distribution must be 0 or higher\n");
		exit(-1);
	}

	clock_t t_ini, t_fin;
	t_ini = clock();


	/********** Get the data base ***********/

	float dataBase[N_INSTANCES * N_FEATURES];
	readDataBase(dataBase, dataBaseName, N_INSTANCES, N_FEATURES);

	// Data base normalization
	normDataBase(dataBase, N_INSTANCES, N_FEATURES);


	/********** Initialize the population and the individuals ***********/

	srand((unsigned int) time(NULL));
	const int totalIndividuals = POPULATION_SIZE << 1;

	// Population will have the parents and children (left half and right half respectively
	// This way is better for the performance
	individual *population = initPopulation(totalIndividuals, N_OBJECTIVES, N_FEATURES, maxFeatures);


	/********** Multiobjective individual evaluation ***********/

	// Get the initial "KMEANS" centroids ***********/
	int selInstances[KMEANS];
	getCentroids(selInstances, N_INSTANCES);

	evaluation(population, 0, POPULATION_SIZE, dataBase, N_INSTANCES, N_FEATURES, N_OBJECTIVES, selInstances);


	/********** Sort the population with the "Non-Domination-Sort" method ***********/

	int nIndFront0 = nonDominationSort(population, POPULATION_SIZE, N_OBJECTIVES, N_INSTANCES, N_FEATURES);


	/********** Get the population quality (calculating the hypervolume) ***********/

	// The reference point will be (X_1 = 1.0, X_2 = 1.0, .... X_N_OBJECTIVES = 1.0)
	double referencePoint[N_OBJECTIVES];
	for (int i = 0; i < N_OBJECTIVES; ++i) {
		referencePoint[i] = 1.0;
	}

	float popHypervolume = getHypervolume(population, nIndFront0, N_OBJECTIVES, referencePoint);
	float auxHypervolume;


	/********** Start the evolution process ***********/

	const int poolSize = POPULATION_SIZE >> 1;
	int pool[poolSize];
	//for (int g = 0; g < nGenerations && popHypervolume > auxHypervolume; ++g) {
	for (int g = 0; g < nGenerations; ++g) {

		// Fill the mating pool
		fillPool(pool, poolSize, tourSize, POPULATION_SIZE);

		// Perform crossover
		int nChildren = crossover(population, POPULATION_SIZE, pool, poolSize, N_OBJECTIVES, N_FEATURES, crossDistribution, mutDistribution);

		// Multiobjective individual evaluation
		int lastChild = POPULATION_SIZE + nChildren;
		evaluation(population, POPULATION_SIZE, lastChild, dataBase, N_INSTANCES, N_FEATURES, N_OBJECTIVES, selInstances);
		
		// The crowding distance of the parents is initialized again for the next nonDominationSort
		for (int i = 0;  i < POPULATION_SIZE; ++i) {
			population[i].crowding = 0.0f;
		}

		// Replace population
		// Parents and children are sorted by rank and crowding distance.
		// The first "populationSize" individuals will advance the next generation
		nIndFront0 = nonDominationSort(population, POPULATION_SIZE + nChildren, N_OBJECTIVES, N_INSTANCES, N_FEATURES);
		
		// Get the population quality (calculating the hypervolume)
		auxHypervolume = getHypervolume(population, nIndFront0, N_OBJECTIVES, referencePoint);
	}

	popHypervolume = auxHypervolume;

	// Finish the time measure
	t_fin = clock();
	double ms = ((double) (t_fin - t_ini) / CLOCKS_PER_SEC) * 1000.0;
	fprintf(stdout, "%.16g\n", ms);
	fprintf(stdout, "%f\n", popHypervolume);

	// Generation of the data file and Gnuplot file for display the Pareto front
	generateGnuplot(dataName, plotName, imageName, population, nIndFront0, N_OBJECTIVES, referencePoint);


	/********** Resources used are released ***********/

	// The individuals (parents and children)
	delete[] population;
}