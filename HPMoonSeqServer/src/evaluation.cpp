/**
 * @file evaluation.cpp
 * @author Juan José Escobar Pérez
 * @date 20/06/2015
 * @brief File with the necessary implementation for the evaluation of the individuals
 *
 */

/********************************** Includes **********************************/

#include "evaluation.h"
#include "hv.h"
#include <stdio.h> // fprintf...
#include <stdlib.h> // malloc, rand...
#include <string.h> // memset...
#include <math.h> // sqrt, INFINITY...

/********************************* Methods ********************************/

/**
 * @brief K-means algorithm which minimize the within-cluster and maximize Inter-cluster sum of squares (WCSS and ICSS)
 * @param pop Current population
 * @param begin The first individual to evaluate
 * @param end The "end-1" position is the last individual to evaluate
 * @param selInstances The instances choosen as initial centroids
 * @param dataBase The database which will contain the instances and the features
 */
void kmeans(individual *pop, const int begin, const int end, const int *const selInstances, const float *const dataBase) {

	bool *mapping = (bool*) malloc(KMEANS * N_INSTANCES * sizeof(bool));
	bool *newMapping = (bool*) malloc(KMEANS * N_INSTANCES * sizeof(bool));

	// Evaluate all the individuals
	for (int ind = begin; ind < end; ++ind) {
		const int totalCoord = KMEANS * N_FEATURES;
		float centroids[KMEANS * N_FEATURES];

		// The centroids will have the selected features of the individual
		for (int k = 0; k < KMEANS; ++k) {
			int posDataBase = selInstances[k] * N_FEATURES;
			int posCentr = k * N_FEATURES;

			for (int f = 0; f < N_FEATURES; ++f) {
				if (pop[ind].chromosome[f] & 1) {
					centroids[posCentr + f] = dataBase[posDataBase + f];
				}
			}
		}



		/******************** Convergence process *********************/

		// Initialize the array of minimum distances and the mapping table
		const int totalDist = KMEANS * N_INSTANCES;
		float distCentroids[KMEANS * N_INSTANCES];
		int samples_in_k[KMEANS];

		// Initialize the mapping table
		for (int i = 0; i < totalDist; ++i) {
			mapping[i] = false;
		}

		// To avoid poor performance, at most "MAX_ITER_KMEANS" iterations are executed
		bool converged = false;
		for (int maxIter = 0; maxIter < MAX_ITER_KMEANS && !converged; ++maxIter) {

			// The mapping table is cleaned in each iteration
			for (int i = 0; i < totalDist; ++i) {
				newMapping[i] = false;
			}
			for (int i = 0; i < KMEANS; ++i) {
				samples_in_k[i] = 0;
			}

			// Calculate all distances (Euclidean distance) between each instance and the centroids
			for (int i = 0; i < N_INSTANCES; ++i) {
				float minDist = INFINITY;
				int selectCentroid = -1;
				int pos = N_FEATURES * i;
				for (int k = 0; k < KMEANS; ++k) {
					float sum = 0.0f;
					int posCentr = k * N_FEATURES;
					int posDistCentr = k * N_INSTANCES;
					for (int f = 0; f < N_FEATURES; ++f) {
						if (pop[ind].chromosome[f] & 1) {
							sum += (dataBase[pos + f] - centroids[posCentr + f]) * (dataBase[pos + f] - centroids[posCentr + f]);
						}
					}

					float euclidean = sqrt(sum);
					distCentroids[posDistCentr + i] = euclidean;
					if (euclidean < minDist) {
						minDist = euclidean;
						selectCentroid = k;
					}
				}

				newMapping[(selectCentroid * N_INSTANCES) + i] = true;
				samples_in_k[selectCentroid]++;
			}

			// Has the algorithm converged?
			converged = true;
			for (int k = 0; k < KMEANS && converged; ++k) { 
				int posMap = k * N_INSTANCES;
				for (int i = 0; i < N_INSTANCES && converged; ++i) {
					if (newMapping[posMap + i] != mapping[posMap + i]) {
						converged = false;
					}
				}
			}

			if (!converged) {

				// Update the position of the centroids
				for (int k = 0; k < KMEANS; ++k) {
					int posCentr = k * N_FEATURES;
					int posMap = k * N_INSTANCES;
					for (int f = 0; f < N_FEATURES; ++f) {
						float sum = 0.0f;
						if (pop[ind].chromosome[f] & 1) {
							for (int i = 0; i < N_INSTANCES; ++i) {
								if (newMapping[posMap + i]) {
									sum += dataBase[(N_FEATURES * i) + f];
								}
							}
							centroids[posCentr + f] = sum / samples_in_k[k];
						}
					}
				}

				// Swap mapping tables
				bool *aux = newMapping;
				newMapping = mapping;
				mapping = aux;
			}
	}


		/************ Minimize the within-cluster and maximize Inter-cluster sum of squares (WCSS and ICSS) *************/

		float sumWithin = 0.0f;
		float sumInter = 0.0f;
		for (int k = 0; k < KMEANS; ++k) {
			int posCentr = k * N_FEATURES;
			int posDistCentr = k * N_INSTANCES;

			// Within-cluster
			for (int i = 0; i < N_INSTANCES; ++i) {
				if (mapping[posDistCentr + i]) {
					sumWithin += distCentroids[posDistCentr + i];
				}
			}

			// Inter-cluster
			for (int i = posCentr + N_FEATURES; i < totalCoord; i += N_FEATURES) {
				float sum = 0.0f;
				for (int f = 0; f < N_FEATURES; ++f) {
					if (pop[ind].chromosome[f] & 1) {
						sum += (centroids[posCentr + f] - centroids[i + f]) * (centroids[posCentr + f] - centroids[i + f]);
					}
				}
				sumInter += sqrt(sum);
			}
		}

		// First objective function (Within-cluster sum of squares (WCSS))
		pop[ind].fitness[0] = sumWithin;

		// Second objective function (Inter-cluster sum of squares (ICSS))
		pop[ind].fitness[1] = sumInter;

		// Third objective function (Number of selected features)
		//pop[ind].fitness[2] = (float) nSelFeatures;
	}

	// Resources used are released
	free(mapping);
	free(newMapping);
}


/**
 * @brief Evaluation of each individual
 * @param pop Current population
 * @param begin The first individual to evaluate
 * @param end The "end-1" position is the last individual to evaluate
 * @param dataBase The database which will contain the instances and the features
 * @param nInstances The number of instances (rows) of the database
 * @param nFeatures The number of features (columns) of the database
 * @param nObjectives The number of objectives
 * @param selInstances The instances choosen as initial centroids
 */
void evaluation(individual *pop, const int begin, const int end, const float *const dataBase, const int nInstances, const int nFeatures, const unsigned char nObjectives, const int *const selInstances) {


	/************ Kmeans algorithm ***********/

	// Evaluate all the individuals and get the first and second objective for them
	kmeans(pop, begin, end, selInstances, dataBase);


	/******************** Fitness normalization *********************/

	int totalInd = end - begin;
	for (unsigned char obj = 0; obj < nObjectives; ++obj) {

		// Fitness vector average
		float average = 0;
		for (int i = begin; i < end; ++i) {
			average += pop[i].fitness[obj];
		}

		average /= totalInd;

		// Fitness vector variance
		float variance = 0;
		for (int i = begin; i < end; ++i) {
			variance += (pop[i].fitness[obj] - average) * (pop[i].fitness[obj] - average);
		}
		variance /= (totalInd - 1);

		// Fitness vector standard deviation
		float std_deviation = sqrt(variance);

		// The second objective is a maximization problem. x_new must be negative
		if (obj == 1) {

			// Normalize a set of continuous values using SoftMax (based on the logistic function)
			for (int i = begin; i < end; ++i) {
				float x_scaled = (pop[i].fitness[obj] - average) / std_deviation;
				float x_new = 1.0f / (1.0f + exp(-x_scaled));
				pop[i].fitness[obj] = -x_new;
			}
		}
		else {

			// Normalize a set of continuous values using SoftMax (based on the logistic function)
			for (int i = begin; i < end; ++i) {
				float x_scaled = (pop[i].fitness[obj] - average) / std_deviation;
				float x_new = 1.0f / (1.0f + exp(-x_scaled));
				pop[i].fitness[obj] = x_new;
			}
		}
	}
}


/**
 * @brief Gets the hypervolume measure of the population
 * @param pop Current population
 * @param nIndFront0 The number of individuals in the front 0
 * @param nObjectives The number of objectives
 * @param referencePoint The necessary reference point for calculation
 * @return The value of the hypervolume
 */
float getHypervolume(const individual *const pop, const int nIndFront0, const unsigned char nObjectives, const double *const referencePoint) {

	// Generation the points for the calculation of the hypervolume
	double *points = new double[nObjectives * nIndFront0];
	for (int i = 0; i < nIndFront0; ++i) {
		for (unsigned char obj = 0; obj < nObjectives; ++obj) {
			points[(i * nObjectives) + obj] = pop[i].fitness[obj];
		}
	}

	float hypervolume = fpli_hv(points, nObjectives, nIndFront0, referencePoint);
	delete[] points;

	return hypervolume;
}


/**
 * @brief Gets the initial centroids (instances choosen randomly)
 * @param selInstances Where the instances choosen as initial centroids will be stored
 * @param nInstances The number of instances (rows) of the database
 */
void getCentroids(int *selInstances, const int nInstances) {

	// The init centroids will be instances choosen randomly (Forgy's Method)
	for (int k = 0; k < KMEANS; ++k) {
		bool exists = false;
		int randomInstance;

		// Avoid repeat centroids
		do {
			randomInstance = rand() % nInstances;
			exists = false;

			// Look if the generated index already exists
			for (int kk = 0; kk < k && !exists; ++kk) {
				exists = (randomInstance == selInstances[kk]);
			}
		} while (exists);

		selInstances[k] = randomInstance;
	}
}


/**
 * @brief Generates gnuplot code for data display
 * @param dataName The name of the file which will contain the fitness of the individuals in the first Pareto front
 * @param plotName The name of the file which will contain the gnuplot code for data display
 * @param imageName The name of the file which will contain the image with the data (graphic)
 * @param pop Current population
 * @param nIndFront0 The number of individuals in the front 0
 * @param nObjectives The number of objectives
 * @param referencePoint The reference point used for the hypervolume calculation
 */
void generateGnuplot(const char *dataName, const char *plotName, const char *imageName, const individual *const pop, const int nIndFront0, const unsigned char nObjectives, const double *const referencePoint) {

	// Open the data file
	FILE *f_data;
	f_data = fopen(dataName, "w");
	if (!f_data) {
		fprintf(stderr, "Error: An error ocurred opening or writting the data file\n");
		exit(-1);
	}

	// Write the data
	fprintf(f_data, "#Objective0");
	for (unsigned char obj = 1; obj < nObjectives; ++obj) {
		fprintf(f_data, "\tObjective%d", obj);
	}
	for (int i = 0; i < nIndFront0; ++i) {
		fprintf(f_data, "\n%f", pop[i].fitness[0]);
		for (unsigned char obj = 1; obj < nObjectives; ++obj) {
			fprintf(f_data, "\t%f", pop[i].fitness[obj]);
		}
	}
	fclose(f_data);

	// Gnuplot is only available for two objectives
	if (nObjectives == 2) {

		// Open the gnuplot script file
		FILE *f_plot;
		f_plot = fopen(plotName, "w");
		if (!f_data) {
			fprintf(stderr, "Error: An error ocurred opening or writting the plot file\n");
			exit(-1);
		}

		// Write the code
		fprintf(f_plot, "#!/usr/bin/gnuplot\n");
		fprintf(f_plot, "set terminal png size 1024,600\n");
		fprintf(f_plot, "set output '%s.png'\n", imageName);
		fprintf(f_plot, "set multiplot\n");
		fprintf(f_plot, "set xlabel \"Objective 0\"\n");
		fprintf(f_plot, "set grid\n");
		fprintf(f_plot, "set title \"Pareto front\"\n");
		fprintf(f_plot, "set ylabel \"Objective 1\"\n");
		fprintf(f_plot, "set size 0.9,0.9\n");
		fprintf(f_plot, "set origin 0.00,0.05\n");
		fprintf(f_plot, "set key center top\n");
		fprintf(f_plot, "plot [0:1][-1:1] '< sort %s' using 1:%d title \"Front 0\" with lp,\\\n", dataName, nObjectives);
		fprintf(f_plot, "\t\"<echo '%f %f'\" title \"Reference point\" with points,\\\n", referencePoint[0], referencePoint[1]);
		fprintf(f_plot, "\t0 title \"Top pareto limit\" with lp;\n");
		fprintf(f_plot, "set nomultiplot\n");
		fprintf(f_plot, "reset\n");
		fclose(f_plot);
	}
	else {
		fprintf(stdout, "Gnuplot is only available for two objectives. Not generated gnuplot file\n");
	}
}