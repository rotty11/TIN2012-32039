/**
 * @file bd.cpp
 * @author Juan José Escobar Pérez
 * @date 15/06/2015
 * @brief File with the necessary implementation to read and store the database
 *
 */

/********************************* Includes *******************************/

#include "bd.h"
#include <stdio.h> // FILE...
#include <stdlib.h> // exit...
#include <math.h> // sqrt...

/********************************* Methods ********************************/

/**
 * @brief Reading the database
 * @param dataBase The database which will contain the instances
 * @param fileName The file name of the database
 * @param nInstances The number of instances (rows) of the database
 * @param nFeatures The number of features (columns) of the database
 */
void readDataBase(float *dataBase, const char *fileName, const int nInstances, const int nFeatures) {


	/********** Open the database ***********/

	FILE *fData = fopen(fileName, "r");
	if (!fData) {
		fprintf(stderr, "An error ocurred opening the data base file\n");
		exit(-1);
	}


	/********** Reading and database storage ***********/

	for(int i = 0; i < nInstances; ++i) {
		for(int j = 0; j < nFeatures; ++j)  {
			fscanf(fData, "%f", &dataBase[(nFeatures * i) + j]);
		}
	}

	// Close the data base and return it
	fclose(fData);
}


/**
 * @brief The database is normalized between 0 and 1
 * @param dataBase Database to be normalized
 * @param nInstances The number of instances (rows) of the database
 * @param nFeatures The number of features (columns) of the database
 */
void normDataBase(float *dataBase, const int nInstances, const int nFeatures) {


	/********** Database normalization ***********/

	for(int j = 0; j < nFeatures; ++j) {

		// Average of the features vector
		float average = 0;
		for(int i = 0; i < nInstances; ++i) {
			int pos = (nFeatures * i) + j;
			average += dataBase[pos];
		}

		average /= nInstances;

		// Variance of the features vector
		float variance = 0;
		for(int i = 0; i < nInstances; ++i) {
			int pos = (nFeatures * i) + j;
			variance += (dataBase[pos] - average) * (dataBase[pos] - average);
		}
		variance /= (nInstances - 1);

		// Standard deviation of the features vector
		float std_deviation = sqrt(variance);

		// Normalize a set of continuous values using SoftMax (based on the logistic function)
		for(int i = 0; i < nInstances; ++i) {
			int pos = (nFeatures * i) + j;
			float x_scaled = (dataBase[pos] - average) / std_deviation;
			float x_new = 1.0f / (1.0f + exp(-x_scaled));
			dataBase[pos] = x_new;
		}
	}
}