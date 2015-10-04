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
#include <stdlib.h> // rand...
#include <math.h> // sqrt...

/********************************* Methods ********************************/

/**
 * @brief Evaluation of each individual
 * @param pop Current population
 * @param begin The first individual to evaluate
 * @param end The "end-1" position is the last individual to evaluate
 * @param totalIndividuals The total number of individuals (parents + children)
 * @param maxIndividualsOnGPUKernel The maximum number of individuals to be processed in a single execution of the kernel
 * @param nInstances The number of instances (rows) of the database
 * @param nObjectives The number of objectives
 * @param deviceType The type of device that will run the kernel
 * @param kernel The OpenCL kernel with the implementation of k-means
 * @param command_queue The command queue which contains the tasks (reads/writes on the device...)
 * @param wiGlobal the total number of work-items (threads) that will run the kernel
 * @param wiLocal The number of work-items (threads) per compute unit that will run the kernel
 * @param objDataBase OpenCL object which contains the database
 * @param objSelInstances OpenCL object which contains the instances choosen as initial centroids
 * @param objPopulation OpenCL object which contains the current population
 */
void evaluation(cl_device_id *device, individual *pop, const int begin, const int end, const int totalIndividuals, const int maxIndividualsOnGPUKernel, const int nInstances, const unsigned char nObjectives, const cl_device_type *deviceType, cl_context *context, cl_kernel *kernel, cl_command_queue *command_queue, size_t wiGlobal, size_t wiLocal, const cl_mem *const objDataBase, const cl_mem *const objSelInstances, cl_mem *objPopulation) {


	/************ Kmeans algorithm in OpenCL ***********/

	cl_int status;
	cl_event kernelEvent;

	// Wait until all the memory objects are copied to the device
	clFinish(command_queue[0]);

	// Sets kernel arguments
	if (clSetKernelArg(kernel[0], 0, sizeof(cl_mem), (void *)&objPopulation[0]) != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not set the first kernel argument \n");
		exit(-1);
	}
	if (clSetKernelArg(kernel[0], 1, sizeof(cl_mem), (void *)&objSelInstances[0]) != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not set the second kernel argument \n");
		exit(-1);
	}
	if (clSetKernelArg(kernel[0], 2, sizeof(cl_mem), (void *)&objDataBase[0]) != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not set the third kernel argument \n");
		exit(-1);
	}

	if (deviceType[0] == CL_DEVICE_TYPE_CPU) {
		if (clSetKernelArg(kernel[0], 3, sizeof(int), &begin) != CL_SUCCESS) {
			fprintf(stderr, "Error: Could not set the fourth kernel argument \n");
			exit(-1);
		}
		if (clSetKernelArg(kernel[0], 4, sizeof(int), &end) != CL_SUCCESS) {
			fprintf(stderr, "Error: Could not set the fifth kernel argument \n");
			exit(-1);
		}
		if (clSetKernelArg(kernel[0], 5, KMEANS * nInstances * sizeof(cl_uchar), NULL) != CL_SUCCESS) {
			fprintf(stderr, "Error: Could not set the sixth kernel argument \n");
			exit(-1);
		}
		if (clSetKernelArg(kernel[0], 6, KMEANS * nInstances * sizeof(cl_uchar), NULL) != CL_SUCCESS) {
			fprintf(stderr, "Error: Could not set the seventh kernel argument \n");
			exit(-1);
		}

		// Running the kernel
		if ((status = clEnqueueNDRangeKernel(command_queue[0], kernel[0], 1, NULL, &wiGlobal, &wiLocal, 0, NULL, &kernelEvent)) != CL_SUCCESS) {
			fprintf(stderr, "Error: Could not run the kernel\n");
			exit(-1);
		}
	}
	else {
		int newBegin = begin;
		int newEnd = (begin + maxIndividualsOnGPUKernel >= end) ? end : begin + maxIndividualsOnGPUKernel;
		bool finished = false;
		while (!finished) {
			if (newEnd == end) {
				finished = true;
			}

			// Sets new kernel arguments
			if (clSetKernelArg(kernel[0], 3, sizeof(int), &newBegin) != CL_SUCCESS) {
				fprintf(stderr, "Error: Could not set the fourth kernel argument \n");
				exit(-1);
			}
			if (clSetKernelArg(kernel[0], 4, sizeof(int), &newEnd) != CL_SUCCESS) {
				fprintf(stderr, "Error: Could not set the fifth kernel argument \n");
				exit(-1);
			}

			// Running the kernel
			if ((status = clEnqueueNDRangeKernel(command_queue[0], kernel[0], 1, NULL, &wiGlobal, &wiLocal, 0, NULL, &kernelEvent)) != CL_SUCCESS) {
				fprintf(stderr, "Error: Could not run the kernel\n");
				exit(-1);
			}

			newBegin = newEnd;
			if (newEnd + maxIndividualsOnGPUKernel >= end) {
				newEnd = end;
			}
			else {
				newEnd += maxIndividualsOnGPUKernel;
			}
		}
	}

	// Read the data from the device
	if ((status = clEnqueueReadBuffer(command_queue[0], objPopulation[0], CL_TRUE, 0, totalIndividuals * sizeof(individual), pop, 1, &kernelEvent, NULL)) != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not read the data from the device\n");
		if (status == CL_OUT_OF_RESOURCES) {
			fprintf(stderr, "Warning: It's likely to have to reduce the value of the \"MaxIndividualsOnGPUKernel\" parameter\n");
		}
		exit(-1);
	}


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