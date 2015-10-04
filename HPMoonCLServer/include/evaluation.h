/**
 * @file evaluation.h
 * @author Juan José Escobar Pérez
 * @date 20/06/2015
 * @brief Header file for the evaluation of the individuals
 *
 */

#ifndef EVALUATION_H
#define EVALUATION_H

/********************************* Includes *******************************/

#include "individual.h" // Individual
#include <CL/cl.h> // OpenCL

/********************************* Defines ********************************/

/**
 * @brief Number of clusters (and centroids) used for the kmeans algorithm
 */
#define KMEANS 3


/**
 * @brief Max iterations for the convergence of centroids
 */
#define MAX_ITER_KMEANS 100

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
void evaluation(cl_device_id *device, individual *pop, const int begin, const int end, const int totalIndividuals, const int maxIndividualsOnGPUKernel, const int nInstances, const unsigned char nObjectives, const cl_device_type *deviceType, cl_context *context, cl_kernel *kernel, cl_command_queue *command_queue, size_t wiGlobal, size_t wiLocal, const cl_mem *const objDataBase, const cl_mem *const objSelInstances, cl_mem *objPopulation);


/**
 * @brief Gets the hypervolume measure of the population
 * @param pop Current population
 * @param nIndFront0 The number of individuals in the front 0
 * @param nObjectives The number of objectives
 * @param referencePoint The necessary reference point for calculation
 * @return The value of the hypervolume
 */
float getHypervolume(const individual *const pop, const int nIndFront0, const unsigned char nObjectives, const double *referencePoint);


/**
 * @brief Gets the initial centroids (instances choosen randomly)
 * @param selInstances Where the instances choosen as initial centroids will be stored
 * @param nInstances The number of instances (rows) of the database
 */
void getCentroids(int *selInstances, const int nInstances);


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
void generateGnuplot(const char *dataName, const char *plotName, const char *imageName, const individual *const pop, const int nIndFront0, const unsigned char nObjectives, const double *const referencePoint);

#endif
