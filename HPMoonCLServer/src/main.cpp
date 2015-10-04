/**
 * @file main.cpp
 * @author Juan José Escobar Pérez
 * @date 08/07/2015
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
#include <omp.h> // OpenMP
#include <stdio.h> // fprintf...
#include <stdlib.h> // malloc, atoi...
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
	const char *kernelName = getKernelName(&configDoc);
	const char *dataName = getDataName(&configDoc);
	const char *plotName = getPlotName(&configDoc);
	const char *imageName = getImageName(&configDoc);
	const cl_device_type deviceType = getDeviceType(&configDoc);
	const char *platformVendor = getPlatformVendor(&configDoc);
	const char *deviceName = getDeviceName(&configDoc);
	const size_t wiGlobal = getWiGlobal(&configDoc);
	size_t wiLocal = getWiLocal(&configDoc);
	const int maxIndividualsOnGPUKernel = getMaxIndividualsOnGPUKernel(&configDoc);

	// OpenCL variables
	cl_platform_id* platforms;
	cl_uint numPlatforms, numDevicesPlatform;
	cl_device_id *devices, device = NULL;
	cl_context context;
	cl_program program;
	cl_command_queue command_queue;
	cl_kernel kernel;
	cl_mem objDataBase, objSelInstances, objPopulation;
	cl_int status;
	cl_float *data;
	cl_int *map;
	cl_ulong kernelStartTime, kernelEndTime;
	FILE *kernels;

	// Time Measure
	double timeStart, timeEnd;


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

	if (maxIndividualsOnGPUKernel < 1) {
		fprintf(stderr, "Error: The maximum of individuals on the GPU kernel must be 1 or higher\n");
		exit(-1);
	}


	/********** OpenCL init ***********/

	// Get the number of platforms
	if (clGetPlatformIDs(0, NULL, &numPlatforms) != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not get the number of platforms\n");
		exit(-1);
	}

	if (numPlatforms == 0) {
		fprintf(stderr, "Error: Platforms not found\n");
		exit(-1);
	}

	// Get the platforms
	platforms = new cl_platform_id[numPlatforms];
	if (clGetPlatformIDs(numPlatforms, platforms, NULL) != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not get the platform id\n");
		exit(-1);
	}

	// Find the platform specified in configuration
	bool found = false;
	for (int i = 0; i < numPlatforms && !found; ++i) {
		char pbuff[120];
		if (clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, sizeof(pbuff), pbuff, NULL) != CL_SUCCESS) {
			fprintf(stderr, "Error: Could not get the information of the platform\n");
			exit(-1);
		}

		// If the platform exists...
		if (strcmp(pbuff, platformVendor) == 0) {

			// Get the number of devices for the specified device type
			if ((status = clGetDeviceIDs(platforms[i], deviceType, 0, 0, (cl_uint*) &numDevicesPlatform)) != CL_SUCCESS) {
				if (status != CL_DEVICE_NOT_FOUND) {
					fprintf(stderr, "Error: Could not get the number of devices for the specified device type\n");
					exit(-1);
				}
				else {
					fprintf(stderr, "Error: Not exists any device on this platform for the specified device type\n");
					exit(-1);
				}
			}
			else {
				devices = new cl_device_id[numDevicesPlatform];

				// Get all devices for the specified device type
				if (clGetDeviceIDs(platforms[i], deviceType, numDevicesPlatform, devices, NULL) != CL_SUCCESS) {
					fprintf(stderr, "Error: Could not get the device id\n");
					exit(-1);
				}

				// Get the CPU or GPU device specified in configuration
				for (int dev = 0; dev < numDevicesPlatform && !found; ++dev) {
					char dbuff[120];
					if (clGetDeviceInfo(devices[dev], CL_DEVICE_NAME, sizeof(dbuff), dbuff, NULL) != CL_SUCCESS) {
						fprintf(stderr, "Error: Could not get the information of the device\n");
						exit(-1);
					}

					// If the device exists...
					if (strcmp(dbuff, deviceName) == 0) {
						device = devices[dev];
						found = true;
					}
				}

				if (!found) {
					fprintf(stderr, "Error: Not exists the device specified in configuration\n");
					exit(-1);
				}
			}
		}
	}

	if (device == NULL) {
		fprintf(stderr, "Error: Not exists the platform specified in configuration\n");
		exit(-1);
	}

	delete[] devices;
    
	// Open the file containing the kernels
	kernels = fopen(kernelName, "r");
	if (!kernels) {
		fprintf(stderr, "Error: An error ocurred opening the kernel file\n");
		exit(-1);
	}

	char *kernelSource = (char *)malloc(0x100000);
	const size_t source_size = fread(kernelSource, 1, 0x100000, kernels);
	fclose(kernels);

	// Create context
	context = clCreateContext(NULL, 1, &device, 0, 0, &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not get the context\n");
		exit(-1);
	}

	// Create Command queue
	command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not get the command queue\n");
		exit(-1);
	}

	// Create program
	program = clCreateProgramWithSource(context, 1, (const char **) &kernelSource, &source_size, &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not create the program\n");
		exit(-1);
	}

	// Build program for the device in the context
	char buildOptions[100];
	sprintf(buildOptions, "-D N_INSTANCES=%d -D N_FEATURES=%d -D N_OBJECTIVES=%d", N_INSTANCES, N_FEATURES, N_OBJECTIVES);
	if (clBuildProgram(program, 1, &device, buildOptions, 0, 0) != CL_SUCCESS) {
		char buffer[4096];
		fprintf(stderr, "Error: Could not build the program\n");
		if (clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, NULL) != CL_SUCCESS) {
			fprintf(stderr, "Error: Could not get the compilation errors\n");
			exit(-1);
		}
		fprintf(stderr, "%s\n", buffer);
		exit(-1);
	}

	// Create kernel
	long int usedMemory;
	if (deviceType == CL_DEVICE_TYPE_CPU) {
		kernel = clCreateKernel(program, "kmeansCPU", &status);

		// If the device is the CPU, the local work-size must be "1". On the contrary the program will fail
		if (wiLocal != 1) {
			wiLocal = 1;
			fprintf(stderr, "Warning: If the device is the CPU, the local work-size must be 1. Local work-size has been set to 1\n");
		}

		// Memory usage
		usedMemory = (KMEANS * N_INSTANCES * sizeof(cl_uchar)) << 1; // Mapping and NewMapping buffers
		usedMemory += sizeof(cl_uchar*); // Swap pointer
	}
	else {
		kernel = clCreateKernel(program, "kmeansGPU", &status);

		// Memory usage
		usedMemory = sizeof(individual); // Individual buffer
		usedMemory += (KMEANS * N_INSTANCES * sizeof(cl_uchar)) * 3; // Mapping, NewMapping and Aux buffers
		usedMemory += KMEANS * N_FEATURES * sizeof(cl_float); // Centroids buffer
		usedMemory += sizeof(cl_bool); // Converged variable
		usedMemory += KMEANS * N_INSTANCES * sizeof(cl_float); // DistCentroids buffer
		usedMemory += KMEANS * sizeof(cl_int); // Samples_in_k buffer
		usedMemory += sizeof(cl_float) * 3; // Local sums
	}

	if (status != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not create the kernel\n");
		exit(-1);
	}

	// Get the maximum local memory size
	long int maxMemory;
	if (clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(long int), &maxMemory, NULL) != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not get the maximum local memory of the device\n");
		exit(-1);
	}

	// Avoid exceeding the maximum local memory available
	if (usedMemory > maxMemory - 1024) { // 1024 bytes of margin
		fprintf(stderr, "Error: Local memory exceeded:\n\tMax memory: %ld bytes\n\tAllow memory: %ld bytes\n\tUsed memory: %ld bytes\n", maxMemory, maxMemory - 1024, usedMemory);
		exit(-1);
	}

	// The global work-size must be a multiple of the local work-size
	if (wiGlobal % wiLocal > 0) {
		fprintf(stderr, "Error: The global work-size must be a multiple of the local work-size\n");
		exit(-1);
	}

	// Start the time measure
	timeStart = omp_get_wtime();


	/********** Get the data base ***********/

	float dataBase[N_INSTANCES * N_FEATURES];
	readDataBase(dataBase, dataBaseName, N_INSTANCES, N_FEATURES);

	// Data base normalization
	normDataBase(dataBase, N_INSTANCES, N_FEATURES);

	// Create the dataBase buffer and write it in the device
	objDataBase = clCreateBuffer(context, CL_MEM_READ_ONLY, N_INSTANCES * N_FEATURES * sizeof(cl_float), 0, &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not create the OpenCL object containing the data base\n");
		exit(-1);
	}
	if (clEnqueueWriteBuffer(command_queue, objDataBase, CL_FALSE, 0, N_INSTANCES * N_FEATURES * sizeof(cl_float), dataBase, 0, NULL, NULL) != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not enqueue the OpenCL object containing the data base\n");
		exit(-1);
	}


	/********** Initialize the population and the individuals ***********/

	srand((unsigned int) time(NULL));
	const int totalIndividuals = POPULATION_SIZE << 1;

	// Population will have the parents and children (left half and right half respectively
	// This way is better for the performance
	individual *population = initPopulation(totalIndividuals, N_OBJECTIVES, N_FEATURES, maxFeatures);

	// Create the population buffer and write it in the device
	objPopulation = clCreateBuffer(context, CL_MEM_READ_WRITE, totalIndividuals * sizeof(individual), 0, &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not create the OpenCL object containing the population\n");
		exit(-1);
	}
	if (clEnqueueWriteBuffer(command_queue, objPopulation, CL_FALSE, 0, totalIndividuals * sizeof(individual), population, 0, NULL, NULL) != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not enqueue the OpenCL object containing the population\n");
		exit(-1);
	}


	/********** Multiobjective individual evaluation ***********/

	// Get the initial "KMEANS" centroids ***********/
	int selInstances[KMEANS];
	getCentroids(selInstances, N_INSTANCES);

	// Create the centroids buffer and write it in the device
	objSelInstances = clCreateBuffer(context, CL_MEM_READ_ONLY, KMEANS * sizeof(cl_int), 0, &status);
	if (status != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not create the OpenCL object containing the init centroids\n");
		exit(-1);
	}
	if (clEnqueueWriteBuffer(command_queue, objSelInstances, CL_FALSE, 0, KMEANS * sizeof(cl_int), selInstances, 0, NULL, NULL) != CL_SUCCESS) {
		fprintf(stderr, "Error: Could not enqueue the OpenCL object containing the init centroids\n");
		exit(-1);
	}

	evaluation(&device, population, 0, POPULATION_SIZE, totalIndividuals, maxIndividualsOnGPUKernel, N_INSTANCES, N_OBJECTIVES, &deviceType, &context, &kernel, &command_queue, wiGlobal, wiLocal, &objDataBase, &objSelInstances, &objPopulation);


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
		if (clEnqueueWriteBuffer(command_queue, objPopulation, CL_FALSE, 0, totalIndividuals * sizeof(individual), population, 0, NULL, NULL) != CL_SUCCESS) {
			fprintf(stderr, "Error: Could not enqueue the OpenCL object containing the population\n");
			exit(-1);
		}

		// Multiobjective individual evaluation
		int lastChild = POPULATION_SIZE + nChildren;
		evaluation(&device, population, POPULATION_SIZE, lastChild, totalIndividuals, maxIndividualsOnGPUKernel, N_INSTANCES, N_OBJECTIVES, &deviceType, &context, &kernel, &command_queue, wiGlobal, wiLocal, &objDataBase, &objSelInstances, &objPopulation);

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
	timeEnd = omp_get_wtime();
	fprintf(stdout, "%.10g ms\n", (timeEnd - timeStart) * 1000.0);
	fprintf(stdout, "%f\n", popHypervolume);

	// Generation of the data file and Gnuplot file for display the Pareto front
	generateGnuplot(dataName, plotName, imageName, population, nIndFront0, N_OBJECTIVES, referencePoint);


	/********** Resources used are released ***********/

	// The individuals (parents and children)
	delete[] population;
	delete[] platforms;

	// OpenCL resources
	status = clReleaseMemObject(objDataBase);
	status = clReleaseMemObject(objSelInstances);
	status = clReleaseMemObject(objPopulation);
	status = clReleaseKernel(kernel);
	status = clReleaseProgram(program);
	status = clReleaseCommandQueue(command_queue);
	status = clReleaseContext(context);
}