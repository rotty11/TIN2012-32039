/**
 * @file evaluation.cl
 * @author Juan José Escobar Pérez
 * @date 12/07/2015
 * @brief File with the necessary implementation for the k-means algorithm in OpenCL
 *
 */

/*********************************** Defines *********************************/

/**
 * @brief Number of clusters (and centroids) used for the kmeans algorithm
 */
#define KMEANS 3


/**
 * @brief Max iterations for the convergence of centroids
 */
#define MAX_ITER_KMEANS 100


/********************************* Structures ********************************/

/**
 * @brief Structure containing the population and their attributes
 *
 * Have five attributes: chromosome, fitness, nSelFeatures, rank and crowding
 */
typedef struct individual {


	/**
	 * @brief Vector denoting the selected features
	 *
	 * Values: Zeros or ones
	 */
	unsigned char chromosome[N_FEATURES];


	/**
	 * @brief Individual fitness for the multiobjective functions
	 *
	 * Values: Each position contains an objective function
	 */
	float fitness[N_OBJECTIVES];


	/**
	 * @brief Number of selected features
	 */
	int nSelFeatures;


	/**
	 * @brief Range of the individual (Pareto front)
	 */
	int rank;


	/**
	 * @brief Crowding distance of the individual
	 *
	 * The values are positives or infinites
	 */
	float crowding;
};

/********************************* Functions *********************************/

/**
 * @brief Computes (*source + operand) and store result at location pointed by source
 * @param source Variable where the atomic sum is perform
 * @param operand Float point operand
 */
inline void AtomicAddL(volatile __local float *source, const float operand) {

	union {
		unsigned int intVal;
		float floatVal;
	} newVal;

	union {
		unsigned int intVal;
		float floatVal;
	} prevVal;

	do {
		prevVal.floatVal = *source;
		newVal.floatVal = prevVal.floatVal + operand;
	} while (atomic_cmpxchg((volatile __local unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}


/**
 * @brief Computes (*source + operand) and store result at location pointed by source
 * @param source Variable where the atomic sum is perform
 * @param operand Unsigned char operand
 */
inline void AtomicUAddL(volatile __local unsigned char *source, const unsigned char operand) {

	union {
		unsigned int intVal;
		unsigned char floatVal;
	} newVal;

	union {
		unsigned int intVal;
		unsigned char floatVal;
	} prevVal;

	do {
		prevVal.floatVal = *source;
		newVal.floatVal = prevVal.floatVal + operand;
	} while (atomic_cmpxchg((volatile __local unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}


/********************************* Kernels ********************************/

/**
 * @brief Computes the k-means algorithm in a OpenCL CPU device
 * @param pop OpenCL object which contains the current population. The object is stored in global memory
 * @param selInstances OpenCL object which contains the instances choosen as initial centroids. The object is stored in constant memory
 * @param dataBase OpenCL object which contains the database. The object is stored in constant memory
 * @param begin The first individual to evaluate
 * @param end The "end-1" position is the last individual to evaluate
 * @param mapping Auxiliar buffer for swapping with the table of minimum distances (newMapping). The reason: avoid copies. The buffer is stored in local memory
 * @param newMapping Auxiliar buffer which contains the table of minimum distances. The buffer is stored in local memory
 */
__kernel void kmeansCPU(__global struct individual *pop, __constant int *selInstances, __constant float *dataBase, const int begin, const int end, __local uchar *mapping, __local uchar *newMapping) {

	uint id = get_global_id(0);
	uint nWorkItems = get_global_size(0);

	// Each work-item compute an individual (master-slave as a deck algorithm)
	for (int ind = begin + id; ind < end; ind += nWorkItems) {
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
				__local uchar *aux = newMapping;
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
		//fitness[2] = (float) nSelFeatures;
	}
}


/**
 * @brief Computes the k-means algorithm in a OpenCL CPU device
 * @param pop OpenCL object which contains the current population. The object is stored in global memory
 * @param selInstances OpenCL object which contains the instances choosen as initial centroids. The object is stored in constant memory
 * @param dataBase OpenCL object which contains the database. The object is stored in constant memory
 * @param begin The first individual to evaluate
 * @param end The "end-1" position is the last individual to evaluate
 */
__kernel void kmeansGPU(__global struct individual *pop, __constant int *selInstances, __global float *dataBase, const int begin, const int end) {

	uint localId = get_local_id(0);
	uint localSize = get_local_size(0);
	uint groupId = get_group_id(0);
	uint numGroups = get_num_groups(0);

	const int totalCoord = KMEANS * N_FEATURES;
	const int totalDist = KMEANS * N_INSTANCES;

	// The individual is cached into local memory to improve performance
	__local struct individual indiv;
	__local uchar mapping[KMEANS * N_INSTANCES];
	__local uchar newMapping[KMEANS * N_INSTANCES];
	__local uchar aux[KMEANS * N_INSTANCES];
	__local float centroids[KMEANS * N_FEATURES];
	__local bool converged;
	__local float distCentroids[KMEANS * N_INSTANCES];
	__local int samples_in_k[KMEANS];
	__local float sumWithin;
	__local float sumInter;
	__local float sumLocal;

	// Each work-group compute an individual (master-slave as a deck algorithm)
	for (int ind = begin + groupId; ind < end; ind += numGroups) {

		// The individual is cached to local memory for improve performance
		if (localId == 0) {
			indiv = pop[ind];
			converged = false;
			sumWithin = 0.0f;
			sumInter = 0.0f;
			sumLocal = 0.0f;
		}

		// Initialize the mapping table
		for (int i = localId; i < totalDist; i += localSize) {
			mapping[i] = 0;
		}

		// Syncpoint
		barrier(CLK_LOCAL_MEM_FENCE);

		// The centroids will have the selected features of the individual
		for (int k = 0; k < KMEANS; ++k) {
			int posDataBase = selInstances[k] * N_FEATURES;
			int posCentr = k * N_FEATURES;
			for (int f = localId; f < N_FEATURES; f += localSize) {
				if (indiv.chromosome[f] & 1) {
					centroids[posCentr + f] = dataBase[posDataBase + f];
				}
			}
		}


		/******************** Convergence process *********************/

		// To avoid poor performance, at most "MAX_ITER_KMEANS" iterations are executed
		for (int maxIter = 0; maxIter < MAX_ITER_KMEANS && !converged; ++maxIter) {

			// The mapping table is cleaned in each iteration
			for (int i = localId; i < totalDist; i += localSize) {
				newMapping[i] = 0;
			}
			for (int i = localId; i < KMEANS; i += localSize) {
				samples_in_k[i] = 0;
			}

			// Syncpoint
			barrier(CLK_LOCAL_MEM_FENCE);

			// Calculate all distances (Euclidean distance) between each instance and the centroids
			for (int i = localId; i < N_INSTANCES; i += localSize) {
				float minDist = INFINITY;
				int selectCentroid = -1;
				int pos = N_FEATURES * i;
				for (int k = 0; k < KMEANS; ++k) {
					float sum = 0.0f;
					int posCentr = k * N_FEATURES;
					int posDistCentr = k * N_INSTANCES;
					for (int f = 0; f < N_FEATURES; ++f) {
						if (indiv.chromosome[f] & 1) {
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

				atomic_inc(&samples_in_k[selectCentroid]);
				newMapping[(selectCentroid * N_INSTANCES) + i] = 1;
			}

			if (localId == 0) {
				converged = true;
			}

			// Syncpoint
			barrier(CLK_LOCAL_MEM_FENCE);

			// Has the algorithm converged?
			for (int k = 0; k < KMEANS && converged; ++k) { 
				int posMap = k * N_INSTANCES;
				for (int i = localId; i < N_INSTANCES && converged; i += localSize) {
					if (newMapping[posMap + i] != mapping[posMap + i]) {
						converged = false;
					}
				}
			}

			// Syncpoint
			barrier(CLK_LOCAL_MEM_FENCE);

			if (!converged) {

				// Update the position of the centroids
				for (int k = 0; k < KMEANS; ++k) {
					int posCentr = k * N_FEATURES;
					int posMap = k * N_INSTANCES;
					for (int f = localId; f < N_FEATURES; f += localSize) {
						float sum = 0.0f;
						if (indiv.chromosome[f] & 1) {
							for (int i = 0; i < N_INSTANCES; ++i) {
								if (newMapping[posMap + i] == 1) {
									sum += dataBase[(N_FEATURES * i) + f];
								}
							}
							centroids[posCentr + f] = sum / samples_in_k[k];
						}
					}
				}

				// Syncpoint
				barrier(CLK_LOCAL_MEM_FENCE);

				// Swap mapping tables
				for (int i = localId; i < totalDist; i += localSize) {
					aux[i] = newMapping[i];
					newMapping[i] = mapping[i];
					mapping[i] = aux[i];
				}
			}

			// Syncpoint
			barrier(CLK_LOCAL_MEM_FENCE);
		}


		/************ Minimize the within-cluster and maximize Inter-cluster sum of squares (WCSS and ICSS) *************/

		for (int k = 0; k < KMEANS; ++k) {
			int posCentr = k * N_FEATURES;
			int posDistCentr = k * N_INSTANCES;
			float privateWithinSum = 0.0f;

			// Within-cluster
			for (int i = localId; i < N_INSTANCES; i += localSize) {
				if (mapping[posDistCentr + i] == 1) {
					privateWithinSum += distCentroids[posDistCentr + i];
				}
			}

			AtomicAddL(&sumWithin, privateWithinSum);

			for (int i = posCentr + N_FEATURES; i < totalCoord; i += N_FEATURES) {

				float sum = 0.0f;
				for (int f = localId; f < N_FEATURES; f += localSize) {
					if (indiv.chromosome[f] & 1) {
						sum += (centroids[posCentr + f] - centroids[i + f]) * (centroids[posCentr + f] - centroids[i + f]);
					}
				}
				AtomicAddL(&sumLocal, sum);

				// Syncpoint
				barrier(CLK_LOCAL_MEM_FENCE);

				if (localId == 0) {
					sumInter += sqrt(sumLocal);
					sumLocal = 0.0f;
				}

				// Syncpoint
				barrier(CLK_LOCAL_MEM_FENCE);
			}
		}

		// Syncpoint
		barrier(CLK_LOCAL_MEM_FENCE);

		if (localId == 0) {

			// First objective function (Within-cluster sum of squares (WCSS))
			pop[ind].fitness[0] = sumWithin;

			// Second objective function (Inter-cluster sum of squares (ICSS))
			pop[ind].fitness[1] = sumInter;

			// Third objective function (Number of selected features)
			//fitness[2] = (float) nSelFeatures;
		}
	}
}