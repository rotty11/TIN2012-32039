/**
 * @file individual.h
 * @author Juan José Escobar Pérez
 * @date 17/06/2015
 * @brief Header file containing the definition of the individuals
 *
 */

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

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
} individual;

#endif