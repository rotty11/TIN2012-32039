/**
 * @file xml.h
 * @author Juan José Escobar Pérez
 * @date 25/06/2015
 * @brief Header file for reading the config file
 *
 */

#ifndef XML_H
#define XML_H

/********************************* Includes *******************************/

#include "tinyxml2.h"

using namespace tinyxml2;

/********************************* Methods ********************************/

/**
 * @brief Gets from an XML file the parameter indicating the name of the database
 * @param configDoc XML type document
 * @return The database filename
 */
const char* getDataBaseName(XMLDocument *configDoc);


/**
 * @brief Gets from an XML file the parameter indicating the number of generations to generate (iterations of the program)
 * @param configDoc XML type document
 * @return The number of generations
 */
const int getNGenerations(XMLDocument *configDoc);


/**
 * @brief Gets from an XML file the parameter indicating the maximum number of features initially set to "1"
 * @param configDoc XML type document
 * @return The maximum number of features initially set to "1"
 */
const int getMaxFeatures(XMLDocument *configDoc);


/**
 * @brief Gets from an XML file the parameter indicating the number of individuals competing in the tournament
 * @param configDoc XML type document
 * @return The number of individuals competing in the tournament
 */
const int getTourSize(XMLDocument *configDoc);


/**
 * @brief Gets from an XML file the parameter indicating the value of the crossover distribution
 * @param configDoc XML type document
 * @return The value of the crossover distribution
 */
const int getCrossDistribution(XMLDocument *configDoc);


/**
 * @brief Gets from an XML file the parameter indicating the value of the mutation distribution
 * @param configDoc XML type document
 * @return The value of the mutation distribution
 */
const int getMutDistribution(XMLDocument *configDoc);


/**
 * @brief Gets from an XML file the parameter indicating the name of the file containing the fitness of the individuals in the first Pareto front
 * @param configDoc XML type document
 * @return The datafile name
 */
const char* getDataName(XMLDocument *configDoc);


/**
 * @brief Gets from an XML file the parameter indicating the name of the file containing the gnuplot code for data display
 * @param configDoc XML type document
 * @return The gnuplot filename
 */
const char* getPlotName(XMLDocument *configDoc);


/**
 * @brief Gets from an XML file the parameter indicating the name of the file containing the image with the data (graphic)
 * @param configDoc XML type document
 * @return The image filename
 */
const char* getImageName(XMLDocument *configDoc);

#endif