/**
 * @file bd.h
 * @author Juan José Escobar Pérez
 * @date 15/06/2015
 * @brief Header file for data base reading
 *
 */

#ifndef BD_H
#define BD_H

/********************************* Methods ********************************/

/**
 * @brief Reading the database
 * @param dataBase The database which will contain the instances
 * @param fileName The file name of the database
 * @param nInstances The number of instances (rows) of the database
 * @param nFeatures The number of features (columns) of the database
 */
void readDataBase(float *dataBase, const char *fileName, const int nInstances, const int nFeatures);


/**
 * @brief The database is normalized between 0 and 1
 * @param dataBase Database to be normalized
 * @param nInstances The number of instances (rows) of the database
 * @param nFeatures The number of features (columns) of the database
 */
void normDataBase(float *dataBase, const int nInstances, const int nFeatures);

#endif