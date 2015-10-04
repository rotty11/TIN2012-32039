/**
 * @file xml.cpp
 * @author Juan José Escobar Pérez
 * @date 25/06/2015
 * @brief File with the necessary implementation for reading the config file
 *
 */

/********************************* Includes *******************************/

#include "xml.h"

/********************************* Methods ********************************/

/**
 * @brief Gets from an XML file the parameter indicating the name of the database
 * @param configDoc XML type document
 * @return The database filename
 */
const char* getDataBaseName(XMLDocument *configDoc) {

	return (*configDoc).FirstChildElement("Config") -> FirstChildElement("DataBaseName") -> GetText();
}


/**
 * @brief Gets from an XML file the parameter indicating the number of generations to generate (iterations of the program)
 * @param configDoc XML type document
 * @return The number of generations
 */
const int getNGenerations(XMLDocument *configDoc) {

	return atoi((*configDoc).FirstChildElement("Config") -> FirstChildElement("NGenerations") -> GetText());
}


/**
 * @brief Gets from an XML file the parameter indicating the maximum number of features initially set to "1"
 * @param configDoc XML type document
 * @return The maximum number of features initially set to "1"
 */
const int getMaxFeatures(XMLDocument *configDoc) {

	return atoi((*configDoc).FirstChildElement("Config") -> FirstChildElement("MaxFeatures") -> GetText());
}


/**
 * @brief Gets from an XML file the parameter indicating the number of individuals competing in the tournament
 * @param configDoc XML type document
 * @return The number of individuals competing in the tournament
 */
const int getTourSize(XMLDocument *configDoc) {

	return atoi((*configDoc).FirstChildElement("Config") -> FirstChildElement("TournamentSize") -> GetText());
}


/**
 * @brief Gets from an XML file the parameter indicating the value of the crossover distribution
 * @param configDoc XML type document
 * @return The value of the crossover distribution
 */
const int getCrossDistribution(XMLDocument *configDoc) {

	return atoi((*configDoc).FirstChildElement("Config") -> FirstChildElement("CrossDistribution") -> GetText());
}


/**
 * @brief Gets from an XML file the parameter indicating the value of the mutation distribution
 * @param configDoc XML type document
 * @return The value of the mutation distribution
 */
const int getMutDistribution(XMLDocument *configDoc) {

	return atoi((*configDoc).FirstChildElement("Config") -> FirstChildElement("MutDistribution") -> GetText());
}


/**
 * @brief Gets from an XML file the parameter indicating the name of the file containing the kernels with the OpenCL code
 * @param configDoc XML type document
 * @return The OpenCL filename
 */
const char* getKernelName(XMLDocument *configDoc) {

	return (*configDoc).FirstChildElement("Config") -> FirstChildElement("KernelName") -> GetText();
}


/**
 * @brief Gets from an XML file the parameter indicating the name of the file containing the fitness of the individuals in the first Pareto front
 * @param configDoc XML type document
 * @return The datafile name
 */
const char* getDataName(XMLDocument *configDoc) {

	return (*configDoc).FirstChildElement("Config") -> FirstChildElement("DataName") -> GetText();
}


/**
 * @brief Gets from an XML file the parameter indicating the name of the file containing the gnuplot code for data display
 * @param configDoc XML type document
 * @return The gnuplot filename
 */
const char* getPlotName(XMLDocument *configDoc) {

	return (*configDoc).FirstChildElement("Config") -> FirstChildElement("PlotName") -> GetText();
}


/**
 * @brief Gets from an XML file the parameter indicating the name of the file containing the image with the data (graphic)
 * @param configDoc XML type document
 * @return The image filename
 */
const char* getImageName(XMLDocument *configDoc) {

	return (*configDoc).FirstChildElement("Config") -> FirstChildElement("ImageName") -> GetText();
}


/**
 * @brief Gets from an XML file the parameter indicating the type of device that will run the program
 * @param configDoc XML type document
 * @return CPU or GPU
 */
const cl_device_type getDeviceType(XMLDocument *configDoc) {

	const char *deviceTypeText = (*configDoc).FirstChildElement("Config") -> FirstChildElement("OpenCL") -> FirstChildElement("DeviceType") -> GetText();
	if (strcmp(deviceTypeText, "CPU") == 0) {
		return CL_DEVICE_TYPE_CPU;
	}
	else if (strcmp(deviceTypeText, "GPU") == 0) {
		return CL_DEVICE_TYPE_GPU;
	}
	else {
		fprintf(stderr, "Error: The device type specified in configuration must be CPU or GPU\n");
		exit(-1);
	}
}


/**
 * @brief Gets from an XML file the parameter indicating the name of the device vendor that will run the program
 * @param configDoc XML type document
 * @return The device vendor name
 */
const char* getPlatformVendor(XMLDocument *configDoc) {

	return (*configDoc).FirstChildElement("Config") -> FirstChildElement("OpenCL") -> FirstChildElement("PlatformVendor") -> GetText();
}


/**
 * @brief Gets from an XML file the parameter indicating the name of the device model that will run the program
 * @param configDoc XML type document
 * @return The device model name
 */
const char* getDeviceName(XMLDocument *configDoc) {

	return (*configDoc).FirstChildElement("Config") -> FirstChildElement("OpenCL") -> FirstChildElement("DeviceName") -> GetText();
}


/**
 * @brief Gets from an XML file the parameter indicating the total number of work-items (threads) that will run the program
 * @param configDoc XML type document
 * @return The total number of work-items
 */
const size_t getWiGlobal(XMLDocument *configDoc) {

	return ((size_t) atoi((*configDoc).FirstChildElement("Config") -> FirstChildElement("OpenCL") -> FirstChildElement("WiGlobal") -> GetText()));
}


/**
 * @brief Gets from an XML file the parameter indicating the number of work-items (threads) per compute unit that will run the program
 * @param configDoc XML type document
 * @return The number of work-items (threads) per compute unit
 */
size_t getWiLocal(XMLDocument *configDoc) {

	return ((size_t) atoi((*configDoc).FirstChildElement("Config") -> FirstChildElement("OpenCL") -> FirstChildElement("WiLocal") -> GetText()));
}


/**
 * @brief Gets from an XML file the parameter indicating the maximum number of individuals to be processed in a single execution of the kernel
 * @param configDoc XML type document
 * @return The maximum number of individuals to be processed in a single execution of the kernel
 */
const int getMaxIndividualsOnGPUKernel(XMLDocument *configDoc) {

	return atoi((*configDoc).FirstChildElement("Config") -> FirstChildElement("OpenCL") -> FirstChildElement("MaxIndividualsOnGPUKernel") -> GetText());
}