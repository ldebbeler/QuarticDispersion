#ifndef WRITERESULTS_H
#define WRITERESULTS_H
#include <H5Cpp.h>
#include <vector>
#include <string>
#include "datastructs.h"

class writeResults {
    scalingValues m_scaling;
    selfEnergyScaling m_seScale;
    seValuesImag m_selfEnergy;

public:
    writeResults(const scalingValues& scaling, const selfEnergyScaling& seScale, const seValuesImag& sev);

private:
    void write1dvector(H5::H5File file, const std::string& datasetname, const std::vector<double>& vector);  // create dataset with 1d vector

    void write2dvector(H5::H5File file, const std::string& datasetname, const std::vector<std::vector<double>>& vector); // create dataset with 2d vector

    void createGroup(H5::H5File file, const std::string& groupname);    // create group 

public:
    void writeMainResults(H5::H5File file);    
};    

#endif //WRITERESULTS_H
