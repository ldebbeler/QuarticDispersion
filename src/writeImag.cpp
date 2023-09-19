#include "writeImag.h"
#include <boost/multi_array.hpp>
#include <iostream>
#include "constants.h"

writeh5::writeh5(const scalingValues& scaling, const selfEnergyScaling& seScale, const seValuesImag& sev):
    m_scaling(scaling), m_seScale(seScale), m_selfEnergy(sev) {}

writeh5::writeh5() {}

void writeh5::write1dvector(H5::H5File file, const std::string& datasetname, const std::vector<double>& vector){
    auto datatype = H5::PredType::NATIVE_DOUBLE;
    datatype.setOrder( H5T_ORDER_LE );
    hsize_t dimv[1];
    dimv[0] = vector.size();
    H5::DataSpace vspace(1, dimv); 

    H5::DataSet dataset = file.createDataSet(datasetname, datatype, vspace);
    dataset.write(&vector[0], datatype);
}

// this probably does not work for an empty vector
void writeh5::write2dvector(H5::H5File file, const std::string& datasetname, const std::vector<std::vector<double>>& vector){
    // replace 2d vector by boost multi_array
    boost::multi_array<double, 2> array( boost::extents[vector.size()][vector[0].size()] );
    for(std::size_t i=0; i<vector.size(); i++){
        for(std::size_t j=0; j<vector[0].size(); j++){
            array[i][j] = vector[i][j];
        }
    }
    auto datatype = H5::PredType::NATIVE_DOUBLE;
    datatype.setOrder( H5T_ORDER_LE );

    std::vector<hsize_t> dimensions(array.shape(), array.shape() +2);
    H5::DataSpace dataspace(2, dimensions.data()); 

    std::cout << "Here3\n";
    H5::DataSet dataset = file.createDataSet(datasetname, datatype, dataspace);
    dataset.write(array.data(), datatype);
}

void writeh5::createGroup(H5::H5File file, const std::string& groupname){
    file.createGroup(groupname.c_str());
}

void writeh5::writeMainResults(H5::H5File file){
    H5::FloatType double_type(H5::PredType::NATIVE_DOUBLE);  
    H5::FloatType int_type(H5::PredType::NATIVE_INT);  
    H5::DataSpace att_space(H5S_SCALAR);

// group and data for scaling functions
    H5::Group scaling = file.createGroup("/BubbleScale");

    write1dvector(file, "/BubbleScale/Variable", m_scaling.m_x);
    write1dvector(file, "/BubbleScale/Real", m_scaling.m_real);
    write1dvector(file, "/BubbleScale/Imag", m_scaling.m_imag);

    H5::Attribute upper_cutoff = scaling.createAttribute( "Upper Cutoff Bubble", double_type, att_space );
    upper_cutoff.write( double_type, &UVCutoff);

    H5::Attribute lower_cutoff = scaling.createAttribute( "Lower Cutoff Bubble", double_type, att_space );
    lower_cutoff.write( double_type, &IRCutoff);

    H5::Attribute precision = scaling.createAttribute( "Precision", double_type, att_space );
    precision.write( double_type, &precision);

    H5::Attribute steps = scaling.createAttribute( "Steps", int_type, att_space );
    steps.write( int_type, &steps);

    H5::Attribute vF_const = scaling.createAttribute( "vF", double_type, att_space );
    vF_const.write( double_type, &vF);

    H5::Attribute b_const = scaling.createAttribute( "b", double_type, att_space );
    b_const.write( double_type, &b);

    H5::Attribute extrap = scaling.createAttribute( "extraPolate", double_type, att_space );
    extrap.write( double_type, &extraPolate);

    H5::Group seScale = file.createGroup("/SelfEnergyScaling");
    write1dvector(file, "/SelfEnergyScaling/krtilde", m_seScale.m_krtilde);
    write1dvector(file, "/SelfEnergyScaling/radPos", m_seScale.m_radPos);
    write1dvector(file, "/SelfEnergyScaling/radNeg", m_seScale.m_radNeg);
    write1dvector(file, "/SelfEnergyScaling/kttilde", m_seScale.m_kttilde);
    write1dvector(file, "/SelfEnergyScaling/tangPos", m_seScale.m_tangPos);
    write1dvector(file, "/SelfEnergyScaling/tangNeg", m_seScale.m_tangNeg);


    H5::Group selfEnergy = file.createGroup("/imaginarySelfEnergy");
    write1dvector(file, "/imaginarySelfEnergy/freqsRad", m_selfEnergy.m_freqsRad);
    write1dvector(file, "/imaginarySelfEnergy/radSe1", m_selfEnergy.m_radSe1);
    write1dvector(file, "/imaginarySelfEnergy/radSe2", m_selfEnergy.m_radSe2);
    write1dvector(file, "/imaginarySelfEnergy/radSe3", m_selfEnergy.m_radSe3);

    H5::Attribute krad1 = selfEnergy.createAttribute( "kr1", double_type, att_space );
    krad1.write( double_type, &m_selfEnergy.m_kr1);
    H5::Attribute krad2 = selfEnergy.createAttribute( "kr2", double_type, att_space );
    krad2.write( double_type, &m_selfEnergy.m_kr2);
    H5::Attribute krad3 = selfEnergy.createAttribute( "kr3", double_type, att_space );
    krad3.write( double_type, &m_selfEnergy.m_kr3);

    write1dvector(file, "/imaginarySelfEnergy/freqsTang", m_selfEnergy.m_freqsTang);
    write1dvector(file, "/imaginarySelfEnergy/tangSe1", m_selfEnergy.m_tangSe1);
    write1dvector(file, "/imaginarySelfEnergy/tangSe2", m_selfEnergy.m_tangSe2);
    write1dvector(file, "/imaginarySelfEnergy/tangSe3", m_selfEnergy.m_tangSe3);

    H5::Attribute ktang1 = selfEnergy.createAttribute( "kt1", double_type, att_space );
    ktang1.write( double_type, &m_selfEnergy.m_kt1);
    H5::Attribute ktang2 = selfEnergy.createAttribute( "kt2", double_type, att_space );
    ktang2.write( double_type, &m_selfEnergy.m_kt2);
    H5::Attribute ktang3 = selfEnergy.createAttribute( "kt3", double_type, att_space );
    ktang3.write( double_type, &m_selfEnergy.m_kt3);


}

