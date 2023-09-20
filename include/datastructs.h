#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H
#include <vector>

struct scalingValues{
    std::vector<double> m_x;        //access via bubblePos::getxVector()
    std::vector<double> m_real;     //access via bubblePos::getRealsPlus()
    std::vector<double> m_imag;     //access via bubblePos::getImagsPlus()

}; 

struct selfEnergyScaling{
    double m_Aplus;
    double m_Aminus;

    std::vector<double> m_krtilde;
    std::vector<double> m_radPos;
    std::vector<double> m_radNeg;

    std::vector<double> m_krFit;
    std::vector<double> m_radPosFit;
    std::vector<double> m_radNegFit;
    
    std::vector<double> m_kttilde;
    std::vector<double> m_tangPos;
    std::vector<double> m_tangNeg;

    std::vector<double> m_ktFit;
    std::vector<double> m_deltaTangPos;
    std::vector<double> m_deltaTangNeg;
};

struct seValuesImag{
    std::vector<double> m_freqsRad;
    std::vector<double> m_radSe1;
    double m_kr1;
    std::vector<double> m_radSe2;
    double m_kr2;
    std::vector<double> m_radSe3;
    double m_kr3;

    std::vector<double> m_freqsTang;
    std::vector<double> m_tangSe1;
    double m_kt1;
    std::vector<double> m_tangSe2;
    double m_kt2;
    std::vector<double> m_tangSe3;
    double m_kt3;

};

#endif //DATASTRUCTS_H
