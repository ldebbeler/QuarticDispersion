#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H
#include <vector>

struct scalingValues{
    std::vector<double> m_x;        //access via bubblePos::getxVector()
    std::vector<double> m_real;     //access via bubblePos::getRealsPlus()
    std::vector<double> m_imag;     //access via bubblePos::getImagsPlus()

}; 

struct bubbleValues{
    std::vector<double> m_omega;
    std::vector<double> m_omegReal;
    std::vector<double> m_omegImag;
    double m_omegqr;
    double m_omegqt;
    std::vector<double> m_qr;
    std::vector<double> m_radReal;
    std::vector<double> m_radImag;
    double m_radOmega;
    double m_radqt;
    std::vector<double> m_qt;
    std::vector<double> m_tangReal;
    std::vector<double> m_tangImag;
    double m_tangOmega;
    double m_tangqr;

};


struct seValuesImag{
    std::vector<double> m_freqs;
    std::vector<double> m_rads;
    //std::vector<double> m_SE;
    std::vector<std::vector<double>> m_SE;

    double m_kr;
    double m_kt;

};
/*
struct seValuesReal{
    std::vector<double> m_val;
    std::vector<double> m_SEPos;
    std::vector<double> m_SENeg;
};
*/
#endif //DATASTRUCTS_H
