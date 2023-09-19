#include <iostream>
#include "grid.h"
#include "bubble.h"
#include "interpolate1d.h"
#include "datastructs.h"
#include "selfenergy.h"
#include "writeImag.h"
#include "timer.hpp"
#include "constants.h"

int main() {

    Timer t;

    bubble bub;

    scalingValues scale;
    scale.m_x = bub.getxVector();
    scale.m_real= bub.getReals();
    scale.m_imag= bub.getImags();

    selfenergy SE(bub); 

    std::vector<double> x(201);
    std::vector<double> radPos(x.size());
    std::vector<double> radNeg(x.size());
    std::vector<double> y(x.size());
    std::vector<double> tangPos(x.size());
    std::vector<double> tangNeg(x.size());
    
    for(std::size_t i=0; i<x.size(); i++){
        double rval{ -5.0 + 0.05*i };
        x[i] = rval;
        radPos[i] = SE.integralPosx(rval);
        radNeg[i] = SE.integralNegx(rval);
        double tval{ 0.025*i };
        y[i] = tval;
        tangPos[i] = SE.integralPosTang(tval);
        tangNeg[i] = SE.integralNegTang(tval);
    }

    selfEnergyScaling seScale;
    seScale.m_krtilde = x;
    seScale.m_radPos = radPos;
    seScale.m_radNeg = radNeg;
    seScale.m_kttilde = y;
    seScale.m_tangPos = tangPos;
    seScale.m_tangNeg = tangNeg;

    std::vector<double> freqsRad(201);
    std::vector<double> radSe1(201);
    double kr1{-0.2 };
    std::vector<double> radSe2(201);
    double kr2{ 0.0 };
    std::vector<double> radSe3(201);
    double kr3{ 0.2 };
    for(std::size_t i=0; i<freqsRad.size(); i++){
        double omega{ -1.0 + 0.01*i };
        freqsRad[i] = omega;
        radSe1[i] = SE.integralLog(omega, kr1, 0.0).imag();
        radSe2[i] = SE.integralLog(omega, kr2, 0.0).imag();
        radSe3[i] = SE.integralLog(omega, kr3, 0.0).imag();
    }

    seValuesImag sev;
    sev.m_freqsRad = freqsRad;
    sev.m_radSe1 = radSe1;
    sev.m_radSe2 = radSe2;
    sev.m_radSe3 = radSe3;
    sev.m_kr1 = kr1;
    sev.m_kr2 = kr2;
    sev.m_kr3 = kr3;

    std::vector<double> freqsTang(201);
    std::vector<double> tangSe1(201);
    double kt1{ 0.0 };
    std::vector<double> tangSe2(201);
    double kt2{ std::pow(0.5, 0.25) };
    std::vector<double> tangSe3(201);
    double kt3{ std::pow(2.5, 0.25) };
    for(std::size_t i=0; i<freqsTang.size(); i++){
        double omega{ -10.0 + 0.1*i };
        freqsTang[i] = omega;
        tangSe1[i] = SE.integralLog(omega, 0.0, kt1).imag();
        tangSe2[i] = SE.integralLog(omega, 0.0, kt2).imag();
        tangSe3[i] = SE.integralLog(omega, 0.0, kt3).imag();
    }

    sev.m_freqsTang = freqsTang;
    sev.m_tangSe1 = tangSe1;
    sev.m_tangSe2 = tangSe2;
    sev.m_tangSe3 = tangSe3;
    sev.m_kt1 = kt1;
    sev.m_kt2 = kt2;
    sev.m_kt3 = kt3;


    H5::H5File file(filenameImag, H5F_ACC_TRUNC );
    writeh5 writeFile(scale, seScale, sev);
    writeFile.writeMainResults(file);


    std::cout << "Results written to: \n " << filenameImag << '\n';

    double tseconds{ t.elapsed() };
    int minutes = (int) tseconds/60.0;
    int seconds = (int) tseconds - 60.0*minutes;
    
    std::cout << "Time taken: " << minutes << "minutes and " << seconds << " seconds\n";
 
    return 0;
}
