//File: selfenergy.h
#ifndef SELFENERGY_H
#define SELFENERGY_H

#include "bubble.h"

class selfenergy{
    bubble m_bubble;
public:
    double dispersion(double qr, double qt);
    double dispersion(double kr, double qr, double kt, double qt);

    selfenergy(bubble foo);

    std::complex<double> integrandPos(double omega, double kr, double kt, double qr, double qt);
    std::complex<double> integrandPosLog(double omega, double kr, double kt, double qr, double Q);
    std::complex<double> integrandNeg(double omega, double kr, double kt, double qr, double qt);
    std::complex<double> integrandNegLog(double omega, double kr, double kt, double qr, double Q);
    std::complex<double> integralLog(double omega, double kr, double kt);

    std::complex<double> integrandPosx(double qr, double qt, double x);
    std::complex<double> integrandPosLogx(double qr, double Qt, double x);
    double integralPosx(double x);
    std::complex<double> integrandNegx(double qr, double qt, double x);
    std::complex<double> integrandNegLogx(double qr, double Qt, double x);
    double integralNegx(double x);

    std::complex<double> integrandPosTang(double qr, double qt, double x);
    std::complex<double> integrandPosLogTang(double qr, double Qt, double x);
    double integralPosTang(double x);
    std::complex<double> integrandNegTang(double qr, double qt, double x);
    std::complex<double> integrandNegLogTang(double qr, double Qt, double x);
    double integralNegTang(double x);


    // real part
    // for now, only kr dependence
    /*
    std::complex<double> fluctuationPropagator(double nu, double qr, double qt);
    std::complex<double> BoseIntegrand1(double nu, double qr, double qt, double kr);
    std::complex<double> BoseIntegrand2(double nu, double qr, double qt, double kr);
    std::complex<double> FermiIntegrand(double nu, double qr, double qt, double kr);
    std::complex<double> RadialIntegrand(double nu, double qr, double qt, double kr);
    std::complex<double> RadialLogIntegrand(double Qn, double Qr, double Qt, double kr);
    double radialValue(double kr);
    */

    /*
    std::complex<double> greenFunction(double omega, double kr, double kt);
    std::complex<double> realIntegrand1(double omega, double kr, double kt, double nu, double qr, double qt);
    std::complex<double> realIntegrand2(double omega, double kr, double kt, double qr, double qt);

    std::complex<double> realLogIntegrand1(double omega, double kr, double kt, double Qn, double Qr, double Qt);
    std::complex<double> realLogIntegrand2(double omega, double kr, double kt, double Qr, double Qt);

    std::complex<double> realLogIntegral1(double omega, double kr, double kt);
    std::complex<double> realLogIntegral2(double omega, double kr, double kt);
    */

    //double realValue(double omega, double kr, double kt);
    
};

#endif //SELFENERGY_H
