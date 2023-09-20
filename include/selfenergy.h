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

    // imaginary part self energy
    std::complex<double> integrandPos(double omega, double kr, double kt, double qr, double qt);
    std::complex<double> integrandPosLog(double omega, double kr, double kt, double qr, double Q);
    std::complex<double> integrandNeg(double omega, double kr, double kt, double qr, double qt);
    std::complex<double> integrandNegLog(double omega, double kr, double kt, double qr, double Q);
    std::complex<double> integralLog(double omega, double kr, double kt);

    // radial scaling function. correct scaling function requires vF=1=b
    std::complex<double> integrandPosx(double qr, double qt, double x);
    std::complex<double> integrandPosLogx(double qr, double Qt, double x);
    double integralPosx(double x);
    std::complex<double> integrandNegx(double qr, double qt, double x);
    std::complex<double> integrandNegLogx(double qr, double Qt, double x);
    double integralNegx(double x);

    // tangential scaling function
    std::complex<double> integrandPosTang(double qr, double qt, double x);
    std::complex<double> integrandPosLogTang(double qr, double Qt, double x);
    double integralPosTang(double x);
    double subtractTangPos(double x);
    std::complex<double> integrandNegTang(double qr, double qt, double x);
    std::complex<double> integrandNegLogTang(double qr, double Qt, double x);
    double integralNegTang(double x);
    double subtractTangNeg(double x);

};

#endif //SELFENERGY_H
