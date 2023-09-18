//File: scaling.h
#ifndef SCALING_H
#define SCALING_H

#include <complex>

class scaling{
    std::complex<double> ktPole(double s, double k0, double x);

    std::complex<double> ktPole(double k0);

    std::complex<double> rootFunction(double k0, double x);

    std::complex<double> rootFunction(double k0);

    std::complex<double> termx(double k0, double x);

    std::complex<double> termZero(double k0);

    std::complex<double> integrand(double k0, double x);

public:
    double realLogIntegral(double x);
    
    double imagLogIntegral(double x);

};
#endif // SCALING_H
