#include "scaling.h"
#include "integrategsl.hpp"
#include "constants.h"


std::complex<double> scaling::ktPole(double s, double k0, double x){
    return s*std::sqrt(-0.75+s*std::sqrt(0.5*(1.0+(2.0*I*k0+x))));
}

std::complex<double> scaling::ktPole(double k0){
    return std::sqrt(std::sqrt(I*k0));
} 

std::complex<double> scaling::rootFunction(double k0, double x){
    return std::sqrt(2.0*(1.0+(2.0*I*k0+x)));
}

std::complex<double> scaling::rootFunction(double k0){
    return std::sqrt(4.0*I*k0);
}

std::complex<double> scaling::termx(double k0, double x){
    return 1.0/rootFunction(k0,x)*(1.0/ktPole(-1.0,k0,x)-1.0/ktPole(1.0,k0,x));
}

std::complex<double> scaling::termZero(double k0){
    return 1.0/rootFunction(k0)*(-I-1.0)/ktPole(k0);
}

std::complex<double> scaling::integrand(double k0, double x){
    double prefactor{ 1.0/(2*M_PI) };
    return prefactor*(termx(k0,x)-termZero(k0));
}

double scaling::realLogIntegral(double x){
    auto realIntegral = [&] (double K) { return std::exp(K)*integrand(std::exp(K),x).real(); };
    return integrater().integrate(realIntegral, std::log(IRCutoff), std::log(UVCutoff));
}

double scaling::imagLogIntegral(double x){
    auto imagIntegral = [&] (double K) { return std::exp(K)*integrand(std::exp(K),x).imag(); };
    return integrater().integrate(imagIntegral, std::log(IRCutoff), std::log(UVCutoff));
}

