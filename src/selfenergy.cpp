#include "selfenergy.h"
#include "constants.h"
#include "d2.hpp"
#include "integrate2d.hpp"

double selfenergy::dispersion(double qr, double qt){
    return vF*qr + b*std::pow(std::abs(qt),4);
}

double selfenergy::dispersion(double kr, double qr, double kt, double qt){
    return -vF*(kr-qr) + b*std::pow(std::abs(kt-qt),4);
}

selfenergy::selfenergy(bubble foo): m_bubble(foo) {}

std::complex<double> selfenergy::integrandPos(double omega, double kr, double kt, double qr, double qt){
    double om{ std::abs(omega) };
    return 1.0/m_bubble.evaluateCombo(vF*qr-om,qr-b/vF*std::pow(qt-kt,4)+kr,qt);
}

std::complex<double> selfenergy::integrandPosLog(double omega, double kr, double kt, double qr, double Q){
    return std::exp(Q)*(integrandPos(omega,kr,kt,qr,std::exp(Q)) + integrandPos(omega,kr,kt,qr,-std::exp(Q)));
}

std::complex<double> selfenergy::integrandNeg(double omega, double kr, double kt, double qr, double qt){
    double om{ std::abs(omega) };
    return -1.0/m_bubble.evaluateCombo(vF*qr+om,qr-b/vF*std::pow(qt-kt,4)+kr,qt);
}

std::complex<double> selfenergy::integrandNegLog(double omega, double kr, double kt, double qr, double Q){
    return std::exp(Q)*(integrandNeg(omega,kr,kt,qr,std::exp(Q)) + integrandNeg(omega,kr,kt,qr,-std::exp(Q)));
}

// no factor 2 here because integral w.r.t. qt not symmetric
std::complex<double> selfenergy::integralLog(double omega, double kr, double kt){
    if(omega>0.0){
        auto integ = [&] (double2 q) -> std::complex<double>{
            double qr = q.x;
            double Q = q.y;
            return integrandPosLog(omega, kr, kt, qr, Q);
        };
        return integrate2Dmeasure2pi(0.0, omega/vF, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d);
    }
    else{
        auto integ = [&] (double2 q) -> std::complex<double>{
            double qr = q.x;
            double Q = q.y;
            return integrandNegLog(omega, kr, kt, qr, Q);
        };
        return integrate2Dmeasure2pi(omega/vF, 0.0, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d);
    }
}

// radial scaling function
std::complex<double> selfenergy::integrandPosx(double qr, double qt, double x){
    return -1.0/m_bubble.evaluateCombo(vF*qr-1.0,qr+x-b/vF*std::pow(qt,4),qt);
}

//factor of 2 because qt integral symmetric
std::complex<double> selfenergy::integrandPosLogx(double qr, double Qt, double x){
    return 2.0*std::exp(Qt)*integrandPosx(qr, std::exp(Qt), x);
}

double selfenergy::integralPosx(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandPosLogx(qr, Qt, x);
    };
    return integrate2Dmeasure2pi(0.0, 1.0/vF, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d).imag();
}

std::complex<double> selfenergy::integrandNegx(double qr, double qt, double x){
    return 1.0/m_bubble.evaluateCombo(vF*qr+1.0,qr+x-b/vF*std::pow(qt,4),qt);
}

std::complex<double> selfenergy::integrandNegLogx(double qr, double Qt, double x){
    return 2.0*std::exp(Qt)*integrandNegx(qr, std::exp(Qt), x);
}

double selfenergy::integralNegx(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandNegLogx(qr, Qt, x);
    };
    return integrate2Dmeasure2pi(-1.0/vF, 0.0, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d).imag();
}


// tangential scaling function
std::complex<double> selfenergy::integrandPosTang(double qr, double qt, double x){
    return -1.0/m_bubble.evaluateCombo(vF*qr-1.0,qr-b/vF*std::pow(qt-x,4),qt);
}

std::complex<double> selfenergy::integrandPosLogTang(double qr, double Qt, double x){
    return std::exp(Qt)*(integrandPosTang(qr, std::exp(Qt), x)+integrandPosTang(qr,-std::exp(Qt), x));
}

double selfenergy::integralPosTang(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandPosLogTang(qr, Qt, x);
    };
    return integrate2Dmeasure2pi(0.0, 1.0/vF, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d).imag();
}

double selfenergy::subtractTangPos(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandPosLogTang(qr, Qt, x) - integrandPosLogTang(qr, Qt, 0.0);
    };
    return integrate2Dmeasure2pi(0.0, 1.0/vF, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d).imag();
}

std::complex<double> selfenergy::integrandNegTang(double qr, double qt, double x){
    return 1.0/m_bubble.evaluateCombo(vF*qr+1.0,qr-b/vF*std::pow(qt-x,4),qt);
}

std::complex<double> selfenergy::integrandNegLogTang(double qr, double Qt, double x){
    return std::exp(Qt)*(integrandNegTang(qr, std::exp(Qt), x) + integrandNegTang(qr,-std::exp(Qt), x));
}

double selfenergy::integralNegTang(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandNegLogTang(qr, Qt, x);
    };
    return integrate2Dmeasure2pi(-1.0/vF, 0.0, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d).imag();
}

double selfenergy::subtractTangNeg(double x){
    auto integ = [&] (double2 q) -> std::complex<double>{
        double qr = q.x;
        double Qt = q.y;
        return integrandNegLogTang(qr, Qt, x) - integrandNegLogTang(qr, Qt, 0.0);
    };
    return integrate2Dmeasure2pi(-1.0/vF, 0.0, std::log(IR), std::log(UV)).integrate(integ, prec2d, steps2d).imag();
}

