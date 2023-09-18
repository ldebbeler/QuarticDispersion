#include "bubble.h"
#include "grid.h"
#include "constants.h"
#include<iostream>

bubble::bubble(): m_extraPos((1.0-I)/(4*M_PI*vF*std::pow(2*b,0.25))),
                        m_extraNeg((std::sqrt(2))/(4*M_PI*vF*std::pow(2*b,0.25))),
                        m_x(gridVector()), m_scalingReal(realVector()), m_scalingImag(imagVector()),
                        m_splineReal(interpolater1d(m_x,m_scalingReal)), 
                        m_splineImag(interpolater1d(m_x,m_scalingImag)) {}

std::vector<double> bubble::gridVector(){
    grid vec(minx, extraPolate, N);
    return vec.m_x;
}

std::vector<double> bubble::realVector(){
    std::vector<double> itx(m_x.size());
    for(std::size_t i=0; i<m_x.size(); i++){
        itx[i] = realLogIntegral(m_x[i]);
    }
    return itx;
}

std::vector<double> bubble::imagVector(){
    std::vector<double> itx(m_x.size());
    for(std::size_t i=0; i<m_x.size(); i++){
        itx[i] = imagLogIntegral(m_x[i]);
    }
    return itx;
}

std::complex<double> bubble::quadraticPos(double omega, double qr, double qt){
    double prefactor{ (3*std::pow(2*b,0.25))/(16*M_PI*vF) };
    return prefactor/std::pow(std::abs(omega - vF*qr),0.25)*std::pow(std::abs(qt),2)*(1.0+I)/2.0;
}

std::complex<double> bubble::quadraticNeg(double omega, double qr, double qt){
    double prefactor{ (3*std::pow(2*b,0.25))/(16*M_PI*vF) };
    return prefactor/std::pow(std::abs(omega - vF*qr),0.25)*std::pow(std::abs(qt),2)*1.0/std::sqrt(2.0);
}

std::complex<double> bubble::evaluateSingle(double omega, double qr, double qt){
    double val{ (omega-vF*qr)/(b*std::pow(std::abs(qt),4)) };
    if(val < -extraPolate+0.1){
        return std::pow(std::abs(omega-vF*qr),0.25)*m_extraNeg + quadraticNeg(omega, qr, qt);
    }
    else if(val > extraPolate-0.1){
        return std::pow(std::abs(omega-vF*qr),0.25)*m_extraPos + quadraticPos(omega, qr, qt);
    }
    else{
        return std::abs(qt)/(4*vF)*(m_splineReal.evaluate(val)+I*m_splineImag.evaluate(val));
    }
}

std::complex<double> bubble::evaluateCombo(double omega, double qr, double qt){
    return evaluateSingle(omega, qr, qt) + std::conj(evaluateSingle(-omega, qr, qt));
}

std::complex<double> bubble::evaluateStatic(double2 q){
    return evaluateSingle(0.0, q.x, q.y) + std::conj(evaluateSingle(0.0, q.x, q.y));
}

std::vector<double> bubble::getxVector(){
    return m_x;
}

std::vector<double> bubble::getReals(){
    return m_scalingReal;
}

std::vector<double> bubble::getImags(){
    return m_scalingImag;
}

