//File: bubble.h
#ifndef BUBBLE_H
#define BUBBLE_H

#include "scaling.h"
#include "interpolate1d.h"
#include "d2.hpp"
#include <vector>
#include <complex>

class bubble : public scaling{
    std::complex<double> m_extraPos;
    std::complex<double> m_extraNeg;
    std::vector<double> m_x;
    std::vector<double> m_scalingReal;
    std::vector<double> m_scalingImag;
    interpolater1d m_splineReal;
    interpolater1d m_splineImag;

    std::vector<double> gridVector();
    std::vector<double> realVector();
    std::vector<double> imagVector();

    std::complex<double> quadraticPos(double omega, double qr, double qt);
    std::complex<double> quadraticNeg(double omega, double qr, double qt);
public:
    bubble();

    std::complex<double> evaluateSingle(double omega, double qr, double qt);

    std::complex<double> evaluateCombo(double omega, double qr, double qt);

    std::complex<double> evaluateStatic(double2 q);

    std::vector<double> getxVector();

    std::vector<double> getReals();
    std::vector<double> getImags();

};

#endif //BUBBLE_H

