//File: constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <complex>


const std::string filenameImag{ "data/quarticDispersion.h5" };


inline constexpr std::complex<double> I{ std::complex<double>(0.0,1.0) };   // imaginary unit
// correct scaling function for self energy requires here vF=1=b
// I keep this version to illustrate the slight deviation in the Paper.
// The correct scaling functions do not depend on vF and b
inline constexpr double vF{ 1.0 };  // Fermi Velocity
inline constexpr double b{ 1.0 };   // quartic coefficient

// bubble calculation via scaling function. With 1d integral (gsl)
inline constexpr double IRCutoff{ 1e-8 };   // lower Cutoff for k Integral Bubble(scaling Function)
inline constexpr double UVCutoff{ 1e5  };   // upper Cutoff for k Integral Bubble(sclaing Function)
inline constexpr int steps{ 5000 };
inline constexpr double precision{ 1e-6 };

inline constexpr bool logSpacing{ true }; //grid points of scaling functions have logarithmic spacing. Otherwise linear spacing
inline constexpr std::size_t N{ 10000 };    //Number of grid points for scaling functions
inline constexpr double minx{ 1e-5 };     //smallest nonzero grid value for logarithmic spacing
inline constexpr double extraPolate{ 100.0 }; // scaling function is extrapolated for values larger than this

// self energy calculation
// limits:
inline constexpr double IR{ 1e-8 };         // lower cutoff for q Integral Self Energy
inline constexpr double UV{ 1e5 };         // upper cutoff for q Integral Self Energy
// 2d precision
inline constexpr int steps2d{ 1000*400 };
inline constexpr double prec2d{ 1e-10 };

#endif  //CONSTANTS_H
