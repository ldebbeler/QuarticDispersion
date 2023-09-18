//File: constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <complex>


const std::string filenameImag{ "data/tangentialDependence.h5" };
const std::string filenameFx{ "data/RadialScaling5.h5" };
const std::string filenameReal{ "data/radialMomentumDependenceIR12.h5" };
const std::string filenameRealParallel{ "data/radialMomentumDependenceParallel12.h5" };
const std::string filenameTang{ "data/staticBubbleTangFlow10.h5" };


inline constexpr std::complex<double> I{ std::complex<double>(0.0,1.0) };   // imaginary unit
inline constexpr double vF{ 1.0 };  // Fermi Velocity
inline constexpr double b{ 0.5 };   // quartic coefficient

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
//inline constexpr bool staticApprox{ true };     // static approximation of the bubble for real part SE calculation
// limits:
inline constexpr double IR{ 1e-12 };         // lower cutoff for q Integral Self Energy
inline constexpr double UV{ 100.0 };         // upper cutoff for q Integral Self Energy
// 2d precision
inline constexpr int steps2d{ 1000*400 };
inline constexpr double prec2d{ 1e-12 };
// 3d precision
inline constexpr double IRt{ IR };
inline constexpr double UVt{ UV };
inline constexpr int steps3d{ 1000*1000*200 };  //attempt to do 1000*1000*1000 on cluster (might take approx. 1 hour per value)
inline constexpr double prec3d{ 1e-8 };

// grid for self energy values (momentum dependence)
inline constexpr std::size_t M { 25 }; //Number of grid points
inline constexpr double mink{ 1e-4 };   //minimal k
inline constexpr double maxk{ 1.0 };    //maximal k

#endif  //CONSTANTS_H
