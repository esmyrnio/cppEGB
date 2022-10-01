#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP
// centimetre–gram–second system of units
namespace CGS{
    extern const double G; //Newton's gravitational constant in cgs
    extern const double C; //speed of light in cgs
    extern const double MSUN; //solar mass in cgs
    extern const double Length; //length dimension in cgs
    extern const double Time; //time dimension in cgs
    extern const double Density; //density dimension in cgs
    extern const double MB; //baryon number density in cgs
};
// code's units
namespace UNITS{

    extern const double KAPPA; // square of length scale = 1e-15*C*C/G
    extern const double KSCALE; // KAPPA*G/(C*C*C*C)
};
// KEH/CST constants
namespace INTEGRATION{
    extern const int SMAX; // maximum value of s-coordinate 
    extern const double DS; // spacing in s-direction
    extern const double RLX_FAC; // relaxation factor
    extern const double TOL; // iteration routine tolerance
    extern const int MAX_ITER; // maximum iterations
};
// EoS related parameters and grid at surface
namespace SURFACE{
    extern const double s_e; // s-coordinate at surface
    extern double p_surface; // pressure at surface for given EoS
    extern double e_surface; // energy density at surface for given EoS
    extern const double enthalpy_min; // min enthalpy
};

namespace EXECUTION{
    extern char* CALL_TYPE;
}
enum ARRAY_SIZE{
    SDIV = 401 // iteration's grid size 
};
#endif