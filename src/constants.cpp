#include "../include/constants.hpp"


char* EXECUTION::call_type = "get_MR";

const double CGS::G = 6.6732e-8;
const double CGS::C = 2.9979e+10;
const double CGS::MSUN = 1.989e+33;

const double UNITS::KSCALE = 1.112668301525780e-36;
const double UNITS::KAPPA = 1.346790806509621e+13;

const int INTEGRATION::RDIV = 901;
const double INTEGRATION::r_max_surface = 10.0;

const int INTEGRATION::SMAX = 1;
const double INTEGRATION::DS = (SMAX/(SDIV-1.0));

const double SURFACE::enthalpy_min = 1.0/(CGS::C*CGS::C);
const double SURFACE::s_e = 0.5;
