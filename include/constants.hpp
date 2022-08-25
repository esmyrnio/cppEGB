#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

namespace CGS{

    extern const double G;
    extern const double C;
    extern const double MSUN;
    extern const double Length;
    extern const double Time;
    extern const double Density;
};

namespace UNITS{

    extern const double KAPPA;
    extern const double KSCALE;
};

namespace INTEGRATION{

    extern const double r_max_surface;
    extern const int SMAX;
    extern const int RDIV;
    extern const double RMIN;
    extern const double DS;
    extern const double RLX_FAC;
    extern const double TOL;
    extern const int MAX_ITER;
};

namespace SURFACE{

    extern const double s_e;
    extern double p_surface;
    extern double e_surface;
    extern const double enthalpy_min;
};

namespace EXECUTION{
    extern char* call_type;
}

enum ARRAY_SIZE{
    SDIV = 401
};

#endif