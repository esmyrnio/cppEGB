#include "../include/keh_cst.hpp"


typedef double* MASS_OUTPUT;
typedef double* RADIUS_OUTPUT;
typedef int* CONVERGENCE_OUTPUT;
typedef double* RLX_FAC_OUTPUT;


void get_MR(std::string _eos_name, double _coupling, double _central_pressure, double _accuracy, 
            int _max_iter, double _relaxation, MASS_OUTPUT M, RADIUS_OUTPUT R, CONVERGENCE_OUTPUT check, RLX_FAC_OUTPUT relaxation){

    KEH_CST model(_eos_name, _coupling, _central_pressure, _accuracy, _relaxation, _max_iter, 0);
    model.compute_MR();
    // model.printModel();
    *M = model.mass;
    *R = model.radius;
    *check = model.convergence_check_output;
    *relaxation = model.relaxation_output;

}
