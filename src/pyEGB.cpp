#include "../include/keh_cst.hpp"

// function to be exposed to Python using SWIG. Computes model's mass and radius.
void get_MR(String _eos_name, double _coupling, double _central_pressure, double _accuracy, 
            int _max_iter, double _relaxation, double* M, double* R, int* check, double* relaxation)
            
    {
        KEH_CST model(_eos_name, _coupling, _central_pressure, _accuracy, _relaxation, _max_iter, 0);
        model.compute_MR();

        *M = model.mass;
        *R = model.radius;
        *check = model.convergence_check_output;
        *relaxation = model.relaxation_output;
    }