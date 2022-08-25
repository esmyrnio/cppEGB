%module pyEGB
%{
    #include "../include/keh_cst.hpp"

    void get_MR(std::string _eos_name, double _coupling, double _central_pressure, double _accuracy, 
            int _max_iter, double _relaxation, double* M, double* R, int* check, double* relaxation);

%}

%include "std_string.i"
%include "typemaps.i"
%apply double *OUTPUT {double *M, double *R, double* relaxation, int* check};
%include <std_vector.i>


void get_MR(std::string _eos_name, double _coupling, double _central_pressure, double _accuracy, 
            int _max_iter, double _relaxation,  double* M, double* R, int* check, double* relaxation);






