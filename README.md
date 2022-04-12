# cpp_KEH_EGB

## Description

Program for obtaining static Neutron Star solutions in regularized 4D Einstein-Gauss-Bonnet gravity using the KEH/CST numerical scheme. Parts of the RNS code (https://gitlab.com/niksterg/rns1-1) have been employed. In addition, a C++ version of the Savitzky-Golay filter (https://github.com/arntanguy/sgsmooth) is used.

## Usage

Compile, using the provided makefile. Simply type "make" in the main directory, where the makefile is also located, to create the executionable file.
   
EGB.exe takes 7 inputs in-turn. The parameters are specified using the following flags:

1. **-f** *eos_name* (The EoS file name).
2. **-c** *coupling* (The EGB coupling constant in km^2).
4. **-e** *central_density* (The central energy density in CGS/10^15).
5. **-t** *relative_error* (The relative error for the iteration scheme).
6. **-m** *maximum_iterations* (The maximum number of iterations).
7. **-l** *relaxation_factor*
8. **-p** *print_option*:
    -  0: Prints gravitational mass M and radius R.
    -  1: Prints (0) along with the distance, metric, scalar, energy density and pressure profiles.


For example,

```
./EGB -f eosSLY.txt -c 10.0 -e 1.2 -t 1e-05 -m 100 -l 0.2 -p 0
```
