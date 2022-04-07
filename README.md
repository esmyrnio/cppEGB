# cpp_KEH_EGB

## Description

Program for obtaining static Neutron Star solutions in regularized 4D Einstein-Gauss-Bonnet gravity using the KEH/CST numerical scheme. Parts of the RNS code (https://gitlab.com/niksterg/rns1-1) have been employed. In addition, a C++ version of the Savitzky-Golay filter (https://github.com/arntanguy/sgsmooth) is used.

## Usage

Compile, using the provided makefile. Simply type "make" in the 
main directory, where the makefile is also located.

run ./EGB
   
EGB.exe takes 7 inputs in-turn:

1. The EoS file name (e.g. eosSLY.txt).
2. The EGB coupling constant in km^2 (e.g. 10.0).
4. The central energy density in CGS/10^15 (e.g 1.2).
5. The relative error for the iteration scheme (e.g 1e-05)
6. The maximum number of iterations (e.g 200)
7. The relaxation factor (e.g 0.2)
8. The print option 0 or 1:
    -  0: Prints gravitational mass M and radius R.
    -  1: Prints (0) along with the distance, metric, scalar, energy density and pressure profiles.
