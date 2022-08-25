# cpp_KEH_EGB

## Description

Program for obtaining static Neutron Star solutions in regularized 4D Einstein-Gauss-Bonnet gravity[^1], for tabulated EoS[^2], using the KEH/CST numerical scheme[^3].
A C++ version of the Savitzky-Golay filter[^3] is used, along with some modules from the BOOST C++ library[^4]. Both are included and implemented as header-only libraries. In addition the C++ code is wrapped using SWIG[^5], into a python library which can be imported and used in any way, as shown in swig/test.py. 

## Usage

Compile, using the provided makefile. Type "make" in the main directory, where the makefile is also located, to create the executionable file.
   
EGB.exe takes 7 inputs in-turn. The parameters are specified using the following flags:

1. **-f** *eos_name* (The EoS file name).
2. **-c** *coupling* (The EGB coupling constant in km^2).
4. **-e** *central_pressure* (The central pressure in dyn/cm^2/10^35).
5. **-t** *relative_error* (The relative error for the iteration scheme).
6. **-m** *maximum_iterations* (The maximum number of iterations).
7. **-l** *relaxation_factor*
8. **-p** *print_option*:
    -  0: Prints gravitational mass M and radius R.
    -  1: Prints (0) along with the distance, metric, scalar, energy density and pressure profiles.


For example,

```
./EGB -f ppsly4.cold.rns1.1.txt -c 20.0 -e 10.0 -t 1e-05 -m 100 -l 0.2 -p 0
```
[^1]:https://iopscience.iop.org/article/10.1088/1475-7516/2022/02/033
[^2]:https://ui.adsabs.harvard.edu/abs/2021PhRvD.103l3004B/abstract (table IX)
[^3]:https://academic.oup.com/mnras/article/237/2/355/976460, https://ui.adsabs.harvard.edu/abs/1992ApJ...398..203C/abstract
[^4]:https://github.com/arntanguy/sgsmooth
[^5]:https://github.com/boostorg/boost
[^6]:https://github.com/swig/swig
