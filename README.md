# cppEGB

## Description

Code for obtaining static neutron star solutions in [regularized 4D Einstein-Gauss-Bonnet gravity](https://iopscience.iop.org/article/10.1088/1475-7516/2022/02/033)[^1], for [tabulated EoS](https://ui.adsabs.harvard.edu/abs/2021PhRvD.103l3004B/abstract)[^2], using the [KEH/CST numerical scheme](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..203C/abstract)[^3].
A [C++ version of the Savitzky-Golay filter](https://github.com/arntanguy/sgsmooth)[^4] is used, along with some modules from the [BOOST C++ library](https://github.com/boostorg/boost)[^5]. In addition the C++ code is wrapped using [SWIG](https://github.com/swig/swig
)[^6], into a python library which can be imported and used in any way, as shown in swig/test.py. 

## Usage

Compile, using the provided makefile.
   
EGB.exe takes 7 inputs in total. The parameters are specified using the following flags:

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

*typical input parameter space:*
   - *central pressure ~ 0.1-30.0 (10^35 dyn/cm^2)*
   - *coupling ~ 0-70 (km^2)*

[^1]:https://iopscience.iop.org/article/10.1088/1475-7516/2022/02/033
[^2]:https://ui.adsabs.harvard.edu/abs/2021PhRvD.103l3004B/abstract (table IX)
[^3]:https://academic.oup.com/mnras/article/237/2/355/976460, https://ui.adsabs.harvard.edu/abs/1992ApJ...398..203C/abstract, https://gitlab.com/niksterg/rns1-1/
[^4]:https://github.com/arntanguy/sgsmooth
[^5]:https://github.com/boostorg/boost
[^6]:https://github.com/swig/swig
