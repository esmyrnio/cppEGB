#!/bin/bash

swig -python -c++ pyEGB.i
g++ -std=c++17 -c -O3 -fpic ../src/main.cpp ../src/pyEGB.cpp ../src/keh_cst.cpp ../src/constants.cpp ../src/eos.cpp ../src/savgol.cpp
g++ -std=c++17 -c -O3 -fpic pyEGB_wrap.cxx -I/usr/include/python2.7
g++ -shared -O3 pyEGB.o pyEGB_wrap.o main.o keh_cst.o constants.o eos.o savgol.o -o _pyEGB.so