#!/bin/bash
# To build the PDE model 
rm *.vtk *.txt *.out *.err *.o test *.stl tmp_data.dat
xlC -q32 -g -c *.cpp
xlC -q32 -g *.o newlib.a -o test -lxlsmp -L /opt/ibmcmp/xlf/bg/11.1/lib/ -lxlf90 
