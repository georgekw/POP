#!/bin/bash
gfortran -c FEMModule.f90
gfortran -o Structure Structure.f90 FEMModule.o 
echo Structure.f90 compilation finished
gcc -lGL -lglut -lm Visual.c -o Visual
echo Visual.f90 compilation finished
echo ====================
./Structure
./Visual




