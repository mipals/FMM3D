# makefile overrides to use Intel ICC compiler, double prec only.

CC=gcc-9 
CXX=g++-9
FC=gfortran

CFLAGS += -I src 

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp 

#MATLAB interface:
MFLAGS += -L/usr/local/Cellar/gcc/9.1.0/lib/gcc/9
MEX = $(shell ls -d /Applications/MATLAB_R201*.app)/bin/mex
LIBS = -lm -lstdc++.6 

