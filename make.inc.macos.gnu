# makefile overrides
# OS:       macOS
# Compiler: gfortran 9.X
# OpenMP:   enabled
#

CC=gcc
CXX=g++
FC=gfortran

ifeq ($(PREFIX),)
    FMMBIE_INSTALL_DIR=/usr/local/lib
endif

ifeq ($(PREFIX_FMM),)
    FMM_INSTALL_DIR=/usr/local/lib
endif

FFLAGS= -fPIC -O3 -funroll-loops -std=legacy -w 

CFLAGS += -I src 

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp

LBLAS=-framework accelerate

# MATLAB interface:
MFLAGS += -L/usr/local/lib/gcc/12
MEX = $(shell ls -d /Applications/MATLAB_R20**.app)/bin/mex
#LIBS = -lm -lstdc++.6
#MEXLIBS= -lm -lstdc++.6 -lgfortran -ldl

test-helm: $(STATICLIB) test/helm
	cd test/helm_wrappers; ./int2-helm
	cd test/helm_wrappers; ./int2-helm-mat

