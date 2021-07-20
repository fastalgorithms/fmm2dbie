# Makefile for fmm3dbie
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 

# compiler, and linking from C, fortran
CC = gcc
CXX = g++
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy 

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 

FMMBIE_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	FMMBIE_INSTALL_DIR = ${HOME}/lib
endif

FMM_INSTALL_DIR=$(PREFIX_FMM)
ifeq ($(PREFIX_FMM),)
	FMM_INSTALL_DIR=${HOME}/lib
endif

LBLAS = -lblas -llapack

LIBS = -lm
DYLIBS = -lm
F2PYDYLIBS = -lm -lblas -llapack

LIBNAME=libfmm2dbie
DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)

LFMMLINKLIB = -lfmm2d 
LLINKLIB = -lfmm2dbie


# flags for MATLAB MEX compilation..
MFLAGS=-largeArrayDims -lgfortran -DMWF77_UNDERSCORE1 -lm  -ldl 
MWFLAGS=-c99complex 

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap/mwrap
MEXLIBS=-lm -lstdc++ -ldl -lgfortran


# For your OS, override the above by placing make variables in make.inc
-include make.inc

# update libs and dynamic libs to include appropriate versions of
# fmm3d
#
# Note: the static library is used for DYLIBS, so that fmm2d 
# does not get bundled in with the fmm3dbie dynamic library
#
LIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB) 
DYLIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB)
F2PYDYLIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB)
MEXLIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB)

# multi-threaded libs & flags needed
ifneq ($(OMP),OFF)
  FFLAGS += $(OMPFLAGS)
  LIBS += $(OMPLIBS)
  DYLIBS += $(OMPLIBS)
  F2PYDYLIBS += $(OMPLIBS)
  MEXLIBS += $(OMPLIBS)
endif

LIBS += $(LBLAS) $(LDBLASINC)
DYLIBS += $(LBLAS) $(LDBLASINC)
MEXLIBS += $(LBLAS) $(LDBLASINC)


# objects to compile
#
# Common objects
COM = src/common
COMOBJS = $(COM)/hkrand.o \
	$(COM)/dlaran.o $(COM)/lapack_wrap.o $(COM)/lapack_f77.o \
	$(COM)/legeexps.o $(COM)/prini_new.o \
	$(COM)/pyplot.o \
	$(COM)/sort.o \
	$(COM)/sparse_reps.o $(COM)/rotmat_gmres.o \
	$(COM)/get_fmm2d_thresh.o \

# Surface wrappers
SURF = src/curve_routs
SOBJS = $(SURF)/chunks.o $(SURF)/curve_routs.o \
	$(SURF)/chunk_near_point.o $(SURF)/chunk_interior.o

# Quadrature wrappers
QUAD = src/quadratures
QOBJS = $(QUAD)/near_field_routs.o $(QUAD)/adap_quads.o \
	$(QUAD)/self_quads2d.o

# Chunk adaptive integration routines
CHUNK = src/chunk_routs
COBJS = $(CHUNK)/dchunkints_main.o $(CHUNK)/zchunkints_main.o

# Helm wrappers
HELM = src/helm_wrappers
HOBJS = $(HELM)/helm_comb_dir2d.o

# Laplace wrappers
LAP = src/lap_wrappers
LOBJS = $(LAP)/lap_comb_dir2d.o

KERN = src/kernels
KOBJS = $(KERN)/helm_kernels.o $(KERN)/lap_kernels.o

OBJS = $(COMOBJS) $(SOBJS) $(COBJS) $(QOBJS) $(HOBJS) $(KOBJS) $(LOBJS)




.PHONY: usage lib install test test-dyn python 

default: usage

usage:
	@echo "-------------------------------------------------------------------------"
	@echo "Makefile for fmm3dbie. Specify what to make:"
	@echo "  make install - compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR) - compile and install the main library at custom location given by PREFIX"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make test - compile and run validation tests (will take around 30 secs)"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo ""
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=ON' for multi-threaded"
	@echo "-------------------------------------------------------------------------"



#
# implicit rules for objects (note -o ensures writes to correct dir)
#
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@



#
# build the library...
#
lib: $(STATICLIB) $(DYNAMICLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif

$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/

$(DYNAMICLIB): $(OBJS) 
	$(FC) -shared -fPIC $(FFLAGS) $(OBJS) -o $(DYNAMICLIB) $(DYLIBS) 
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/

install: $(STATICLIB) $(DYNAMICLIB)
	echo $(FMMBIE_INSTALL_DIR)
	mkdir -p $(FMMBIE_INSTALL_DIR)
	cp -f lib/$(DYNAMICLIB) $(FMMBIE_INSTALL_DIR)/
	cp -f lib-static/$(STATICLIB) $(FMMBIE_INSTALL_DIR)/
	[ ! -f lib/$(LIMPLIB) ] || cp lib/$(LIMPLIB) $(FMMBIE_INSTALL_DIR)/
	@echo "Make sure to include " $(FMMBIE_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(FMMBIE_INSTALL_DIR)  " "$(LLINKLIB) " -L"$(FMM_INSTALL_DIR)  " "$(LFMMLINKLIB)


#
# testing routines
#
test: $(STATICLIB) test/curv test/chunk test/quad test/helm test/near-point test/interior test/lap
	cd test/curve_routs; ./int2-curv
	cd test/curve_routs; ./int2-near-point
	cd test/curve_routs; ./int2-chunk-interior
	cd test/chunk_routs; ./int2-chunk
	cd test/quadratures; ./int2-quad
	cd test/helm_wrappers; ./int2-helm
	cd test/lap_wrappers; ./int2-lap
	cat print_testres.txt
	rm print_testres.txt


test/curv:
	$(FC) $(FFLAGS) test/curve_routs/test_curve_routs.f -o test/curve_routs/int2-curv lib-static/$(STATICLIB) $(LIBS)

NEARPTOBJS = test/curve_routs/test_near_point.o
test/near-point: $(NEARPTOBJS)
	$(FC) $(FFLAGS) test/curve_routs/test_near_point.f -o test/curve_routs/int2-near-point lib-static/$(STATICLIB) $(LIBS)

INTERIOROBJS = test/curve_routs/test_chunk_interior.o
test/interior: $(INTERIOROBJS)
	$(FC) $(FFLAGS) test/curve_routs/test_chunk_interior.f90 -o test/curve_routs/int2-chunk-interior lib-static/$(STATICLIB) $(LIBS)

CTOBJS = test/chunk_routs/test_dchunkints.o test/chunk_routs/test_zchunkints.o test/chunk_routs/test_bary.o

test/chunk: $(CTOBJS)
	$(FC) $(FFLAGS) test/chunk_routs/test_chunk_routs.f -o test/chunk_routs/int2-chunk $(CTOBJS) lib-static/$(STATICLIB) $(LIBS)

QTOBJS = test/quadratures/test_near_field_routs.o test/quadratures/test_self_adap.o

test/quad: $(QTOBJS)
	$(FC) $(FFLAGS) test/quadratures/test_quadratures.f -o test/quadratures/int2-quad $(QTOBJS) lib-static/$(STATICLIB) $(LIBS)

test/helm:
	$(FC) $(FFLAGS) test/helm_wrappers/test_helm_wrappers_qg_lp.f -o test/helm_wrappers/int2-helm lib-static/$(STATICLIB) $(LIBS)

test/lap:
	$(FC) $(FFLAGS) test/lap_wrappers/test_lap_wrappers_qg_lp.f -o test/lap_wrappers/int2-lap lib-static/$(STATICLIB) $(LIBS)

TESTCOMOBJS = test/common/test_rsc_to_csc.o

test/common: $(TESTCOMOBJS)
	$(FC) $(FFLAGS) test/common/test_rsc_to_csc.f -o test/common/int2-rsc lib-static/$(STATICLIB) $(LIBS)


#
# matlab..
#
MWRAPFILE = curve_routs
GATEWAY = $(MWRAPFILE)


matlab:	$(STATICLIB) matlab/$(GATEWAY).c 
	$(MEX) matlab/$(GATEWAY).c lib-static/$(STATICLIB) $(MFLAGS) -output matlab/curve_routs $(MEXLIBS);


mex:  $(STATICLIB)
	cd matlab;  $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw;\
	$(MEX) $(GATEWAY).c ../lib-static/$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE) $(MEXLIBS); \

#
# housekeeping routines
#
clean: objclean
	rm -f lib-static/*.a lib/*.so
	rm -f test/curve_routs/int2-curv
	rm -f test/chunk_routs/int2-chunk
	rm -f test/quadratures/int2-quad

objclean: 
	rm -f $(OBJS) 
	rm -f test/curve_routs/*.o
	rm -f test/chunk_routs/*.o
	rm -f test/quadratures/*.o
