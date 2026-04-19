#$Id: makefile_oneapi  $

##Intel oneAPI (source setvars.sh)
ONEAPI_ROOT = /opt/intel/oneapi
COMPILER_ROOT = $(ONEAPI_ROOT)/compiler/latest
MKLROOT = $(ONEAPI_ROOT)/mkl/latest

##includefile
filein  = -I$(MKLROOT)/include -I/usr/include

##lib
LDLIBS  = -L$(COMPILER_ROOT)/linux/compiler/lib/intel64 -lifcore \
          -L$(MKLROOT)/lib/intel64 \
          -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
          -liomp5 -lpthread

# GSLlib
GSLIB = -L/usr/lib64 -lgsl

##option
CXXAPPFLAGS = -O1 -Wno-deprecated -Dfortran3 -Dnewc
f90appflags = -O1 -fpp

##Intel MPI 
f90 = mpiifx
f77 = mpiifx
CXX = mpiicpx
CC  = mpiicx
CLINKER = mpiicpx


###openMPI
#f90 = mpif90
#f77 = mpif77
#CXX = mpicxx
#CC  = mpicc
#CLINKER = mpicxx
