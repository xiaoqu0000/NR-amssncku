
## filein  = -I/usr/include -I/usr/lib/x86_64-linux-gnu/mpich/include -I/usr/lib/x86_64-linux-gnu/openmpi/lib/ -I/usr/lib/gcc/x86_64-linux-gnu/11/ -I/usr/include/c++/11/

filein  = -I/usr/include/ -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/x86_64-linux-gnu/openmpi/lib/ -I/usr/lib/gcc/x86_64-linux-gnu/11/ -I/usr/include/c++/11/ -I/usr/lib/cuda/include

## LDLIBS  = -L/usr/lib/x86_64-linux-gnu -lmpich -L/usr/lib64 -L/usr/lib/gcc/x86_64-linux-gnu/11 -lgfortran
LDLIBS  = -L/usr/lib/x86_64-linux-gnu -L/usr/lib64 -L/usr/lib/gcc/x86_64-linux-gnu/11 -lgfortran -L/usr/lib/cuda/lib64 -lcudart -lmpi -lgfortran

CXXAPPFLAGS  = -O3 -Wno-deprecated -Dfortran3 -Dnewc
#f90appflags = -O3 -fpp
f90appflags  = -O3 -x f95-cpp-input
f90          = gfortran 
f77          = gfortran 
CXX          = g++
CC           = gcc
CLINKER      = mpic++ 

Cu = nvcc
CUDA_LIB_PATH = -L/usr/lib/cuda/lib64 -I/usr/include -I/usr/lib/cuda/include
#CUDA_APP_FLAGS = -c -g -O3 --ptxas-options=-v -arch compute_13 -code compute_13,sm_13 -Dfortran3 -Dnewc
CUDA_APP_FLAGS = -c -g -O3 --ptxas-options=-v -Dfortran3 -Dnewc
