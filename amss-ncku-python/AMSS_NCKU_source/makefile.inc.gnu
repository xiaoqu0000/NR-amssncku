# GCC + OpenMPI configuration (default open-source build).
# Selected automatically when COMPILER=gnu (the default in ../config.env).
#
# Prerequisites (Ubuntu/Debian):
#   sudo apt-get install build-essential gfortran openmpi-bin libopenmpi-dev
# Recommended optional libraries (auto-detected; build still works without them):
#   sudo apt-get install libhwloc-dev libnuma-dev

# --- 1. MPI compiler wrappers ---
# `mpicxx`/`mpicc`/`mpifort` are provided by libopenmpi-dev (or MPICH).
CXX      = mpicxx
CC       = mpicc
f90      = mpifort
f77      = mpifort
CLINKER  = mpicxx

# --- 2. Include paths ---
filein_CXX = -I. -I/usr/include
filein_C   = -I. -I/usr/include
# Legacy alias still referenced by the top-level makefile rules.
filein     = -I. -I/usr/include

# --- 3. Link libraries ---
# hwloc / numa are optional — only added when pkg-config finds them.
HWLOC_CFLAGS  := $(shell pkg-config --cflags hwloc 2>/dev/null)
HWLOC_LIBS    := $(shell pkg-config --libs   hwloc 2>/dev/null)
NUMA_LIBS     := $(shell pkg-config --libs   numa  2>/dev/null)
filein_CXX += $(HWLOC_CFLAGS)
LDLIBS    = -lm -lgfortran -lquadmath $(HWLOC_LIBS) $(NUMA_LIBS)

# --- 4. Compiler flags ---
EXTRA_CXXFLAGS ?=
EXTRA_F90FLAGS ?=
CXXAPPFLAGS  = -std=c++20 -fopenmp -O3 -g -Wno-deprecated -Dfortran3 -Dnewc -fPIC $(filein_CXX) $(EXTRA_CXXFLAGS)
CFLAGS       = -fopenmp -O3 -g -fPIC -I. $(filein_C) -lm
f90appflags  = -fopenmp -O3 -g -cpp -ffree-line-length-none -funroll-loops $(EXTRA_F90FLAGS)
# Smaller interp loops do not benefit from aggressive fast-math; keep them plain.
f90flags_interp = -fopenmp -O3 -g -cpp -ffree-line-length-none $(EXTRA_F90FLAGS)

# --- 5. CUDA (only used when building ABEGPU) ---
Cu = nvcc
CUDA_LIB_PATH = -L/usr/lib/cuda/lib64 -I/usr/include -I/usr/lib/cuda/include
CUDA_APP_FLAGS = -c -g -O3 --ptxas-options=-v -Dfortran3 -Dnewc
