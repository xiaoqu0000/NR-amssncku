# GCC + OpenMPI configuration (default open-source build).
# Selected automatically when COMPILER=gnu (the default in ../config.env).
# Self-contained for both CPU (ABE) and GPU (ABEGPU) builds.
#
# Prerequisites (Ubuntu/Debian):
#   sudo apt-get install build-essential gfortran openmpi-bin libopenmpi-dev
# Recommended optional libraries (auto-detected; build still works without them):
#   sudo apt-get install libhwloc-dev libnuma-dev
# GPU builds additionally need the CUDA Toolkit (nvcc + libcudart).

# --- 1. MPI compiler wrappers ---
# `mpicxx`/`mpicc`/`mpifort` are provided by libopenmpi-dev (or MPICH).
CXX      = mpicxx
CC       = mpicc
f90      = mpifort
f77      = mpifort
CLINKER  = mpicxx

# --- 2. Include paths ---
filein_CXX = -I. -IAHF -IAHF_GPU -I/usr/include
filein_C   = -I. -IAHF -IAHF_GPU -I/usr/include
# Legacy alias still referenced by the top-level makefile rules.
filein     = -I. -IAHF -IAHF_GPU -I/usr/include
HWLOC_CFLAGS  := $(shell pkg-config --cflags hwloc 2>/dev/null)
filein_CXX += $(HWLOC_CFLAGS)

# --- 3. CUDA (only used when building ABEGPU) ---
CUDA_CANDIDATE_DIRS = $(CUDA_HOME) $(CUDA_PATH) /usr/local/cuda \
                      /usr/local/cuda-13.1 /usr/local/cuda-13.0 /usr/local/cuda-12.9 \
                      /usr/local/cuda-12.8 /usr/local/cuda-12.6 /usr/local/cuda-12.4 \
                      $(patsubst %/bin/nvcc,%,$(shell which nvcc 2>/dev/null))
CUDA_DIR ?= $(firstword $(foreach dir,$(CUDA_CANDIDATE_DIRS),$(if $(wildcard $(dir)/bin/nvcc),$(dir),)))
# ============================================================================
# RTX5090(Blackwell)=120, A100=80, V100=70, RTX3090=86, RTX4090=89
# ============================================================================
CUDA_COMPUTE_CAPABILITY ?= 80

# The .cu files include <mpi.h>, so nvcc needs the MPI include dirs as well.
# OpenMPI wrappers answer --showme:incdirs; MPICH / Intel MPI wrappers use -show.
MPI_INC_DIRS := $(shell $(CXX) --showme:incdirs 2>/dev/null)
ifeq ($(strip $(MPI_INC_DIRS)),)
  MPI_INC_DIRS := $(shell $(CXX) -show 2>/dev/null | tr ' ' '\n' | sed -n 's/^-I//p' | sort -u)
endif
MPI_CUDA_INCLUDES := $(foreach dir,$(MPI_INC_DIRS),-I$(dir))
Cu = $(CUDA_DIR)/bin/nvcc
## Device-side (LLVM) optimization level for nvcc.
## The ABE_GPU_pro fused kernels exhaust host memory under nvcc -O3; its
## original makefile.inc.example uses -g (no LLVM opt) plus -Xptxas -O3.
## Compile that variant with `make CUDA_DEV_OPT=-g ABEGPU`.
CUDA_DEV_OPT ?= -O3
CUDA_APP_FLAGS = $(CUDA_DEV_OPT) -Xptxas -O3 --ptxas-options=-v -Dfortran3 -Dnewc -arch=sm_$(CUDA_COMPUTE_CAPABILITY) -I. -IAHF -IAHF_GPU -I$(CUDA_DIR)/include $(MPI_CUDA_INCLUDES)
CUDA_LIB_PATH = -L$(CUDA_DIR)/lib64 -I$(CUDA_DIR)/include $(MPI_CUDA_INCLUDES)
# Only pull in the CUDA runtime when it is actually installed, so CPU-only
# machines can still link ABE without a CUDA toolkit.
CUDA_LDFLAGS = $(if $(wildcard $(CUDA_DIR)/lib64/libcudart.so*),-L$(CUDA_DIR)/lib64 -lcudart)

# --- 4. Link libraries ---
# hwloc / numa are optional — only added when pkg-config finds them.
HWLOC_LIBS    := $(shell pkg-config --libs   hwloc 2>/dev/null)
NUMA_LIBS     := $(shell pkg-config --libs   numa  2>/dev/null)
LDLIBS    = -lm -lgfortran -lquadmath $(HWLOC_LIBS) $(NUMA_LIBS) $(CUDA_LDFLAGS)

# --- 5. Compiler flags ---
EXTRA_CXXFLAGS ?=
EXTRA_F90FLAGS ?=
CXXAPPFLAGS  = -std=c++20 -fopenmp -O3 -g -Wno-deprecated -Dfortran3 -Dnewc -fPIC $(filein_CXX) $(EXTRA_CXXFLAGS)
CFLAGS       = -fopenmp -O3 -g -fPIC -I. $(filein_C) -lm
f90appflags  = -fopenmp -O3 -g -cpp -ffree-line-length-none -funroll-loops $(EXTRA_F90FLAGS)
# Smaller interp loops do not benefit from aggressive fast-math; keep them plain.
f90flags_interp = -fopenmp -O3 -g -cpp -ffree-line-length-none $(EXTRA_F90FLAGS)
