#ifndef MPI_CUDA_AWARE_RUNTIME_H
#define MPI_CUDA_AWARE_RUNTIME_H

#include "macrodef.h"

#if MPI_CUDA_AWARE
extern "C" int MPIX_Query_cuda_support(void) __attribute__((weak));
#endif

inline bool amss_mpi_cuda_aware_runtime_enabled() {
#if MPI_CUDA_AWARE
    return MPIX_Query_cuda_support && MPIX_Query_cuda_support() == 1;
#else
    return false;
#endif
}

inline const char* amss_mpi_cuda_aware_runtime_status() {
#if MPI_CUDA_AWARE
    if (!MPIX_Query_cuda_support) {
        return "compiled_with_cuda_aware_mpi_but_runtime_query_unavailable";
    }

    return MPIX_Query_cuda_support() == 1
        ? "cuda_aware_mpi_runtime_enabled"
        : "compiled_with_cuda_aware_mpi_but_runtime_reports_no_cuda_support";
#else
    return "cuda_aware_mpi_disabled_at_compile_time";
#endif
}

#endif
