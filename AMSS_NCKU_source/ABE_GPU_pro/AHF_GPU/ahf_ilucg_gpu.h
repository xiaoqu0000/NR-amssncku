#ifndef AHF_ILUCG_GPU_H
#define AHF_ILUCG_GPU_H

#include <cuda_runtime.h>

struct GpuAhfLinearSolver;

enum GpuAhfLinearStatus {
    GPU_AHF_LINEAR_SUCCESS = 1,
    GPU_AHF_LINEAR_INVALID_ARGUMENT = -1,
    GPU_AHF_LINEAR_CUDA_FAILURE = -2,
    GPU_AHF_LINEAR_NONFINITE_INPUT = -3,
    GPU_AHF_LINEAR_NONFINITE_OUTPUT = -4,
    GPU_AHF_LINEAR_GETRF_FAILURE = -5,
    GPU_AHF_LINEAR_GETRS_FAILURE = -6,
    GPU_AHF_LINEAR_GETRF_ARGUMENT = -7,
    GPU_AHF_LINEAR_GETRS_INFO = -8,
    GPU_AHF_LINEAR_NEWTON_STATE = -9,
    GPU_AHF_LINEAR_SINGULAR_BASE = -32
};

// Persistent sparse-input FP64 partial-pivot dense LU (cuSOLVER Dgetrf/Dgetrs).
// n and nnz describe the validated, zero-based, strictly sorted static CSR.
// Device rows/columns are borrowed, not revalidated; they and stream must outlive
// the solver. Creation checks allocation sizes and throws on invalid dimensions,
// allocation or library errors, releasing any partially constructed resources.
GpuAhfLinearSolver* gpu_ahf_linear_create(cudaStream_t stream, int n, int nnz,
                                        const int* d_rows, const int* d_columns);
void gpu_ahf_linear_destroy(GpuAhfLinearSolver* solver) noexcept;

// Device arrays: values[nnz], rhs[n], solution[n]. Output must not overlap any
// input. Calls on one solver must be serialized on its creation device/stream.
// No per-solve allocation, CPU fallback, mixed precision or refinement. Only
// scalar status leaves the device; input and output finiteness is checked there.
// Returns 1 on success (not an iteration count), a named negative failure, or
// GPU_AHF_LINEAR_SINGULAR_BASE-k for a zero U(k,k), using cuSOLVER's 1-based k.
// Failure leaves solution unspecified. Success is not a residual/error bound.
int gpu_ahf_linear_solve(GpuAhfLinearSolver* solver, const double* d_values,
                         const double* d_rhs, double* d_solution) noexcept;

#endif
