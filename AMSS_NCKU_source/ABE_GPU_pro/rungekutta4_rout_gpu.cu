#include "macrodef.h"
#include "rungekutta4_rout.h"
#include "gpu_common.cuh"
#include <cuda_runtime.h>

// =============================================================================
// RK4 array kernel (single variable)
// Reference: rungekutta4_rout.f90 subroutine rungekutta4_rout (lines 94-139)
// =============================================================================
__global__ void rungekutta4_rout_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    const double* __restrict__ f0,
    double* __restrict__ f1,
    double* __restrict__ f_rhs,
    double dT,
    int step
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = i + pex0 * (j + pex1 * k);

    if (step == 0) {
        f1[idx] = f0[idx] + 0.5 * dT * f_rhs[idx];
    } else if (step == 1) {
        double f1_val = f1[idx];
        f_rhs[idx] = f_rhs[idx] + 2.0 * f1_val;
        f1[idx]    = f0[idx] + 0.5 * dT * f1_val;
    } else if (step == 2) {
        double f1_val = f1[idx];
        f_rhs[idx] = f_rhs[idx] + 2.0 * f1_val;
        f1[idx]    = f0[idx] + dT * f1_val;
    } else if (step == 3) {
        f1[idx] = f0[idx] + (1.0 / 6.0) * dT * (f1[idx] + f_rhs[idx]);
    }
}

// =============================================================================
// RK4 fused kernel (multiple variables in one launch)
// Reference: rungekutta4_rout.f90 subroutine rungekutta4_rout (lines 94-139)
// =============================================================================
__global__ void fused_rungekutta4_rout_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    FusedRK4Vars vars,
    double dT,
    int step
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int idx = i + pex0 * (j + pex1 * k);

    if (step == 0) {
        for (int v = 0; v < vars.num_vars; ++v) {
            vars.f1[v][idx] = vars.f0[v][idx] + 0.5 * dT * vars.f_rhs[v][idx];
        }
    }
    else if (step == 1) {
        for (int v = 0; v < vars.num_vars; ++v) {
            double f1_val = vars.f1[v][idx];
            vars.f_rhs[v][idx] = vars.f_rhs[v][idx] + 2.0 * f1_val;
            vars.f1[v][idx]    = vars.f0[v][idx] + 0.5 * dT * f1_val;
        }
    }
    else if (step == 2) {
        for (int v = 0; v < vars.num_vars; ++v) {
            double f1_val = vars.f1[v][idx];
            vars.f_rhs[v][idx] = vars.f_rhs[v][idx] + 2.0 * f1_val;
            vars.f1[v][idx]    = vars.f0[v][idx] + dT * f1_val;
        }
    }
    else if (step == 3) {
        for (int v = 0; v < vars.num_vars; ++v) {
            vars.f1[v][idx] = vars.f0[v][idx] + (1.0 / 6.0) * dT * (vars.f1[v][idx] + vars.f_rhs[v][idx]);
        }
    }
}

void gpu_fused_rungekutta4_rout_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3], double dT,
    FusedRK4Vars vars, int RK4_step
) {
    if (vars.num_vars == 0) return;

    dim3 block(32, 4, 2);
    dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );

    fused_rungekutta4_rout_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], vars, dT, RK4_step
    );
}

void gpu_rungekutta4_rout_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3], double dT,
    const double* d_f0, double* d_f1, double* d_f_rhs, int RK4_step
) {
    dim3 block(32, 4, 2);
    dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );

    rungekutta4_rout_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], d_f0, d_f1, d_f_rhs, dT, RK4_step
    );
}

// =============================================================================
// ICN (Implicit Crank-Nicolson) array kernel
// Reference: rungekutta4_rout.f90 subroutine icn_rout (lines 200-228)
// =============================================================================
__global__ void icn_rout_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    const double* __restrict__ f0,
    double* __restrict__ f1,
    double* __restrict__ f_rhs,
    double dT,
    int step
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = i + pex0 * (j + pex1 * k);

    if (step == 0) {
        f1[idx] = f0[idx] + dT * f_rhs[idx];
    } else {
        f1[idx] = f0[idx] + 0.5 * dT * (f1[idx] + f_rhs[idx]);
    }
}

void gpu_icn_rout_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3], double dT,
    const double* d_f0, double* d_f1, double* d_f_rhs, int step
) {
    dim3 block(32, 4, 2);
    dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );

    icn_rout_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], d_f0, d_f1, d_f_rhs, dT, step
    );
}

// =============================================================================
// Euler array kernel
// Reference: rungekutta4_rout.f90 subroutine euler_rout (lines 230-246)
// =============================================================================
__global__ void euler_rout_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    const double* __restrict__ f0,
    const double* __restrict__ f_rhs,
    double dT,
    double* __restrict__ f1
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = i + pex0 * (j + pex1 * k);

    f1[idx] = f0[idx] + dT * f_rhs[idx];
}

void gpu_euler_rout_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3], double dT,
    const double* d_f0, double* d_f1, const double* d_f_rhs
) {
    dim3 block(32, 4, 2);
    dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );

    euler_rout_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], d_f0, d_f_rhs, dT, d_f1
    );
}
