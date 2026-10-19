#include "macrodef.h"
// ============================================================
// diff_new_sh_gpu.cu
// CUDA kernels and host launch functions for shell coordinate
// derivatives and KO dissipation.
// ============================================================

#include "diff_new_sh_gpu.cuh"
#include "kodiss_sh_gpu.cuh"
#include "gpu_manager.h"
#include "derivatives.h"
#include "macrodef.fh"

#include <cuda_runtime.h>

// ============================================================
// IDX3D_PAD: 3D padded array index
// ============================================================
#ifndef IDX3D_PAD
#define IDX3D_PAD(i, j, k, pnx, pny) ((i) + (pnx) * ((j) + (pny) * (k)))
#endif

// ============================================================
// 1. fderivs_sh_kernel<GW> -- all first derivatives (shell)
// ============================================================
template<int GW>
__global__ void fderivs_sh_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ f,
    double* __restrict__ fx, double* __restrict__ fy, double* __restrict__ fz,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int ex[3] = {ex0, ex1, ex2};
    int pex[3] = {pex0, pex1, pex2};

    double lfx, lfy, lfz;
    d_fderivs_sh_point_gw<GW>(ex, pex, f, &lfx, &lfy, &lfz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);

    int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    fx[idx] = lfx;
    fy[idx] = lfy;
    fz[idx] = lfz;
}

// ============================================================
// 2. fdderivs_sh_kernel<GW> -- all second derivatives (shell)
// ============================================================
template<int GW>
__global__ void fdderivs_sh_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ f,
    double* __restrict__ fxx, double* __restrict__ fxy, double* __restrict__ fxz,
    double* __restrict__ fyy, double* __restrict__ fyz, double* __restrict__ fzz,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int ex[3] = {ex0, ex1, ex2};
    int pex[3] = {pex0, pex1, pex2};

    double lfxx, lfxy, lfxz, lfyy, lfyz, lfzz;
    d_fdderivs_sh_point_gw<GW>(ex, pex, f,
        &lfxx, &lfxy, &lfxz, &lfyy, &lfyz, &lfzz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);

    int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    fxx[idx] = lfxx;
    fxy[idx] = lfxy;
    fxz[idx] = lfxz;
    fyy[idx] = lfyy;
    fyz[idx] = lfyz;
    fzz[idx] = lfzz;
}

// ============================================================
// 3. kodis_sh_kernel<GW> -- main KO dissipation (all 3 dirs)
// ============================================================
template<int GW>
__global__ void kodis_sh_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ f,
    double* __restrict__ f_rhs,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int ex[3] = {ex0, ex1, ex2};
    int pex[3] = {pex0, pex1, pex2};

    int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    f_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, f, X, Y, Z,
        SYM1, SYM2, SYM3, symmetry, eps, sst, i, j, k);
}

// ============================================================
// 4. kodis_shor_kernel<GW> -- radial-only (z direction) KO
// ============================================================
template<int GW>
__global__ void kodis_shor_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ f,
    double* __restrict__ f_rhs,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int ex[3] = {ex0, ex1, ex2};
    int pex[3] = {pex0, pex1, pex2};

    int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    f_rhs[idx] += d_kodis_shor_point<GW>(ex, pex, f, X, Y, Z,
        SYM1, SYM2, SYM3, symmetry, eps, sst, i, j, k);
}

// ============================================================
// 5. kodis_shcr_kernel<GW> -- constant-r (x,y only) KO
// ============================================================
template<int GW>
__global__ void kodis_shcr_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ f,
    double* __restrict__ f_rhs,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int ex[3] = {ex0, ex1, ex2};
    int pex[3] = {pex0, pex1, pex2};

    int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    f_rhs[idx] += d_kodis_shcr_point<GW>(ex, pex, f, X, Y, Z,
        SYM1, SYM2, SYM3, symmetry, eps, sst, i, j, k);
}

// ============================================================
// 6. kodis_sh_new_kernel<GW> -- near-boundary KO (GW=2 only)
// ============================================================
template<int GW>
__global__ void kodis_sh_new_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ f,
    double* __restrict__ f_rhs,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int ex[3] = {ex0, ex1, ex2};
    int pex[3] = {pex0, pex1, pex2};

    int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    f_rhs[idx] += d_kodis_sh_new_point<GW>(ex, pex, f, X, Y, Z,
        SYM1, SYM2, SYM3, symmetry, eps, sst, i, j, k);
}

// ============================================================
// Host launch functions
// ============================================================

// ----- fderivs_sh -----
template<int GW>
void gpu_fderivs_sh_launch_impl(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* fx, double* fy, double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst)
{
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    fderivs_sh_kernel<GW><<<gridDim, blockDim, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        f, fx, fy, fz, X, Y, Z, SYM1, SYM2, SYM3,
        symmetry, onoff, sst);
}

void gpu_fderivs_sh_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* fx, double* fy, double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst)
{
    gpu_fderivs_sh_launch_impl<ghost_width>(
        stream, ex, pex, f, fx, fy, fz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst);
}

// ----- fdderivs_sh -----
template<int GW>
void gpu_fdderivs_sh_launch_impl(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst)
{
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    fdderivs_sh_kernel<GW><<<gridDim, blockDim, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        f, fxx, fxy, fxz, fyy, fyz, fzz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst);
}

void gpu_fdderivs_sh_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst)
{
    gpu_fdderivs_sh_launch_impl<ghost_width>(
        stream, ex, pex, f, fxx, fxy, fxz, fyy, fyz, fzz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst);
}

// ----- kodis_sh (main KO, all 3 directions) -----
template<int GW>
void gpu_kodis_sh_launch_impl(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* f_rhs,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    kodis_sh_kernel<GW><<<gridDim, blockDim, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        f, f_rhs, X, Y, Z, SYM1, SYM2, SYM3,
        symmetry, eps, sst);
}

void gpu_kodis_sh_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* f_rhs,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    gpu_kodis_sh_launch_impl<ghost_width>(
        stream, ex, pex, f, f_rhs,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, eps, sst);
}

// ----- kodis_shor (radial-only KO) -----
template<int GW>
void gpu_kodis_shor_launch_impl(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* f_rhs,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    kodis_shor_kernel<GW><<<gridDim, blockDim, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        f, f_rhs, X, Y, Z, SYM1, SYM2, SYM3,
        symmetry, eps, sst);
}

void gpu_kodis_shor_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* f_rhs,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    gpu_kodis_shor_launch_impl<ghost_width>(
        stream, ex, pex, f, f_rhs,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, eps, sst);
}

// ----- kodis_shcr (constant-r KO, x,y only) -----
template<int GW>
void gpu_kodis_shcr_launch_impl(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* f_rhs,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    kodis_shcr_kernel<GW><<<gridDim, blockDim, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        f, f_rhs, X, Y, Z, SYM1, SYM2, SYM3,
        symmetry, eps, sst);
}

void gpu_kodis_shcr_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* f_rhs,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    gpu_kodis_shcr_launch_impl<ghost_width>(
        stream, ex, pex, f, f_rhs,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, eps, sst);
}

// ----- kodis_sh_new (near-boundary KO, GW=2 only) -----
template<int GW>
void gpu_kodis_sh_new_launch_impl(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* f_rhs,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    kodis_sh_new_kernel<GW><<<gridDim, blockDim, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        f, f_rhs, X, Y, Z, SYM1, SYM2, SYM3,
        symmetry, eps, sst);
}

void gpu_kodis_sh_new_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* f,
    double* f_rhs,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst)
{
    gpu_kodis_sh_new_launch_impl<ghost_width>(
        stream, ex, pex, f, f_rhs,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, eps, sst);
}
