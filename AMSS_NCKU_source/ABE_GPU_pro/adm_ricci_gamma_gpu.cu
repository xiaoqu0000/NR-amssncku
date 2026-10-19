#include "macrodef.h"
#include "adm_ricci_gamma_gpu.cuh"
#include "gpu_common.cuh"
#include "macrodef.fh"

#include <cuda_runtime.h>

// ============================================================
// Standalone kernel: compute Ricci tensor and Christoffel
// symbols for ADM variables at every grid point.
// ============================================================
template<int GW, bool Shell>
__global__ void adm_ricci_gamma_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ crho,
    const double* __restrict__ sigma,
    const double* __restrict__ Rcoord,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    double* __restrict__ Gamxxx, double* __restrict__ Gamxxy, double* __restrict__ Gamxxz,
    double* __restrict__ Gamxyy, double* __restrict__ Gamxyz, double* __restrict__ Gamxzz,
    double* __restrict__ Gamyxx, double* __restrict__ Gamyxy, double* __restrict__ Gamyxz,
    double* __restrict__ Gamyyy, double* __restrict__ Gamyyz, double* __restrict__ Gamyzz,
    double* __restrict__ Gamzxx, double* __restrict__ Gamzxy, double* __restrict__ Gamzxz,
    double* __restrict__ Gamzyy, double* __restrict__ Gamzyz, double* __restrict__ Gamzzz,
    double* __restrict__ Rxx, double* __restrict__ Rxy, double* __restrict__ Rxz,
    double* __restrict__ Ryy, double* __restrict__ Ryz, double* __restrict__ Rzz,
    int symmetry, int sst)
{
    const int n = ex0 * ex1 * ex2;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;

    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);

    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const size_t id = static_cast<size_t>(i) + static_cast<size_t>(pex0) *
        (static_cast<size_t>(j) + static_cast<size_t>(pex1) * static_cast<size_t>(k));

    double l_Gamxxx, l_Gamxxy, l_Gamxxz;
    double l_Gamxyy, l_Gamxyz, l_Gamxzz;
    double l_Gamyxx, l_Gamyxy, l_Gamyxz;
    double l_Gamyyy, l_Gamyyz, l_Gamyzz;
    double l_Gamzxx, l_Gamzxy, l_Gamzxz;
    double l_Gamzyy, l_Gamzyz, l_Gamzzz;
    double l_Rxx, l_Rxy, l_Rxz, l_Ryy, l_Ryz, l_Rzz;

    if constexpr (Shell) {
        d_adm_ricci_gamma_point_ss<GW>(
            ex, pex, i, j, k, id,
            crho, sigma, Rcoord,
            dxx, gxy, gxz, dyy, gyz, dzz,
            l_Gamxxx, l_Gamxxy, l_Gamxxz,
            l_Gamxyy, l_Gamxyz, l_Gamxzz,
            l_Gamyxx, l_Gamyxy, l_Gamyxz,
            l_Gamyyy, l_Gamyyz, l_Gamyzz,
            l_Gamzxx, l_Gamzxy, l_Gamzxz,
            l_Gamzyy, l_Gamzyz, l_Gamzzz,
            l_Rxx, l_Rxy, l_Rxz, l_Ryy, l_Ryz, l_Rzz,
            symmetry, sst);
    } else {
        d_adm_ricci_gamma_point<GW>(
            ex, pex, i, j, k, id,
            X, Y, Z,
            dxx, gxy, gxz, dyy, gyz, dzz,
            l_Gamxxx, l_Gamxxy, l_Gamxxz,
            l_Gamxyy, l_Gamxyz, l_Gamxzz,
            l_Gamyxx, l_Gamyxy, l_Gamyxz,
            l_Gamyyy, l_Gamyyz, l_Gamyzz,
            l_Gamzxx, l_Gamzxy, l_Gamzxz,
            l_Gamzyy, l_Gamzyz, l_Gamzzz,
            l_Rxx, l_Rxy, l_Rxz, l_Ryy, l_Ryz, l_Rzz,
            symmetry);
    }

    // Store results
    Gamxxx[id] = l_Gamxxx; Gamxxy[id] = l_Gamxxy; Gamxxz[id] = l_Gamxxz;
    Gamxyy[id] = l_Gamxyy; Gamxyz[id] = l_Gamxyz; Gamxzz[id] = l_Gamxzz;
    Gamyxx[id] = l_Gamyxx; Gamyxy[id] = l_Gamyxy; Gamyxz[id] = l_Gamyxz;
    Gamyyy[id] = l_Gamyyy; Gamyyz[id] = l_Gamyyz; Gamyzz[id] = l_Gamyzz;
    Gamzxx[id] = l_Gamzxx; Gamzxy[id] = l_Gamzxy; Gamzxz[id] = l_Gamzxz;
    Gamzyy[id] = l_Gamzyy; Gamzyz[id] = l_Gamzyz; Gamzzz[id] = l_Gamzzz;
    Rxx[id] = l_Rxx; Rxy[id] = l_Rxy; Rxz[id] = l_Rxz;
    Ryy[id] = l_Ryy; Ryz[id] = l_Ryz; Rzz[id] = l_Rzz;
}

// ============================================================
// Host launch helpers
// ============================================================
template<bool Shell>
void launch_adm_ricci_gamma(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* crho, const double* sigma, const double* Rcoord,
    const double* X, const double* Y, const double* Z,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    int symmetry, int sst)
{
    const int n = ex[0] * ex[1] * ex[2];
    const int block = 256;
    const int grid = (n + block - 1) / block;
    adm_ricci_gamma_kernel<ghost_width, Shell><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        crho, sigma, Rcoord, X, Y, Z,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        symmetry, sst);
}

// Explicit template instantiation macro
#ifndef INSTANTIATE_SHELL_VARIANTS
#define INSTANTIATE_ADM_RICCI_GAMMA_LAUNCH(CartShell) \
template void launch_adm_ricci_gamma<CartShell>( \
    cudaStream_t, int[3], int[3], \
    const double*, const double*, const double*, \
    const double*, const double*, const double*, \
    const double*, const double*, const double*, \
    const double*, const double*, const double*, \
    double*, double*, double*, double*, double*, double*, \
    double*, double*, double*, double*, double*, double*, \
    double*, double*, double*, double*, double*, double*, \
    double*, double*, double*, double*, double*, double*, \
    int, int);
#else
#define INSTANTIATE_ADM_RICCI_GAMMA_LAUNCH(CartShell)
#endif

INSTANTIATE_ADM_RICCI_GAMMA_LAUNCH(false)
INSTANTIATE_ADM_RICCI_GAMMA_LAUNCH(true)

// C-linkage wrappers
extern "C" void gpu_adm_ricci_gamma_launch(
    cudaStream_t stream, int ex[3], int pex[3],
    double* X, double* Y, double* Z,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    int symmetry)
{
    launch_adm_ricci_gamma<false>(stream, ex, pex,
        nullptr, nullptr, nullptr, X, Y, Z,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        symmetry, 0);
}

extern "C" void gpu_adm_ricci_gamma_ss_launch(
    cudaStream_t stream, int ex[3], int pex[3],
    double* crho, double* sigma, double* Rcoord,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    int symmetry, int sst)
{
    launch_adm_ricci_gamma<true>(stream, ex, pex,
        crho, sigma, Rcoord, nullptr, nullptr, nullptr,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        symmetry, sst);
}
