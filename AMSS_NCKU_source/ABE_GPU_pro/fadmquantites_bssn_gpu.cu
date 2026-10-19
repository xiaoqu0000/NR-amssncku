#include "macrodef.h"
#include "fadmquantites_bssn.h"

#include <cuda_runtime.h>
#include <math.h>

#include "fmisc.h"
#include "gpu_manager.h"
#include "derivatives.h"
#include "diff_new_gpu.cuh"

// =============================================================================
// admmass_bssn: compute ADM mass density from BSSN variables
// Reference: fadmquantites_bssn.f90 subroutine admmass_bssn (lines 5-65)
// =============================================================================
__global__ void admmass_bssn_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ d_X, const double* __restrict__ d_Y, const double* __restrict__ d_Z,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Axx, const double* __restrict__ Axy, const double* __restrict__ Axz,
    const double* __restrict__ Ayy, const double* __restrict__ Ayz, const double* __restrict__ Azz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    double* __restrict__ massx, double* __restrict__ massy, double* __restrict__ massz,
    int symmetry
) {
    int i0 = blockIdx.x * blockDim.x + threadIdx.x;
    int j0 = blockIdx.y * blockDim.y + threadIdx.y;
    int k0 = blockIdx.z * blockDim.z + threadIdx.z;

    if (i0 >= ex0 || j0 >= ex1 || k0 >= ex2) return;

    int idx = (k0 * pex1 + j0) * pex0 + i0;

    double F1o2pi = 1.0 / (2.0 * M_PI);
    double F1o8 = 1.0 / 8.0;

    int ext[3]  = {ex0, ex1, ex2};
    int pext[3] = {pex0, pex1, pex2};

    double chix = 0.0, chiy = 0.0, chiz = 0.0;
    d_fderivs_point(
        ext, pext, chi,
        &chix, &chiy, &chiz,
        d_X, d_Y, d_Z,
        1.0, 1.0, 1.0,
        symmetry, 0,
        i0, j0, k0
    );

    double gxx_v = dxx[idx] + 1.0;
    double gyy_v = dyy[idx] + 1.0;
    double gzz_v = dzz[idx] + 1.0;
    double gxy_v = gxy[idx];
    double gxz_v = gxz[idx];
    double gyz_v = gyz[idx];

    double det = gxx_v * gyy_v * gzz_v + gxy_v * gyz_v * gxz_v + gxz_v * gxy_v * gyz_v
               - gxz_v * gyy_v * gxz_v - gxy_v * gxy_v * gzz_v - gxx_v * gyz_v * gyz_v;

    double gupxx = (gyy_v * gzz_v - gyz_v * gyz_v) / det;
    double gupxy = -(gxy_v * gzz_v - gyz_v * gxz_v) / det;
    double gupxz = (gxy_v * gyz_v - gyy_v * gxz_v) / det;
    double gupyy = (gxx_v * gzz_v - gxz_v * gxz_v) / det;
    double gupyz = -(gxx_v * gyz_v - gxy_v * gxz_v) / det;
    double gupzz = (gxx_v * gyy_v - gxy_v * gxy_v) / det;

    double f = 0.25 * pow(chi[idx] + 1.0, -1.25);

    massx[idx] = (F1o8 * Gamx[idx] + f * (gupxx * chix + gupxy * chiy + gupxz * chiz)) * F1o2pi;
    massy[idx] = (F1o8 * Gamy[idx] + f * (gupxy * chix + gupyy * chiy + gupyz * chiz)) * F1o2pi;
    massz[idx] = (F1o8 * Gamz[idx] + f * (gupxz * chix + gupyz * chiy + gupzz * chiz)) * F1o2pi;
}

void gpu_admmass_bssn_launch(
    cudaStream_t stream, const int ext[3], const int pext[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* chi, const double* trK,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    const double* Axx, const double* Axy, const double* Axz,
    const double* Ayy, const double* Ayz, const double* Azz,
    const double* Gamx, const double* Gamy, const double* Gamz,
    double* massx, double* massy, double* massz,
    int symmetry
) {
    dim3 block(8, 8, 8);
    dim3 grid((ext[0] + block.x - 1) / block.x,
              (ext[1] + block.y - 1) / block.y,
              (ext[2] + block.z - 1) / block.z);

    admmass_bssn_kernel<<<grid, block, 0, stream>>>(
        ext[0], ext[1], ext[2],
        pext[0], pext[1], pext[2],
        d_X, d_Y, d_Z,
        chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz,
        massx, massy, massz,
        symmetry
    );
}

// =============================================================================
// admmomentum_bssn: compute ADM momentum density from BSSN variables
// Reference: fadmquantites_bssn.f90 subroutine admmomentum_bssn (lines 69-115)
// =============================================================================
__global__ void admmomentum_bssn_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Axx, const double* __restrict__ Axy, const double* __restrict__ Axz,
    const double* __restrict__ Ayy, const double* __restrict__ Ayz, const double* __restrict__ Azz,
    double* __restrict__ pxx, double* __restrict__ pxy, double* __restrict__ pxz,
    double* __restrict__ pyy, double* __restrict__ pyz, double* __restrict__ pzz
) {
    constexpr double ONE = 1.0;
    constexpr double F1o3 = 1.0 / 3.0;
    constexpr double F1o8pi = 1.0 / (8.0 * M_PI);

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int idx = i + pex0 * (j + pex1 * k);

    // Build conformal metric
    double gxx = dxx[idx] + ONE;
    double gyy = dyy[idx] + ONE;
    double gzz = dzz[idx] + ONE;
    double gxyv = gxy[idx];
    double gxzv = gxz[idx];
    double gyzv = gyz[idx];

    // chim4 = 1/(chi+1)^4, used to reconstruct Kij from Aij
    double inv_chi_plus_one = 1.0 / (chi[idx] + ONE);
    double chim4 = inv_chi_plus_one * inv_chi_plus_one;
    chim4 *= chim4;
    double trK_val = trK[idx];

    // Reconstruct extrinsic curvature: Kij = chim4 * (Aij + 1/3 * gij * trK)
    double Kxx = chim4 * (Axx[idx] + F1o3 * gxx * trK_val);
    double Kxy = chim4 * (Axy[idx] + F1o3 * gxyv * trK_val);
    double Kxz = chim4 * (Axz[idx] + F1o3 * gxzv * trK_val);
    double Kyy = chim4 * (Ayy[idx] + F1o3 * gyy * trK_val);
    double Kyz = chim4 * (Ayz[idx] + F1o3 * gyzv * trK_val);
    double Kzz = chim4 * (Azz[idx] + F1o3 * gzz * trK_val);

    // Momentum density: pij = (Kij - delta_ij * trK) / (8*PI)
    pxx[idx] = (Kxx - trK_val) * F1o8pi;
    pxy[idx] = Kxy * F1o8pi;
    pxz[idx] = Kxz * F1o8pi;
    pyy[idx] = (Kyy - trK_val) * F1o8pi;
    pyz[idx] = Kyz * F1o8pi;
    pzz[idx] = (Kzz - trK_val) * F1o8pi;
}

void gpu_admmomentum_bssn_launch(
    cudaStream_t stream, const int ext[3], const int pext[3],
    const double* chi, const double* trK,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    const double* Axx, const double* Axy, const double* Axz,
    const double* Ayy, const double* Ayz, const double* Azz,
    double* pxx, double* pxy, double* pxz,
    double* pyy, double* pyz, double* pzz
) {
    dim3 block(8, 8, 8);
    dim3 grid((ext[0] + block.x - 1) / block.x,
              (ext[1] + block.y - 1) / block.y,
              (ext[2] + block.z - 1) / block.z);

    admmomentum_bssn_kernel<<<grid, block, 0, stream>>>(
        ext[0], ext[1], ext[2],
        pext[0], pext[1],
        chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        pxx, pxy, pxz, pyy, pyz, pzz
    );
}
