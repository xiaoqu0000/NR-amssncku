#include "macrodef.h"
#include "bssn_constraint_gpu.cuh"
#include "gpu_common.cuh"

#include <cuda_runtime.h>
#include <math_constants.h>

using namespace bssn_constraint_gpu;

namespace {

constexpr double F8 = 8.0;
constexpr double THREE = 3.0;
constexpr double ONE_THIRD = 1.0 / 3.0;

__device__ __forceinline__ double d_get4ricciscalar(
    double chi, double trK, double rho,
    double dxx, double gxy, double gxz, double dyy, double gyz, double dzz,
    double Axx, double Axy, double Axz, double Ayy, double Ayz, double Azz,
    double Rxx, double Rxy, double Rxz, double Ryy, double Ryz, double Rzz,
    double Sxx, double Sxy, double Sxz, double Syy, double Syz, double Szz)
{
    // a) Physical metric: g_ij = d_ij + delta_ij
    const double l_gxx = dxx + ONE;
    const double l_gyy = dyy + ONE;
    const double l_gzz = dzz + ONE;
    // gxy, gxz, gyz are the physical off-diagonal components as-is

    // b) chipn1 = chi + 1
    const double chipn1 = chi + ONE;

    // c) Invert the physical 3-metric
    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    invert_metric(l_gxx, gxy, gxz, l_gyy, gyz, l_gzz,
                  gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);

    // d) trR = gup^{ij} R_{ij}
    const double trR = gupxx * Rxx + gupyy * Ryy + gupzz * Rzz
                     + TWO * (gupxy * Rxy + gupxz * Rxz + gupyz * Ryz);

    // e) trS = gup^{ij} S_{ij}
    const double trS = gupxx * Sxx + gupyy * Syy + gupzz * Szz
                     + TWO * (gupxy * Sxy + gupxz * Sxz + gupyz * Syz);

    // f) A_ij A^{ij} — uses input A_ij directly (matches Fortran fourdcurvature.f90)
    const double Asq = contract_symmetric_square(
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz);

    // g) 4D Ricci scalar:
    //    RR = A_{ij}A^{ij} - (2/3)trK^2 - (1+chi)trR - 8*pi*(3*rho - (1+chi)trS)
    const double RR = Asq
                    - F2O3 * trK * trK
                    - chipn1 * trR
                    - F8 * CUDART_PI * (THREE * rho - chipn1 * trS);

    return RR;
}

} // namespace

__global__ void fourdcurvature_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ rho,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Axx, const double* __restrict__ Axy, const double* __restrict__ Axz,
    const double* __restrict__ Ayy, const double* __restrict__ Ayz, const double* __restrict__ Azz,
    const double* __restrict__ Rxx, const double* __restrict__ Rxy, const double* __restrict__ Rxz,
    const double* __restrict__ Ryy, const double* __restrict__ Ryz, const double* __restrict__ Rzz,
    const double* __restrict__ Sxx, const double* __restrict__ Sxy, const double* __restrict__ Sxz,
    const double* __restrict__ Syy, const double* __restrict__ Syz, const double* __restrict__ Szz,
    double* __restrict__ RR)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int n = ex0 * ex1 * ex2;
    if (tid >= n) return;

    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);

    const int pex[3] = {pex0, pex1, pex2};
    const size_t id = idx3(i, j, k, pex);

    RR[id] = d_get4ricciscalar(
        chi[id], trK[id], rho[id],
        dxx[id], gxy[id], gxz[id], dyy[id], gyz[id], dzz[id],
        Axx[id], Axy[id], Axz[id], Ayy[id], Ayz[id], Azz[id],
        Rxx[id], Rxy[id], Rxz[id], Ryy[id], Ryz[id], Rzz[id],
        Sxx[id], Sxy[id], Sxz[id], Syy[id], Syz[id], Szz[id]);
}

extern "C" void gpu_get4ricciscalar_launch(
    cudaStream_t stream, int ex[3], int pex[3],
    double* chi, double* trK, double* rho,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* RR)
{
    const int n = ex[0] * ex[1] * ex[2];
    const int block = 256;
    const int grid = (n + block - 1) / block;

    fourdcurvature_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        chi, trK, rho,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        Sxx, Sxy, Sxz, Syy, Syz, Szz,
        RR);
}
