#include "macrodef.h"
#include "enforce_algebra.h"

#include <cuda_runtime.h>
#include <math.h>

#include "gpu_manager.h"

// =============================================================================
// enforce_ga: enforce det(gamma~) = 1 first, then tr(A~) = 0
// Reference: enforce_algebra.f90 subroutine enforce_ga (#if 1 block, lines 72-133)
// =============================================================================
__global__ void enforce_ga_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz,
    double* Ayy, double* Ayz, double* Azz
) {
    constexpr double ONE = 1.0;
    constexpr double F1o3 = 1.0 / 3.0;
    constexpr double TWO = 2.0;

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

    // Enforce det(gamma~) = 1: scale metric by detg^(-1/3)
    double detg = gxx * gyy * gzz + gxyv * gyzv * gxzv + gxzv * gxyv * gyzv
                - gxzv * gyy * gxzv - gxyv * gxyv * gzz - gxx * gyzv * gyzv;

    double scale = ONE / cbrt(detg);

    gxx *= scale; gxyv *= scale; gxzv *= scale;
    gyy *= scale; gyzv *= scale; gzz *= scale;

    dxx[idx] = gxx - ONE;
    dyy[idx] = gyy - ONE;
    dzz[idx] = gzz - ONE;
    gxy[idx] = gxyv;
    gxz[idx] = gxzv;
    gyz[idx] = gyzv;

    // Enforce tr(A~) = 0: compute trace and subtract
    double gupxx =   ( gyy * gzz - gyzv * gyzv );
    double gupxy = - ( gxyv * gzz - gyzv * gxzv );
    double gupxz =   ( gxyv * gyzv - gyy * gxzv );
    double gupyy =   ( gxx * gzz - gxzv * gxzv );
    double gupyz = - ( gxx * gyzv - gxyv * gxzv );
    double gupzz =   ( gxx * gyy - gxyv * gxyv );

    double trA = gupxx * Axx[idx] + gupyy * Ayy[idx] + gupzz * Azz[idx]
                + TWO * (gupxy * Axy[idx] + gupxz * Axz[idx] + gupyz * Ayz[idx]);

    Axx[idx] -= F1o3 * gxx * trA;
    Axy[idx] -= F1o3 * gxyv * trA;
    Axz[idx] -= F1o3 * gxzv * trA;
    Ayy[idx] -= F1o3 * gyy * trA;
    Ayz[idx] -= F1o3 * gyzv * trA;
    Azz[idx] -= F1o3 * gzz * trA;
}

void gpu_enforce_ga_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    double* d_dxx, double* d_gxy, double* d_gxz,
    double* d_dyy, double* d_gyz, double* d_dzz,
    double* d_Axx, double* d_Axy, double* d_Axz,
    double* d_Ayy, double* d_Ayz, double* d_Azz
) {
    dim3 block(32, 4, 2);
    dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );

    enforce_ga_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2],
        pex[0], pex[1],
        d_dxx, d_gxy, d_gxz,
        d_dyy, d_gyz, d_dzz,
        d_Axx, d_Axy, d_Axz,
        d_Ayy, d_Ayz, d_Azz
    );
}

// =============================================================================
// enforce_ag: enforce tr(A~) = 0 first, then det(gamma~) = 1
// Reference: enforce_algebra.f90 subroutine enforce_ag (lines 8-67)
// =============================================================================
__global__ void enforce_ag_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz,
    double* Ayy, double* Ayz, double* Azz
) {
    constexpr double ONE = 1.0;
    constexpr double F1o3 = 1.0 / 3.0;
    constexpr double TWO = 2.0;

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

    // Step 1: Compute detg and inverse metric, then enforce tr(A~) = 0
    double detg = gxx * gyy * gzz + gxyv * gyzv * gxzv + gxzv * gxyv * gyzv
                - gxzv * gyy * gxzv - gxyv * gxyv * gzz - gxx * gyzv * gyzv;

    double gupxx =   ( gyy * gzz - gyzv * gyzv ) / detg;
    double gupxy = - ( gxyv * gzz - gyzv * gxzv ) / detg;
    double gupxz =   ( gxyv * gyzv - gyy * gxzv ) / detg;
    double gupyy =   ( gxx * gzz - gxzv * gxzv ) / detg;
    double gupyz = - ( gxx * gyzv - gxyv * gxzv ) / detg;
    double gupzz =   ( gxx * gyy - gxyv * gxyv ) / detg;

    double trA = gupxx * Axx[idx] + gupyy * Ayy[idx] + gupzz * Azz[idx]
                + TWO * (gupxy * Axy[idx] + gupxz * Axz[idx] + gupyz * Ayz[idx]);

    Axx[idx] -= F1o3 * gxx * trA;
    Axy[idx] -= F1o3 * gxyv * trA;
    Axz[idx] -= F1o3 * gxzv * trA;
    Ayy[idx] -= F1o3 * gyy * trA;
    Ayz[idx] -= F1o3 * gyzv * trA;
    Azz[idx] -= F1o3 * gzz * trA;

    // Step 2: Enforce det(gamma~) = 1
    double scale = ONE / cbrt(detg);

    gxx *= scale; gxyv *= scale; gxzv *= scale;
    gyy *= scale; gyzv *= scale; gzz *= scale;

    dxx[idx] = gxx - ONE;
    dyy[idx] = gyy - ONE;
    dzz[idx] = gzz - ONE;
    gxy[idx] = gxyv;
    gxz[idx] = gxzv;
    gyz[idx] = gyzv;
}

void gpu_enforce_ag_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    double* d_dxx, double* d_gxy, double* d_gxz,
    double* d_dyy, double* d_gyz, double* d_dzz,
    double* d_Axx, double* d_Axy, double* d_Axz,
    double* d_Ayy, double* d_Ayz, double* d_Azz
) {
    dim3 block(32, 4, 2);
    dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );

    enforce_ag_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2],
        pex[0], pex[1],
        d_dxx, d_gxy, d_gxz,
        d_dyy, d_gyz, d_dzz,
        d_Axx, d_Axy, d_Axz,
        d_Ayy, d_Ayz, d_Azz
    );
}
