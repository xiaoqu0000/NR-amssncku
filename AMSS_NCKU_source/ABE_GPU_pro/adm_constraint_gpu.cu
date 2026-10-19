#include "macrodef.h"
#include "adm_ricci_gamma_gpu.cuh"
#include "bssn_constraint_gpu.cuh"
#include "diff_new_gpu.cuh"
#include "diff_new_sh_gpu.cuh"
#include "gpu_common.cuh"
#include "gpu_manager.h"
#include "macrodef.fh"

#include <cuda_runtime.h>
#include <math_constants.h>

using namespace bssn_constraint_gpu;

namespace {

// ============================================================
// Per-point ADM constraint computation
// ============================================================
template<int GW, bool Shell>
__device__ __forceinline__ void d_adm_constraint_point(
    const int ex[3], const int pex[3],
    int i, int j, int k, size_t idx,
    const double* X, const double* Y, const double* Z,
    const double* crho, const double* sigma, const double* Rcoord,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    const double* Kxx, const double* Kxy, const double* Kxz,
    const double* Kyy, const double* Kyz, const double* Kzz,
    const double* trK_arr,
    const double* rho, const double* Sx, const double* Sy, const double* Sz,
    double& ham_Res, double& movx_Res, double& movy_Res, double& movz_Res,
    int symmetry, int sst)
{
    // ----------------------------------------------------------
    // Step 1: Compute Ricci tensor R_ij and Christoffel Gamma^i_jk
    // ----------------------------------------------------------
    double Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz;
    double Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz;
    double Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz;
    double Rxx, Rxy, Rxz, Ryy, Ryz, Rzz;

    if constexpr (Shell) {
        d_adm_ricci_gamma_point_ss<GW>(
            ex, pex, i, j, k, idx,
            crho, sigma, Rcoord,
            dxx, gxy, gxz, dyy, gyz, dzz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
            symmetry, sst);
    } else {
        d_adm_ricci_gamma_point<GW>(
            ex, pex, i, j, k, idx, X, Y, Z,
            dxx, gxy, gxz, dyy, gyz, dzz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
            symmetry);
    }

    // ----------------------------------------------------------
    // Step 2: Build physical metric g_ij = d_ij + delta_ij, invert
    // ----------------------------------------------------------
    const double l_gxx = dxx[idx] + ONE;
    const double l_gxy = gxy[idx];
    const double l_gxz = gxz[idx];
    const double l_gyy = dyy[idx] + ONE;
    const double l_gyz = gyz[idx];
    const double l_gzz = dzz[idx] + ONE;

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    invert_metric(l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz,
                  gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);

    // ----------------------------------------------------------
    // Step 3: Compute trK = gup^ij K_ij
    // ----------------------------------------------------------
    const double l_Kxx = Kxx[idx];
    const double l_Kxy = Kxy[idx];
    const double l_Kxz = Kxz[idx];
    const double l_Kyy = Kyy[idx];
    const double l_Kyz = Kyz[idx];
    const double l_Kzz = Kzz[idx];

    const double trK_val = gupxx * l_Kxx + gupyy * l_Kyy + gupzz * l_Kzz
                         + TWO * (gupxy * l_Kxy + gupxz * l_Kxz + gupyz * l_Kyz);

    // ----------------------------------------------------------
    // Step 4: Hamiltonian constraint
    //   ham = trR + trK^2 - K_ij K^ij - 16*pi*rho
    // ----------------------------------------------------------
    const double trR = gupxx * Rxx + gupyy * Ryy + gupzz * Rzz
                     + TWO * (gupxy * Rxy + gupxz * Rxz + gupyz * Ryz);

    ham_Res = trR
            + trK_val * trK_val
            - contract_symmetric_square(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
                                        l_Kxx, l_Kxy, l_Kxz, l_Kyy, l_Kyz, l_Kzz)
            - F16 * CUDART_PI * rho[idx];

    // ----------------------------------------------------------
    // Step 5: First derivatives of K_ij (7 calls, 6 components + trK)
    // ----------------------------------------------------------
    // Naming: Kijk = partial_k K_ij  (e.g. Kxxx = d_x K_xx, Kxyy = d_y K_xy)
    // Symmetries:
    //   Kxx: SYM, SYM, SYM      Kxy: ANTI, ANTI, SYM
    //   Kxz: ANTI, SYM, ANTI    Kyy: SYM, SYM, SYM
    //   Kyz: SYM, ANTI, ANTI    Kzz: SYM, SYM, SYM
    //   trK: SYM, SYM, SYM

    double Kxxx, Kxxy, Kxxz;
    double Kxyx, Kxyy, Kxyz;
    double Kxzx, Kxzy, Kxzz;
    double Kyyx, Kyyy, Kyyz;
    double Kyzx, Kyzy, Kyzz;
    double Kzzx, Kzzy, Kzzz;
    double trK_x, trK_y, trK_z;

    if constexpr (Shell) {
        d_fderivs_sh_point_gw<GW>(ex, pex, Kxx, &Kxxx, &Kxxy, &Kxxz,
            crho, sigma, Rcoord, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);
        d_fderivs_sh_point_gw<GW>(ex, pex, Kxy, &Kxyx, &Kxyy, &Kxyz,
            crho, sigma, Rcoord, ANTI, ANTI, SYM, symmetry, 0, sst, i, j, k);
        d_fderivs_sh_point_gw<GW>(ex, pex, Kxz, &Kxzx, &Kxzy, &Kxzz,
            crho, sigma, Rcoord, ANTI, SYM, ANTI, symmetry, 0, sst, i, j, k);
        d_fderivs_sh_point_gw<GW>(ex, pex, Kyy, &Kyyx, &Kyyy, &Kyyz,
            crho, sigma, Rcoord, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);
        d_fderivs_sh_point_gw<GW>(ex, pex, Kyz, &Kyzx, &Kyzy, &Kyzz,
            crho, sigma, Rcoord, SYM, ANTI, ANTI, symmetry, 0, sst, i, j, k);
        d_fderivs_sh_point_gw<GW>(ex, pex, Kzz, &Kzzx, &Kzzy, &Kzzz,
            crho, sigma, Rcoord, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);
        d_fderivs_sh_point_gw<GW>(ex, pex, trK_arr, &trK_x, &trK_y, &trK_z,
            crho, sigma, Rcoord, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);
    } else {
        d_fderivs_point_gw<GW>(ex, pex, Kxx, &Kxxx, &Kxxy, &Kxxz,
            X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
        d_fderivs_point_gw<GW>(ex, pex, Kxy, &Kxyx, &Kxyy, &Kxyz,
            X, Y, Z, ANTI, ANTI, SYM, symmetry, 0, i, j, k);
        d_fderivs_point_gw<GW>(ex, pex, Kxz, &Kxzx, &Kxzy, &Kxzz,
            X, Y, Z, ANTI, SYM, ANTI, symmetry, 0, i, j, k);
        d_fderivs_point_gw<GW>(ex, pex, Kyy, &Kyyx, &Kyyy, &Kyyz,
            X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
        d_fderivs_point_gw<GW>(ex, pex, Kyz, &Kyzx, &Kyzy, &Kyzz,
            X, Y, Z, SYM, ANTI, ANTI, symmetry, 0, i, j, k);
        d_fderivs_point_gw<GW>(ex, pex, Kzz, &Kzzx, &Kzzy, &Kzzz,
            X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
        d_fderivs_point_gw<GW>(ex, pex, trK_arr, &trK_x, &trK_y, &trK_z,
            X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
    }

    // ----------------------------------------------------------
    // Steps 6-7: Convert partial to covariant derivatives and contract
    //   D_i K_jk = partial_i K_jk - Gamma^l_ij K_lk - Gamma^l_ik K_jl
    //   mov_Res_j = gup^kl D_k K_jl
    // This is handled by covariant_momentum with chi terms disabled.
    // ----------------------------------------------------------
    covariant_momentum(
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        l_Kxx, l_Kxy, l_Kxz, l_Kyy, l_Kyz, l_Kzz,
        Kxxx, Kxyx, Kxzx, Kyyx, Kyzx, Kzzx,
        Kxxy, Kxyy, Kxzy, Kyyy, Kyzy, Kzzy,
        Kxxz, Kxyz, Kxzz, Kyyz, Kyzz, Kzzz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        0.0, 0.0, 0.0, 0.0,
        movx_Res, movy_Res, movz_Res);

    // ----------------------------------------------------------
    // Step 8: Subtract d_j(trK) and source terms
    //   mov_Res_j -= d_j(trK) + 8*pi*S_j
    // ----------------------------------------------------------
    movx_Res -= trK_x + F8 * CUDART_PI * Sx[idx];
    movy_Res -= trK_y + F8 * CUDART_PI * Sy[idx];
    movz_Res -= trK_z + F8 * CUDART_PI * Sz[idx];
}

// Same as d_adm_constraint_point but also writes Gamma and Ricci outputs.
// Used when the caller (Fortran constraint_adm_ss) needs them as intent(out).
template<int GW, bool Shell>
__device__ __forceinline__ void d_adm_constraint_point_with_outputs(
    const int ex[3], const int pex[3],
    int i, int j, int k, size_t idx,
    const double* X, const double* Y, const double* Z,
    const double* crho, const double* sigma, const double* Rcoord,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    const double* Kxx, const double* Kxy, const double* Kxz,
    const double* Kyy, const double* Kyz, const double* Kzz,
    const double* trK_arr,
    const double* rho, const double* Sx, const double* Sy, const double* Sz,
    double& ham_Res, double& movx_Res, double& movy_Res, double& movz_Res,
    // Gamma outputs (may be null)
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    // Ricci outputs (may be null)
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    int symmetry, int sst)
{
    double l_Gamxxx, l_Gamxxy, l_Gamxxz, l_Gamxyy, l_Gamxyz, l_Gamxzz;
    double l_Gamyxx, l_Gamyxy, l_Gamyxz, l_Gamyyy, l_Gamyyz, l_Gamyzz;
    double l_Gamzxx, l_Gamzxy, l_Gamzxz, l_Gamzyy, l_Gamzyz, l_Gamzzz;
    double l_Rxx, l_Rxy, l_Rxz, l_Ryy, l_Ryz, l_Rzz;

    if constexpr (Shell) {
        d_adm_ricci_gamma_point_ss<GW>(
            ex, pex, i, j, k, idx, crho, sigma, Rcoord,
            dxx, gxy, gxz, dyy, gyz, dzz,
            l_Gamxxx, l_Gamxxy, l_Gamxxz, l_Gamxyy, l_Gamxyz, l_Gamxzz,
            l_Gamyxx, l_Gamyxy, l_Gamyxz, l_Gamyyy, l_Gamyyz, l_Gamyzz,
            l_Gamzxx, l_Gamzxy, l_Gamzxz, l_Gamzyy, l_Gamzyz, l_Gamzzz,
            l_Rxx, l_Rxy, l_Rxz, l_Ryy, l_Ryz, l_Rzz, symmetry, sst);
    } else {
        d_adm_ricci_gamma_point<GW>(
            ex, pex, i, j, k, idx, X, Y, Z,
            dxx, gxy, gxz, dyy, gyz, dzz,
            l_Gamxxx, l_Gamxxy, l_Gamxxz, l_Gamxyy, l_Gamxyz, l_Gamxzz,
            l_Gamyxx, l_Gamyxy, l_Gamyxz, l_Gamyyy, l_Gamyyz, l_Gamyzz,
            l_Gamzxx, l_Gamzxy, l_Gamzxz, l_Gamzyy, l_Gamzyz, l_Gamzzz,
            l_Rxx, l_Rxy, l_Rxz, l_Ryy, l_Ryz, l_Rzz, symmetry);
    }

    // Store Gamma outputs if requested
    if (Gamxxx) {
        *Gamxxx=l_Gamxxx;*Gamxxy=l_Gamxxy;*Gamxxz=l_Gamxxz;
        *Gamxyy=l_Gamxyy;*Gamxyz=l_Gamxyz;*Gamxzz=l_Gamxzz;
        *Gamyxx=l_Gamyxx;*Gamyxy=l_Gamyxy;*Gamyxz=l_Gamyxz;
        *Gamyyy=l_Gamyyy;*Gamyyz=l_Gamyyz;*Gamyzz=l_Gamyzz;
        *Gamzxx=l_Gamzxx;*Gamzxy=l_Gamzxy;*Gamzxz=l_Gamzxz;
        *Gamzyy=l_Gamzyy;*Gamzyz=l_Gamzyz;*Gamzzz=l_Gamzzz;
    }
    if (Rxx) {
        *Rxx=l_Rxx;*Rxy=l_Rxy;*Rxz=l_Rxz;
        *Ryy=l_Ryy;*Ryz=l_Ryz;*Rzz=l_Rzz;
    }

    // Continue with constraint computation using the same Ricci/Gamma
    const double l_gxx = dxx[idx] + ONE;
    const double l_gxy = gxy[idx], l_gxz = gxz[idx];
    const double l_gyy = dyy[idx] + ONE, l_gyz = gyz[idx];
    const double l_gzz = dzz[idx] + ONE;

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    invert_metric(l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz,
                  gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);

    const double l_Kxx=Kxx[idx], l_Kxy=Kxy[idx], l_Kxz=Kxz[idx];
    const double l_Kyy=Kyy[idx], l_Kyz=Kyz[idx], l_Kzz=Kzz[idx];
    const double trK_val = gupxx*l_Kxx+gupyy*l_Kyy+gupzz*l_Kzz
                         + TWO*(gupxy*l_Kxy+gupxz*l_Kxz+gupyz*l_Kyz);
    const double trR = gupxx*l_Rxx+gupyy*l_Ryy+gupzz*l_Rzz
                     + TWO*(gupxy*l_Rxy+gupxz*l_Rxz+gupyz*l_Ryz);
    ham_Res = trR + trK_val*trK_val
            - contract_symmetric_square(gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
                l_Kxx,l_Kxy,l_Kxz,l_Kyy,l_Kyz,l_Kzz)
            - F16*CUDART_PI*rho[idx];

    double Kxxx,Kxxy,Kxxz,Kxyx,Kxyy,Kxyz,Kxzx,Kxzy,Kxzz;
    double Kyyx,Kyyy,Kyyz,Kyzx,Kyzy,Kyzz,Kzzx,Kzzy,Kzzz;
    double trK_x,trK_y,trK_z;
    if constexpr (Shell) {
        d_fderivs_sh_point_gw<GW>(ex,pex,Kxx,&Kxxx,&Kxxy,&Kxxz,crho,sigma,Rcoord,SYM,SYM,SYM,symmetry,0,sst,i,j,k);
        d_fderivs_sh_point_gw<GW>(ex,pex,Kxy,&Kxyx,&Kxyy,&Kxyz,crho,sigma,Rcoord,ANTI,ANTI,SYM,symmetry,0,sst,i,j,k);
        d_fderivs_sh_point_gw<GW>(ex,pex,Kxz,&Kxzx,&Kxzy,&Kxzz,crho,sigma,Rcoord,ANTI,SYM,ANTI,symmetry,0,sst,i,j,k);
        d_fderivs_sh_point_gw<GW>(ex,pex,Kyy,&Kyyx,&Kyyy,&Kyyz,crho,sigma,Rcoord,SYM,SYM,SYM,symmetry,0,sst,i,j,k);
        d_fderivs_sh_point_gw<GW>(ex,pex,Kyz,&Kyzx,&Kyzy,&Kyzz,crho,sigma,Rcoord,SYM,ANTI,ANTI,symmetry,0,sst,i,j,k);
        d_fderivs_sh_point_gw<GW>(ex,pex,Kzz,&Kzzx,&Kzzy,&Kzzz,crho,sigma,Rcoord,SYM,SYM,SYM,symmetry,0,sst,i,j,k);
        d_fderivs_sh_point_gw<GW>(ex,pex,trK_arr,&trK_x,&trK_y,&trK_z,crho,sigma,Rcoord,SYM,SYM,SYM,symmetry,0,sst,i,j,k);
    } else {
        d_fderivs_point_gw<GW>(ex,pex,Kxx,&Kxxx,&Kxxy,&Kxxz,X,Y,Z,SYM,SYM,SYM,symmetry,0,i,j,k);
        d_fderivs_point_gw<GW>(ex,pex,Kxy,&Kxyx,&Kxyy,&Kxyz,X,Y,Z,ANTI,ANTI,SYM,symmetry,0,i,j,k);
        d_fderivs_point_gw<GW>(ex,pex,Kxz,&Kxzx,&Kxzy,&Kxzz,X,Y,Z,ANTI,SYM,ANTI,symmetry,0,i,j,k);
        d_fderivs_point_gw<GW>(ex,pex,Kyy,&Kyyx,&Kyyy,&Kyyz,X,Y,Z,SYM,SYM,SYM,symmetry,0,i,j,k);
        d_fderivs_point_gw<GW>(ex,pex,Kyz,&Kyzx,&Kyzy,&Kyzz,X,Y,Z,SYM,ANTI,ANTI,symmetry,0,i,j,k);
        d_fderivs_point_gw<GW>(ex,pex,Kzz,&Kzzx,&Kzzy,&Kzzz,X,Y,Z,SYM,SYM,SYM,symmetry,0,i,j,k);
        d_fderivs_point_gw<GW>(ex,pex,trK_arr,&trK_x,&trK_y,&trK_z,X,Y,Z,SYM,SYM,SYM,symmetry,0,i,j,k);
    }
    covariant_momentum(gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
        l_Kxx,l_Kxy,l_Kxz,l_Kyy,l_Kyz,l_Kzz,
        Kxxx,Kxyx,Kxzx,Kyyx,Kyzx,Kzzx,Kxxy,Kxyy,Kxzy,Kyyy,Kyzy,Kzzy,Kxxz,Kxyz,Kxzz,Kyyz,Kyzz,Kzzz,
        l_Gamxxx,l_Gamxxy,l_Gamxxz,l_Gamxyy,l_Gamxyz,l_Gamxzz,
        l_Gamyxx,l_Gamyxy,l_Gamyxz,l_Gamyyy,l_Gamyyz,l_Gamyzz,
        l_Gamzxx,l_Gamzxy,l_Gamzxz,l_Gamzyy,l_Gamzyz,l_Gamzzz,
        0.0,0.0,0.0,0.0, movx_Res,movy_Res,movz_Res);
    movx_Res -= trK_x + F8*CUDART_PI*Sx[idx];
    movy_Res -= trK_y + F8*CUDART_PI*Sy[idx];
    movz_Res -= trK_z + F8*CUDART_PI*Sz[idx];
}

// ============================================================
// Helper kernel: precompute trK = gup^ij K_ij at all grid points
// ============================================================
template<int GW>
__global__ void compute_trK_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ dxx,
    const double* __restrict__ gxy,
    const double* __restrict__ gxz,
    const double* __restrict__ dyy,
    const double* __restrict__ gyz,
    const double* __restrict__ dzz,
    const double* __restrict__ Kxx,
    const double* __restrict__ Kxy,
    const double* __restrict__ Kxz,
    const double* __restrict__ Kyy,
    const double* __restrict__ Kyz,
    const double* __restrict__ Kzz,
    double* __restrict__ trK_out)
{
    const int n = ex0 * ex1 * ex2;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;

    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);

    const int pex[3] = {pex0, pex1, pex2};
    const size_t id = idx3(i, j, k, pex);

    const double l_gxx = dxx[id] + ONE;
    const double l_gxy = gxy[id];
    const double l_gxz = gxz[id];
    const double l_gyy = dyy[id] + ONE;
    const double l_gyz = gyz[id];
    const double l_gzz = dzz[id] + ONE;

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    invert_metric(l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz,
                  gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);

    trK_out[id] = gupxx * Kxx[id] + gupyy * Kyy[id] + gupzz * Kzz[id]
                + TWO * (gupxy * Kxy[id] + gupxz * Kxz[id] + gupyz * Kyz[id]);

    (void)GW;
}

// ============================================================
// Main ADM constraint kernel: iterate over all grid points
// ============================================================
template<int GW, bool Shell>
__global__ void constraint_adm_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ X,
    const double* __restrict__ Y,
    const double* __restrict__ Z,
    const double* __restrict__ crho,
    const double* __restrict__ sigma,
    const double* __restrict__ Rcoord,
    const double* __restrict__ dxx,
    const double* __restrict__ gxy,
    const double* __restrict__ gxz,
    const double* __restrict__ dyy,
    const double* __restrict__ gyz,
    const double* __restrict__ dzz,
    const double* __restrict__ Kxx,
    const double* __restrict__ Kxy,
    const double* __restrict__ Kxz,
    const double* __restrict__ Kyy,
    const double* __restrict__ Kyz,
    const double* __restrict__ Kzz,
    const double* __restrict__ trK_arr,
    const double* __restrict__ rho,
    const double* __restrict__ Sx,
    const double* __restrict__ Sy,
    const double* __restrict__ Sz,
    double* __restrict__ ham_Res,
    double* __restrict__ movx_Res,
    double* __restrict__ movy_Res,
    double* __restrict__ movz_Res,
    int symmetry,
    int sst)
{
    const int n = ex0 * ex1 * ex2;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;

    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);

    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const size_t idx = idx3(i, j, k, pex);

    d_adm_constraint_point<GW, Shell>(
        ex, pex, i, j, k, idx,
        X, Y, Z, crho, sigma, Rcoord,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Kxx, Kxy, Kxz, Kyy, Kyz, Kzz,
        trK_arr,
        rho, Sx, Sy, Sz,
        ham_Res[idx], movx_Res[idx], movy_Res[idx], movz_Res[idx],
        symmetry, sst);
}

// ============================================================
// Kernel to populate Gamma and Ricci output arrays (shell)
// ============================================================
template<int GW>
__global__ void fill_gamma_ricci_ss_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* crho, const double* sigma, const double* Rcoord,
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
    const int n = ex0 * ex1 * ex2;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;
    int i, j, k; linear_to_ijk(tid, ex0, ex1, i, j, k);
    const int ex[3] = {ex0, ex1, ex2}, pex[3] = {pex0, pex1, pex2};
    size_t idx = idx3(i, j, k, pex);
    d_adm_ricci_gamma_point_ss<GW>(ex, pex, i, j, k, idx,
        crho, sigma, Rcoord,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Gamxxx[idx], Gamxxy[idx], Gamxxz[idx],
        Gamxyy[idx], Gamxyz[idx], Gamxzz[idx],
        Gamyxx[idx], Gamyxy[idx], Gamyxz[idx],
        Gamyyy[idx], Gamyyz[idx], Gamyzz[idx],
        Gamzxx[idx], Gamzxy[idx], Gamzxz[idx],
        Gamzyy[idx], Gamzyz[idx], Gamzzz[idx],
        Rxx[idx], Rxy[idx], Rxz[idx],
        Ryy[idx], Ryz[idx], Rzz[idx], symmetry, sst);
}

} // namespace

// ============================================================
// Host launch functions
// ============================================================

extern "C" void gpu_constraint_adm_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    double* X, double* Y, double* Z,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Kxx, double* Kxy, double* Kxz,
    double* Kyy, double* Kyz, double* Kzz,
    double* Lap, double* Sfx, double* Sfy, double* Sfz,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    int symmetry)
{
    (void)Lap; (void)Sfx; (void)Sfy; (void)Sfz;

    const int n = ex[0] * ex[1] * ex[2];
    const int block = 256;
    const int grid = (n + block - 1) / block;
    const size_t num_elements = static_cast<size_t>(pex[0]) * pex[1] * pex[2];

    GPUManager& gpu_manager = GPUManager::getInstance();
    double* trK_temp = gpu_manager.allocate_device_memory<double>(num_elements, stream);

    // Kernel 1: precompute trK at all grid points
    compute_trK_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        dxx, gxy, gxz, dyy, gyz, dzz,
        Kxx, Kxy, Kxz, Kyy, Kyz, Kzz,
        trK_temp);

    // Kernel 2: constraint computation (Cartesian)
    constraint_adm_kernel<ghost_width, false><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        X, Y, Z,
        nullptr, nullptr, nullptr,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Kxx, Kxy, Kxz, Kyy, Kyz, Kzz,
        trK_temp,
        rho, Sx, Sy, Sz,
        ham_Res, movx_Res, movy_Res, movz_Res,
        symmetry, 0);

    gpu_manager.free_device_memory(trK_temp, num_elements, stream);
}

extern "C" void gpu_constraint_adm_ss_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    double* crho, double* sigma, double* Rcoord,
    double* X, double* Y, double* Z,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* drhodxx, double* drhodxy, double* drhodxz,
    double* drhodyy, double* drhodyz, double* drhodzz,
    double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
    double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
    double* dRdxx, double* dRdxy, double* dRdxz,
    double* dRdyy, double* dRdyz, double* dRdzz,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Kxx, double* Kxy, double* Kxz,
    double* Kyy, double* Kyz, double* Kzz,
    double* Lap, double* Sfx, double* Sfy, double* Sfz,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    int symmetry, int lev, int sst)
{
    // Shell derivative stencils compute directly in shell coordinates;
    // Jacobian arrays (drhod*, dsigmad*, dRd*, drhod??, ...) are unused.
    (void)X; (void)Y; (void)Z;
    (void)Lap; (void)Sfx; (void)Sfy; (void)Sfz;
    (void)drhodxx; (void)drhodxy; (void)drhodxz;
    (void)drhodyy; (void)drhodyz; (void)drhodzz;
    (void)dsigmadxx; (void)dsigmadxy; (void)dsigmadxz;
    (void)dsigmadyy; (void)dsigmadyz; (void)dsigmadzz;
    (void)dRdxx; (void)dRdxy; (void)dRdxz;
    (void)dRdyy; (void)dRdyz; (void)dRdzz;
    (void)lev;
    (void)drhodx; (void)drhody; (void)drhodz;
    (void)dsigmadx; (void)dsigmady; (void)dsigmadz;
    (void)dRdx; (void)dRdy; (void)dRdz;

    const int n = ex[0] * ex[1] * ex[2];
    const int block = 256;
    const int grid = (n + block - 1) / block;
    const size_t num_elements = static_cast<size_t>(pex[0]) * pex[1] * pex[2];

    GPUManager& gpu_manager = GPUManager::getInstance();

    // Kernel 1: populate Gamma and Ricci output arrays (matches Fortran intent(out))
    fill_gamma_ricci_ss_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        crho, sigma, Rcoord,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, symmetry, sst);

    // Kernel 2: precompute trK at all grid points
    double* trK_temp = gpu_manager.allocate_device_memory<double>(num_elements, stream);
    compute_trK_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        dxx, gxy, gxz, dyy, gyz, dzz,
        Kxx, Kxy, Kxz, Kyy, Kyz, Kzz,
        trK_temp);

    // Kernel 3: constraint computation (shell)
    constraint_adm_kernel<ghost_width, true><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        crho, sigma, Rcoord,     // X, Y, Z (unused by Shell path)
        crho, sigma, Rcoord,     // crho, sigma, Rcoord (used by Shell path)
        dxx, gxy, gxz, dyy, gyz, dzz,
        Kxx, Kxy, Kxz, Kyy, Kyz, Kzz,
        trK_temp,
        rho, Sx, Sy, Sz,
        ham_Res, movx_Res, movy_Res, movz_Res,
        symmetry, sst);

    gpu_manager.free_device_memory(trK_temp, num_elements, stream);
}
