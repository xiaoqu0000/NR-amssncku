#include "macrodef.h"
// Z4c shell-coordinate RHS GPU implementation (Issue #20)
// Translates compute_rhs_z4c_ss from Z4c_rhs_ss.f90 to CUDA.
// Uses existing shell derivative device functions from point_diff_new_sh_gpu.cuh
// and KO dissipation from kodiss_sh_gpu.cuh.

#include "Z4c_rhs_ss_gpu.cuh"
#include "Z4c_rhs_gpu.cuh"
#include "point_diff_new_sh_gpu.cuh"
#include "kodiss_sh_gpu.cuh"
#include "diff_new_sh_gpu.cuh"
#include "gpu_info.h"
#include "macrodef.fh"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math_constants.h>
#include <math.h>

#if (ABV != 0)
#error "Z4c_rhs_ss GPU only supports ABV==0. ABV==1 (ricci_gamma_ss + constraint_bssn_ss) is not yet ported."
#endif

#ifndef Z4C_KAPPA1
#define Z4C_KAPPA1 0.02
#endif
#ifndef Z4C_KAPPA2
#define Z4C_KAPPA2 0.0
#endif

namespace {

// Physical constants matching Fortran get_Z4cparameters
constexpr double ZEO = 0.0;
constexpr double ONE = 1.0;
constexpr double TWO = 2.0;
constexpr double FOUR = 4.0;
constexpr double HALF = 0.5;
constexpr double EIGHT = 8.0;
constexpr double F1o3 = 1.0 / 3.0;
constexpr double F2o3 = 2.0 / 3.0;
constexpr double F3o2 = 1.5;
constexpr double F8 = 8.0;
constexpr double F16 = 16.0;
constexpr double SYM_VAL = 1.0;
constexpr double ANTI_VAL = -1.0;
constexpr double PI_VAL = CUDART_PI;
constexpr double kappa1 = Z4C_KAPPA1;
constexpr double kappa2 = Z4C_KAPPA2;
constexpr double FF = 0.75;
constexpr double eta_val = 2.0;

// ==========================================================================
// Shell derivative wrappers (file scope, replace per-kernel macros)
// ==========================================================================
template<int GW>
__device__ __forceinline__ void deriv1_sh(
    const int ex[3], const int pex[3], const double* f,
    double &fx, double &fy, double &fz,
    double s1, double s2, double s3,
    const double* crho, const double* sigma, const double* R,
    int symmetry, int onoff, int sst,
    const double* drhodx, const double* drhody, const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx, const double* dRdy, const double* dRdz,
    int i, int j, int k)
{
    d_point_fderivs_shc_gw<GW>(ex, pex, f, &fx, &fy, &fz,
        crho, sigma, R, s1, s2, s3, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, i, j, k);
}

template<int GW>
__device__ __forceinline__ void deriv2_sh(
    const int ex[3], const int pex[3], const double* f,
    double &fxx, double &fxy, double &fxz, double &fyy, double &fyz, double &fzz,
    double s1, double s2, double s3,
    const double* crho, const double* sigma, const double* R,
    int symmetry, int onoff, int sst,
    const double* drhodx, const double* drhody, const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx, const double* dRdy, const double* dRdz,
    const double* drhodxx, const double* drhodxy, const double* drhodxz,
    const double* drhodyy, const double* drhodyz, const double* drhodzz,
    const double* dsigmadxx, const double* dsigmadxy, const double* dsigmadxz,
    const double* dsigmadyy, const double* dsigmadyz, const double* dsigmadzz,
    const double* dRdxx, const double* dRdxy, const double* dRdxz,
    const double* dRdyy, const double* dRdyz, const double* dRdzz,
    int i, int j, int k)
{
    d_point_fdderivs_shc_gw<GW>(ex, pex, f, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz,
        crho, sigma, R, s1, s2, s3, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);
}

// ==========================================================================
// Kernel 1+2 merged: Ricci + RHS eval + constraints
// Translates lines 5-1184 of Z4c_rhs_ss.f90
// ==========================================================================
template<int GW>
__launch_bounds__(GPU_DIM_X * GPU_DIM_Y * GPU_DIM_Z, 1)
__global__ void z4c_ss_main_kernel(
    const int ex0, const int ex1, const int ex2,
    const int pex0, const int pex1, const int pex2,
    const int symmetry, const int lev, const double T, const int sst,
    const double* __restrict__ crho,
    const double* __restrict__ sigma,
    const double* __restrict__ R,
    const double* __restrict__ drhodx,
    const double* __restrict__ drhody,
    const double* __restrict__ drhodz,
    const double* __restrict__ dsigmadx,
    const double* __restrict__ dsigmady,
    const double* __restrict__ dsigmadz,
    const double* __restrict__ dRdx,
    const double* __restrict__ dRdy,
    const double* __restrict__ dRdz,
    const double* __restrict__ drhodxx,
    const double* __restrict__ drhodxy,
    const double* __restrict__ drhodxz,
    const double* __restrict__ drhodyy,
    const double* __restrict__ drhodyz,
    const double* __restrict__ drhodzz,
    const double* __restrict__ dsigmadxx,
    const double* __restrict__ dsigmadxy,
    const double* __restrict__ dsigmadxz,
    const double* __restrict__ dsigmadyy,
    const double* __restrict__ dsigmadyz,
    const double* __restrict__ dsigmadzz,
    const double* __restrict__ dRdxx,
    const double* __restrict__ dRdxy,
    const double* __restrict__ dRdxz,
    const double* __restrict__ dRdyy,
    const double* __restrict__ dRdyz,
    const double* __restrict__ dRdzz,
    const double* __restrict__ chi,
    const double* __restrict__ trK,
    const double* __restrict__ dxx,
    const double* __restrict__ gxy,
    const double* __restrict__ gxz,
    const double* __restrict__ dyy,
    const double* __restrict__ gyz,
    const double* __restrict__ dzz,
    const double* __restrict__ Axx,
    const double* __restrict__ Axy,
    const double* __restrict__ Axz,
    const double* __restrict__ Ayy,
    const double* __restrict__ Ayz,
    const double* __restrict__ Azz,
    const double* __restrict__ Gamx,
    const double* __restrict__ Gamy,
    const double* __restrict__ Gamz,
    const double* __restrict__ Lap,
    const double* __restrict__ betax,
    const double* __restrict__ betay,
    const double* __restrict__ betaz,
    const double* __restrict__ dtSfx,
    const double* __restrict__ dtSfy,
    const double* __restrict__ dtSfz,
    const double* __restrict__ TZ,
    double* __restrict__ chi_rhs,
    double* __restrict__ trK_rhs,
    double* __restrict__ gxx_rhs,
    double* __restrict__ gxy_rhs,
    double* __restrict__ gxz_rhs,
    double* __restrict__ gyy_rhs,
    double* __restrict__ gyz_rhs,
    double* __restrict__ gzz_rhs,
    double* __restrict__ Axx_rhs,
    double* __restrict__ Axy_rhs,
    double* __restrict__ Axz_rhs,
    double* __restrict__ Ayy_rhs,
    double* __restrict__ Ayz_rhs,
    double* __restrict__ Azz_rhs,
    double* __restrict__ Gamx_rhs,
    double* __restrict__ Gamy_rhs,
    double* __restrict__ Gamz_rhs,
    double* __restrict__ Lap_rhs,
    double* __restrict__ betax_rhs,
    double* __restrict__ betay_rhs,
    double* __restrict__ betaz_rhs,
    double* __restrict__ dtSfx_rhs,
    double* __restrict__ dtSfy_rhs,
    double* __restrict__ dtSfz_rhs,
    double* __restrict__ TZ_rhs,
    const double* __restrict__ rho,
    const double* __restrict__ Sx,
    const double* __restrict__ Sy,
    const double* __restrict__ Sz,
    const double* __restrict__ Sxx,
    const double* __restrict__ Sxy,
    const double* __restrict__ Sxz,
    const double* __restrict__ Syy,
    const double* __restrict__ Syz,
    const double* __restrict__ Szz,
    double* __restrict__ Gamxxx_out,
    double* __restrict__ Gamxxy_out,
    double* __restrict__ Gamxxz_out,
    double* __restrict__ Gamxyy_out,
    double* __restrict__ Gamxyz_out,
    double* __restrict__ Gamxzz_out,
    double* __restrict__ Gamyxx_out,
    double* __restrict__ Gamyxy_out,
    double* __restrict__ Gamyxz_out,
    double* __restrict__ Gamyyy_out,
    double* __restrict__ Gamyyz_out,
    double* __restrict__ Gamyzz_out,
    double* __restrict__ Gamzxx_out,
    double* __restrict__ Gamzxy_out,
    double* __restrict__ Gamzxz_out,
    double* __restrict__ Gamzyy_out,
    double* __restrict__ Gamzyz_out,
    double* __restrict__ Gamzzz_out,
    double* __restrict__ Rxx_out,
    double* __restrict__ Rxy_out,
    double* __restrict__ Rxz_out,
    double* __restrict__ Ryy_out,
    double* __restrict__ Ryz_out,
    double* __restrict__ Rzz_out,
    double* __restrict__ Hcon,
    double* __restrict__ Mxcon,
    double* __restrict__ Mycon,
    double* __restrict__ Mzcon,
    double* __restrict__ Gmxcon,
    double* __restrict__ Gmycon,
    double* __restrict__ Gmzcon
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    (void)T;   // unused in current implementation
    (void)lev; // passed through but onoff=0 used instead

    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    int idx = z4c_idx3d(i, j, k, pex0, pex1);

    // =====================================================================
    // Phase 0: Basic quantities
    // =====================================================================
    double l_chi = chi[idx];
    double l_trK = trK[idx];
    double l_dxx = dxx[idx];
    double l_gxy = gxy[idx];
    double l_gxz = gxz[idx];
    double l_dyy = dyy[idx];
    double l_gyz = gyz[idx];
    double l_dzz = dzz[idx];
    double l_Axx = Axx[idx];
    double l_Axy = Axy[idx];
    double l_Axz = Axz[idx];
    double l_Ayy = Ayy[idx];
    double l_Ayz = Ayz[idx];
    double l_Azz = Azz[idx];
    double l_Gamx = Gamx[idx];
    double l_Gamy = Gamy[idx];
    double l_Gamz = Gamz[idx];
    double l_Lap  = Lap[idx];
    double l_dtSfx = dtSfx[idx];
    double l_dtSfy = dtSfy[idx];
    double l_dtSfz = dtSfz[idx];
    double l_TZ   = TZ[idx];
    double l_betax = betax[idx];
    double l_betay = betay[idx];
    double l_betaz = betaz[idx];
    double l_rho  = rho[idx];
    double l_Sx   = Sx[idx];
    double l_Sy   = Sy[idx];
    double l_Sz   = Sz[idx];
    double l_Sxx  = Sxx[idx];
    double l_Sxy  = Sxy[idx];
    double l_Sxz  = Sxz[idx];
    double l_Syy  = Syy[idx];
    double l_Syz  = Syz[idx];
    double l_Szz  = Szz[idx];

    double alpn1 = l_Lap + ONE;
    double chin1 = l_chi + ONE;
    double l_gxx = l_dxx + ONE;
    double l_gyy = l_dyy + ONE;
    double l_gzz = l_dzz + ONE;
    double trKd  = l_trK + TWO * l_TZ;

    // =====================================================================
    // Phase 1: Shift derivatives, div_beta, chi_rhs
    // =====================================================================
    double betaxx, betaxy, betaxz;
    double betayx, betayy, betayz;
    double betazx, betazy, betazz;
    deriv1_sh<GW>(ex, pex, betax, betaxx, betaxy, betaxz, ANTI_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    deriv1_sh<GW>(ex, pex, betay, betayx, betayy, betayz, SYM_VAL, ANTI_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    deriv1_sh<GW>(ex, pex, betaz, betazx, betazy, betazz, SYM_VAL, SYM_VAL, ANTI_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double div_beta = betaxx + betayy + betazz;

    double chix, chiy, chiz;
    deriv1_sh<GW>(ex, pex, chi, chix, chiy, chiz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double chi_rhs_val = F2o3 * chin1 * (alpn1 * trKd - div_beta);

    // =====================================================================
    // Phase 2: Metric derivatives, gij_rhs
    // =====================================================================
    double gxxx, gxxy, gxxz;
    double gxyx, gxyy, gxyz;
    double gxzx, gxzy, gxzz;
    double gyyx, gyyy, gyyz;
    double gyzx, gyzy, gyzz;
    double gzzx, gzzy, gzzz;

    deriv1_sh<GW>(ex, pex, dxx, gxxx, gxxy, gxxz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    deriv1_sh<GW>(ex, pex, gxy, gxyx, gxyy, gxyz, ANTI_VAL, ANTI_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    deriv1_sh<GW>(ex, pex, gxz, gxzx, gxzy, gxzz, ANTI_VAL, SYM_VAL, ANTI_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    deriv1_sh<GW>(ex, pex, dyy, gyyx, gyyy, gyyz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    deriv1_sh<GW>(ex, pex, gyz, gyzx, gyzy, gyzz, SYM_VAL, ANTI_VAL, ANTI_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    deriv1_sh<GW>(ex, pex, dzz, gzzx, gzzy, gzzz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double gxx_rhs_val = -TWO * alpn1 * l_Axx - F2o3 * l_gxx * div_beta
        + TWO * (l_gxx * betaxx + l_gxy * betayx + l_gxz * betazx);
    double gyy_rhs_val = -TWO * alpn1 * l_Ayy - F2o3 * l_gyy * div_beta
        + TWO * (l_gxy * betaxy + l_gyy * betayy + l_gyz * betazy);
    double gzz_rhs_val = -TWO * alpn1 * l_Azz - F2o3 * l_gzz * div_beta
        + TWO * (l_gxz * betaxz + l_gyz * betayz + l_gzz * betazz);
    double gxy_rhs_val = -TWO * alpn1 * l_Axy + F1o3 * l_gxy * div_beta
        + l_gxx * betaxy + l_gxz * betazy + l_gyy * betayx + l_gyz * betazx - l_gxy * betazz;
    double gyz_rhs_val = -TWO * alpn1 * l_Ayz + F1o3 * l_gyz * div_beta
        + l_gxy * betaxz + l_gyy * betayz + l_gxz * betaxy + l_gzz * betazy - l_gyz * betaxx;
    double gxz_rhs_val = -TWO * alpn1 * l_Axz + F1o3 * l_gxz * div_beta
        + l_gxx * betaxz + l_gxy * betayz + l_gyz * betayx + l_gzz * betazx - l_gxz * betayy;

    // =====================================================================
    // Phase 3: Inverse metric
    // =====================================================================
    double gupzz_det = l_gxx * l_gyy * l_gzz + l_gxy * l_gyz * l_gxz
        + l_gxz * l_gxy * l_gyz - l_gxz * l_gyy * l_gxz
        - l_gxy * l_gxy * l_gzz - l_gxx * l_gyz * l_gyz;
    double inv_det = ONE / gupzz_det;
    double gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) * inv_det;
    double gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) * inv_det;
    double gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) * inv_det;
    double gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) * inv_det;
    double gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) * inv_det;
    gupzz_det = (l_gxx * l_gyy - l_gxy * l_gxy) * inv_det;
    double gupzz = gupzz_det;

    // =====================================================================
    // Phase 4: Second derivatives of metric (gij_,kl)
    // =====================================================================
    double gxxxx, gxxxy, gxxxz, gxxyy, gxxyz, gxxzz;
    double gyyxx, gyyxy, gyyxz, gyyyy, gyyyz, gyyzz;
    double gzzxx, gzzxy, gzzxz, gzzyy, gzzyz, gzzzz;
    double gxyxx, gxyxy, gxyxz, gxyyy, gxyyz, gxyzz;
    double gxzxx, gxzxy, gxzxz, gxzyy, gxzyz, gxzzz;
    double gyzxx, gyzxy, gyzxz, gyzyy, gyzyz, gyzzz;

    deriv2_sh<GW>(ex, pex, dxx, gxxxx, gxxxy, gxxxz, gxxyy, gxxyz, gxxzz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);
    deriv2_sh<GW>(ex, pex, dyy, gyyxx, gyyxy, gyyxz, gyyyy, gyyyz, gyyzz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);
    deriv2_sh<GW>(ex, pex, dzz, gzzxx, gzzxy, gzzxz, gzzyy, gzzyz, gzzzz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);
    deriv2_sh<GW>(ex, pex, gxy, gxyxx, gxyxy, gxyxz, gxyyy, gxyyz, gxyzz, ANTI_VAL, ANTI_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);
    deriv2_sh<GW>(ex, pex, gxz, gxzxx, gxzxy, gxzxz, gxzyy, gxzyz, gxzzz, ANTI_VAL, SYM_VAL, ANTI_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);
    deriv2_sh<GW>(ex, pex, gyz, gyzxx, gyzxy, gyzxz, gyzyy, gyzyz, gyzzz, SYM_VAL, ANTI_VAL, ANTI_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);

    // =====================================================================
    // Phase 5: Christoffel symbols (second kind)
    // =====================================================================
    double Gamxxx = HALF * (gupxx * gxxx + gupxy * (TWO * gxyx - gxxy) + gupxz * (TWO * gxzx - gxxz));
    double Gamyxx = HALF * (gupxy * gxxx + gupyy * (TWO * gxyx - gxxy) + gupyz * (TWO * gxzx - gxxz));
    double Gamzxx = HALF * (gupxz * gxxx + gupyz * (TWO * gxyx - gxxy) + gupzz * (TWO * gxzx - gxxz));

    double Gamxyy = HALF * (gupxx * (TWO * gxyy - gyyx) + gupxy * gyyy + gupxz * (TWO * gyzy - gyyz));
    double Gamyyy = HALF * (gupxy * (TWO * gxyy - gyyx) + gupyy * gyyy + gupyz * (TWO * gyzy - gyyz));
    double Gamzyy = HALF * (gupxz * (TWO * gxyy - gyyx) + gupyz * gyyy + gupzz * (TWO * gyzy - gyyz));

    double Gamxzz = HALF * (gupxx * (TWO * gxzz - gzzx) + gupxy * (TWO * gyzz - gzzy) + gupxz * gzzz);
    double Gamyzz = HALF * (gupxy * (TWO * gxzz - gzzx) + gupyy * (TWO * gyzz - gzzy) + gupyz * gzzz);
    double Gamzzz = HALF * (gupxz * (TWO * gxzz - gzzx) + gupyz * (TWO * gyzz - gzzy) + gupzz * gzzz);

    double Gamxxy = HALF * (gupxx * gxxy + gupxy * gyyx + gupxz * (gxzy + gyzx - gxyz));
    double Gamyxy = HALF * (gupxy * gxxy + gupyy * gyyx + gupyz * (gxzy + gyzx - gxyz));
    double Gamzxy = HALF * (gupxz * gxxy + gupyz * gyyx + gupzz * (gxzy + gyzx - gxyz));

    double Gamxxz = HALF * (gupxx * gxxz + gupxy * (gxyz + gyzx - gxzy) + gupxz * gzzx);
    double Gamyxz = HALF * (gupxy * gxxz + gupyy * (gxyz + gyzx - gxzy) + gupyz * gzzx);
    double Gamzxz = HALF * (gupxz * gxxz + gupyz * (gxyz + gyzx - gxzy) + gupzz * gzzx);

    double Gamxyz = HALF * (gupxx * (gxyz + gxzy - gyzx) + gupxy * gyyz + gupxz * gzzy);
    double Gamyyz = HALF * (gupxy * (gxyz + gxzy - gyzx) + gupyy * gyyz + gupyz * gzzy);
    double Gamzyz = HALF * (gupxz * (gxyz + gxzy - gyzx) + gupyz * gyyz + gupzz * gzzy);

    // Gam^a = g^{ki} Gam^a_{ki}
    double Gamxa = gupxx * Gamxxx + gupyy * Gamxyy + gupzz * Gamxzz
        + TWO * (gupxy * Gamxxy + gupxz * Gamxxz + gupyz * Gamxyz);
    double Gamya = gupxx * Gamyxx + gupyy * Gamyyy + gupzz * Gamyzz
        + TWO * (gupxy * Gamyxy + gupxz * Gamyxz + gupyz * Gamyyz);
    double Gamza = gupxx * Gamzxx + gupyy * Gamzyy + gupzz * Gamzzz
        + TWO * (gupxy * Gamzxy + gupxz * Gamzxz + gupyz * Gamzyz);

    // Gmxcon = Gam^i - Gam^a
    double l_Gmxcon = l_Gamx - Gamxa;
    double l_Gmycon = l_Gamy - Gamya;
    double l_Gmzcon = l_Gamz - Gamza;

    // =====================================================================
    // Phase 6: Maple contraction for TZ_rhs (g^ki g^jm g^ln g_mn,k g_ij,l)
    // =====================================================================
    // Reuse Gamxx, Gamyy, Gamxz, Gamxy, Gamyx, Gamyz, Gamzx, Gamzy, Gamzz as temp
    double tmp_Gamyy = 3*gupxz*gupxz*gupzz*gxzz*gxzz + gupxx*gupxz*gupxz*gxxz*gxxz
        + 2*gxyx*gupxy*gupxy*gupxy*gxyy + 2*gxyx*gupxy*gupxy*gupxy*gyyx
        + gupxx*gupxx*gupzz*gxzx*gxzx + 3*gupxx*gupxy*gupxy*gxyx*gxyx
        + 6*gxyx*gupxy*gupxz*gupyy*gyzy + gupxx*gupxx*gupyy*gxyx*gxyx
        + 2*gxyz*gupxy*gupyz*gupyz*gyyz + 2*gxxz*gupxx*gupxx*gupyz*gxyx
        + gupxz*gupxz*gupyy*gyzx*gyzx + 2*gxxy*gupxx*gupxy*gupxz*gxxz
        + 2*gyzx*gupxy*gupxz*gupzz*gzzx + 3*gupyy*gupyz*gupyz*gyzy*gyzy
        + 2*gyyy*gupyz*gupyz*gupyz*gzzz + 2*gxxz*gupxz*gupxz*gupxz*gxzz
        + 4*gxzy*gupxx*gupxz*gupyy*gxyx + gupyy*gupyz*gupyz*gyyz*gyyz;

    double tmp_Gamxz = tmp_Gamyy + 2*gxxz*gupxy*gupxy*gupzz*gyzy
        + 4*gxyz*gupxx*gupxy*gupxz*gxxx + 6*gxzz*gupxy*gupyz*gupzz*gyzy
        + 2*gxxy*gupxx*gupxz*gupyz*gxzz + 3*gupxy*gupxy*gupyy*gxyy*gxyy
        + 2*gxyz*gupxx*gupyy*gupzz*gyzx + 4*gxyy*gupxx*gupyy*gupyz*gyzx
        + 6*gxyy*gupxy*gupxz*gupyz*gxzz + 4*gxzz*gupxx*gupyz*gupzz*gyzx
        + 3*gupxx*gupxz*gupxz*gxzx*gxzx + 4*gxyz*gupxx*gupxy*gupyz*gxyx
        + 2*gxxz*gupxx*gupxz*gupyz*gxyz + 2*gxxy*gupxy*gupxz*gupyz*gyzz
        + 2*gxzx*gupxy*gupxz*gupyz*gyyz + gupyz*gupyz*gupzz*gzzy*gzzy
        + gupxz*gupxz*gupzz*gzzx*gzzx + gupyy*gupzz*gupzz*gyzz*gyzz
        + 2*gyzy*gupyz*gupyz*gupyz*gzzy + gupxx*gupzz*gupzz*gxzz*gxzz;

    tmp_Gamyy = tmp_Gamxz + gupxx*gupyz*gupyz*gxzy*gxzy + 2*gxzx*gupxz*gupxz*gupxz*gzzx
        + 3*gupyz*gupyz*gupzz*gyzz*gyzz + 2*gyzy*gupyz*gupyz*gupyz*gyzz
        + gupyy*gupyy*gupzz*gyzy*gyzy + gupxy*gupxy*gupzz*gyzx*gyzx
        + 2*gyyz*gupyz*gupyz*gupyz*gyzz + gupxy*gupxy*gupyy*gyyx*gyyx
        + gupxx*gupyz*gupyz*gxyz*gxyz + gupxx*gupyy*gupyy*gxyy*gxyy
        + gupxy*gupxy*gupzz*gxzy*gxzy + 2*gxzx*gupxz*gupxz*gupxz*gxzz
        + 2*gyyx*gupxy*gupxz*gupyy*gyzx + gupxx*gupxy*gupxy*gxxy*gxxy
        + 2*gxxx*gupxz*gupxz*gupxz*gzzz + 2*gxxx*gupxy*gupxy*gupxy*gyyy
        + gupxz*gupxz*gupyy*gxyz*gxyz + 2*gxyy*gupxy*gupxy*gupxy*gxxy;

    tmp_Gamxz = tmp_Gamyy + 2*gxyy*gupxz*gupyy*gupyy*gyzy
        + 6*gxyy*gupxx*gupxy*gupyz*gxzx + 4*gxyy*gupxy*gupxz*gupyy*gxyz
        + 2*gyzx*gupxz*gupyy*gupzz*gzzy + 2*gxzy*gupxy*gupxz*gupyy*gxyy
        + 4*gxzy*gupxy*gupxz*gupzz*gxzz + 2*gyyx*gupxz*gupyy*gupyz*gyzz
        + 6*gxyx*gupxx*gupxz*gupyz*gxzz + 2*gxyz*gupxy*gupxy*gupzz*gxzy
        + 2*gxyz*gupxy*gupxy*gupyz*gxyy + 2*gxyz*gupxy*gupxy*gupxz*gxxy
        + 2*gupxy*gupxz*gupyz*gxyz*gxyz + 4*gxyy*gupxz*gupyz*gupyz*gzzz
        + 2*gxyy*gupxy*gupyz*gupyz*gzzy + 4*gxyy*gupxy*gupxy*gupyz*gxzy
        + 2*gxyy*gupxy*gupxy*gupxz*gxxz + 4*gxyy*gupxx*gupxy*gupxy*gxxx
        + 2*gxyx*gupxy*gupxy*gupxz*gxzy + 2*gxyx*gupxy*gupxy*gupyz*gyzy;

    tmp_Gamyy = tmp_Gamxz + 2*gxyx*gupxx*gupxy*gupxy*gxxy
        + 4*gyzz*gupyz*gupzz*gupzz*gzzz + 4*gxzy*gupxx*gupxz*gupyz*gxzx
        + 2*gxzy*gupxx*gupyy*gupzz*gyzx + 4*gxxx*gupxx*gupxy*gupxz*gyzx
        + 2*gxyx*gupxx*gupxx*gupyz*gxzx + 2*gxyx*gupxy*gupxy*gupxz*gxyz
        + 2*gxzy*gupxz*gupyy*gupyz*gyyz + 4*gxzy*gupxy*gupyy*gupyz*gyyy
        + 2*gxzy*gupxx*gupyy*gupyz*gyyx + 2*gyyx*gupxy*gupxz*gupyy*gxzy
        + 2*gyyx*gupxy*gupyy*gupyz*gyyz + 2*gyyx*gupxy*gupyy*gupyz*gyzy
        + 4*gxzy*gupxz*gupyy*gupzz*gyzz + 2*gyyx*gupxy*gupxz*gupyz*gxzz
        + 2*gxyz*gupxx*gupyy*gupzz*gxzy + 2*gxyy*gupxz*gupyy*gupyz*gzzy;

    tmp_Gamxz = tmp_Gamyy + 2*gxyy*gupxy*gupxz*gupyz*gzzx
        + 2*gxyy*gupxy*gupxz*gupyy*gyzx + 2*gxyy*gupxy*gupyy*gupyz*gyyz
        + 2*gxyy*gupxx*gupyy*gupyz*gxzy + 2*gxxy*gupxy*gupxy*gupxz*gxzy
        + 2*gxxy*gupxy*gupxy*gupyz*gyzy + 2*gxxy*gupxy*gupxy*gupyy*gyyy
        + 2*gxxy*gupxx*gupxx*gupyz*gxzx + 2*gxxy*gupxx*gupxx*gupyy*gxyx
        + 2*gxxx*gupxx*gupxz*gupxz*gzzx + 4*gxxx*gupxy*gupxz*gupxz*gyzz
        + 4*gxxx*gupxy*gupxy*gupxz*gyzy + 2*gxxx*gupxx*gupxy*gupxy*gyyx
        + 4*gxxx*gupxx*gupxz*gupxz*gxzz + 4*gxxx*gupxx*gupxx*gupxz*gxzx
        + 2*gxxx*gupxx*gupxx*gupxz*gxxz + 4*gxyz*gupxz*gupyz*gupyz*gyzz
        + 2*gxyz*gupxy*gupyz*gupyz*gyzy + 2*gxzy*gupxy*gupyy*gupzz*gyzy;

    tmp_Gamyy = tmp_Gamxz + 2*gxyy*gupxx*gupyy*gupyz*gxyz
        + 6*gxzz*gupxz*gupyz*gupzz*gyzz + 4*gxzy*gupxz*gupyz*gupzz*gzzz
        + gupyy*gupyy*gupyy*gyyy*gyyy + 2*gxzy*gupxy*gupyz*gupzz*gzzy
        + 2*gxzy*gupxx*gupyz*gupzz*gzzx + 2*gxyz*gupxx*gupyz*gupzz*gxzz
        + 2*gxzy*gupxx*gupyz*gupzz*gxzz + 2*gyzy*gupxy*gupyz*gupzz*gzzx
        + 2*gyzy*gupxz*gupyy*gupyz*gxzy + 6*gyzy*gupyy*gupyz*gupzz*gyzz
        + 4*gyzx*gupxz*gupyy*gupyz*gyzy + 4*gyzx*gupxy*gupyz*gupzz*gyzz
        + 2*gxxy*gupxx*gupxy*gupyy*gxyy + 4*gyzx*gupxz*gupyz*gupzz*gzzz
        + 2*gyzx*gupxy*gupyy*gupzz*gyzy + 2*gyyz*gupyy*gupyz*gupzz*gzzy
        + 2*gyyz*gupxy*gupyz*gupzz*gzzx;

    tmp_Gamxz = tmp_Gamyy + 2*gyyz*gupyy*gupyz*gupzz*gyzz
        + 2*gyyz*gupxy*gupyy*gupzz*gyzx + 2*gyyz*gupxy*gupyz*gupzz*gxzz
        + 2*gxxy*gupxx*gupxy*gupyz*gyzx + 4*gyyy*gupxy*gupyy*gupyz*gyzx
        + 2*gyyx*gupxy*gupxz*gupyz*gzzx + 2*gxyz*gupxy*gupyz*gupzz*gyzz
        + 2*gxxz*gupxz*gupxz*gupzz*gzzz + 2*gxxz*gupxz*gupxz*gupyz*gyzz
        + 2*gxxz*gupxy*gupxz*gupxz*gxzy + 2*gxxz*gupxx*gupxz*gupxz*gxzx
        + 2*gxxz*gupxy*gupxy*gupyz*gyyy + 2*gxxz*gupxx*gupxx*gupzz*gxzx
        + 2*gxxy*gupxz*gupxz*gupyz*gzzz + 2*gxxy*gupxz*gupxz*gupyy*gyzz
        + 2*gxxy*gupxy*gupxz*gupxz*gxzz + 2*gzzx*gupxz*gupyz*gupzz*gzzy
        + 2*gyzz*gupxz*gupyz*gupzz*gzzx + 2*gxzx*gupxx*gupxz*gupzz*gzzx
        + 2*gyzx*gupxz*gupyy*gupzz*gyzz;

    tmp_Gamyy = tmp_Gamxz + gupzz*gupzz*gupzz*gzzz*gzzz
        + 2*gxzz*gupxy*gupxz*gupzz*gyzx + 6*gxzx*gupxy*gupxz*gupyz*gyzy
        + 2*gxxy*gupxy*gupxz*gupyz*gzzy + 4*gxzz*gupxy*gupyz*gupyz*gyyy
        + 2*gxzy*gupxz*gupyz*gupyz*gyzz + 2*gxzy*gupxz*gupxz*gupyz*gxzz
        + 2*gxzy*gupxz*gupxz*gupyy*gxyz + 2*gupxy*gupxz*gupyz*gxzy*gxzy
        + 4*gxzx*gupxz*gupxz*gupzz*gzzz + 2*gxzx*gupxz*gupxz*gupyz*gyzz
        + 2*gxyz*gupxy*gupxz*gupzz*gzzx + 2*gxyz*gupxz*gupyy*gupzz*gzzy
        + 2*gxyx*gupxx*gupxz*gupyy*gxyz + 2*gxzz*gupxz*gupyz*gupyz*gyyz
        + 2*gxxy*gupxx*gupxy*gupxz*gxzx + 2*gyyx*gupxy*gupxy*gupxz*gxzx;

    tmp_Gamxz = tmp_Gamyy + 2*gxyx*gupxy*gupxz*gupyz*gzzy
        + 2*gyzy*gupyy*gupyz*gupzz*gzzy + 2*gxyx*gupxx*gupxz*gupyy*gyzx
        + 2*gyyx*gupxy*gupyz*gupyz*gyzz + 2*gyyx*gupxy*gupxy*gupyz*gyzx
        + 2*gyyx*gupxz*gupyz*gupyz*gzzz + 2*gyyx*gupxy*gupyy*gupyy*gyyy
        + 2*gxyz*gupxy*gupxy*gupzz*gyzx + 2*gxyz*gupxy*gupxy*gupyz*gyyx
        + 2*gxyy*gupxy*gupyz*gupyz*gyzz + 2*gxyy*gupxy*gupxy*gupyz*gyzx
        + 2*gxyy*gupxy*gupxy*gupyy*gyyx + 2*gxyx*gupxy*gupxz*gupxz*gzzx
        + 2*gxyx*gupxy*gupxy*gupyz*gyyz + 4*gxzz*gupxz*gupzz*gupzz*gzzz
        + 2*gxzz*gupxy*gupzz*gupzz*gzzy + 2*gxzz*gupxx*gupzz*gupzz*gzzx
        + 6*gxyx*gupxx*gupxy*gupxz*gxzx + 2*gxyz*gupxy*gupxz*gupyz*gyzx;

    tmp_Gamyy = tmp_Gamxz + 2*gyyx*gupxz*gupyy*gupyy*gyzy
        + 2*gyyx*gupxz*gupyy*gupyz*gzzy + 2*gxxz*gupxx*gupxy*gupyz*gxyy
        + 2*gyzx*gupxz*gupxz*gupyy*gxzy + 4*gyzx*gupxy*gupxz*gupxz*gxzx
        + 2*gyzx*gupxz*gupyz*gupyz*gyzz + 2*gyzx*gupxz*gupxz*gupyz*gxzz
        + 2*gupxy*gupxz*gupyz*gyzx*gyzx + 2*gyyz*gupyz*gupyz*gupzz*gzzz
        + 2*gyyz*gupyy*gupyz*gupyz*gyzy + 2*gyyz*gupxy*gupyz*gupyz*gyzx
        + 2*gyyz*gupyy*gupyy*gupzz*gyzy + 2*gyyz*gupxy*gupxy*gupzz*gxzx
        + 2*gyyy*gupyy*gupyz*gupyz*gzzy + 2*gyyy*gupxy*gupyz*gupyz*gzzx
        + 4*gyyy*gupyy*gupyz*gupyz*gyzz + 4*gyyy*gupyy*gupyy*gupyz*gyzy
        + 2*gyyy*gupyy*gupyy*gupyz*gyyz;

    tmp_Gamxz = tmp_Gamyy + 2*gxyz*gupxz*gupyy*gupyz*gyzy
        + 2*gxyz*gupxx*gupyy*gupyz*gyyx + 2*gzzx*gupxz*gupzz*gupzz*gzzz
        + 2*gxzy*gupxy*gupxz*gupyz*gyzx + 2*gyzz*gupyz*gupyz*gupzz*gzzy
        + 2*gyzy*gupxz*gupyz*gupyz*gzzx + 2*gyzx*gupxz*gupyz*gupyz*gzzy
        + 2*gyzx*gupxz*gupxz*gupyz*gzzx + 2*gxzz*gupxz*gupxz*gupzz*gzzx
        + 2*gxzz*gupxy*gupzz*gupzz*gyzz + 2*gxzy*gupxz*gupyz*gupyz*gzzy
        + 2*gxzy*gupxz*gupxz*gupyz*gzzx + 2*gxzx*gupxz*gupxz*gupyz*gzzy
        + 2*gyzz*gupyy*gupzz*gupzz*gzzy + 2*gyzz*gupxy*gupzz*gupzz*gzzx
        + 4*gyzy*gupyz*gupyz*gupzz*gzzz + 2*gyzy*gupxy*gupyz*gupyz*gyzx
        + 2*gyzy*gupxz*gupyz*gupyz*gxzz + 2*gxzy*gupxy*gupyz*gupzz*gyzz
        + 2*gxyx*gupxx*gupxy*gupyz*gxzy;

    tmp_Gamyy = tmp_Gamxz + gupxx*gupxx*gupxx*gxxx*gxxx
        + 2*gzzy*gupyz*gupzz*gupzz*gzzz + 6*gxyx*gupxx*gupxy*gupyy*gxyy
        + 2*gxzz*gupxz*gupyz*gupzz*gzzy + 6*gxyx*gupxy*gupxz*gupyz*gyzz
        + 2*gxzx*gupxy*gupxz*gupxz*gxzy + 2*gxyx*gupxx*gupxy*gupyy*gyyx
        + 2*gxyx*gupxx*gupxz*gupyz*gzzx + 2*gxyx*gupxx*gupxy*gupxz*gxxz
        + 4*gxyx*gupxx*gupxx*gupxy*gxxx + 2*gxyx*gupxy*gupxz*gupyy*gyyz
        + 6*gxyy*gupxy*gupyy*gupyz*gyzy + 2*gxyx*gupxx*gupxy*gupyz*gyzx
        + 6*gxyy*gupxz*gupyy*gupyz*gyzz + 4*gxyz*gupxx*gupxy*gupzz*gxzx
        + 2*gxyz*gupxy*gupxz*gupzz*gxzz + 4*gxyx*gupxy*gupxy*gupyy*gyyy
        + 2*gxyz*gupxz*gupyy*gupyz*gyyz;

    tmp_Gamxz = tmp_Gamyy + 4*gxyz*gupxy*gupyy*gupyz*gyyy
        + 2*gxyx*gupxz*gupxz*gupyy*gyzz + 2*gxyz*gupxz*gupyy*gupzz*gyzz
        + 2*gxyx*gupxy*gupxz*gupxz*gxzz + 4*gxyz*gupxy*gupyy*gupzz*gyzy
        + 2*gxzx*gupxy*gupxy*gupzz*gyzy + 2*gxyz*gupxx*gupxz*gupyz*gxzx
        + 4*gxyx*gupxz*gupxz*gupyz*gzzz + 4*gxzx*gupxy*gupxy*gupyz*gyyy
        + 2*gyyz*gupxy*gupyy*gupzz*gxzy + 2*gxyz*gupxy*gupxz*gupyz*gxzy
        + 2*gxyz*gupxx*gupyz*gupzz*gzzx + 4*gxyy*gupxy*gupyy*gupyy*gyyy
        + 2*gxyy*gupxx*gupyy*gupyy*gyyx + 2*gxyy*gupxx*gupyz*gupyz*gzzx
        + 2*gxyz*gupxy*gupyz*gupzz*gzzy + 2*gxyy*gupxz*gupyy*gupyy*gyyz
        + 4*gxyz*gupxz*gupyz*gupzz*gzzz + 2*gxxy*gupxx*gupxz*gupyy*gxyz;

    tmp_Gamyy = tmp_Gamxz + 2*gxzx*gupxy*gupxz*gupxz*gxyz
        + 2*gxxy*gupxy*gupxz*gupyy*gyzy + 4*gxxx*gupxx*gupxy*gupxz*gxzy
        + 2*gxxy*gupxy*gupxz*gupyy*gyyz + 2*gxxy*gupxx*gupxz*gupyy*gyzx
        + 2*gxxy*gupxx*gupxz*gupyz*gzzx + 2*gxzx*gupxy*gupxy*gupxz*gxyy
        + 2*gxxy*gupxx*gupxy*gupyz*gxzy + 2*gxyz*gupxy*gupxz*gupxz*gxxz
        + 2*gxxy*gupxx*gupxy*gupyy*gyyx + 2*gxyz*gupxx*gupyz*gupyz*gyzx
        + 4*gxyz*gupxz*gupxz*gupyz*gxzz + 2*gxxz*gupxx*gupxy*gupzz*gxzy
        + 2*gxxz*gupxx*gupxz*gupzz*gxzz + 2*gxxx*gupxx*gupxx*gupxy*gxxy
        + 2*gxxz*gupxx*gupxy*gupyz*gyyx + 2*gxxz*gupxy*gupxz*gupzz*gyzz
        + 2*gxxz*gupxx*gupxy*gupzz*gyzx;

    double TZ_rhs_val = tmp_Gamyy + 2*gxxz*gupxy*gupxz*gupyz*gyyz
        + 2*gxxz*gupxx*gupxz*gupyz*gyzx + 2*gxxz*gupxx*gupxz*gupzz*gzzx
        + 2*gxxz*gupxy*gupxz*gupyz*gyzy + 2*gxzx*gupxx*gupxy*gupzz*gxzy
        + 2*gxxz*gupxy*gupxz*gupzz*gzzy + 6*gxzx*gupxy*gupxz*gupzz*gyzz
        + 2*gxzx*gupxx*gupxy*gupzz*gyzx + 2*gxzx*gupxx*gupxy*gupyz*gyyx
        + 6*gxzx*gupxx*gupxz*gupzz*gxzz + 2*gxxx*gupxy*gupxy*gupxz*gyyz
        + 2*gxzx*gupxy*gupxz*gupzz*gzzy + 2*gxzx*gupxx*gupxz*gupyz*gyzx
        + 2*gxxx*gupxy*gupxz*gupxz*gzzy + 4*gxzy*gupxy*gupyz*gupyz*gyzy
        + 2*gxzy*gupxx*gupyz*gupyz*gyzx + 2*gxzz*gupxx*gupyz*gupyz*gyyx
        + 4*gxyx*gupxy*gupxy*gupxz*gyzx + 2*gxyx*gupxz*gupxz*gupyy*gzzy
        + 2*gxyy*gupxx*gupyz*gupyz*gxzz;

    // =====================================================================
    // Phase 7: Gam^i_,j derivatives, TZ_rhs continued
    // =====================================================================
    double Gamxx, Gamxy, Gamxz;
    double Gamyx, Gamyy_v, Gamyz;
    double Gamzx, Gamzy, Gamzz;
    deriv1_sh<GW>(ex, pex, Gamx, Gamxx, Gamxy, Gamxz, ANTI_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    deriv1_sh<GW>(ex, pex, Gamy, Gamyx, Gamyy_v, Gamyz, SYM_VAL, ANTI_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    deriv1_sh<GW>(ex, pex, Gamz, Gamzx, Gamzy, Gamzz, SYM_VAL, SYM_VAL, ANTI_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    TZ_rhs_val = chix * l_Gmxcon + chiy * l_Gmycon + chiz * l_Gmzcon
        + chin1 * (Gamxx + Gamyy_v + Gamzz -
            (TWO * (gupxz*gupyz*gyzxz + gupxx*gupyy*gxyxy + gupxy*gupyz*gxzyy
                + gupxx*gupxy*gxxxy + gupxx*gupxz*gxxxz + gupxx*gupxy*gxyxx
                + gupxx*gupyz*gxyxz + gupxx*gupxz*gxzxx + gupxx*gupyz*gxzxy
                + gupxx*gupzz*gxzxz + gupxy*gupxz*gxxyz + gupxy*gupyy*gxyyy
                + gupxy*gupyz*gxyyz + gupxy*gupxz*gxzxy + gupxy*gupzz*gxzyz
                + gupxy*gupxz*gxyxz + gupxz*gupyy*gxyyz + gupxz*gupyz*gxyzz
                + gupxz*gupyz*gxzyz + gupxz*gupzz*gxzzz + gupxy*gupyy*gyyxy
                + gupxy*gupyz*gyyxz + gupxy*gupxz*gyzxx + gupxy*gupyz*gyzxy
                + gupxy*gupzz*gyzxz + gupyy*gupyz*gyyyz + gupxz*gupyy*gyzxy
                + gupyy*gupyz*gyzyy + gupyy*gupzz*gyzyz + gupyz*gupzz*gyzzz
                + gupxz*gupyz*gzzxy + gupxz*gupzz*gzzxz + gupyz*gupzz*gzzyz
                + gupxy*gupxy*gxyxy + gupxz*gupxz*gxzxz + gupyz*gupyz*gyzyz)
            + gupxx*gupxx*gxxxx + gupxy*gupxy*gxxyy + gupxz*gupxz*gxxzz
            + gupxy*gupxy*gyyxx + gupyy*gupyy*gyyyy + gupyz*gupyz*gyyzz
            + gupxz*gupxz*gzzxx + gupyz*gupyz*gzzyy + gupzz*gupzz*gzzzz)
        + (l_gxx*Gamxa*Gamxa + l_gyy*Gamya*Gamya + l_gzz*Gamza*Gamza
            + TWO*(l_gxy*Gamxa*Gamya + l_gxz*Gamxa*Gamza + l_gyz*Gamya*Gamza))
        + TZ_rhs_val);

    // =====================================================================
    // Phase 8: Raise Aij indices (store in Rij)
    // =====================================================================
    double Rxx_v = gupxx*gupxx*l_Axx + gupxy*gupxy*l_Ayy + gupxz*gupxz*l_Azz
        + TWO*(gupxx*gupxy*l_Axy + gupxx*gupxz*l_Axz + gupxy*gupxz*l_Ayz);
    double Ryy_v = gupxy*gupxy*l_Axx + gupyy*gupyy*l_Ayy + gupyz*gupyz*l_Azz
        + TWO*(gupxy*gupyy*l_Axy + gupxy*gupyz*l_Axz + gupyy*gupyz*l_Ayz);
    double Rzz_v = gupxz*gupxz*l_Axx + gupyz*gupyz*l_Ayy + gupzz*gupzz*l_Azz
        + TWO*(gupxz*gupyz*l_Axy + gupxz*gupzz*l_Axz + gupyz*gupzz*l_Ayz);
    double Rxy_v = gupxx*gupxy*l_Axx + gupxy*gupyy*l_Ayy + gupxz*gupyz*l_Azz
        + (gupxx*gupyy + gupxy*gupxy)*l_Axy
        + (gupxx*gupyz + gupxz*gupxy)*l_Axz
        + (gupxy*gupyz + gupxz*gupyy)*l_Ayz;
    double Rxz_v = gupxx*gupxz*l_Axx + gupxy*gupyz*l_Ayy + gupxz*gupzz*l_Azz
        + (gupxx*gupyz + gupxy*gupxz)*l_Axy
        + (gupxx*gupzz + gupxz*gupxz)*l_Axz
        + (gupxy*gupzz + gupxz*gupyz)*l_Ayz;
    double Ryz_v = gupxy*gupxz*l_Axx + gupyy*gupyz*l_Ayy + gupyz*gupzz*l_Azz
        + (gupxy*gupyz + gupyy*gupxz)*l_Axy
        + (gupxy*gupzz + gupyz*gupxz)*l_Axz
        + (gupyy*gupzz + gupyz*gupyz)*l_Ayz;

    // =====================================================================
    // Phase 9: Gam^i_rhs (without shift terms)
    // =====================================================================
    double Lapx, Lapy, Lapz;
    deriv1_sh<GW>(ex, pex, Lap, Lapx, Lapy, Lapz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double Kx, Ky, Kz;
    deriv1_sh<GW>(ex, pex, trK, Kx, Ky, Kz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double fxx, fxy, fxz;
    deriv1_sh<GW>(ex, pex, TZ, fxx, fxy, fxz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    Kx = Kx + fxx / TWO;
    Ky = Ky + fxy / TWO;
    Kz = Kz + fxz / TWO;

    double Gamx_rhs_val = -TWO * (Lapx * Rxx_v + Lapy * Rxy_v + Lapz * Rxz_v)
        + TWO * alpn1 * (
            -F3o2/chin1 * (chix * Rxx_v + chiy * Rxy_v + chiz * Rxz_v)
            - gupxx * (F2o3 * Kx + EIGHT * PI_VAL * l_Sx)
            - gupxy * (F2o3 * Ky + EIGHT * PI_VAL * l_Sy)
            - gupxz * (F2o3 * Kz + EIGHT * PI_VAL * l_Sz)
            + Gamxxx * Rxx_v + Gamxyy * Ryy_v + Gamxzz * Rzz_v
            + TWO * (Gamxxy * Rxy_v + Gamxxz * Rxz_v + Gamxyz * Ryz_v));

    double Gamy_rhs_val = -TWO * (Lapx * Rxy_v + Lapy * Ryy_v + Lapz * Ryz_v)
        + TWO * alpn1 * (
            -F3o2/chin1 * (chix * Rxy_v + chiy * Ryy_v + chiz * Ryz_v)
            - gupxy * (F2o3 * Kx + EIGHT * PI_VAL * l_Sx)
            - gupyy * (F2o3 * Ky + EIGHT * PI_VAL * l_Sy)
            - gupyz * (F2o3 * Kz + EIGHT * PI_VAL * l_Sz)
            + Gamyxx * Rxx_v + Gamyyy * Ryy_v + Gamyzz * Rzz_v
            + TWO * (Gamyxy * Rxy_v + Gamyxz * Rxz_v + Gamyyz * Ryz_v));

    double Gamz_rhs_val = -TWO * (Lapx * Rxz_v + Lapy * Ryz_v + Lapz * Rzz_v)
        + TWO * alpn1 * (
            -F3o2/chin1 * (chix * Rxz_v + chiy * Ryz_v + chiz * Rzz_v)
            - gupxz * (F2o3 * Kx + EIGHT * PI_VAL * l_Sx)
            - gupyz * (F2o3 * Ky + EIGHT * PI_VAL * l_Sy)
            - gupzz * (F2o3 * Kz + EIGHT * PI_VAL * l_Sz)
            + Gamzxx * Rxx_v + Gamzyy * Ryy_v + Gamzzz * Rzz_v
            + TWO * (Gamzxy * Rxy_v + Gamzxz * Rxz_v + Gamzyz * Ryz_v));

    // =====================================================================
    // Phase 10: Second derivatives of shift, Gam^i_rhs continued
    // =====================================================================
    deriv2_sh<GW>(ex, pex, betax, gxxx, gxyx, gxzx, gyyx, gyzx, gzzx, ANTI_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);
    deriv2_sh<GW>(ex, pex, betay, gxxy, gxyy, gxzy, gyyy, gyzy, gzzy, SYM_VAL, ANTI_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);
    deriv2_sh<GW>(ex, pex, betaz, gxxz, gxyz, gxzz, gyyz, gyzz, gzzz, SYM_VAL, SYM_VAL, ANTI_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);

    fxx = gxxx + gxyy + gxzz;
    fxy = gxyx + gyyy + gyzz;
    fxz = gxzx + gyzy + gzzz;

    Gamx_rhs_val += F2o3 * Gamxa * div_beta
        - Gamxa * betaxx - Gamya * betaxy - Gamza * betaxz
        + F1o3 * (gupxx * fxx + gupxy * fxy + gupxz * fxz)
        + gupxx * gxxx + gupyy * gyyx + gupzz * gzzx
        + TWO * (gupxy * gxyx + gupxz * gxzx + gupyz * gyzx);

    Gamy_rhs_val += F2o3 * Gamya * div_beta
        - Gamxa * betayx - Gamya * betayy - Gamza * betayz
        + F1o3 * (gupxy * fxx + gupyy * fxy + gupyz * fxz)
        + gupxx * gxxy + gupyy * gyyy + gupzz * gzzy
        + TWO * (gupxy * gxyy + gupxz * gxzy + gupyz * gyzy);

    Gamz_rhs_val += F2o3 * Gamza * div_beta
        - Gamxa * betazx - Gamya * betazy - Gamza * betazz
        + F1o3 * (gupxz * fxx + gupyz * fxy + gupzz * fxz)
        + gupxx * gxxz + gupyy * gyyz + gupzz * gzzz
        + TWO * (gupxy * gxyz + gupxz * gxzz + gupyz * gyzz);

    // =====================================================================
    // Phase 11: First kind of connection, Ricci tensor
    // =====================================================================
    // Overwrite gij,k with first-kind connection: g_ij Gam^k_ij
    gxxx = l_gxx * Gamxxx + l_gxy * Gamyxx + l_gxz * Gamzxx;
    gxyx = l_gxx * Gamxxy + l_gxy * Gamyxy + l_gxz * Gamzxy;
    gxzx = l_gxx * Gamxxz + l_gxy * Gamyxz + l_gxz * Gamzxz;
    gyyx = l_gxx * Gamxyy + l_gxy * Gamyyy + l_gxz * Gamzyy;
    gyzx = l_gxx * Gamxyz + l_gxy * Gamyyz + l_gxz * Gamzyz;
    gzzx = l_gxx * Gamxzz + l_gxy * Gamyzz + l_gxz * Gamzzz;

    gxxy = l_gxy * Gamxxx + l_gyy * Gamyxx + l_gyz * Gamzxx;
    gxyy = l_gxy * Gamxxy + l_gyy * Gamyxy + l_gyz * Gamzxy;
    gxzy = l_gxy * Gamxxz + l_gyy * Gamyxz + l_gyz * Gamzxz;
    gyyy = l_gxy * Gamxyy + l_gyy * Gamyyy + l_gyz * Gamzyy;
    gyzy = l_gxy * Gamxyz + l_gyy * Gamyyz + l_gyz * Gamzyz;
    gzzy = l_gxy * Gamxzz + l_gyy * Gamyzz + l_gyz * Gamzzz;

    gxxz = l_gxz * Gamxxx + l_gyz * Gamyxx + l_gzz * Gamzxx;
    gxyz = l_gxz * Gamxxy + l_gyz * Gamyxy + l_gzz * Gamzxy;
    gxzz = l_gxz * Gamxxz + l_gyz * Gamyxz + l_gzz * Gamzxz;
    gyyz = l_gxz * Gamxyy + l_gyz * Gamyyy + l_gzz * Gamzyy;
    gyzz = l_gxz * Gamxyz + l_gyz * Gamyyz + l_gzz * Gamzyz;
    gzzz = l_gxz * Gamxzz + l_gyz * Gamyzz + l_gzz * Gamzzz;

    // Ricci tensor
    double Rxx_tmp = gupxx*gxxxx + gupyy*gxxyy + gupzz*gxxzz
        + (gupxy*gxxxy + gupxz*gxxxz + gupyz*gxxyz) * TWO;
    double Ryy_tmp = gupxx*gyyxx + gupyy*gyyyy + gupzz*gyyzz
        + (gupxy*gyyxy + gupxz*gyyxz + gupyz*gyyyz) * TWO;
    double Rzz_tmp = gupxx*gzzxx + gupyy*gzzyy + gupzz*gzzzz
        + (gupxy*gzzxy + gupxz*gzzxz + gupyz*gzzyz) * TWO;
    double Rxy_tmp = gupxx*gxyxx + gupyy*gxyyy + gupzz*gxyzz
        + (gupxy*gxyxy + gupxz*gxyxz + gupyz*gxyyz) * TWO;
    double Rxz_tmp = gupxx*gxzxx + gupyy*gxzyy + gupzz*gxzzz
        + (gupxy*gxzxy + gupxz*gxzxz + gupyz*gxzyz) * TWO;
    double Ryz_tmp = gupxx*gyzxx + gupyy*gyzyy + gupzz*gyzzz
        + (gupxy*gyzxy + gupxz*gyzxz + gupyz*gyzyz) * TWO;

    // Rxx
    Rxx_v = -HALF * Rxx_tmp
        + l_gxx * Gamxx + l_gxy * Gamyx + l_gxz * Gamzx
        + Gamxa * gxxx + Gamya * gxyx + Gamza * gxzx
        + gupxx * (TWO*(Gamxxx*gxxx + Gamyxx*gxyx + Gamzxx*gxzx) + Gamxxx*gxxx + Gamyxx*gxxy + Gamzxx*gxxz)
        + gupxy * (TWO*(Gamxxx*gxyx + Gamyxx*gyyx + Gamzxx*gyzx + Gamxxy*gxxx + Gamyxy*gxyx + Gamzxy*gxzx)
            + Gamxxy*gxxx + Gamyxy*gxxy + Gamzxy*gxxz + Gamxxx*gxyx + Gamyxx*gxyy + Gamzxx*gxyz)
        + gupxz * (TWO*(Gamxxx*gxzx + Gamyxx*gyzx + Gamzxx*gzzx + Gamxxz*gxxx + Gamyxz*gxyx + Gamzxz*gxzx)
            + Gamxxz*gxxx + Gamyxz*gxxy + Gamzxz*gxxz + Gamxxx*gxzx + Gamyxx*gxzy + Gamzxx*gxzz)
        + gupyy * (TWO*(Gamxxy*gxyx + Gamyxy*gyyx + Gamzxy*gyzx) + Gamxxy*gxyx + Gamyxy*gxyy + Gamzxy*gxyz)
        + gupyz * (TWO*(Gamxxy*gxzx + Gamyxy*gyzx + Gamzxy*gzzx + Gamxxz*gxyx + Gamyxz*gyyx + Gamzxz*gyzx)
            + Gamxxz*gxyx + Gamyxz*gxyy + Gamzxz*gxyz + Gamxxy*gxzx + Gamyxy*gxzy + Gamzxy*gxzz)
        + gupzz * (TWO*(Gamxxz*gxzx + Gamyxz*gyzx + Gamzxz*gzzx) + Gamxxz*gxzx + Gamyxz*gxzy + Gamzxz*gxzz);

    // Ryy
    Ryy_v = -HALF * Ryy_tmp
        + l_gxy * Gamxy + l_gyy * Gamyy_v + l_gyz * Gamzy
        + Gamxa * gxyy + Gamya * gyyy + Gamza * gyzy
        + gupxx * (TWO*(Gamxxy*gxxy + Gamyxy*gxyy + Gamzxy*gxzy) + Gamxxy*gxyx + Gamyxy*gxyy + Gamzxy*gxyz)
        + gupxy * (TWO*(Gamxxy*gxyy + Gamyxy*gyyy + Gamzxy*gyzy + Gamxyy*gxxy + Gamyyy*gxyy + Gamzyy*gxzy)
            + Gamxyy*gxyx + Gamyyy*gxyy + Gamzyy*gxyz + Gamxxy*gyyx + Gamyxy*gyyy + Gamzxy*gyyz)
        + gupxz * (TWO*(Gamxxy*gxzy + Gamyxy*gyzy + Gamzxy*gzzy + Gamxyz*gxxy + Gamyyz*gxyy + Gamzyz*gxzy)
            + Gamxyz*gxyx + Gamyyz*gxyy + Gamzyz*gxyz + Gamxxy*gyzx + Gamyxy*gyzy + Gamzxy*gyzz)
        + gupyy * (TWO*(Gamxyy*gxyy + Gamyyy*gyyy + Gamzyy*gyzy) + Gamxyy*gyyx + Gamyyy*gyyy + Gamzyy*gyyz)
        + gupyz * (TWO*(Gamxyy*gxzy + Gamyyy*gyzy + Gamzyy*gzzy + Gamxyz*gxyy + Gamyyz*gyyy + Gamzyz*gyzy)
            + Gamxyz*gyyx + Gamyyz*gyyy + Gamzyz*gyyz + Gamxyy*gyzx + Gamyyy*gyzy + Gamzyy*gyzz)
        + gupzz * (TWO*(Gamxyz*gxzy + Gamyyz*gyzy + Gamzyz*gzzy) + Gamxyz*gyzx + Gamyyz*gyzy + Gamzyz*gyzz);

    // Rzz
    Rzz_v = -HALF * Rzz_tmp
        + l_gxz * Gamxz + l_gyz * Gamyz + l_gzz * Gamzz
        + Gamxa * gxzz + Gamya * gyzz + Gamza * gzzz
        + gupxx * (TWO*(Gamxxz*gxxz + Gamyxz*gxyz + Gamzxz*gxzz) + Gamxxz*gxzx + Gamyxz*gxzy + Gamzxz*gxzz)
        + gupxy * (TWO*(Gamxxz*gxyz + Gamyxz*gyyz + Gamzxz*gyzz + Gamxyz*gxxz + Gamyyz*gxyz + Gamzyz*gxzz)
            + Gamxyz*gxzx + Gamyyz*gxzy + Gamzyz*gxzz + Gamxxz*gyzx + Gamyxz*gyzy + Gamzxz*gyzz)
        + gupxz * (TWO*(Gamxxz*gxzz + Gamyxz*gyzz + Gamzxz*gzzz + Gamxzz*gxxz + Gamyzz*gxyz + Gamzzz*gxzz)
            + Gamxzz*gxzx + Gamyzz*gxzy + Gamzzz*gxzz + Gamxxz*gzzx + Gamyxz*gzzy + Gamzxz*gzzz)
        + gupyy * (TWO*(Gamxyz*gxyz + Gamyyz*gyyz + Gamzyz*gyzz) + Gamxyz*gyzx + Gamyyz*gyzy + Gamzyz*gyzz)
        + gupyz * (TWO*(Gamxyz*gxzz + Gamyyz*gyzz + Gamzyz*gzzz + Gamxzz*gxyz + Gamyzz*gyyz + Gamzzz*gyzz)
            + Gamxzz*gyzx + Gamyzz*gyzy + Gamzzz*gyzz + Gamxyz*gzzx + Gamyyz*gzzy + Gamzyz*gzzz)
        + gupzz * (TWO*(Gamxzz*gxzz + Gamyzz*gyzz + Gamzzz*gzzz) + Gamxzz*gzzx + Gamyzz*gzzy + Gamzzz*gzzz);

    // Rxy
    Rxy_v = HALF * (-Rxy_tmp
        + l_gxx * Gamxy + l_gxy * Gamyy_v + l_gxz * Gamzy
        + l_gxy * Gamxx + l_gyy * Gamyx + l_gyz * Gamzx
        + Gamxa * gxyx + Gamya * gyyx + Gamza * gyzx
        + Gamxa * gxxy + Gamya * gxyy + Gamza * gxzy)
        + gupxx * (Gamxxx*gxxy + Gamyxx*gxyy + Gamzxx*gxzy + Gamxxy*gxxx + Gamyxy*gxyx + Gamzxy*gxzx + Gamxxx*gxyx + Gamyxx*gxyy + Gamzxx*gxyz)
        + gupxy * (Gamxxx*gxyy + Gamyxx*gyyy + Gamzxx*gyzy + Gamxxy*gxyx + Gamyxy*gyyx + Gamzxy*gyzx + Gamxxy*gxyx + Gamyxy*gxyy + Gamzxy*gxyz + Gamxxy*gxxy + Gamyxy*gxyy + Gamzxy*gxzy + Gamxyy*gxxx + Gamyyy*gxyx + Gamzyy*gxzx + Gamxxx*gyyx + Gamyxx*gyyy + Gamzxx*gyyz)
        + gupxz * (Gamxxx*gxzy + Gamyxx*gyzy + Gamzxx*gzzy + Gamxxy*gxzx + Gamyxy*gyzx + Gamzxy*gzzx + Gamxxz*gxyx + Gamyxz*gxyy + Gamzxz*gxyz + Gamxxz*gxxy + Gamyxz*gxyy + Gamzxz*gxzy + Gamxyz*gxxx + Gamyyz*gxyx + Gamzyz*gxzx + Gamxxx*gyzx + Gamyxx*gyzy + Gamzxx*gyzz)
        + gupyy * (Gamxxy*gxyy + Gamyxy*gyyy + Gamzxy*gyzy + Gamxyy*gxyx + Gamyyy*gyyx + Gamzyy*gyzx + Gamxxy*gyyx + Gamyxy*gyyy + Gamzxy*gyyz)
        + gupyz * (Gamxxy*gxzy + Gamyxy*gyzy + Gamzxy*gzzy + Gamxyy*gxzx + Gamyyy*gyzx + Gamzyy*gzzx + Gamxxz*gyyx + Gamyxz*gyyy + Gamzxz*gyyz + Gamxxz*gxyy + Gamyxz*gyyy + Gamzxz*gyzy + Gamxyz*gxyx + Gamyyz*gyyx + Gamzyz*gyzx + Gamxxy*gyzx + Gamyxy*gyzy + Gamzxy*gyzz)
        + gupzz * (Gamxxz*gxzy + Gamyxz*gyzy + Gamzxz*gzzy + Gamxyz*gxzx + Gamyyz*gyzx + Gamzyz*gzzx + Gamxxz*gyzx + Gamyxz*gyzy + Gamzxz*gyzz);

    // Rxz
    Rxz_v = HALF * (-Rxz_tmp
        + l_gxx * Gamxz + l_gxy * Gamyz + l_gxz * Gamzz
        + l_gxz * Gamxx + l_gyz * Gamyx + l_gzz * Gamzx
        + Gamxa * gxzx + Gamya * gyzx + Gamza * gzzx
        + Gamxa * gxxz + Gamya * gxyz + Gamza * gxzz)
        + gupxx * (Gamxxx*gxxz + Gamyxx*gxyz + Gamzxx*gxzz + Gamxxz*gxxx + Gamyxz*gxyx + Gamzxz*gxzx + Gamxxx*gxzx + Gamyxx*gxzy + Gamzxx*gxzz)
        + gupxy * (Gamxxx*gxyz + Gamyxx*gyyz + Gamzxx*gyzz + Gamxxz*gxyx + Gamyxz*gyyx + Gamzxz*gyzx + Gamxxy*gxzx + Gamyxy*gxzy + Gamzxy*gxzz + Gamxxy*gxxz + Gamyxy*gxyz + Gamzxy*gxzz + Gamxyz*gxxx + Gamyyz*gxyx + Gamzyz*gxzx + Gamxxx*gyzx + Gamyxx*gyzy + Gamzxx*gyzz)
        + gupxz * (Gamxxx*gxzz + Gamyxx*gyzz + Gamzxx*gzzz + Gamxxz*gxzx + Gamyxz*gyzx + Gamzxz*gzzx + Gamxxz*gxzx + Gamyxz*gxzy + Gamzxz*gxzz + Gamxxz*gxxz + Gamyxz*gxyz + Gamzxz*gxzz + Gamxzz*gxxx + Gamyzz*gxyx + Gamzzz*gxzx + Gamxxx*gzzx + Gamyxx*gzzy + Gamzxx*gzzz)
        + gupyy * (Gamxxy*gxyz + Gamyxy*gyyz + Gamzxy*gyzz + Gamxyz*gxyx + Gamyyz*gyyx + Gamzyz*gyzx + Gamxxy*gyzx + Gamyxy*gyzy + Gamzxy*gyzz)
        + gupyz * (Gamxxy*gxzz + Gamyxy*gyzz + Gamzxy*gzzz + Gamxyz*gxzx + Gamyyz*gyzx + Gamzyz*gzzx + Gamxxz*gyzx + Gamyxz*gyzy + Gamzxz*gyzz + Gamxxz*gxyz + Gamyxz*gyyz + Gamzxz*gyzz + Gamxzz*gxyx + Gamyzz*gyyx + Gamzzz*gyzx + Gamxxy*gzzx + Gamyxy*gzzy + Gamzxy*gzzz)
        + gupzz * (Gamxxz*gxzz + Gamyxz*gyzz + Gamzxz*gzzz + Gamxzz*gxzx + Gamyzz*gyzx + Gamzzz*gzzx + Gamxxz*gzzx + Gamyxz*gzzy + Gamzxz*gzzz);

    // Ryz
    Ryz_v = HALF * (-Ryz_tmp
        + l_gxy * Gamxz + l_gyy * Gamyz + l_gyz * Gamzz
        + l_gxz * Gamxy + l_gyz * Gamyy_v + l_gzz * Gamzy
        + Gamxa * gxzy + Gamya * gyzy + Gamza * gzzy
        + Gamxa * gxyz + Gamya * gyyz + Gamza * gyzz)
        + gupxx * (Gamxxy*gxxz + Gamyxy*gxyz + Gamzxy*gxzz + Gamxxz*gxxy + Gamyxz*gxyy + Gamzxz*gxzy + Gamxxy*gxzx + Gamyxy*gxzy + Gamzxy*gxzz)
        + gupxy * (Gamxxy*gxyz + Gamyxy*gyyz + Gamzxy*gyzz + Gamxxz*gxyy + Gamyxz*gyyy + Gamzxz*gyzy + Gamxyy*gxzx + Gamyyy*gxzy + Gamzyy*gxzz + Gamxyy*gxxz + Gamyyy*gxyz + Gamzyy*gxzz + Gamxyz*gxxy + Gamyyz*gxyy + Gamzyz*gxzy + Gamxxy*gyzx + Gamyxy*gyzy + Gamzxy*gyzz)
        + gupxz * (Gamxxy*gxzz + Gamyxy*gyzz + Gamzxy*gzzz + Gamxxz*gxzy + Gamyxz*gyzy + Gamzxz*gzzy + Gamxyz*gxzx + Gamyyz*gxzy + Gamzyz*gxzz + Gamxyz*gxxz + Gamyyz*gxyz + Gamzyz*gxzz + Gamxzz*gxxy + Gamyzz*gxyy + Gamzzz*gxzy + Gamxxy*gzzx + Gamyxy*gzzy + Gamzxy*gzzz)
        + gupyy * (Gamxyy*gxyz + Gamyyy*gyyz + Gamzyy*gyzz + Gamxyz*gxyy + Gamyyz*gyyy + Gamzyz*gyzy + Gamxyy*gyzx + Gamyyy*gyzy + Gamzyy*gyzz)
        + gupyz * (Gamxyy*gxzz + Gamyyy*gyzz + Gamzyy*gzzz + Gamxyz*gxzy + Gamyyz*gyzy + Gamzyz*gzzy + Gamxyz*gyzx + Gamyyz*gyzy + Gamzyz*gyzz + Gamxyz*gxyz + Gamyyz*gyyz + Gamzyz*gyzz + Gamxzz*gxyy + Gamyzz*gyyy + Gamzzz*gyzy + Gamxyy*gzzx + Gamyyy*gzzy + Gamzyy*gzzz)
        + gupzz * (Gamxyz*gxzz + Gamyyz*gyzz + Gamzyz*gzzz + Gamxzz*gxzy + Gamyzz*gyzy + Gamzzz*gzzy + Gamxyz*gzzx + Gamyyz*gzzy + Gamzyz*gzzz);

    // =====================================================================
    // Phase 12: Covariant D^2 chi, chi correction to Ricci
    // =====================================================================
    double fyy, fyz, fzz;
    deriv2_sh<GW>(ex, pex, chi, fxx, fxy, fxz, fyy, fyz, fzz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);

    fxx = fxx - Gamxxx * chix - Gamyxx * chiy - Gamzxx * chiz;
    fxy = fxy - Gamxxy * chix - Gamyxy * chiy - Gamzxy * chiz;
    fxz = fxz - Gamxxz * chix - Gamyxz * chiy - Gamzxz * chiz;
    fyy = fyy - Gamxyy * chix - Gamyyy * chiy - Gamzyy * chiz;
    fyz = fyz - Gamxyz * chix - Gamyyz * chiy - Gamzyz * chiz;
    fzz = fzz - Gamxzz * chix - Gamyzz * chiy - Gamzzz * chiz;

    double Dchi_Dchi = gupxx * (fxx - F3o2/chin1 * chix * chix)
        + gupyy * (fyy - F3o2/chin1 * chiy * chiy)
        + gupzz * (fzz - F3o2/chin1 * chiz * chiz)
        + TWO * gupxy * (fxy - F3o2/chin1 * chix * chiy)
        + TWO * gupxz * (fxz - F3o2/chin1 * chix * chiz)
        + TWO * gupyz * (fyz - F3o2/chin1 * chiy * chiz);

    fxx = Rxx_v + (fxx - chix*chix/chin1/TWO + l_gxx * Dchi_Dchi) / chin1 / TWO;
    fyy = Ryy_v + (fyy - chiy*chiy/chin1/TWO + l_gyy * Dchi_Dchi) / chin1 / TWO;
    fzz = Rzz_v + (fzz - chiz*chiz/chin1/TWO + l_gzz * Dchi_Dchi) / chin1 / TWO;
    fxy = Rxy_v + (fxy - chix*chiy/chin1/TWO + l_gxy * Dchi_Dchi) / chin1 / TWO;
    fxz = Rxz_v + (fxz - chix*chiz/chin1/TWO + l_gxz * Dchi_Dchi) / chin1 / TWO;
    fyz = Ryz_v + (fyz - chiy*chiz/chin1/TWO + l_gyz * Dchi_Dchi) / chin1 / TWO;

    // R/chi in Hcon
    double Hcon_val = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    Rxx_v = fxx; Ryy_v = fyy; Rzz_v = fzz;
    Rxy_v = fxy; Rxz_v = fxz; Ryz_v = fyz;

    // =====================================================================
    // Phase 13: Physical connection (chi correction)
    // =====================================================================
    double dchidx = (gupxx * chix + gupxy * chiy + gupxz * chiz) / chin1;
    double dchidy = (gupxy * chix + gupyy * chiy + gupyz * chiz) / chin1;
    double dchidz = (gupxz * chix + gupyz * chiy + gupzz * chiz) / chin1;

    Gamxxx = Gamxxx - ((chix + chix)/chin1 - l_gxx * dchidx) * HALF;
    Gamyxx = Gamyxx - (                - l_gxx * dchidy) * HALF;
    Gamzxx = Gamzxx - (                - l_gxx * dchidz) * HALF;
    Gamxyy = Gamxyy - (                - l_gyy * dchidx) * HALF;
    Gamyyy = Gamyyy - ((chiy + chiy)/chin1 - l_gyy * dchidy) * HALF;
    Gamzyy = Gamzyy - (                - l_gyy * dchidz) * HALF;
    Gamxzz = Gamxzz - (                - l_gzz * dchidx) * HALF;
    Gamyzz = Gamyzz - (                - l_gzz * dchidy) * HALF;
    Gamzzz = Gamzzz - ((chiz + chiz)/chin1 - l_gzz * dchidz) * HALF;
    Gamxxy = Gamxxy - (chiy/chin1     - l_gxy * dchidx) * HALF;
    Gamyxy = Gamyxy - (chix/chin1     - l_gxy * dchidy) * HALF;
    Gamzxy = Gamzxy - (                - l_gxy * dchidz) * HALF;
    Gamxxz = Gamxxz - (chiz/chin1     - l_gxz * dchidx) * HALF;
    Gamyxz = Gamyxz - (                - l_gxz * dchidy) * HALF;
    Gamzxz = Gamzxz - (chix/chin1     - l_gxz * dchidz) * HALF;
    Gamxyz = Gamxyz - (                - l_gyz * dchidx) * HALF;
    Gamyyz = Gamyyz - (chiz/chin1     - l_gyz * dchidy) * HALF;
    Gamzyz = Gamzyz - (chiy/chin1     - l_gyz * dchidz) * HALF;

    // =====================================================================
    // Phase 14: Covariant D^2 Lap, Aij_rhs, trK_rhs
    // =====================================================================
    deriv2_sh<GW>(ex, pex, Lap, fxx, fxy, fxz, fyy, fyz, fzz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);

    fxx = fxx - Gamxxx * Lapx - Gamyxx * Lapy - Gamzxx * Lapz;
    fxy = fxy - Gamxxy * Lapx - Gamyxy * Lapy - Gamzxy * Lapz;
    fxz = fxz - Gamxxz * Lapx - Gamyxz * Lapy - Gamzxz * Lapz;
    fyy = fyy - Gamxyy * Lapx - Gamyyy * Lapy - Gamzyy * Lapz;
    fyz = fyz - Gamxyz * Lapx - Gamyyz * Lapy - Gamzyz * Lapz;
    fzz = fzz - Gamxzz * Lapx - Gamyzz * Lapy - Gamzzz * Lapz;

    double D2Lap = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    fxx = EIGHT * PI_VAL * alpn1 * l_Sxx + fxx;
    fxy = EIGHT * PI_VAL * alpn1 * l_Sxy + fxy;
    fxz = EIGHT * PI_VAL * alpn1 * l_Sxz + fxz;
    fyy = EIGHT * PI_VAL * alpn1 * l_Syy + fyy;
    fyz = EIGHT * PI_VAL * alpn1 * l_Syz + fyz;
    fzz = EIGHT * PI_VAL * alpn1 * l_Szz + fzz;

    double f_trace = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);
    f_trace = F1o3 * (Hcon_val * alpn1 - f_trace);

    fxx = alpn1 * Rxx_v - fxx;
    fxy = alpn1 * Rxy_v - fxy;
    fxz = alpn1 * Rxz_v - fxz;
    fyy = alpn1 * Ryy_v - fyy;
    fyz = alpn1 * Ryz_v - fyz;
    fzz = alpn1 * Rzz_v - fzz;

    double Axx_rhs_val = fxx - l_gxx * f_trace;
    double Ayy_rhs_val = fyy - l_gyy * f_trace;
    double Azz_rhs_val = fzz - l_gzz * f_trace;
    double Axy_rhs_val = fxy - l_gxy * f_trace;
    double Axz_rhs_val = fxz - l_gxz * f_trace;
    double Ayz_rhs_val = fyz - l_gyz * f_trace;

    // A_il A^l_j
    double fA_xx = gupxx*l_Axx*l_Axx + gupyy*l_Axy*l_Axy + gupzz*l_Axz*l_Axz
        + TWO*(gupxy*l_Axx*l_Axy + gupxz*l_Axx*l_Axz + gupyz*l_Axy*l_Axz);
    double fA_yy = gupxx*l_Axy*l_Axy + gupyy*l_Ayy*l_Ayy + gupzz*l_Ayz*l_Ayz
        + TWO*(gupxy*l_Axy*l_Ayy + gupxz*l_Axy*l_Ayz + gupyz*l_Ayy*l_Ayz);
    double fA_zz = gupxx*l_Axz*l_Axz + gupyy*l_Ayz*l_Ayz + gupzz*l_Azz*l_Azz
        + TWO*(gupxy*l_Axz*l_Ayz + gupxz*l_Axz*l_Azz + gupyz*l_Ayz*l_Azz);
    double fA_xy = gupxx*l_Axx*l_Axy + gupyy*l_Axy*l_Ayy + gupzz*l_Axz*l_Ayz
        + gupxy*(l_Axx*l_Ayy + l_Axy*l_Axy)
        + gupxz*(l_Axx*l_Ayz + l_Axz*l_Axy)
        + gupyz*(l_Axy*l_Ayz + l_Axz*l_Ayy);
    double fA_xz = gupxx*l_Axx*l_Axz + gupyy*l_Axy*l_Ayz + gupzz*l_Axz*l_Azz
        + gupxy*(l_Axx*l_Ayz + l_Axy*l_Axz)
        + gupxz*(l_Axx*l_Azz + l_Axz*l_Axz)
        + gupyz*(l_Axy*l_Azz + l_Axz*l_Ayz);
    double fA_yz = gupxx*l_Axy*l_Axz + gupyy*l_Ayy*l_Ayz + gupzz*l_Ayz*l_Azz
        + gupxy*(l_Axy*l_Ayz + l_Ayy*l_Axz)
        + gupxz*(l_Axy*l_Azz + l_Ayz*l_Axz)
        + gupyz*(l_Ayy*l_Azz + l_Ayz*l_Ayz);

    double chin1_val = chin1;
    D2Lap = chin1_val * D2Lap;

    Axx_rhs_val = chin1_val * Axx_rhs_val + alpn1 * (trKd * l_Axx - TWO * fA_xx)
        + TWO * (l_Axx * betaxx + l_Axy * betayx + l_Axz * betazx) - F2o3 * l_Axx * div_beta;
    Ayy_rhs_val = chin1_val * Ayy_rhs_val + alpn1 * (trKd * l_Ayy - TWO * fA_yy)
        + TWO * (l_Axy * betaxy + l_Ayy * betayy + l_Ayz * betazy) - F2o3 * l_Ayy * div_beta;
    Azz_rhs_val = chin1_val * Azz_rhs_val + alpn1 * (trKd * l_Azz - TWO * fA_zz)
        + TWO * (l_Axz * betaxz + l_Ayz * betayz + l_Azz * betazz) - F2o3 * l_Azz * div_beta;
    Axy_rhs_val = chin1_val * Axy_rhs_val + alpn1 * (trKd * l_Axy - TWO * fA_xy)
        + l_Axx * betaxy + l_Axz * betazy + l_Ayy * betayx + l_Ayz * betazx
        + F1o3 * l_Axy * div_beta - l_Axy * betazz;
    Ayz_rhs_val = chin1_val * Ayz_rhs_val + alpn1 * (trKd * l_Ayz - TWO * fA_yz)
        + l_Axy * betaxz + l_Ayy * betayz + l_Axz * betaxy + l_Azz * betazy
        + F1o3 * l_Ayz * div_beta - l_Ayz * betaxx;
    Axz_rhs_val = chin1_val * Axz_rhs_val + alpn1 * (trKd * l_Axz - TWO * fA_xz)
        + l_Axx * betaxz + l_Axy * betayz + l_Ayz * betayx + l_Azz * betazx
        + F1o3 * l_Axz * div_beta - l_Axz * betayy;

    // trK_rhs
    double S_trace = chin1_val * (gupxx*l_Sxx + gupyy*l_Syy + gupzz*l_Szz
        + TWO*(gupxy*l_Sxy + gupxz*l_Sxz + gupyz*l_Syz));

    double trK_rhs_val = -D2Lap + alpn1 * (F1o3 * trKd * trKd
        + gupxx*fA_xx + gupyy*fA_yy + gupzz*fA_zz
        + TWO*(gupxy*fA_xy + gupxz*fA_xz + gupyz*fA_yz)
        + FOUR * PI_VAL * (l_rho + S_trace));

    // =====================================================================
    // Phase 15: Gauge RHS
    // =====================================================================
    double Lap_rhs_val = -TWO * alpn1 * l_trK;
    double betax_rhs_val, betay_rhs_val, betaz_rhs_val;
    double dtSfx_rhs_val, dtSfy_rhs_val, dtSfz_rhs_val;

#if (GAUGE == 0)
    betax_rhs_val = FF * l_dtSfx;
    betay_rhs_val = FF * l_dtSfy;
    betaz_rhs_val = FF * l_dtSfz;
    dtSfx_rhs_val = Gamx_rhs_val - eta_val * l_dtSfx;
    dtSfy_rhs_val = Gamy_rhs_val - eta_val * l_dtSfy;
    dtSfz_rhs_val = Gamz_rhs_val - eta_val * l_dtSfz;
#elif (GAUGE == 1)
    betax_rhs_val = l_Gamx - eta_val * l_betax;
    betay_rhs_val = l_Gamy - eta_val * l_betay;
    betaz_rhs_val = l_Gamz - eta_val * l_betaz;
    dtSfx_rhs_val = ZEO;
    dtSfy_rhs_val = ZEO;
    dtSfz_rhs_val = ZEO;
#else
    // GAUGE>=2: Z4c_rhs_ss does not define beta/dtSf RHS for these gauges.
    // Match CPU Fortran behavior: variables remain at their zero-initialized state
    // (Fortran intent(out) arrays are zero-initialized by most compilers).
    betax_rhs_val = ZEO;
    betay_rhs_val = ZEO;
    betaz_rhs_val = ZEO;
    dtSfx_rhs_val = ZEO;
    dtSfy_rhs_val = ZEO;
    dtSfz_rhs_val = ZEO;
#endif

    // =====================================================================
    // Phase 16: Hamiltonian constraint Hcon
    // =====================================================================
    // A_ij A^ij (computed from the raised A^i_j stored in R*_v)
    // Actually recompute: R*_v was overwritten. Use fA_** which stores A_il A^l_j
    double AijAij = gupxx * fA_xx + gupyy * fA_yy + gupzz * fA_zz
        + TWO * (gupxy * fA_xy + gupxz * fA_xz + gupyz * fA_yz);

    Hcon_val = chin1_val * Hcon_val + F2o3 * trKd * trKd - AijAij - F16 * PI_VAL * l_rho;

    // =====================================================================
    // Phase 17: Momentum constraint Mxcon, Mycon, Mzcon
    // =====================================================================
    // Need D_i A_jk derivatives of Aij
    double gxxx_A, gxyx_A, gxzx_A, gyyx_A, gyzx_A, gzzx_A;
    double gxxy_A, gxyy_A, gxzy_A, gyyy_A, gyzy_A, gzzy_A;
    double gxxz_A, gxyz_A, gxzz_A, gyyz_A, gyzz_A, gzzz_A;

    double dAxx_dx, dAxx_dy, dAxx_dz;
    deriv1_sh<GW>(ex, pex, Axx, dAxx_dx, dAxx_dy, dAxx_dz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double dAxy_dx, dAxy_dy, dAxy_dz;
    deriv1_sh<GW>(ex, pex, Axy, dAxy_dx, dAxy_dy, dAxy_dz, ANTI_VAL, ANTI_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double dAxz_dx, dAxz_dy, dAxz_dz;
    deriv1_sh<GW>(ex, pex, Axz, dAxz_dx, dAxz_dy, dAxz_dz, ANTI_VAL, SYM_VAL, ANTI_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double dAyy_dx, dAyy_dy, dAyy_dz;
    deriv1_sh<GW>(ex, pex, Ayy, dAyy_dx, dAyy_dy, dAyy_dz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double dAyz_dx, dAyz_dy, dAyz_dz;
    deriv1_sh<GW>(ex, pex, Ayz, dAyz_dx, dAyz_dy, dAyz_dz, SYM_VAL, ANTI_VAL, ANTI_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double dAzz_dx, dAzz_dy, dAzz_dz;
    deriv1_sh<GW>(ex, pex, Azz, dAzz_dx, dAzz_dy, dAzz_dz, SYM_VAL, SYM_VAL, SYM_VAL, crho, sigma, R, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    // D_i A_jk - 1/chi d_i chi * A_jk
    // x-direction
    gxxx_A = dAxx_dx - (Gamxxx*l_Axx + Gamyxx*l_Axy + Gamzxx*l_Axz
                       + Gamxxx*l_Axx + Gamyxx*l_Axy + Gamzxx*l_Axz) - chix*l_Axx/chin1;
    gxyx_A = dAxy_dx - (Gamxxy*l_Axx + Gamyxy*l_Axy + Gamzxy*l_Axz
                       + Gamxxx*l_Axy + Gamyxx*l_Ayy + Gamzxx*l_Ayz) - chix*l_Axy/chin1;
    gxzx_A = dAxz_dx - (Gamxxz*l_Axx + Gamyxz*l_Axy + Gamzxz*l_Axz
                       + Gamxxx*l_Axz + Gamyxx*l_Ayz + Gamzxx*l_Azz) - chix*l_Axz/chin1;
    gyyx_A = dAyy_dx - (Gamxxy*l_Axy + Gamyxy*l_Ayy + Gamzxy*l_Ayz
                       + Gamxxy*l_Axy + Gamyxy*l_Ayy + Gamzxy*l_Ayz) - chix*l_Ayy/chin1;
    gyzx_A = dAyz_dx - (Gamxxz*l_Axy + Gamyxz*l_Ayy + Gamzxz*l_Ayz
                       + Gamxxy*l_Axz + Gamyxy*l_Ayz + Gamzxy*l_Azz) - chix*l_Ayz/chin1;
    gzzx_A = dAzz_dx - (Gamxxz*l_Axz + Gamyxz*l_Ayz + Gamzxz*l_Azz
                       + Gamxxz*l_Axz + Gamyxz*l_Ayz + Gamzxz*l_Azz) - chix*l_Azz/chin1;

    // y-direction
    gxxy_A = dAxx_dy - (Gamxxy*l_Axx + Gamyxy*l_Axy + Gamzxy*l_Axz
                       + Gamxxy*l_Axx + Gamyxy*l_Axy + Gamzxy*l_Axz) - chiy*l_Axx/chin1;
    gxyy_A = dAxy_dy - (Gamxyy*l_Axx + Gamyyy*l_Axy + Gamzyy*l_Axz
                       + Gamxxy*l_Axy + Gamyxy*l_Ayy + Gamzxy*l_Ayz) - chiy*l_Axy/chin1;
    gxzy_A = dAxz_dy - (Gamxyz*l_Axx + Gamyyz*l_Axy + Gamzyz*l_Axz
                       + Gamxxy*l_Axz + Gamyxy*l_Ayz + Gamzxy*l_Azz) - chiy*l_Axz/chin1;
    gyyy_A = dAyy_dy - (Gamxyy*l_Axy + Gamyyy*l_Ayy + Gamzyy*l_Ayz
                       + Gamxyy*l_Axy + Gamyyy*l_Ayy + Gamzyy*l_Ayz) - chiy*l_Ayy/chin1;
    gyzy_A = dAyz_dy - (Gamxyz*l_Axy + Gamyyz*l_Ayy + Gamzyz*l_Ayz
                       + Gamxyy*l_Axz + Gamyyy*l_Ayz + Gamzyy*l_Azz) - chiy*l_Ayz/chin1;
    gzzy_A = dAzz_dy - (Gamxyz*l_Axz + Gamyyz*l_Ayz + Gamzyz*l_Azz
                       + Gamxyz*l_Axz + Gamyyz*l_Ayz + Gamzyz*l_Azz) - chiy*l_Azz/chin1;

    // z-direction
    gxxz_A = dAxx_dz - (Gamxxz*l_Axx + Gamyxz*l_Axy + Gamzxz*l_Axz
                       + Gamxxz*l_Axx + Gamyxz*l_Axy + Gamzxz*l_Axz) - chiz*l_Axx/chin1;
    gxyz_A = dAxy_dz - (Gamxyz*l_Axx + Gamyyz*l_Axy + Gamzyz*l_Axz
                       + Gamxxz*l_Axy + Gamyxz*l_Ayy + Gamzxz*l_Ayz) - chiz*l_Axy/chin1;
    gxzz_A = dAxz_dz - (Gamxzz*l_Axx + Gamyzz*l_Axy + Gamzzz*l_Axz
                       + Gamxxz*l_Axz + Gamyxz*l_Ayz + Gamzxz*l_Azz) - chiz*l_Axz/chin1;
    gyyz_A = dAyy_dz - (Gamxyz*l_Axy + Gamyyz*l_Ayy + Gamzyz*l_Ayz
                       + Gamxyz*l_Axy + Gamyyz*l_Ayy + Gamzyz*l_Ayz) - chiz*l_Ayy/chin1;
    gyzz_A = dAyz_dz - (Gamxzz*l_Axy + Gamyzz*l_Ayy + Gamzzz*l_Ayz
                       + Gamxyz*l_Axz + Gamyyz*l_Ayz + Gamzyz*l_Azz) - chiz*l_Ayz/chin1;
    gzzz_A = dAzz_dz - (Gamxzz*l_Axz + Gamyzz*l_Ayz + Gamzzz*l_Azz
                       + Gamxzz*l_Axz + Gamyzz*l_Ayz + Gamzzz*l_Azz) - chiz*l_Azz/chin1;

    double l_Mxcon = gupxx*gxxx_A + gupyy*gxyy_A + gupzz*gxzz_A
        + gupxy*gxyx_A + gupxz*gxzx_A + gupyz*gxzy_A
        + gupxy*gxxy_A + gupxz*gxxz_A + gupyz*gxyz_A;
    double l_Mycon = gupxx*gxyx_A + gupyy*gyyy_A + gupzz*gyzz_A
        + gupxy*gyyx_A + gupxz*gyzx_A + gupyz*gyzy_A
        + gupxy*gxyy_A + gupxz*gxyz_A + gupyz*gyyz_A;
    double l_Mzcon = gupxx*gxzx_A + gupyy*gyzy_A + gupzz*gzzz_A
        + gupxy*gyzx_A + gupxz*gzzx_A + gupyz*gzzy_A
        + gupxy*gxzy_A + gupxz*gxzz_A + gupyz*gyzz_A;

    l_Mxcon = l_Mxcon - F2o3 * Kx - F8 * PI_VAL * l_Sx;
    l_Mycon = l_Mycon - F2o3 * Ky - F8 * PI_VAL * l_Sy;
    l_Mzcon = l_Mzcon - F2o3 * Kz - F8 * PI_VAL * l_Sz;

    // =====================================================================
    // Phase 18: TZ_rhs finalization
    // =====================================================================
    double TZ_rhs_final = alpn1 * Hcon_val / TWO;

    // =====================================================================
    // Write all outputs
    // =====================================================================
    chi_rhs[idx] = chi_rhs_val;
    trK_rhs[idx] = trK_rhs_val;
    gxx_rhs[idx] = gxx_rhs_val;
    gxy_rhs[idx] = gxy_rhs_val;
    gxz_rhs[idx] = gxz_rhs_val;
    gyy_rhs[idx] = gyy_rhs_val;
    gyz_rhs[idx] = gyz_rhs_val;
    gzz_rhs[idx] = gzz_rhs_val;
    Axx_rhs[idx] = Axx_rhs_val;
    Axy_rhs[idx] = Axy_rhs_val;
    Axz_rhs[idx] = Axz_rhs_val;
    Ayy_rhs[idx] = Ayy_rhs_val;
    Ayz_rhs[idx] = Ayz_rhs_val;
    Azz_rhs[idx] = Azz_rhs_val;
    Gamx_rhs[idx] = Gamx_rhs_val;
    Gamy_rhs[idx] = Gamy_rhs_val;
    Gamz_rhs[idx] = Gamz_rhs_val;
    Lap_rhs[idx]  = Lap_rhs_val;
    betax_rhs[idx] = betax_rhs_val;
    betay_rhs[idx] = betay_rhs_val;
    betaz_rhs[idx] = betaz_rhs_val;
    dtSfx_rhs[idx] = dtSfx_rhs_val;
    dtSfy_rhs[idx] = dtSfy_rhs_val;
    dtSfz_rhs[idx] = dtSfz_rhs_val;
    TZ_rhs[idx]   = TZ_rhs_final;

    Gamxxx_out[idx] = Gamxxx; Gamxxy_out[idx] = Gamxxy; Gamxxz_out[idx] = Gamxxz;
    Gamxyy_out[idx] = Gamxyy; Gamxyz_out[idx] = Gamxyz; Gamxzz_out[idx] = Gamxzz;
    Gamyxx_out[idx] = Gamyxx; Gamyxy_out[idx] = Gamyxy; Gamyxz_out[idx] = Gamyxz;
    Gamyyy_out[idx] = Gamyyy; Gamyyz_out[idx] = Gamyyz; Gamyzz_out[idx] = Gamyzz;
    Gamzxx_out[idx] = Gamzxx; Gamzxy_out[idx] = Gamzxy; Gamzxz_out[idx] = Gamzxz;
    Gamzyy_out[idx] = Gamzyy; Gamzyz_out[idx] = Gamzyz; Gamzzz_out[idx] = Gamzzz;

    Rxx_out[idx] = Rxx_v; Rxy_out[idx] = Rxy_v; Rxz_out[idx] = Rxz_v;
    Ryy_out[idx] = Ryy_v; Ryz_out[idx] = Ryz_v; Rzz_out[idx] = Rzz_v;

    Hcon[idx]   = Hcon_val;
    Mxcon[idx]  = l_Mxcon;
    Mycon[idx]  = l_Mycon;
    Mzcon[idx]  = l_Mzcon;
    Gmxcon[idx] = l_Gmxcon;
    Gmycon[idx] = l_Gmycon;
    Gmzcon[idx] = l_Gmzcon;
}

// ==========================================================================
// Kernel 3: Advection + constraint damping + KO dissipation
// Translates lines 1186-1389 of Z4c_rhs_ss.f90
// ==========================================================================
template<int GW>
__launch_bounds__(GPU_DIM_X * GPU_DIM_Y * GPU_DIM_Z, 1)
__global__ void z4c_ss_advection_dissipation_kernel(
    const int ex0, const int ex1, const int ex2,
    const int pex0, const int pex1, const int pex2,
    const int symmetry, const int onoff, const double eps, const int sst,
    const double* __restrict__ crho,
    const double* __restrict__ sigma,
    const double* __restrict__ R,
    const double* __restrict__ drhodx,
    const double* __restrict__ drhody,
    const double* __restrict__ drhodz,
    const double* __restrict__ dsigmadx,
    const double* __restrict__ dsigmady,
    const double* __restrict__ dsigmadz,
    const double* __restrict__ dRdx,
    const double* __restrict__ dRdy,
    const double* __restrict__ dRdz,
    const double* __restrict__ chi,
    const double* __restrict__ trK,
    const double* __restrict__ dxx,
    const double* __restrict__ gxy,
    const double* __restrict__ gxz,
    const double* __restrict__ dyy,
    const double* __restrict__ gyz,
    const double* __restrict__ dzz,
    const double* __restrict__ Axx,
    const double* __restrict__ Axy,
    const double* __restrict__ Axz,
    const double* __restrict__ Ayy,
    const double* __restrict__ Ayz,
    const double* __restrict__ Azz,
    const double* __restrict__ Gamx,
    const double* __restrict__ Gamy,
    const double* __restrict__ Gamz,
    const double* __restrict__ Gmxcon,
    const double* __restrict__ Gmycon,
    const double* __restrict__ Gmzcon,
    const double* __restrict__ Lap,
    const double* __restrict__ betax,
    const double* __restrict__ betay,
    const double* __restrict__ betaz,
    const double* __restrict__ dtSfx,
    const double* __restrict__ dtSfy,
    const double* __restrict__ dtSfz,
    const double* __restrict__ TZ,
    double* __restrict__ chi_rhs,
    double* __restrict__ trK_rhs,
    double* __restrict__ gxx_rhs,
    double* __restrict__ gxy_rhs,
    double* __restrict__ gxz_rhs,
    double* __restrict__ gyy_rhs,
    double* __restrict__ gyz_rhs,
    double* __restrict__ gzz_rhs,
    double* __restrict__ Axx_rhs,
    double* __restrict__ Axy_rhs,
    double* __restrict__ Axz_rhs,
    double* __restrict__ Ayy_rhs,
    double* __restrict__ Ayz_rhs,
    double* __restrict__ Azz_rhs,
    double* __restrict__ Gamx_rhs,
    double* __restrict__ Gamy_rhs,
    double* __restrict__ Gamz_rhs,
    double* __restrict__ Lap_rhs,
    double* __restrict__ betax_rhs,
    double* __restrict__ betay_rhs,
    double* __restrict__ betaz_rhs,
    double* __restrict__ dtSfx_rhs,
    double* __restrict__ dtSfy_rhs,
    double* __restrict__ dtSfz_rhs,
    double* __restrict__ TZ_rhs
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    int idx = z4c_idx3d(i, j, k, pex0, pex1);

    // Macro: advection term for a field
    #define ADVECT_SH(field, rhs, s1, s2, s3) do { \
        double _fx, _fy, _fz; \
        d_point_fderivs_shc_gw<GW>(ex, pex, field, &_fx, &_fy, &_fz, \
            crho, sigma, R, s1, s2, s3, symmetry, onoff, sst, \
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, \
            dRdx, dRdy, dRdz, i, j, k); \
        rhs[idx] += l_betax*_fx + l_betay*_fy + l_betaz*_fz; \
    } while(0)

    double l_betax = betax[idx];
    double l_betay = betay[idx];
    double l_betaz = betaz[idx];

    // Advection for all 24 evolution variables
    ADVECT_SH(dxx,  gxx_rhs,  SYM_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(gxy,  gxy_rhs,  ANTI_VAL, ANTI_VAL, SYM_VAL);
    ADVECT_SH(gxz,  gxz_rhs,  ANTI_VAL, SYM_VAL, ANTI_VAL);
    ADVECT_SH(dyy,  gyy_rhs,  SYM_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(gyz,  gyz_rhs,  SYM_VAL, ANTI_VAL, ANTI_VAL);
    ADVECT_SH(dzz,  gzz_rhs,  SYM_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(Axx,  Axx_rhs,  SYM_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(Axy,  Axy_rhs,  ANTI_VAL, ANTI_VAL, SYM_VAL);
    ADVECT_SH(Axz,  Axz_rhs,  ANTI_VAL, SYM_VAL, ANTI_VAL);
    ADVECT_SH(Ayy,  Ayy_rhs,  SYM_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(Ayz,  Ayz_rhs,  SYM_VAL, ANTI_VAL, ANTI_VAL);
    ADVECT_SH(Azz,  Azz_rhs,  SYM_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(chi,  chi_rhs,  SYM_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(trK,  trK_rhs,  SYM_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(Gamx, Gamx_rhs, ANTI_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(Gamy, Gamy_rhs, SYM_VAL, ANTI_VAL, SYM_VAL);
    ADVECT_SH(Gamz, Gamz_rhs, SYM_VAL, SYM_VAL, ANTI_VAL);
    ADVECT_SH(Lap,  Lap_rhs,  SYM_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(betax, betax_rhs, ANTI_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(betay, betay_rhs, SYM_VAL, ANTI_VAL, SYM_VAL);
    ADVECT_SH(betaz, betaz_rhs, SYM_VAL, SYM_VAL, ANTI_VAL);
    ADVECT_SH(dtSfx, dtSfx_rhs, ANTI_VAL, SYM_VAL, SYM_VAL);
    ADVECT_SH(dtSfy, dtSfy_rhs, SYM_VAL, ANTI_VAL, SYM_VAL);
    ADVECT_SH(dtSfz, dtSfz_rhs, SYM_VAL, SYM_VAL, ANTI_VAL);
    ADVECT_SH(TZ,    TZ_rhs,    SYM_VAL, SYM_VAL, SYM_VAL);

    // Constraint damping: use precomputed Gmxcon/Gmycon/Gmzcon from main kernel
    // Gmxcon = Gamx - Gamx_a, Gmycon = Gamy - Gamy_a, Gmzcon = Gamz - Gamz_a
    double l_TZ = TZ[idx];
    double l_alpn1 = Lap[idx] + ONE;
    TZ_rhs[idx]   -= l_alpn1 * (TWO + kappa2) * kappa1 * l_TZ;
    trK_rhs[idx]   += l_alpn1 * kappa1 * (ONE - kappa2) * l_TZ;
    Gamx_rhs[idx]  -= TWO * l_alpn1 * kappa1 * Gmxcon[idx];
    Gamy_rhs[idx]  -= TWO * l_alpn1 * kappa1 * Gmycon[idx];
    Gamz_rhs[idx]  -= TWO * l_alpn1 * kappa1 * Gmzcon[idx];

    // KO dissipation
    if (eps > 0) {
        chi_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, chi, crho, sigma, R, SYM_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        trK_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, trK, crho, sigma, R, SYM_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        gxx_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, dxx, crho, sigma, R, SYM_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        gxy_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, gxy, crho, sigma, R, ANTI_VAL, ANTI_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        gxz_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, gxz, crho, sigma, R, ANTI_VAL, SYM_VAL, ANTI_VAL, symmetry, eps, sst, i, j, k);
        gyy_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, dyy, crho, sigma, R, SYM_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        gyz_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, gyz, crho, sigma, R, SYM_VAL, ANTI_VAL, ANTI_VAL, symmetry, eps, sst, i, j, k);
        gzz_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, dzz, crho, sigma, R, SYM_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        Axx_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, Axx, crho, sigma, R, SYM_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        Axy_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, Axy, crho, sigma, R, ANTI_VAL, ANTI_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        Axz_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, Axz, crho, sigma, R, ANTI_VAL, SYM_VAL, ANTI_VAL, symmetry, eps, sst, i, j, k);
        Ayy_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, Ayy, crho, sigma, R, SYM_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        Ayz_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, Ayz, crho, sigma, R, SYM_VAL, ANTI_VAL, ANTI_VAL, symmetry, eps, sst, i, j, k);
        Azz_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, Azz, crho, sigma, R, SYM_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        Gamx_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Gamx, crho, sigma, R, ANTI_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        Gamy_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Gamy, crho, sigma, R, SYM_VAL, ANTI_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        Gamz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Gamz, crho, sigma, R, SYM_VAL, SYM_VAL, ANTI_VAL, symmetry, eps, sst, i, j, k);
        Lap_rhs[idx]  += d_kodis_sh_point<GW>(ex, pex, Lap, crho, sigma, R, SYM_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        betax_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, betax, crho, sigma, R, ANTI_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        betay_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, betay, crho, sigma, R, SYM_VAL, ANTI_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        betaz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, betaz, crho, sigma, R, SYM_VAL, SYM_VAL, ANTI_VAL, symmetry, eps, sst, i, j, k);
        dtSfx_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, dtSfx, crho, sigma, R, ANTI_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        dtSfy_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, dtSfy, crho, sigma, R, SYM_VAL, ANTI_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
        dtSfz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, dtSfz, crho, sigma, R, SYM_VAL, SYM_VAL, ANTI_VAL, symmetry, eps, sst, i, j, k);
        TZ_rhs[idx]    += d_kodis_sh_point<GW>(ex, pex, TZ, crho, sigma, R, SYM_VAL, SYM_VAL, SYM_VAL, symmetry, eps, sst, i, j, k);
    }
}

#undef ADVECT_SH

// ==========================================================================
// Launch helper: compute grid/block dims
// ==========================================================================
static void compute_grid_block(const int* ex, dim3& grid, dim3& block)
{
    block = dim3(GPU_DIM_X, GPU_DIM_Y, GPU_DIM_Z);
    grid = dim3(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z);
}

} // anonymous namespace

// ==========================================================================
// Public launch function
// ==========================================================================
void gpu_compute_rhs_z4c_ss_launch(
    cudaStream_t &stream, int* ex, int* pex, double T,
    double* crho, double* sigma, double* R,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* drhodxx, double* drhodxy, double* drhodxz,
    double* drhodyy, double* drhodyz, double* drhodzz,
    double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
    double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
    double* dRdxx, double* dRdxy, double* dRdxz,
    double* dRdyy, double* dRdyz, double* dRdzz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz, double* TZ,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs, double* TZ_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* Hcon, double* Mxcon, double* Mycon, double* Mzcon,
    double* Gmxcon, double* Gmycon, double* Gmzcon,
    int symmetry, int lev, double eps, int sst, int co)
{
    (void)co;  // unused in current implementation
    dim3 grid, block;
    compute_grid_block(ex, grid, block);

    // Kernel 1+2: Ricci + RHS + constraints
    z4c_ss_main_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        symmetry, lev, T, sst,
        crho, sigma, R,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap, betax, betay, betaz, dtSfx, dtSfy, dtSfz, TZ,
        chi_rhs, trK_rhs, gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs, betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs, TZ_rhs,
        rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        Hcon, Mxcon, Mycon, Mzcon, Gmxcon, Gmycon, Gmzcon);

    // Kernel 3: Advection + damping + KO dissipation
    z4c_ss_advection_dissipation_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        symmetry, 0, eps, sst,
        crho, sigma, R,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Gmxcon, Gmycon, Gmzcon, Lap, betax, betay, betaz, dtSfx, dtSfy, dtSfz, TZ,
        chi_rhs, trK_rhs, gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs, betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs, TZ_rhs);
}
