#include "macrodef.h"
// =========================================================================
// cpbc_gpu.cu -- CUDA kernels and host launch functions for the
//               compactified periodic boundary condition (CPBC) module.
//
// Mirrors the Fortran implementations:
//   - david_milton_cpbc_ss      -> gpu_david_milton_cpbc_ss_launch
//   - david_milton_extroplate_ss -> gpu_david_milton_extroplate_ss_launch
// =========================================================================

#include "cpbc_gpu.cuh"
#include "diff_new_sh_gpu.cuh"
#include "point_diff_new_sh_gpu.cuh"
#include "derivatives.h"

#include "gpu_manager.h"
#include "fmisc_gpu.cuh"

#include <cuda_runtime.h>
#include <math.h>
#include <cmath>
#include <stdio.h>

// ---------------------------------------------------------------------
// Compile-time constants matching macrodef.fh / Fortran parameters
// ---------------------------------------------------------------------
constexpr double ZEO = 0.0;
constexpr double ONE = 1.0;
constexpr double SYM = 1.0;
constexpr double ANTI = -1.0;

// Z4c parameters (must match get_Z4cparameters in cpbc.f90)
static __constant__ double d_kappa1 = 0.02;
static __constant__ double d_kappa2 = 0.0;
static __constant__ double d_eta    = 2.0;

// =====================================================================
// NOTE: Shell-to-Cartesian chain rule conversion is now handled by the
// proper device functions d_point_fderivs_shc_gw and
// d_point_fdderivs_shc_gw from point_diff_new_sh_gpu.cuh.
// These combine shell differentiation and coordinate transformation
// into a single call, replacing the former manual two-step approach.
// =====================================================================

// =====================================================================
// Kernel 1: cpbc_main_kernel<GW>
//   Template on ghost_width GW.
//   Each thread processes one (i,j,k) point in the boundary layer.
//   Computes derivatives via d_point_fderivs_shc_gw /
//   d_point_fdderivs_shc_gw which handle shell-to-Cartesian conversion
//   internally, then calls d_cpbc_point to produce RHS outputs.
// =====================================================================
template<int GW>
__global__ void cpbc_main_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ x,
    const double* __restrict__ y,
    const double* __restrict__ z,
    const double* __restrict__ crho,
    const double* __restrict__ sigma,
    const double* __restrict__ R,
    const double* __restrict__ TZ,
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
    // RHS arrays (output)
    double* __restrict__ trK_rhs,
    double* __restrict__ TZ_rhs,
    double* __restrict__ Axx_rhs,
    double* __restrict__ Axy_rhs,
    double* __restrict__ Axz_rhs,
    double* __restrict__ Ayy_rhs,
    double* __restrict__ Ayz_rhs,
    double* __restrict__ Azz_rhs,
    double* __restrict__ Gamx_rhs,
    double* __restrict__ Gamy_rhs,
    double* __restrict__ Gamz_rhs,
    int Symmetry, int sst,
    int k_start, int k_end,
    // Jacobian arrays for shell-to-Cartesian chain rule
    const double* __restrict__ drhodx, const double* __restrict__ drhody, const double* __restrict__ drhodz,
    const double* __restrict__ dsigmadx, const double* __restrict__ dsigmady, const double* __restrict__ dsigmadz,
    const double* __restrict__ dRdx, const double* __restrict__ dRdy, const double* __restrict__ dRdz,
    const double* __restrict__ drhodxx, const double* __restrict__ drhodxy, const double* __restrict__ drhodxz,
    const double* __restrict__ drhodyy, const double* __restrict__ drhodyz, const double* __restrict__ drhodzz,
    const double* __restrict__ dsigmadxx, const double* __restrict__ dsigmadxy, const double* __restrict__ dsigmadxz,
    const double* __restrict__ dsigmadyy, const double* __restrict__ dsigmadyz, const double* __restrict__ dsigmadzz,
    const double* __restrict__ dRdxx, const double* __restrict__ dRdxy, const double* __restrict__ dRdxz,
    const double* __restrict__ dRdyy, const double* __restrict__ dRdyz, const double* __restrict__ dRdzz)
{
    int gx = blockIdx.x * blockDim.x + threadIdx.x;
    int gy = blockIdx.y * blockDim.y + threadIdx.y;
    int gz = k_start + blockIdx.z * blockDim.z + threadIdx.z;

    if (gx >= ex0 || gy >= ex1) return;
    if (gz > k_end) return;

    size_t idx = (size_t)gx + (size_t)pex0 * ((size_t)gy + (size_t)pex1 * gz);

    int ex[3]  = {ex0, ex1, ex2};
    int pex[3] = {pex0, pex1, pex2};

    // Derived variables
    double alpn1 = Lap[idx] + ONE;
    double chin1 = chi[idx] + ONE;
    double gxx   = dxx[idx] + ONE;
    double gyy   = dyy[idx] + ONE;
    double gzz   = dzz[idx] + ONE;

    // -----------------------------------------------------------------
    // First derivatives of 23 BSSN fields
    // Each call computes shell derivatives and converts to Cartesian
    // in one step via d_point_fderivs_shc_gw.
    // -----------------------------------------------------------------

    // trK: SYM,SYM,SYM
    double Kx, Ky, Kz;
    d_point_fderivs_shc_gw<GW>(ex, pex, trK, &Kx, &Ky, &Kz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // chi: SYM,SYM,SYM
    double chix, chiy, chiz;
    d_point_fderivs_shc_gw<GW>(ex, pex, chi, &chix, &chiy, &chiz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // Lap: SYM,SYM,SYM
    double Lapx, Lapy, Lapz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Lap, &Lapx, &Lapy, &Lapz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // TZ: SYM,SYM,SYM
    double TZx, TZy, TZz;
    d_point_fderivs_shc_gw<GW>(ex, pex, TZ, &TZx, &TZy, &TZz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // Gamx: ANTI,SYM,SYM
    double Gamxx, Gamxy, Gamxz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Gamx, &Gamxx, &Gamxy, &Gamxz,
        crho, sigma, R, ANTI, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // Gamy: SYM,ANTI,SYM
    double Gamyx, Gamyy, Gamyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Gamy, &Gamyx, &Gamyy, &Gamyz,
        crho, sigma, R, SYM, ANTI, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // Gamz: SYM,SYM,ANTI
    double Gamzx, Gamzy, Gamzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Gamz, &Gamzx, &Gamzy, &Gamzz,
        crho, sigma, R, SYM, SYM, ANTI, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // betax: ANTI,SYM,SYM
    double sfxx, sfxy, sfxz;
    d_point_fderivs_shc_gw<GW>(ex, pex, betax, &sfxx, &sfxy, &sfxz,
        crho, sigma, R, ANTI, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // betay: SYM,ANTI,SYM
    double sfyx, sfyy, sfyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, betay, &sfyx, &sfyy, &sfyz,
        crho, sigma, R, SYM, ANTI, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // betaz: SYM,SYM,ANTI
    double sfzx, sfzy, sfzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, betaz, &sfzx, &sfzy, &sfzz,
        crho, sigma, R, SYM, SYM, ANTI, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // Axx: SYM,SYM,SYM
    double Axxx, Axxy, Axxz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Axx, &Axxx, &Axxy, &Axxz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // Axy: ANTI,ANTI,SYM
    double Axyx, Axyy, Axyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Axy, &Axyx, &Axyy, &Axyz,
        crho, sigma, R, ANTI, ANTI, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // Axz: ANTI,SYM,ANTI
    double Axzx, Axzy, Axzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Axz, &Axzx, &Axzy, &Axzz,
        crho, sigma, R, ANTI, SYM, ANTI, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // Ayy: SYM,SYM,SYM
    double Ayyx, Ayyy, Ayyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Ayy, &Ayyx, &Ayyy, &Ayyz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // Ayz: SYM,ANTI,ANTI
    double Ayzx, Ayzy, Ayzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Ayz, &Ayzx, &Ayzy, &Ayzz,
        crho, sigma, R, SYM, ANTI, ANTI, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // Azz: SYM,SYM,SYM
    double Azzx, Azzy, Azzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Azz, &Azzx, &Azzy, &Azzz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // dxx (= gxx - 1): SYM,SYM,SYM
    double gxxx, gxxy, gxxz;
    d_point_fderivs_shc_gw<GW>(ex, pex, dxx, &gxxx, &gxxy, &gxxz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // gxy: ANTI,ANTI,SYM
    double gxyx, gxyy, gxyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, gxy, &gxyx, &gxyy, &gxyz,
        crho, sigma, R, ANTI, ANTI, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // gxz: ANTI,SYM,ANTI
    double gxzx, gxzy, gxzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, gxz, &gxzx, &gxzy, &gxzz,
        crho, sigma, R, ANTI, SYM, ANTI, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // dyy (= gyy - 1): SYM,SYM,SYM
    double gyyx, gyyy, gyyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, dyy, &gyyx, &gyyy, &gyyz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // gyz: SYM,ANTI,ANTI
    double gyzx, gyzy, gyzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, gyz, &gyzx, &gyzy, &gyzz,
        crho, sigma, R, SYM, ANTI, ANTI, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // dzz (= gzz - 1): SYM,SYM,SYM
    double gzzx, gzzy, gzzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, dzz, &gzzx, &gzzy, &gzzz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        gx, gy, gz);

    // -----------------------------------------------------------------
    // Second derivatives (11 fields: Lap, chi, dxx, dyy, dzz, gxy, gxz,
    //                         gyz, betax, betay, betaz)
    // Each call computes shell derivatives and converts to Cartesian
    // in one step via d_point_fdderivs_shc_gw.
    // -----------------------------------------------------------------

    // Lap: SYM,SYM,SYM
    double Lapxx, Lapxy, Lapxz, Lapyy, Lapyz, Lapzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, Lap,
        &Lapxx, &Lapxy, &Lapxz, &Lapyy, &Lapyz, &Lapzz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        gx, gy, gz);

    // chi: SYM,SYM,SYM
    double chixx, chixy, chixz, chiyy, chiyz, chizz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, chi,
        &chixx, &chixy, &chixz, &chiyy, &chiyz, &chizz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        gx, gy, gz);

    // dxx: SYM,SYM,SYM
    double gxxxx, gxxxy, gxxxz, gxxyy, gxxyz, gxxzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, dxx,
        &gxxxx, &gxxxy, &gxxxz, &gxxyy, &gxxyz, &gxxzz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        gx, gy, gz);

    // dyy: SYM,SYM,SYM
    double gyyxx, gyyxy, gyyxz, gyyyy, gyyyz, gyyzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, dyy,
        &gyyxx, &gyyxy, &gyyxz, &gyyyy, &gyyyz, &gyyzz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        gx, gy, gz);

    // dzz: SYM,SYM,SYM
    double gzzxx, gzzxy, gzzxz, gzzyy, gzzyz, gzzzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, dzz,
        &gzzxx, &gzzxy, &gzzxz, &gzzyy, &gzzyz, &gzzzz,
        crho, sigma, R, SYM, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        gx, gy, gz);

    // gxy: ANTI,ANTI,SYM
    double gxyxx, gxyxy, gxyxz, gxyyy, gxyyz, gxyzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, gxy,
        &gxyxx, &gxyxy, &gxyxz, &gxyyy, &gxyyz, &gxyzz,
        crho, sigma, R, ANTI, ANTI, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        gx, gy, gz);

    // gxz: ANTI,SYM,ANTI
    double gxzxx, gxzxy, gxzxz, gxzyy, gxzyz, gxzzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, gxz,
        &gxzxx, &gxzxy, &gxzxz, &gxzyy, &gxzyz, &gxzzz,
        crho, sigma, R, ANTI, SYM, ANTI, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        gx, gy, gz);

    // gyz: SYM,ANTI,ANTI
    double gyzxx, gyzxy, gyzxz, gyzyy, gyzyz, gyzzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, gyz,
        &gyzxx, &gyzxy, &gyzxz, &gyzyy, &gyzyz, &gyzzz,
        crho, sigma, R, SYM, ANTI, ANTI, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        gx, gy, gz);

    // betax: ANTI,SYM,SYM
    double sfxxx, sfxxy, sfxxz, sfxyy, sfxyz, sfxzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, betax,
        &sfxxx, &sfxxy, &sfxxz, &sfxyy, &sfxyz, &sfxzz,
        crho, sigma, R, ANTI, SYM, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        gx, gy, gz);

    // betay: SYM,ANTI,SYM
    double sfyxx, sfyxy, sfyxz, sfyyy, sfyyz, sfyzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, betay,
        &sfyxx, &sfyxy, &sfyxz, &sfyyy, &sfyyz, &sfyzz,
        crho, sigma, R, SYM, ANTI, SYM, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        gx, gy, gz);

    // betaz: SYM,SYM,ANTI
    double sfzxx, sfzxy, sfzxz, sfzyy, sfzyz, sfzzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, betaz,
        &sfzxx, &sfzxy, &sfzxz, &sfzyy, &sfzyz, &sfzzz,
        crho, sigma, R, SYM, SYM, ANTI, Symmetry, 0, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        gx, gy, gz);

    // -----------------------------------------------------------------
    // Call d_cpbc_point to compute the RHS outputs
    // -----------------------------------------------------------------
    double rKhat, rTheta;
    double rA11, rA12, rA13, rA22, rA23, rA33;
    double rG1, rG2, rG3;

    d_cpbc_point(
        R[gz], x[idx], y[idx], z[idx],
        TZ[idx], chin1, trK[idx],
        gxx, gxy[idx], gxz[idx], gyy, gyz[idx], gzz,
        Axx[idx], Axy[idx], Axz[idx], Ayy[idx], Ayz[idx], Azz[idx],
        Gamx[idx], Gamy[idx], Gamz[idx],
        alpn1, betax[idx], betay[idx], betaz[idx],
        // Lap first derivatives
        Lapx, Lapy, Lapz,
        // Lap second derivatives
        Lapxx, Lapxy, Lapxz, Lapyy, Lapyz, Lapzz,
        // Shift first derivatives (sfx = betax derivatives, etc.)
        sfxx, sfxy, sfxz,
        sfyx, sfyy, sfyz,
        sfzx, sfzy, sfzz,
        // Shift second derivatives
        sfxxx, sfxxy, sfxxz, sfxyy, sfxyz, sfxzz,
        sfyxx, sfyxy, sfyxz, sfyyy, sfyyz, sfyzz,
        sfzxx, sfzxy, sfzxz, sfzyy, sfzyz, sfzzz,
        // chi first derivatives
        chix, chiy, chiz,
        // chi second derivatives
        chixx, chixy, chixz, chiyy, chiyz, chizz,
        // gij first derivatives (column-major order: dg_{ij,k})
        //   Fortran: gxxx,gxyx,gxzx, gyyx,gyzx,gzzx  -> k=1 (x)
        gxxx, gxyx, gxzx, gyyx, gyzx, gzzx,
        //   Fortran: gxxy,gxyy,gxzy, gyyy,gyzy,gzzy  -> k=2 (y)
        gxxy, gxyy, gxzy, gyyy, gyzy, gzzy,
        //   Fortran: gxxz,gxyz,gxzz, gyyz,gyzz,gzzz  -> k=3 (z)
        gxxz, gxyz, gxzz, gyyz, gyzz, gzzz,
        // gij second derivatives (ddg_{ij,kl})
        gxxxx, gxxxy, gxxxz, gxxyy, gxxyz, gxxzz,
        gxyxx, gxyxy, gxyxz, gxyyy, gxyyz, gxyzz,
        gxzxx, gxzxy, gxzxz, gxzyy, gxzyz, gxzzz,
        gyyxx, gyyxy, gyyxz, gyyyy, gyyyz, gyyzz,
        gyzxx, gyzxy, gyzxz, gyzyy, gyzyz, gyzzz,
        gzzxx, gzzxy, gzzxz, gzzyy, gzzyz, gzzzz,
        // trK first derivatives
        Kx, Ky, Kz,
        // Aij first derivatives
        Axxx, Axyx, Axzx, Ayyx, Ayzx, Azzx,
        Axxy, Axyy, Axzy, Ayyy, Ayzy, Azzy,
        Axxz, Axyz, Axzz, Ayyz, Ayzz, Azzz,
        // Gam first derivatives
        Gamxx, Gamxy, Gamxz,
        Gamyx, Gamyy, Gamyz,
        Gamzx, Gamzy, Gamzz,
        // TZ first derivatives
        TZx, TZy, TZz,
        // Outputs
        rKhat, rTheta,
        rA11, rA12, rA13, rA22, rA23, rA33,
        rG1, rG2, rG3,
        d_kappa1, d_kappa2, d_eta);

    // -----------------------------------------------------------------
    // Write outputs to RHS arrays
    // -----------------------------------------------------------------
    trK_rhs[idx]  = rKhat;
    TZ_rhs[idx]   = rTheta;
    Axx_rhs[idx]  = rA11;
    Axy_rhs[idx]  = rA12;
    Axz_rhs[idx]  = rA13;
    Ayy_rhs[idx]  = rA22;
    Ayz_rhs[idx]  = rA23;
    Azz_rhs[idx]  = rA33;
    Gamx_rhs[idx] = rG1;
    Gamy_rhs[idx] = rG2;
    Gamz_rhs[idx] = rG3;
}

// =====================================================================
// Kernel 2: cpbc_zero_rhs_kernel
//   Zeros all 25 RHS arrays over a given k-range.
//   Used for the inner ghost zone (gp=6) before applying CPBC.
// =====================================================================
__global__ void cpbc_zero_rhs_kernel(
    int ex0, int ex1,
    int pex0, int pex1,
    double* __restrict__ TZ_rhs,
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
    int k_start, int k_end)
{
    int gx = blockIdx.x * blockDim.x + threadIdx.x;
    int gy = blockIdx.y * blockDim.y + threadIdx.y;
    int gz = k_start + blockIdx.z * blockDim.z + threadIdx.z;

    if (gx >= ex0 || gy >= ex1) return;
    if (gz > k_end) return;

    size_t idx = (size_t)gx + (size_t)pex0 * ((size_t)gy + (size_t)pex1 * gz);

    TZ_rhs[idx]    = ZEO;
    chi_rhs[idx]   = ZEO;
    trK_rhs[idx]   = ZEO;
    gxx_rhs[idx]   = ZEO;
    gxy_rhs[idx]   = ZEO;
    gxz_rhs[idx]   = ZEO;
    gyy_rhs[idx]   = ZEO;
    gyz_rhs[idx]   = ZEO;
    gzz_rhs[idx]   = ZEO;
    Axx_rhs[idx]   = ZEO;
    Axy_rhs[idx]   = ZEO;
    Axz_rhs[idx]   = ZEO;
    Ayy_rhs[idx]   = ZEO;
    Ayz_rhs[idx]   = ZEO;
    Azz_rhs[idx]   = ZEO;
    Gamx_rhs[idx]  = ZEO;
    Gamy_rhs[idx]  = ZEO;
    Gamz_rhs[idx]  = ZEO;
    Lap_rhs[idx]   = ZEO;
    betax_rhs[idx] = ZEO;
    betay_rhs[idx] = ZEO;
    betaz_rhs[idx] = ZEO;
    dtSfx_rhs[idx] = ZEO;
    dtSfy_rhs[idx] = ZEO;
    dtSfz_rhs[idx] = ZEO;
}

// =====================================================================
// Kernel 3: repo_extro_kernel
//   2D grid over (i,j). Each thread extrapolates a full z-column in
//   the ghost zone using one of three extrapolation methods (tpp).
// =====================================================================
__global__ void repo_extro_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    double* __restrict__ f,
    int k_base,      // 0-based start of ghost zone (ex2 - ghost_width)
    int tpp,         // 0=Lagrange, 1=CG, 2=extended Lagrange
    int gw)          // ghost_width
{
    int gi = blockIdx.x * blockDim.x + threadIdx.x;
    int gj = blockIdx.y * blockDim.y + threadIdx.y;

    if (gi >= ex0 || gj >= ex1) return;

    const size_t zstride = (size_t)pex0 * pex1;

    if (tpp == 0) {
        // Lagrange: NP = ghost_width*2+1
        int NP = gw * 2 + 1;
        for (int k = k_base; k < ex2; ++k) {
            size_t idx = (size_t)gi + (size_t)pex0 * ((size_t)gj + (size_t)pex1 * k);
            double f_local[16];
            size_t base = (size_t)gi + (size_t)pex0 * ((size_t)gj + (size_t)pex1 * (k - NP));
            for (int n = 0; n < NP; ++n)
                f_local[n] = f[base + n * zstride];
            f[idx] = d_extroplate_lag(NP, f_local);
        }
    } else if (tpp == 1) {
        // Chebyshev-Gauss: NP = ghost_width*2
        int NP = gw * 2;
        for (int k = k_base; k < ex2; ++k) {
            size_t idx = (size_t)gi + (size_t)pex0 * ((size_t)gj + (size_t)pex1 * k);
            double f_local[16];
            size_t base = (size_t)gi + (size_t)pex0 * ((size_t)gj + (size_t)pex1 * (k - NP));
            for (int n = 0; n < NP; ++n)
                f_local[n] = f[base + n * zstride];
            f[idx] = d_extroplate_cg(NP, f_local);
        }
    } else if (tpp == 2) {
        // Extended Lagrange: NP = ghost_width*2+1 + CPBC_ghost_width
        int NP = gw * 2 + 1 + gw;
        int NI = NP - gw;
        double f_local[32];
        // Gather the tail NP values from (ex2-NP) to (ex2-1)
        size_t base = (size_t)gi + (size_t)pex0 * ((size_t)gj + (size_t)pex1 * (ex2 - NP));
        for (int n = 0; n < NP; ++n)
            f_local[n] = f[base + n * zstride];
        // In-place modification of ghost zone points via d_extroplate_lag2
        d_extroplate_lag2(NP, f_local, gw);
        // Write back the modified ghost zone points (last gw elements)
        for (int n = NI; n < NP; ++n) {
            size_t dst = base + n * zstride;
            f[dst] = f_local[n];
        }
    }
}

// =====================================================================
// Host launch function: gpu_david_milton_cpbc_ss_launch (templated impl)
// =====================================================================
template<int GW>
void gpu_david_milton_cpbc_ss_launch_impl(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    double *crho, double *sigma, double *R,
    double *x, double *y, double *z,
    double *drhodx, double *drhody, double *drhodz,
    double *dsigmadx, double *dsigmady, double *dsigmadz,
    double *dRdx, double *dRdy, double *dRdz,
    double *drhodxx, double *drhodxy, double *drhodxz,
    double *drhodyy, double *drhodyz, double *drhodzz,
    double *dsigmadxx, double *dsigmadxy, double *dsigmadxz,
    double *dsigmadyy, double *dsigmadyz, double *dsigmadzz,
    double *dRdxx, double *dRdxy, double *dRdxz,
    double *dRdyy, double *dRdyz, double *dRdzz,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double *TZ, double *chi, double *trK,
    double *dxx, double *gxy, double *gxz,
    double *dyy, double *gyz, double *dzz,
    double *Axx, double *Axy, double *Axz,
    double *Ayy, double *Ayz, double *Azz,
    double *Gamx, double *Gamy, double *Gamz,
    double *Lap, double *betax, double *betay, double *betaz,
    double *dtSfx, double *dtSfy, double *dtSfz,
    double *TZ_rhs, double *chi_rhs, double *trK_rhs,
    double *gxx_rhs, double *gxy_rhs, double *gxz_rhs,
    double *gyy_rhs, double *gyz_rhs, double *gzz_rhs,
    double *Axx_rhs, double *Axy_rhs, double *Axz_rhs,
    double *Ayy_rhs, double *Ayz_rhs, double *Azz_rhs,
    double *Gamx_rhs, double *Gamy_rhs, double *Gamz_rhs,
    double *Lap_rhs, double *betax_rhs, double *betay_rhs, double *betaz_rhs,
    double *dtSfx_rhs, double *dtSfy_rhs, double *dtSfz_rhs,
    double *pGamxxx, double *pGamxxy, double *pGamxxz,
    double *pGamxyy, double *pGamxyz, double *pGamxzz,
    double *pGamyxx, double *pGamyxy, double *pGamyxz,
    double *pGamyyy, double *pGamyyz, double *pGamyzz,
    double *pGamzxx, double *pGamzxy, double *pGamzxz,
    double *pGamzyy, double *pGamzyz, double *pGamzzz,
    double *Rxx, double *Rxy, double *Rxz,
    double *Ryy, double *Ryz, double *Rzz,
    double *Gmxcon, double *Gmycon, double *Gmzcon,
    int Symmetry, double eps, int sst,
    bool has_outer)
{
    (void)dtSfx; (void)dtSfy; (void)dtSfz;
    (void)pGamxxx; (void)pGamxxy; (void)pGamxxz;
    (void)pGamxyy; (void)pGamxyz; (void)pGamxzz;
    (void)pGamyxx; (void)pGamyxy; (void)pGamyxz;
    (void)pGamyyy; (void)pGamyyz; (void)pGamyzz;
    (void)pGamzxx; (void)pGamzxy; (void)pGamzxz;
    (void)pGamzyy; (void)pGamzyz; (void)pGamzzz;
    (void)Rxx; (void)Rxy; (void)Rxz;
    (void)Ryy; (void)Ryz; (void)Rzz;
    (void)Gmxcon; (void)Gmycon; (void)Gmzcon;

    // -----------------------------------------------------------------
    // Compute boundary layer indices on host.
    // has_outer is passed from the caller, who determines whether this
    // patch touches the outer boundary (Fortran: abs(R(ex3)-zmax) < dR).
    // -----------------------------------------------------------------

    // 0-based k indices:
    // CPBC ghost zone: from (ex[2] - CPBC_ghost_width) to (ex[2] - 1)
    // The boundary layer (gp=3) is a single layer at k = ex[2] - CPBC_ghost_width - 1
    // In 1-based Fortran: k from ex3-CPBC_ghost_width to ex3-CPBC_ghost_width
    // In 0-based: k from ex[2]-CPBC_ghost_width-1 to ex[2]-CPBC_ghost_width-1
    int cpbc_k = ex[2] - GW - 1;  // 0-based boundary layer index
    int ghost_k_start = ex[2] - GW; // 0-based first ghost zone k index
    int ghost_k_end = ex[2] - 1;    // 0-based last ghost zone k index

    dim3 block3d(8, 8, 4);

    // -----------------------------------------------------------------
    // Step 1: Zero RHS in ghost zone (gp=6) if outer boundary exists
    // -----------------------------------------------------------------
    if (has_outer) {
        dim3 grid_zero(
            (ex[0] + block3d.x - 1) / block3d.x,
            (ex[1] + block3d.y - 1) / block3d.y,
            (ghost_k_end - ghost_k_start + 1 + block3d.z - 1) / block3d.z);

        cpbc_zero_rhs_kernel<<<grid_zero, block3d, 0, stream>>>(
            ex[0], ex[1],
            pex[0], pex[1],
            TZ_rhs, chi_rhs, trK_rhs,
            gxx_rhs, gxy_rhs, gxz_rhs,
            gyy_rhs, gyz_rhs, gzz_rhs,
            Axx_rhs, Axy_rhs, Axz_rhs,
            Ayy_rhs, Ayz_rhs, Azz_rhs,
            Gamx_rhs, Gamy_rhs, Gamz_rhs,
            Lap_rhs, betax_rhs, betay_rhs, betaz_rhs,
            dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
            ghost_k_start, ghost_k_end);
    }

    // -----------------------------------------------------------------
    // Step 2: Main CPBC kernel on boundary layer (gp=3)
    //   Fortran iterates k from layer(3,3)=ex3-CPBC_ghost_width to
    //   layer(6,3)=ex3-CPBC_ghost_width (a single k-plane).
    //   0-based: k from cpbc_k to cpbc_k.
    // -----------------------------------------------------------------
    if (has_outer) {
        // Single k-layer: grid z-dimension = 1
        dim3 grid_cpbc(
            (ex[0] + block3d.x - 1) / block3d.x,
            (ex[1] + block3d.y - 1) / block3d.y,
            1);

        cpbc_main_kernel<GW><<<grid_cpbc, block3d, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1], pex[2],
            x, y, z,
            crho, sigma, R,
            TZ, chi, trK,
            dxx, gxy, gxz, dyy, gyz, dzz,
            Axx, Axy, Axz, Ayy, Ayz, Azz,
            Gamx, Gamy, Gamz,
            Lap, betax, betay, betaz,
            trK_rhs, TZ_rhs,
            Axx_rhs, Axy_rhs, Axz_rhs,
            Ayy_rhs, Ayz_rhs, Azz_rhs,
            Gamx_rhs, Gamy_rhs, Gamz_rhs,
            Symmetry, sst,
            cpbc_k, cpbc_k,
            drhodx, drhody, drhodz,
            dsigmadx, dsigmady, dsigmadz,
            dRdx, dRdy, dRdz,
            drhodxx, drhodxy, drhodxz,
            drhodyy, drhodyz, drhodzz,
            dsigmadxx, dsigmadxy, dsigmadxz,
            dsigmadyy, dsigmadyz, dsigmadzz,
            dRdxx, dRdxy, dRdxz,
            dRdyy, dRdyz, dRdzz);
    }

    // -----------------------------------------------------------------
    // Step 3: KO dissipation (if eps > 0)
    //   Call gpu_kodis_sh_launch for each of the 25 BSSN variables
    //   with the appropriate SoA symmetry flags.
    // -----------------------------------------------------------------
    if (eps > 0.0) {
        // SoA symmetry arrays

        // chi: SSS
        gpu_kodis_sh_launch(stream, ex, pex, chi, chi_rhs,
            crho, sigma, R, 1.0, 1.0, 1.0, Symmetry, eps, sst);
        // trK: SSS
        gpu_kodis_sh_launch(stream, ex, pex, trK, trK_rhs,
            crho, sigma, R, 1.0, 1.0, 1.0, Symmetry, eps, sst);
        // dxx -> gxx_rhs: SSS
        gpu_kodis_sh_launch(stream, ex, pex, dxx, gxx_rhs,
            crho, sigma, R, 1.0, 1.0, 1.0, Symmetry, eps, sst);
        // gxy: AAS
        gpu_kodis_sh_launch(stream, ex, pex, gxy, gxy_rhs,
            crho, sigma, R, -1.0, -1.0, 1.0, Symmetry, eps, sst);
        // gxz: ASA
        gpu_kodis_sh_launch(stream, ex, pex, gxz, gxz_rhs,
            crho, sigma, R, -1.0, 1.0, -1.0, Symmetry, eps, sst);
        // dyy -> gyy_rhs: SSS
        gpu_kodis_sh_launch(stream, ex, pex, dyy, gyy_rhs,
            crho, sigma, R, 1.0, 1.0, 1.0, Symmetry, eps, sst);
        // gyz: SAA
        gpu_kodis_sh_launch(stream, ex, pex, gyz, gyz_rhs,
            crho, sigma, R, 1.0, -1.0, -1.0, Symmetry, eps, sst);
        // dzz -> gzz_rhs: SSS
        gpu_kodis_sh_launch(stream, ex, pex, dzz, gzz_rhs,
            crho, sigma, R, 1.0, 1.0, 1.0, Symmetry, eps, sst);
        // Axx: SSS
        gpu_kodis_sh_launch(stream, ex, pex, Axx, Axx_rhs,
            crho, sigma, R, 1.0, 1.0, 1.0, Symmetry, eps, sst);
        // Axy: AAS
        gpu_kodis_sh_launch(stream, ex, pex, Axy, Axy_rhs,
            crho, sigma, R, -1.0, -1.0, 1.0, Symmetry, eps, sst);
        // Axz: ASA
        gpu_kodis_sh_launch(stream, ex, pex, Axz, Axz_rhs,
            crho, sigma, R, -1.0, 1.0, -1.0, Symmetry, eps, sst);
        // Ayy: SSS
        gpu_kodis_sh_launch(stream, ex, pex, Ayy, Ayy_rhs,
            crho, sigma, R, 1.0, 1.0, 1.0, Symmetry, eps, sst);
        // Ayz: SAA
        gpu_kodis_sh_launch(stream, ex, pex, Ayz, Ayz_rhs,
            crho, sigma, R, 1.0, -1.0, -1.0, Symmetry, eps, sst);
        // Azz: SSS
        gpu_kodis_sh_launch(stream, ex, pex, Azz, Azz_rhs,
            crho, sigma, R, 1.0, 1.0, 1.0, Symmetry, eps, sst);
        // Gamx: ASS
        gpu_kodis_sh_launch(stream, ex, pex, Gamx, Gamx_rhs,
            crho, sigma, R, -1.0, 1.0, 1.0, Symmetry, eps, sst);
        // Gamy: SAS
        gpu_kodis_sh_launch(stream, ex, pex, Gamy, Gamy_rhs,
            crho, sigma, R, 1.0, -1.0, 1.0, Symmetry, eps, sst);
        // Gamz: SSA
        gpu_kodis_sh_launch(stream, ex, pex, Gamz, Gamz_rhs,
            crho, sigma, R, 1.0, 1.0, -1.0, Symmetry, eps, sst);
        // Lap: SSS
        gpu_kodis_sh_launch(stream, ex, pex, Lap, Lap_rhs,
            crho, sigma, R, 1.0, 1.0, 1.0, Symmetry, eps, sst);
        // betax: ASS
        gpu_kodis_sh_launch(stream, ex, pex, betax, betax_rhs,
            crho, sigma, R, -1.0, 1.0, 1.0, Symmetry, eps, sst);
        // betay: SAS
        gpu_kodis_sh_launch(stream, ex, pex, betay, betay_rhs,
            crho, sigma, R, 1.0, -1.0, 1.0, Symmetry, eps, sst);
        // betaz: SSA
        gpu_kodis_sh_launch(stream, ex, pex, betaz, betaz_rhs,
            crho, sigma, R, 1.0, 1.0, -1.0, Symmetry, eps, sst);
        // TZ: SSS
        gpu_kodis_sh_launch(stream, ex, pex, TZ, TZ_rhs,
            crho, sigma, R, 1.0, 1.0, 1.0, Symmetry, eps, sst);
        // dtSfx: ASS
        gpu_kodis_sh_launch(stream, ex, pex, dtSfx, dtSfx_rhs,
            crho, sigma, R, -1.0, 1.0, 1.0, Symmetry, eps, sst);
        // dtSfy: SAS
        gpu_kodis_sh_launch(stream, ex, pex, dtSfy, dtSfy_rhs,
            crho, sigma, R, 1.0, -1.0, 1.0, Symmetry, eps, sst);
        // dtSfz: SSA
        gpu_kodis_sh_launch(stream, ex, pex, dtSfz, dtSfz_rhs,
            crho, sigma, R, 1.0, 1.0, -1.0, Symmetry, eps, sst);
    }
}

// =====================================================================
// Non-templated wrapper: gpu_david_milton_cpbc_ss_launch
//   Instantiates with the compile-time ghost_width from macrodef.fh.
// =====================================================================
void gpu_david_milton_cpbc_ss_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    double *crho, double *sigma, double *R,
    double *x, double *y, double *z,
    double *drhodx, double *drhody, double *drhodz,
    double *dsigmadx, double *dsigmady, double *dsigmadz,
    double *dRdx, double *dRdy, double *dRdz,
    double *drhodxx, double *drhodxy, double *drhodxz,
    double *drhodyy, double *drhodyz, double *drhodzz,
    double *dsigmadxx, double *dsigmadxy, double *dsigmadxz,
    double *dsigmadyy, double *dsigmadyz, double *dsigmadzz,
    double *dRdxx, double *dRdxy, double *dRdxz,
    double *dRdyy, double *dRdyz, double *dRdzz,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double *TZ, double *chi, double *trK,
    double *dxx, double *gxy, double *gxz,
    double *dyy, double *gyz, double *dzz,
    double *Axx, double *Axy, double *Axz,
    double *Ayy, double *Ayz, double *Azz,
    double *Gamx, double *Gamy, double *Gamz,
    double *Lap, double *betax, double *betay, double *betaz,
    double *dtSfx, double *dtSfy, double *dtSfz,
    double *TZ_rhs, double *chi_rhs, double *trK_rhs,
    double *gxx_rhs, double *gxy_rhs, double *gxz_rhs,
    double *gyy_rhs, double *gyz_rhs, double *gzz_rhs,
    double *Axx_rhs, double *Axy_rhs, double *Axz_rhs,
    double *Ayy_rhs, double *Ayz_rhs, double *Azz_rhs,
    double *Gamx_rhs, double *Gamy_rhs, double *Gamz_rhs,
    double *Lap_rhs, double *betax_rhs, double *betay_rhs, double *betaz_rhs,
    double *dtSfx_rhs, double *dtSfy_rhs, double *dtSfz_rhs,
    double *pGamxxx, double *pGamxxy, double *pGamxxz,
    double *pGamxyy, double *pGamxyz, double *pGamxzz,
    double *pGamyxx, double *pGamyxy, double *pGamyxz,
    double *pGamyyy, double *pGamyyz, double *pGamyzz,
    double *pGamzxx, double *pGamzxy, double *pGamzxz,
    double *pGamzyy, double *pGamzyz, double *pGamzzz,
    double *Rxx, double *Rxy, double *Rxz,
    double *Ryy, double *Ryz, double *Rzz,
    double *Gmxcon, double *Gmycon, double *Gmzcon,
    int Symmetry, double eps, int sst,
    bool has_outer)
{
    gpu_david_milton_cpbc_ss_launch_impl<ghost_width>(
        stream, ex, pex,
        crho, sigma, R,
        x, y, z,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        xmin, ymin, zmin, xmax, ymax, zmax,
        TZ, chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz,
        Lap, betax, betay, betaz,
        dtSfx, dtSfy, dtSfz,
        TZ_rhs, chi_rhs, trK_rhs,
        gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs,
        Lap_rhs, betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        pGamxxx, pGamxxy, pGamxxz,
        pGamxyy, pGamxyz, pGamxzz,
        pGamyxx, pGamyxy, pGamyxz,
        pGamyyy, pGamyyz, pGamyzz,
        pGamzxx, pGamzxy, pGamzxz,
        pGamzyy, pGamzyz, pGamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        Gmxcon, Gmycon, Gmzcon,
        Symmetry, eps, sst,
        has_outer);
}

// =====================================================================
// Host launch function: gpu_david_milton_extroplate_ss_launch
//
// Mirrors f_david_milton_extroplate_ss.  Launches repo_extro_kernel
// 25 times (once per BSSN variable), each with tpp=1 (Chebyshev-Gauss).
// =====================================================================
void gpu_david_milton_extroplate_ss_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    double *crho, double *sigma, double *R,
    double *TZ, double *chi, double *trK,
    double *dxx, double *gxy, double *gxz,
    double *dyy, double *gyz, double *dzz,
    double *Axx, double *Axy, double *Axz,
    double *Ayy, double *Ayz, double *Azz,
    double *Gmx, double *Gmy, double *Gmz,
    double *Lap, double *betax, double *betay, double *betaz,
    double *dtSfx, double *dtSfy, double *dtSfz,
    double zmin, double zmax,
    bool has_outer)
{
    (void)crho; (void)sigma; (void)R;
    (void)zmin; (void)zmax;

    // Only extrapolate for patches that touch the outer boundary
    if (!has_outer) return;

    // Fortran bounds check: ex(3) > CPBC_ghost_width + (ghost_width*2+1)
    if (ex[2] <= ghost_width + (ghost_width * 2 + 1)) {
        fprintf(stderr, "repo_extro: ex[2]=%d too small for ghost_width=%d\n",
                ex[2], ghost_width);
        return;
    }

    // Use tpp=1 (Chebyshev-Gauss extrapolation) matching Fortran default
    const int tpp = 1;

    // 0-based ghost zone start: ex[2] - ghost_width (0-based)
    // Fortran: k from layer(3,gp)+1 to ex(3), where layer(3,gp) = ex(3)-CPBC_ghost_width
    //   => 1-based k_start = ex(3)-CPBC_ghost_width+1
    //   => 0-based k_start = ex[2]-ghost_width
    int k_base = ex[2] - ghost_width;

    dim3 block2d(16, 16, 1);
    dim3 grid2d(
        (ex[0] + block2d.x - 1) / block2d.x,
        (ex[1] + block2d.y - 1) / block2d.y,
        1);

    // List of 25 field pointers to extrapolate
    double* fields[25] = {
        TZ, chi, dxx, dyy, dzz,
        gxy, gxz, gyz,
        trK,
        Axx, Ayy, Azz, Axy, Axz, Ayz,
        Gmx, Gmy, Gmz,
        Lap, betax, betay, betaz,
        dtSfx, dtSfy, dtSfz
    };

    for (int v = 0; v < 25; ++v) {
        repo_extro_kernel<<<grid2d, block2d, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1],
            fields[v],
            k_base, tpp, ghost_width);
    }
}
