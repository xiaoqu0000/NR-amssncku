#include "macrodef.h"
#include "bssn_rhs.h"
#include "fmisc_gpu.cuh"
#include "point_diff_new_sh_gpu.cuh"
#include "kodiss_sh_gpu.cuh"
#include <cuda_runtime.h>
#include <math.h>

extern int BHN;
extern double Mass[3];
extern double PBH[9];

constexpr double ONE = 1.0;
constexpr double TWO = 2.0;
constexpr double HALF = 0.5;
constexpr double F1o3 = 1.0 / 3.0;
constexpr double F2o3 = 2.0 / 3.0;
constexpr double F3o2 = 1.5;
constexpr double FF = 0.75;
constexpr double ETA = 2.0;
constexpr double SYM = 1.0;
constexpr double ANTI = -1.0;
constexpr double EIGHT = 8.0;
constexpr double FOUR = 4.0;
constexpr double F16 = 16.0;
constexpr double F8 = 8.0;
constexpr double PIO = M_PI;

template<int GW, int GAUGE_PARAM>
__global__ void bssn_rhs_ss_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    int symmetry, int lev, double eps, int co, int sst,
    int BHN, double Porg1, double Porg2, double Porg3,
    double Porg4, double Porg5, double Porg6,
    double Mass1, double Mass2,
    const double* __restrict__ crho, const double* __restrict__ sigma, const double* __restrict__ Rcrd,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ drhodx, const double* __restrict__ drhody, const double* __restrict__ drhodz,
    const double* __restrict__ dsigmadx, const double* __restrict__ dsigmady, const double* __restrict__ dsigmadz,
    const double* __restrict__ dRdx, const double* __restrict__ dRdy, const double* __restrict__ dRdz,
    const double* __restrict__ drhodxx, const double* __restrict__ drhodxy, const double* __restrict__ drhodxz,
    const double* __restrict__ drhodyy, const double* __restrict__ drhodyz, const double* __restrict__ drhodzz,
    const double* __restrict__ dsigmadxx, const double* __restrict__ dsigmadxy, const double* __restrict__ dsigmadxz,
    const double* __restrict__ dsigmadyy, const double* __restrict__ dsigmadyz, const double* __restrict__ dsigmadzz,
    const double* __restrict__ dRdxx, const double* __restrict__ dRdxy, const double* __restrict__ dRdxz,
    const double* __restrict__ dRdyy, const double* __restrict__ dRdyz, const double* __restrict__ dRdzz,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ dxx, const double* __restrict__ gxy_f, const double* __restrict__ gxz_f,
    const double* __restrict__ dyy, const double* __restrict__ gyz_f, const double* __restrict__ dzz,
    const double* __restrict__ Axx, const double* __restrict__ Axy, const double* __restrict__ Axz,
    const double* __restrict__ Ayy, const double* __restrict__ Ayz, const double* __restrict__ Azz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    const double* __restrict__ Lap, const double* __restrict__ betax, const double* __restrict__ betay, const double* __restrict__ betaz,
    const double* __restrict__ dtSfx, const double* __restrict__ dtSfy, const double* __restrict__ dtSfz,
    const double* __restrict__ rho, const double* __restrict__ Sx, const double* __restrict__ Sy, const double* __restrict__ Sz,
    const double* __restrict__ Sxx, const double* __restrict__ Sxy_f, const double* __restrict__ Sxz_f,
    const double* __restrict__ Syy, const double* __restrict__ Syz_f, const double* __restrict__ Szz,
    double* __restrict__ chi_rhs, double* __restrict__ trK_rhs,
    double* __restrict__ gxx_rhs, double* __restrict__ gxy_rhs, double* __restrict__ gxz_rhs,
    double* __restrict__ gyy_rhs, double* __restrict__ gyz_rhs, double* __restrict__ gzz_rhs,
    double* __restrict__ Axx_rhs, double* __restrict__ Axy_rhs, double* __restrict__ Axz_rhs,
    double* __restrict__ Ayy_rhs, double* __restrict__ Ayz_rhs, double* __restrict__ Azz_rhs,
    double* __restrict__ Gamx_rhs, double* __restrict__ Gamy_rhs, double* __restrict__ Gamz_rhs,
    double* __restrict__ Lap_rhs, double* __restrict__ betax_rhs, double* __restrict__ betay_rhs, double* __restrict__ betaz_rhs,
    double* __restrict__ dtSfx_rhs, double* __restrict__ dtSfy_rhs, double* __restrict__ dtSfz_rhs,
    double* __restrict__ Gamxxx, double* __restrict__ Gamxxy, double* __restrict__ Gamxxz,
    double* __restrict__ Gamxyy, double* __restrict__ Gamxyz, double* __restrict__ Gamxzz,
    double* __restrict__ Gamyxx, double* __restrict__ Gamyxy, double* __restrict__ Gamyxz,
    double* __restrict__ Gamyyy, double* __restrict__ Gamyyz, double* __restrict__ Gamyzz,
    double* __restrict__ Gamzxx, double* __restrict__ Gamzxy, double* __restrict__ Gamzxz,
    double* __restrict__ Gamzyy, double* __restrict__ Gamzyz, double* __restrict__ Gamzzz,
    double* __restrict__ Ricxx, double* __restrict__ Ricxy, double* __restrict__ Ricxz,
    double* __restrict__ Ricyy, double* __restrict__ Ricyz, double* __restrict__ Riczz,
    double* __restrict__ ham_Res, double* __restrict__ movx_Res, double* __restrict__ movy_Res, double* __restrict__ movz_Res,
    double* __restrict__ Gmx_Res, double* __restrict__ Gmy_Res, double* __restrict__ Gmz_Res
)

{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int ex[3] = {ex0, ex1, ex2};
    int pex[3] = {pex0, pex1, pex2};
    int idx = (k * pex1 + j) * pex0 + i;
    int onoff = 0;

    // ===================================================================
    // Load local field values and compute constants
    // ===================================================================

    double l_alpn1 = Lap[idx] + ONE;
    double l_chin1 = chi[idx] + ONE;
    double inv_chin1 = ONE / l_chin1;

    double l_gxx = dxx[idx] + ONE;
    double l_gxy = gxy_f[idx];
    double l_gxz = gxz_f[idx];
    double l_gyy = dyy[idx] + ONE;
    double l_gyz = gyz_f[idx];
    double l_gzz = dzz[idx] + ONE;

    double l_Axx = Axx[idx]; double l_Axy = Axy[idx]; double l_Axz = Axz[idx];
    double l_Ayy = Ayy[idx]; double l_Ayz = Ayz[idx]; double l_Azz = Azz[idx];

    double l_Gamx = Gamx[idx]; double l_Gamy = Gamy[idx]; double l_Gamz = Gamz[idx];
    double l_Lap = Lap[idx];
    double l_betax = betax[idx]; double l_betay = betay[idx]; double l_betaz = betaz[idx];
    double l_dtSfx = dtSfx[idx]; double l_dtSfy = dtSfy[idx]; double l_dtSfz = dtSfz[idx];
    double l_trK = trK[idx];
    double l_chi = chi[idx];

    double l_rho = rho[idx];
    double l_Sx = Sx[idx]; double l_Sy = Sy[idx]; double l_Sz = Sz[idx];
    double l_Sxx = Sxx[idx]; double l_Sxy = Sxy_f[idx]; double l_Sxz = Sxz_f[idx];
    double l_Syy = Syy[idx]; double l_Syz = Syz_f[idx]; double l_Szz = Szz[idx];

    // ===================================================================
    // beta first derivatives
    // ===================================================================
    double betaxx, betaxy, betaxz;
    d_point_fderivs_shc_gw<GW>(ex, pex, betax, &betaxx, &betaxy, &betaxz,
        crho, sigma, Rcrd, ANTI, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double betayx, betayy, betayz;
    d_point_fderivs_shc_gw<GW>(ex, pex, betay, &betayx, &betayy, &betayz,
        crho, sigma, Rcrd, SYM, ANTI, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double betazx, betazy, betazz;
    d_point_fderivs_shc_gw<GW>(ex, pex, betaz, &betazx, &betazy, &betazz,
        crho, sigma, Rcrd, SYM, SYM, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double div_beta = betaxx + betayy + betazz;

    // ===================================================================
    // chi first derivatives -> chi_rhs
    // ===================================================================
    double chix, chiy, chiz;
    d_point_fderivs_shc_gw<GW>(ex, pex, chi, &chix, &chiy, &chiz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);
    chi_rhs[idx] = F2o3 * l_chin1 * (l_alpn1 * l_trK - div_beta);

    // ===================================================================
    // Metric first derivatives
    // ===================================================================
    double gxxx, gxxy, gxxz, gxyx, gxyy, gxyz, gxzx, gxzy, gxzz;
    double gyyx, gyyy, gyyz, gyzx, gyzy, gyzz, gzzx, gzzy, gzzz;

    d_point_fderivs_shc_gw<GW>(ex, pex, dxx, &gxxx, &gxxy, &gxxz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);
    d_point_fderivs_shc_gw<GW>(ex, pex, gxy_f, &gxyx, &gxyy, &gxyz,
        crho, sigma, Rcrd, ANTI, ANTI, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);
    d_point_fderivs_shc_gw<GW>(ex, pex, gxz_f, &gxzx, &gxzy, &gxzz,
        crho, sigma, Rcrd, ANTI, SYM, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);
    d_point_fderivs_shc_gw<GW>(ex, pex, dyy, &gyyx, &gyyy, &gyyz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);
    d_point_fderivs_shc_gw<GW>(ex, pex, gyz_f, &gyzx, &gyzy, &gyzz,
        crho, sigma, Rcrd, SYM, ANTI, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);
    d_point_fderivs_shc_gw<GW>(ex, pex, dzz, &gzzx, &gzzy, &gzzz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    // ===================================================================
    // Metric inversion
    // ===================================================================
    double det = l_gxx * l_gyy * l_gzz + TWO * l_gxy * l_gyz * l_gxz
               - l_gxz * l_gyy * l_gxz - l_gxx * l_gyz * l_gyz - l_gzz * l_gxy * l_gxy;
    double inv_det = ONE / det;

    double gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) * inv_det;
    double gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) * inv_det;
    double gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) * inv_det;
    double gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) * inv_det;
    double gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) * inv_det;
    double gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) * inv_det;

    // ===================================================================
    // Gmx_Res/Gmy_Res/Gmz_Res (only when co == 0)
    // ===================================================================
    if (co == 0) {
        Gmx_Res[idx] = l_Gamx - (
            gupxx * (gupxx * gxxx + gupxy * gxyx + gupxz * gxzx) +
            gupxy * (gupxx * gxyx + gupxy * gyyx + gupxz * gyzx) +
            gupxz * (gupxx * gxzx + gupxy * gyzx + gupxz * gzzx) +
            gupxx * (gupxy * gxxy + gupyy * gxyy + gupyz * gxzy) +
            gupxy * (gupxy * gxyy + gupyy * gyyy + gupyz * gyzy) +
            gupxz * (gupxy * gxzy + gupyy * gyzy + gupyz * gzzy) +
            gupxx * (gupxz * gxxz + gupyz * gxyz + gupzz * gxzz) +
            gupxy * (gupxz * gxyz + gupyz * gyyz + gupzz * gyzz) +
            gupxz * (gupxz * gxzz + gupyz * gyzz + gupzz * gzzz)
        );
        Gmy_Res[idx] = l_Gamy - (
            gupxx * (gupxy * gxxx + gupyy * gxyx + gupyz * gxzx) +
            gupxy * (gupxy * gxyx + gupyy * gyyx + gupyz * gyzx) +
            gupxz * (gupxy * gxzx + gupyy * gyzx + gupyz * gzzx) +
            gupxy * (gupxy * gxxy + gupyy * gxyy + gupyz * gxzy) +
            gupyy * (gupxy * gxyy + gupyy * gyyy + gupyz * gyzy) +
            gupyz * (gupxy * gxzy + gupyy * gyzy + gupyz * gzzy) +
            gupxy * (gupxz * gxxz + gupyz * gxyz + gupzz * gxzz) +
            gupyy * (gupxz * gxyz + gupyz * gyyz + gupzz * gyzz) +
            gupyz * (gupxz * gxzz + gupyz * gyzz + gupzz * gzzz)
        );
        Gmz_Res[idx] = l_Gamz - (
            gupxx * (gupxz * gxxx + gupyz * gxyx + gupzz * gxzx) +
            gupxy * (gupxz * gxyx + gupyz * gyyx + gupzz * gyzx) +
            gupxz * (gupxz * gxzx + gupyz * gyzx + gupzz * gzzx) +
            gupxy * (gupxz * gxxy + gupyz * gxyy + gupzz * gxzy) +
            gupyy * (gupxz * gxyy + gupyz * gyyy + gupzz * gyzy) +
            gupyz * (gupxz * gxzy + gupyz * gyzy + gupzz * gzzy) +
            gupxz * (gupxz * gxxz + gupyz * gxyz + gupzz * gxzz) +
            gupyz * (gupxz * gxyz + gupyz * gyyz + gupzz * gyzz) +
            gupzz * (gupxz * gxzz + gupyz * gyzz + gupzz * gzzz)
        );
    }

    // ===================================================================
    // gij_rhs
    // ===================================================================
    gxx_rhs[idx] = -TWO * l_alpn1 * l_Axx - F2o3 * l_gxx * div_beta
                 + TWO * (l_gxx * betaxx + l_gxy * betayx + l_gxz * betazx);
    gyy_rhs[idx] = -TWO * l_alpn1 * l_Ayy - F2o3 * l_gyy * div_beta
                 + TWO * (l_gxy * betaxy + l_gyy * betayy + l_gyz * betazy);
    gzz_rhs[idx] = -TWO * l_alpn1 * l_Azz - F2o3 * l_gzz * div_beta
                 + TWO * (l_gxz * betaxz + l_gyz * betayz + l_gzz * betazz);
    gxy_rhs[idx] = -TWO * l_alpn1 * l_Axy + F1o3 * l_gxy * div_beta
                 + l_gxx * betaxy + l_gxz * betazy
                 + l_gyy * betayx + l_gyz * betazx - l_gxy * betazz;
    gyz_rhs[idx] = -TWO * l_alpn1 * l_Ayz + F1o3 * l_gyz * div_beta
                 + l_gxy * betaxz + l_gyy * betayz
                 + l_gxz * betaxy + l_gzz * betazy - l_gyz * betaxx;
    gxz_rhs[idx] = -TWO * l_alpn1 * l_Axz + F1o3 * l_gxz * div_beta
                 + l_gxx * betaxz + l_gxy * betayz
                 + l_gyz * betayx + l_gzz * betazx - l_gxz * betayy;

    // ===================================================================
    // Christoffel symbols (second kind)
    // ===================================================================
    double l_Gamxxx = HALF * (gupxx * gxxx + gupxy * (TWO * gxyx - gxxy) + gupxz * (TWO * gxzx - gxxz));
    double l_Gamyxx = HALF * (gupxy * gxxx + gupyy * (TWO * gxyx - gxxy) + gupyz * (TWO * gxzx - gxxz));
    double l_Gamzxx = HALF * (gupxz * gxxx + gupyz * (TWO * gxyx - gxxy) + gupzz * (TWO * gxzx - gxxz));

    double l_Gamxyy = HALF * (gupxx * (TWO * gxyy - gyyx) + gupxy * gyyy + gupxz * (TWO * gyzy - gyyz));
    double l_Gamyyy = HALF * (gupxy * (TWO * gxyy - gyyx) + gupyy * gyyy + gupyz * (TWO * gyzy - gyyz));
    double l_Gamzyy = HALF * (gupxz * (TWO * gxyy - gyyx) + gupyz * gyyy + gupzz * (TWO * gyzy - gyyz));

    double l_Gamxzz = HALF * (gupxx * (TWO * gxzz - gzzx) + gupxy * (TWO * gyzz - gzzy) + gupxz * gzzz);
    double l_Gamyzz = HALF * (gupxy * (TWO * gxzz - gzzx) + gupyy * (TWO * gyzz - gzzy) + gupyz * gzzz);
    double l_Gamzzz = HALF * (gupxz * (TWO * gxzz - gzzx) + gupyz * (TWO * gyzz - gzzy) + gupzz * gzzz);

    double l_Gamxxy = HALF * (gupxx * gxxy + gupxy * gyyx + gupxz * (gxzy + gyzx - gxyz));
    double l_Gamyxy = HALF * (gupxy * gxxy + gupyy * gyyx + gupyz * (gxzy + gyzx - gxyz));
    double l_Gamzxy = HALF * (gupxz * gxxy + gupyz * gyyx + gupzz * (gxzy + gyzx - gxyz));

    double l_Gamxxz = HALF * (gupxx * gxxz + gupxy * (gxyz + gyzx - gxzy) + gupxz * gzzx);
    double l_Gamyxz = HALF * (gupxy * gxxz + gupyy * (gxyz + gyzx - gxzy) + gupyz * gzzx);
    double l_Gamzxz = HALF * (gupxz * gxxz + gupyz * (gxyz + gyzx - gxzy) + gupzz * gzzx);

    double l_Gamxyz = HALF * (gupxx * (gxyz + gxzy - gyzx) + gupxy * gyyz + gupxz * gzzy);
    double l_Gamyyz = HALF * (gupxy * (gxyz + gxzy - gyzx) + gupyy * gyyz + gupyz * gzzy);
    double l_Gamzyz = HALF * (gupxz * (gxyz + gxzy - gyzx) + gupyz * gyyz + gupzz * gzzy);

    // Store Christoffel symbols
    Gamxxx[idx] = l_Gamxxx; Gamxxy[idx] = l_Gamxxy; Gamxxz[idx] = l_Gamxxz;
    Gamxyy[idx] = l_Gamxyy; Gamxyz[idx] = l_Gamxyz; Gamxzz[idx] = l_Gamxzz;
    Gamyxx[idx] = l_Gamyxx; Gamyxy[idx] = l_Gamyxy; Gamyxz[idx] = l_Gamyxz;
    Gamyyy[idx] = l_Gamyyy; Gamyyz[idx] = l_Gamyyz; Gamyzz[idx] = l_Gamyzz;
    Gamzxx[idx] = l_Gamzxx; Gamzxy[idx] = l_Gamzxy; Gamzxz[idx] = l_Gamzxz;
    Gamzyy[idx] = l_Gamzyy; Gamzyz[idx] = l_Gamzyz; Gamzzz[idx] = l_Gamzzz;

    // ===================================================================
    // Raise indices of A_ij -> Aup^ij
    // ===================================================================
    double Aup_xx = gupxx * gupxx * l_Axx + gupxy * gupxy * l_Ayy + gupxz * gupxz * l_Azz
                  + TWO * (gupxx * gupxy * l_Axy + gupxx * gupxz * l_Axz + gupxy * gupxz * l_Ayz);
    double Aup_yy = gupxy * gupxy * l_Axx + gupyy * gupyy * l_Ayy + gupyz * gupyz * l_Azz
                  + TWO * (gupxy * gupyy * l_Axy + gupxy * gupyz * l_Axz + gupyy * gupyz * l_Ayz);
    double Aup_zz = gupxz * gupxz * l_Axx + gupyz * gupyz * l_Ayy + gupzz * gupzz * l_Azz
                  + TWO * (gupxz * gupyz * l_Axy + gupxz * gupzz * l_Axz + gupyz * gupzz * l_Ayz);
    double Aup_xy = gupxx * gupxy * l_Axx + gupxy * gupyy * l_Ayy + gupxz * gupyz * l_Azz
                  + (gupxx * gupyy + gupxy * gupxy) * l_Axy
                  + (gupxx * gupyz + gupxz * gupxy) * l_Axz
                  + (gupxy * gupyz + gupxz * gupyy) * l_Ayz;
    double Aup_xz = gupxx * gupxz * l_Axx + gupxy * gupyz * l_Ayy + gupxz * gupzz * l_Azz
                  + (gupxx * gupyz + gupxy * gupxz) * l_Axy
                  + (gupxx * gupzz + gupxz * gupxz) * l_Axz
                  + (gupxy * gupzz + gupxz * gupyz) * l_Ayz;
    double Aup_yz = gupxy * gupxz * l_Axx + gupyy * gupyz * l_Ayy + gupyz * gupzz * l_Azz
                  + (gupxy * gupyz + gupyy * gupxz) * l_Axy
                  + (gupxy * gupzz + gupyz * gupxz) * l_Axz
                  + (gupyy * gupzz + gupyz * gupyz) * l_Ayz;

    // ===================================================================
    // Lap and trK first derivatives; Gam^i RHS part 1
    // ===================================================================
    double Lapx, Lapy, Lapz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Lap, &Lapx, &Lapy, &Lapz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double Kx, Ky, Kz;
    d_point_fderivs_shc_gw<GW>(ex, pex, trK, &Kx, &Ky, &Kz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double val_Gamx_rhs = -TWO * (Lapx * Aup_xx + Lapy * Aup_xy + Lapz * Aup_xz)
        + TWO * l_alpn1 * (
            -F3o2 * inv_chin1 * (chix * Aup_xx + chiy * Aup_xy + chiz * Aup_xz)
            - gupxx * (F2o3 * Kx + EIGHT * PIO * l_Sx)
            - gupxy * (F2o3 * Ky + EIGHT * PIO * l_Sy)
            - gupxz * (F2o3 * Kz + EIGHT * PIO * l_Sz)
            + l_Gamxxx * Aup_xx + l_Gamxyy * Aup_yy + l_Gamxzz * Aup_zz
            + TWO * (l_Gamxxy * Aup_xy + l_Gamxxz * Aup_xz + l_Gamxyz * Aup_yz)
        );

    double val_Gamy_rhs = -TWO * (Lapx * Aup_xy + Lapy * Aup_yy + Lapz * Aup_yz)
        + TWO * l_alpn1 * (
            -F3o2 * inv_chin1 * (chix * Aup_xy + chiy * Aup_yy + chiz * Aup_yz)
            - gupxy * (F2o3 * Kx + EIGHT * PIO * l_Sx)
            - gupyy * (F2o3 * Ky + EIGHT * PIO * l_Sy)
            - gupyz * (F2o3 * Kz + EIGHT * PIO * l_Sz)
            + l_Gamyxx * Aup_xx + l_Gamyyy * Aup_yy + l_Gamyzz * Aup_zz
            + TWO * (l_Gamyxy * Aup_xy + l_Gamyxz * Aup_xz + l_Gamyyz * Aup_yz)
        );

    double val_Gamz_rhs = -TWO * (Lapx * Aup_xz + Lapy * Aup_yz + Lapz * Aup_zz)
        + TWO * l_alpn1 * (
            -F3o2 * inv_chin1 * (chix * Aup_xz + chiy * Aup_yz + chiz * Aup_zz)
            - gupxz * (F2o3 * Kx + EIGHT * PIO * l_Sx)
            - gupyz * (F2o3 * Ky + EIGHT * PIO * l_Sy)
            - gupzz * (F2o3 * Kz + EIGHT * PIO * l_Sz)
            + l_Gamzxx * Aup_xx + l_Gamzyy * Aup_yy + l_Gamzzz * Aup_zz
            + TWO * (l_Gamzxy * Aup_xy + l_Gamzxz * Aup_xz + l_Gamzyz * Aup_yz)
        );

    // ===================================================================
    // beta second derivatives
    // ===================================================================
    double bx_fxx, bx_fxy, bx_fxz, bx_fyy, bx_fyz, bx_fzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, betax,
        &bx_fxx, &bx_fxy, &bx_fxz, &bx_fyy, &bx_fyz, &bx_fzz,
        crho, sigma, Rcrd, ANTI, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);

    double by_fxx, by_fxy, by_fxz, by_fyy, by_fyz, by_fzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, betay,
        &by_fxx, &by_fxy, &by_fxz, &by_fyy, &by_fyz, &by_fzz,
        crho, sigma, Rcrd, SYM, ANTI, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);

    double bz_fxx, bz_fxy, bz_fxz, bz_fyy, bz_fyz, bz_fzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, betaz,
        &bz_fxx, &bz_fxy, &bz_fxz, &bz_fyy, &bz_fyz, &bz_fzz,
        crho, sigma, Rcrd, SYM, SYM, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);

    double fxx_beta = bx_fxx + by_fxy + bz_fxz;
    double fxy_beta = bx_fxy + by_fyy + bz_fyz;
    double fxz_beta = bx_fxz + by_fyz + bz_fzz;

    // ===================================================================
    // Gamxa, Gamya, Gamza
    // ===================================================================
    double Gamxa = gupxx * l_Gamxxx + gupyy * l_Gamxyy + gupzz * l_Gamxzz
                 + TWO * (gupxy * l_Gamxxy + gupxz * l_Gamxxz + gupyz * l_Gamxyz);
    double Gamya = gupxx * l_Gamyxx + gupyy * l_Gamyyy + gupzz * l_Gamyzz
                 + TWO * (gupxy * l_Gamyxy + gupxz * l_Gamyxz + gupyz * l_Gamyyz);
    double Gamza = gupxx * l_Gamzxx + gupyy * l_Gamzyy + gupzz * l_Gamzzz
                 + TWO * (gupxy * l_Gamzxy + gupxz * l_Gamzxz + gupyz * l_Gamzyz);

    // ===================================================================
    // Gam first derivatives
    // ===================================================================
    double Gamxx, Gamxy, Gamxz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Gamx, &Gamxx, &Gamxy, &Gamxz,
        crho, sigma, Rcrd, ANTI, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double Gamyx, Gamyy, Gamyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Gamy, &Gamyx, &Gamyy, &Gamyz,
        crho, sigma, Rcrd, SYM, ANTI, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double Gamzx, Gamzy, Gamzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Gamz, &Gamzx, &Gamzy, &Gamzz,
        crho, sigma, Rcrd, SYM, SYM, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    // ===================================================================
    // Gam^i RHS part 2
    // ===================================================================
    val_Gamx_rhs = val_Gamx_rhs + F2o3 * Gamxa * div_beta
                 - Gamxa * betaxx - Gamya * betaxy - Gamza * betaxz
                 + F1o3 * (gupxx * fxx_beta + gupxy * fxy_beta + gupxz * fxz_beta)
                 + gupxx * bx_fxx + gupyy * bx_fyy + gupzz * bx_fzz
                 + TWO * (gupxy * bx_fxy + gupxz * bx_fxz + gupyz * bx_fyz);

    val_Gamy_rhs = val_Gamy_rhs + F2o3 * Gamya * div_beta
                 - Gamxa * betayx - Gamya * betayy - Gamza * betayz
                 + F1o3 * (gupxy * fxx_beta + gupyy * fxy_beta + gupyz * fxz_beta)
                 + gupxx * by_fxx + gupyy * by_fyy + gupzz * by_fzz
                 + TWO * (gupxy * by_fxy + gupxz * by_fxz + gupyz * by_fyz);

    val_Gamz_rhs = val_Gamz_rhs + F2o3 * Gamza * div_beta
                 - Gamxa * betazx - Gamya * betazy - Gamza * betazz
                 + F1o3 * (gupxz * fxx_beta + gupyz * fxy_beta + gupzz * fxz_beta)
                 + gupxx * bz_fxx + gupyy * bz_fyy + gupzz * bz_fzz
                 + TWO * (gupxy * bz_fxy + gupxz * bz_fxz + gupyz * bz_fyz);


    // ===================================================================
    // First kind of connection (covariant derivatives of metric)
    // ===================================================================
    double Cgxxx = l_gxx * l_Gamxxx + l_gxy * l_Gamyxx + l_gxz * l_Gamzxx;
    double Cgxyx = l_gxx * l_Gamxxy + l_gxy * l_Gamyxy + l_gxz * l_Gamzxy;
    double Cgxzx = l_gxx * l_Gamxxz + l_gxy * l_Gamyxz + l_gxz * l_Gamzxz;
    double Cgyyx = l_gxx * l_Gamxyy + l_gxy * l_Gamyyy + l_gxz * l_Gamzyy;
    double Cgyzx = l_gxx * l_Gamxyz + l_gxy * l_Gamyyz + l_gxz * l_Gamzyz;
    double Cgzzx = l_gxx * l_Gamxzz + l_gxy * l_Gamyzz + l_gxz * l_Gamzzz;

    double Cgxxy = l_gxy * l_Gamxxx + l_gyy * l_Gamyxx + l_gyz * l_Gamzxx;
    double Cgxyy = l_gxy * l_Gamxxy + l_gyy * l_Gamyxy + l_gyz * l_Gamzxy;
    double Cgxzy = l_gxy * l_Gamxxz + l_gyy * l_Gamyxz + l_gyz * l_Gamzxz;
    double Cgyyy = l_gxy * l_Gamxyy + l_gyy * l_Gamyyy + l_gyz * l_Gamzyy;
    double Cgyzy = l_gxy * l_Gamxyz + l_gyy * l_Gamyyz + l_gyz * l_Gamzyz;
    double Cgzzy = l_gxy * l_Gamxzz + l_gyy * l_Gamyzz + l_gyz * l_Gamzzz;

    double Cgxxz = l_gxz * l_Gamxxx + l_gyz * l_Gamyxx + l_gzz * l_Gamzxx;
    double Cgxyz = l_gxz * l_Gamxxy + l_gyz * l_Gamyxy + l_gzz * l_Gamzxy;
    double Cgxzz = l_gxz * l_Gamxxz + l_gyz * l_Gamyxz + l_gzz * l_Gamzxz;
    double Cgyyz = l_gxz * l_Gamxyy + l_gyz * l_Gamyyy + l_gzz * l_Gamzyy;
    double Cgyzz = l_gxz * l_Gamxyz + l_gyz * l_Gamyyz + l_gzz * l_Gamzyz;
    double Cgzzz = l_gxz * l_Gamxzz + l_gyz * l_Gamyzz + l_gzz * l_Gamzzz;

    // ===================================================================
    // Ricci tensor for conformal metric
    // ===================================================================
    double d2_fxx, d2_fxy, d2_fxz, d2_fyy, d2_fyz, d2_fzz;

    // dxx second derivatives
    d_point_fdderivs_shc_gw<GW>(ex, pex, dxx,
        &d2_fxx, &d2_fxy, &d2_fxz, &d2_fyy, &d2_fyz, &d2_fzz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);
    double contract_xx = gupxx * d2_fxx + gupyy * d2_fyy + gupzz * d2_fzz
                       + (gupxy * d2_fxy + gupxz * d2_fxz + gupyz * d2_fyz) * TWO;

    // dyy second derivatives
    d_point_fdderivs_shc_gw<GW>(ex, pex, dyy,
        &d2_fxx, &d2_fxy, &d2_fxz, &d2_fyy, &d2_fyz, &d2_fzz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);
    double contract_yy = gupxx * d2_fxx + gupyy * d2_fyy + gupzz * d2_fzz
                       + (gupxy * d2_fxy + gupxz * d2_fxz + gupyz * d2_fyz) * TWO;

    // dzz second derivatives
    d_point_fdderivs_shc_gw<GW>(ex, pex, dzz,
        &d2_fxx, &d2_fxy, &d2_fxz, &d2_fyy, &d2_fyz, &d2_fzz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);
    double contract_zz = gupxx * d2_fxx + gupyy * d2_fyy + gupzz * d2_fzz
                       + (gupxy * d2_fxy + gupxz * d2_fxz + gupyz * d2_fyz) * TWO;

    // gxy second derivatives
    d_point_fdderivs_shc_gw<GW>(ex, pex, gxy_f,
        &d2_fxx, &d2_fxy, &d2_fxz, &d2_fyy, &d2_fyz, &d2_fzz,
        crho, sigma, Rcrd, ANTI, ANTI, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);
    double contract_xy = gupxx * d2_fxx + gupyy * d2_fyy + gupzz * d2_fzz
                       + (gupxy * d2_fxy + gupxz * d2_fxz + gupyz * d2_fyz) * TWO;

    // gxz second derivatives
    d_point_fdderivs_shc_gw<GW>(ex, pex, gxz_f,
        &d2_fxx, &d2_fxy, &d2_fxz, &d2_fyy, &d2_fyz, &d2_fzz,
        crho, sigma, Rcrd, ANTI, SYM, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);
    double contract_xz = gupxx * d2_fxx + gupyy * d2_fyy + gupzz * d2_fzz
                       + (gupxy * d2_fxy + gupxz * d2_fxz + gupyz * d2_fyz) * TWO;

    // gyz second derivatives
    d_point_fdderivs_shc_gw<GW>(ex, pex, gyz_f,
        &d2_fxx, &d2_fxy, &d2_fxz, &d2_fyy, &d2_fyz, &d2_fzz,
        crho, sigma, Rcrd, SYM, ANTI, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);
    double contract_yz = gupxx * d2_fxx + gupyy * d2_fyy + gupzz * d2_fzz
                       + (gupxy * d2_fxy + gupxz * d2_fxz + gupyz * d2_fyz) * TWO;


    // ===================================================================
    // Ricci tensor assembly (conformal part)
    // ===================================================================
    double l_Rxx = -HALF * contract_xx
                 + l_gxx * Gamxx + l_gxy * Gamyx + l_gxz * Gamzx
                 + Gamxa * Cgxxx + Gamya * Cgxyx + Gamza * Cgxzx
                 + gupxx * (
                     TWO * (l_Gamxxx * Cgxxx + l_Gamyxx * Cgxyx + l_Gamzxx * Cgxzx)
                     + l_Gamxxx * Cgxxx + l_Gamyxx * Cgxxy + l_Gamzxx * Cgxxz)
                 + gupxy * (
                     TWO * (l_Gamxxx * Cgxyx + l_Gamyxx * Cgyyx + l_Gamzxx * Cgyzx
                          + l_Gamxxy * Cgxxx + l_Gamyxy * Cgxyx + l_Gamzxy * Cgxzx)
                     + l_Gamxxy * Cgxxx + l_Gamyxy * Cgxxy + l_Gamzxy * Cgxxz
                     + l_Gamxxx * Cgxyx + l_Gamyxx * Cgxyy + l_Gamzxx * Cgxyz)
                 + gupxz * (
                     TWO * (l_Gamxxx * Cgxzx + l_Gamyxx * Cgyzx + l_Gamzxx * Cgzzx
                          + l_Gamxxz * Cgxxx + l_Gamyxz * Cgxyx + l_Gamzxz * Cgxzx)
                     + l_Gamxxz * Cgxxx + l_Gamyxz * Cgxxy + l_Gamzxz * Cgxxz
                     + l_Gamxxx * Cgxzx + l_Gamyxx * Cgxzy + l_Gamzxx * Cgxzz)
                 + gupyy * (
                     TWO * (l_Gamxxy * Cgxyx + l_Gamyxy * Cgyyx + l_Gamzxy * Cgyzx)
                     + l_Gamxxy * Cgxyx + l_Gamyxy * Cgxyy + l_Gamzxy * Cgxyz)
                 + gupyz * (
                     TWO * (l_Gamxxy * Cgxzx + l_Gamyxy * Cgyzx + l_Gamzxy * Cgzzx
                          + l_Gamxxz * Cgxyx + l_Gamyxz * Cgyyx + l_Gamzxz * Cgyzx)
                     + l_Gamxxz * Cgxyx + l_Gamyxz * Cgxyy + l_Gamzxz * Cgxyz
                     + l_Gamxxy * Cgxzx + l_Gamyxy * Cgxzy + l_Gamzxy * Cgxzz)
                 + gupzz * (
                     TWO * (l_Gamxxz * Cgxzx + l_Gamyxz * Cgyzx + l_Gamzxz * Cgzzx)
                     + l_Gamxxz * Cgxzx + l_Gamyxz * Cgxzy + l_Gamzxz * Cgxzz);

    double l_Ryy = -HALF * contract_yy
                 + l_gxy * Gamxy + l_gyy * Gamyy + l_gyz * Gamzy
                 + Gamxa * Cgxyy + Gamya * Cgyyy + Gamza * Cgyzy
                 + gupxx * (
                     TWO * (l_Gamxxy * Cgxxy + l_Gamyxy * Cgxyy + l_Gamzxy * Cgxzy)
                     + l_Gamxxy * Cgxyx + l_Gamyxy * Cgxyy + l_Gamzxy * Cgxyz)
                 + gupxy * (
                     TWO * (l_Gamxxy * Cgxyy + l_Gamyxy * Cgyyy + l_Gamzxy * Cgyzy
                          + l_Gamxyy * Cgxxy + l_Gamyyy * Cgxyy + l_Gamzyy * Cgxzy)
                     + l_Gamxyy * Cgxyx + l_Gamyyy * Cgxyy + l_Gamzyy * Cgxyz
                     + l_Gamxxy * Cgyyx + l_Gamyxy * Cgyyy + l_Gamzxy * Cgyyz)
                 + gupxz * (
                     TWO * (l_Gamxxy * Cgxzy + l_Gamyxy * Cgyzy + l_Gamzxy * Cgzzy
                          + l_Gamxyz * Cgxxy + l_Gamyyz * Cgxyy + l_Gamzyz * Cgxzy)
                     + l_Gamxyz * Cgxyx + l_Gamyyz * Cgxyy + l_Gamzyz * Cgxyz
                     + l_Gamxxy * Cgyzx + l_Gamyxy * Cgyzy + l_Gamzxy * Cgyzz)
                 + gupyy * (
                     TWO * (l_Gamxyy * Cgxyy + l_Gamyyy * Cgyyy + l_Gamzyy * Cgyzy)
                     + l_Gamxyy * Cgyyx + l_Gamyyy * Cgyyy + l_Gamzyy * Cgyyz)
                 + gupyz * (
                     TWO * (l_Gamxyy * Cgxzy + l_Gamyyy * Cgyzy + l_Gamzyy * Cgzzy
                          + l_Gamxyz * Cgxyy + l_Gamyyz * Cgyyy + l_Gamzyz * Cgyzy)
                     + l_Gamxyz * Cgyyx + l_Gamyyz * Cgyyy + l_Gamzyz * Cgyyz
                     + l_Gamxyy * Cgyzx + l_Gamyyy * Cgyzy + l_Gamzyy * Cgyzz)
                 + gupzz * (
                     TWO * (l_Gamxyz * Cgxzy + l_Gamyyz * Cgyzy + l_Gamzyz * Cgzzy)
                     + l_Gamxyz * Cgyzx + l_Gamyyz * Cgyzy + l_Gamzyz * Cgyzz);

    double l_Rzz = -HALF * contract_zz
                 + l_gxz * Gamxz + l_gyz * Gamyz + l_gzz * Gamzz
                 + Gamxa * Cgxzz + Gamya * Cgyzz + Gamza * Cgzzz
                 + gupxx * (
                     TWO * (l_Gamxxz * Cgxxz + l_Gamyxz * Cgxyz + l_Gamzxz * Cgxzz)
                     + l_Gamxxz * Cgxzx + l_Gamyxz * Cgxzy + l_Gamzxz * Cgxzz)
                 + gupxy * (
                     TWO * (l_Gamxxz * Cgxyz + l_Gamyxz * Cgyyz + l_Gamzxz * Cgyzz
                          + l_Gamxyz * Cgxxz + l_Gamyyz * Cgxyz + l_Gamzyz * Cgxzz)
                     + l_Gamxyz * Cgxzx + l_Gamyyz * Cgxzy + l_Gamzyz * Cgxzz
                     + l_Gamxxz * Cgyzx + l_Gamyxz * Cgyzy + l_Gamzxz * Cgyzz)
                 + gupxz * (
                     TWO * (l_Gamxxz * Cgxzz + l_Gamyxz * Cgyzz + l_Gamzxz * Cgzzz
                          + l_Gamxzz * Cgxxz + l_Gamyzz * Cgxyz + l_Gamzzz * Cgxzz)
                     + l_Gamxzz * Cgxzx + l_Gamyzz * Cgxzy + l_Gamzzz * Cgxzz
                     + l_Gamxxz * Cgzzx + l_Gamyxz * Cgzzy + l_Gamzxz * Cgzzz)
                 + gupyy * (
                     TWO * (l_Gamxyz * Cgxyz + l_Gamyyz * Cgyyz + l_Gamzyz * Cgyzz)
                     + l_Gamxyz * Cgyzx + l_Gamyyz * Cgyzy + l_Gamzyz * Cgyzz)
                 + gupyz * (
                     TWO * (l_Gamxyz * Cgxzz + l_Gamyyz * Cgyzz + l_Gamzyz * Cgzzz
                          + l_Gamxzz * Cgxyz + l_Gamyzz * Cgyyz + l_Gamzzz * Cgyzz)
                     + l_Gamxzz * Cgyzx + l_Gamyzz * Cgyzy + l_Gamzzz * Cgyzz
                     + l_Gamxyz * Cgzzx + l_Gamyyz * Cgzzy + l_Gamzyz * Cgzzz)
                 + gupzz * (
                     TWO * (l_Gamxzz * Cgxzz + l_Gamyzz * Cgyzz + l_Gamzzz * Cgzzz)
                     + l_Gamxzz * Cgzzx + l_Gamyzz * Cgzzy + l_Gamzzz * Cgzzz);

    // Rxy
    double l_Rxy = HALF * (-contract_xy
                 + l_gxx * Gamxy + l_gxy * Gamyy + l_gxz * Gamzy
                 + l_gxy * Gamxx + l_gyy * Gamyx + l_gyz * Gamzx
                 + Gamxa * Cgxyx + Gamya * Cgyyx + Gamza * Cgyzx
                 + Gamxa * Cgxxy + Gamya * Cgxyy + Gamza * Cgxzy)
                 + gupxx * (
                     l_Gamxxx * Cgxxy + l_Gamyxx * Cgxyy + l_Gamzxx * Cgxzy
                   + l_Gamxxy * Cgxxx + l_Gamyxy * Cgxyx + l_Gamzxy * Cgxzx
                   + l_Gamxxx * Cgxyx + l_Gamyxx * Cgxyy + l_Gamzxx * Cgxyz)
                 + gupxy * (
                     l_Gamxxx * Cgxyy + l_Gamyxx * Cgyyy + l_Gamzxx * Cgyzy
                   + l_Gamxxy * Cgxyx + l_Gamyxy * Cgyyx + l_Gamzxy * Cgyzx
                   + l_Gamxxy * Cgxyx + l_Gamyxy * Cgxyy + l_Gamzxy * Cgxyz
                   + l_Gamxxy * Cgxxy + l_Gamyxy * Cgxyy + l_Gamzxy * Cgxzy
                   + l_Gamxyy * Cgxxx + l_Gamyyy * Cgxyx + l_Gamzyy * Cgxzx
                   + l_Gamxxx * Cgyyx + l_Gamyxx * Cgyyy + l_Gamzxx * Cgyyz)
                 + gupxz * (
                     l_Gamxxx * Cgxzy + l_Gamyxx * Cgyzy + l_Gamzxx * Cgzzy
                   + l_Gamxxy * Cgxzx + l_Gamyxy * Cgyzx + l_Gamzxy * Cgzzx
                   + l_Gamxxz * Cgxyx + l_Gamyxz * Cgxyy + l_Gamzxz * Cgxyz
                   + l_Gamxxz * Cgxxy + l_Gamyxz * Cgxyy + l_Gamzxz * Cgxzy
                   + l_Gamxyz * Cgxxx + l_Gamyyz * Cgxyx + l_Gamzyz * Cgxzx
                   + l_Gamxxx * Cgyzx + l_Gamyxx * Cgyzy + l_Gamzxx * Cgyzz)
                 + gupyy * (
                     l_Gamxxy * Cgxyy + l_Gamyxy * Cgyyy + l_Gamzxy * Cgyzy
                   + l_Gamxyy * Cgxyx + l_Gamyyy * Cgyyx + l_Gamzyy * Cgyzx
                   + l_Gamxxy * Cgyyx + l_Gamyxy * Cgyyy + l_Gamzxy * Cgyyz)
                 + gupyz * (
                     l_Gamxxy * Cgxzy + l_Gamyxy * Cgyzy + l_Gamzxy * Cgzzy
                   + l_Gamxyy * Cgxzx + l_Gamyyy * Cgyzx + l_Gamzyy * Cgzzx
                   + l_Gamxxz * Cgyyx + l_Gamyxz * Cgyyy + l_Gamzxz * Cgyyz
                   + l_Gamxxz * Cgxyy + l_Gamyxz * Cgyyy + l_Gamzxz * Cgyzy
                   + l_Gamxyz * Cgxyx + l_Gamyyz * Cgyyx + l_Gamzyz * Cgyzx
                   + l_Gamxxy * Cgyzx + l_Gamyxy * Cgyzy + l_Gamzxy * Cgyzz)
                 + gupzz * (
                     l_Gamxxz * Cgxzy + l_Gamyxz * Cgyzy + l_Gamzxz * Cgzzy
                   + l_Gamxyz * Cgxzx + l_Gamyyz * Cgyzx + l_Gamzyz * Cgzzx
                   + l_Gamxxz * Cgyzx + l_Gamyxz * Cgyzy + l_Gamzxz * Cgyzz);

    // Rxz
    double l_Rxz = HALF * (-contract_xz
                 + l_gxx * Gamxz + l_gxy * Gamyz + l_gxz * Gamzz
                 + l_gxz * Gamxx + l_gyz * Gamyx + l_gzz * Gamzx
                 + Gamxa * Cgxzx + Gamya * Cgyzx + Gamza * Cgzzx
                 + Gamxa * Cgxxz + Gamya * Cgxyz + Gamza * Cgxzz)
                 + gupxx * (
                     l_Gamxxx * Cgxxz + l_Gamyxx * Cgxyz + l_Gamzxx * Cgxzz
                   + l_Gamxxz * Cgxxx + l_Gamyxz * Cgxyx + l_Gamzxz * Cgxzx
                   + l_Gamxxx * Cgxzx + l_Gamyxx * Cgxzy + l_Gamzxx * Cgxzz)
                 + gupxy * (
                     l_Gamxxx * Cgxyz + l_Gamyxx * Cgyyz + l_Gamzxx * Cgyzz
                   + l_Gamxxz * Cgxyx + l_Gamyxz * Cgyyx + l_Gamzxz * Cgyzx
                   + l_Gamxxy * Cgxzx + l_Gamyxy * Cgxzy + l_Gamzxy * Cgxzz
                   + l_Gamxxy * Cgxxz + l_Gamyxy * Cgxyz + l_Gamzxy * Cgxzz
                   + l_Gamxyz * Cgxxx + l_Gamyyz * Cgxyx + l_Gamzyz * Cgxzx
                   + l_Gamxxx * Cgyzx + l_Gamyxx * Cgyzy + l_Gamzxx * Cgyzz)
                 + gupxz * (
                     l_Gamxxx * Cgxzz + l_Gamyxx * Cgyzz + l_Gamzxx * Cgzzz
                   + l_Gamxxz * Cgxzx + l_Gamyxz * Cgyzx + l_Gamzxz * Cgzzx
                   + l_Gamxxz * Cgxzx + l_Gamyxz * Cgxzy + l_Gamzxz * Cgxzz
                   + l_Gamxxz * Cgxxz + l_Gamyxz * Cgxyz + l_Gamzxz * Cgxzz
                   + l_Gamxzz * Cgxxx + l_Gamyzz * Cgxyx + l_Gamzzz * Cgxzx
                   + l_Gamxxx * Cgzzx + l_Gamyxx * Cgzzy + l_Gamzxx * Cgzzz)
                 + gupyy * (
                     l_Gamxxy * Cgxyz + l_Gamyxy * Cgyyz + l_Gamzxy * Cgyzz
                   + l_Gamxyz * Cgxyx + l_Gamyyz * Cgyyx + l_Gamzyz * Cgyzx
                   + l_Gamxxy * Cgyzx + l_Gamyxy * Cgyzy + l_Gamzxy * Cgyzz)
                 + gupyz * (
                     l_Gamxxy * Cgxzz + l_Gamyxy * Cgyzz + l_Gamzxy * Cgzzz
                   + l_Gamxyz * Cgxzx + l_Gamyyz * Cgyzx + l_Gamzyz * Cgzzx
                   + l_Gamxxz * Cgyzx + l_Gamyxz * Cgyzy + l_Gamzxz * Cgyzz
                   + l_Gamxxz * Cgxyz + l_Gamyxz * Cgyyz + l_Gamzxz * Cgyzz
                   + l_Gamxzz * Cgxyx + l_Gamyzz * Cgyyx + l_Gamzzz * Cgyzx
                   + l_Gamxxy * Cgzzx + l_Gamyxy * Cgzzy + l_Gamzxy * Cgzzz)
                 + gupzz * (
                     l_Gamxxz * Cgxzz + l_Gamyxz * Cgyzz + l_Gamzxz * Cgzzz
                   + l_Gamxzz * Cgxzx + l_Gamyzz * Cgyzx + l_Gamzzz * Cgzzx
                   + l_Gamxxz * Cgzzx + l_Gamyxz * Cgzzy + l_Gamzxz * Cgzzz);

    // Ryz
    double l_Ryz = HALF * (-contract_yz
                 + l_gxy * Gamxz + l_gyy * Gamyz + l_gyz * Gamzz
                 + l_gxz * Gamxy + l_gyz * Gamyy + l_gzz * Gamzy
                 + Gamxa * Cgxzy + Gamya * Cgyzy + Gamza * Cgzzy
                 + Gamxa * Cgxyz + Gamya * Cgyyz + Gamza * Cgyzz)
                 + gupxx * (
                     l_Gamxxy * Cgxxz + l_Gamyxy * Cgxyz + l_Gamzxy * Cgxzz
                   + l_Gamxxz * Cgxxy + l_Gamyxz * Cgxyy + l_Gamzxz * Cgxzy
                   + l_Gamxxy * Cgxzx + l_Gamyxy * Cgxzy + l_Gamzxy * Cgxzz)
                 + gupxy * (
                     l_Gamxxy * Cgxyz + l_Gamyxy * Cgyyz + l_Gamzxy * Cgyzz
                   + l_Gamxxz * Cgxyy + l_Gamyxz * Cgyyy + l_Gamzxz * Cgyzy
                   + l_Gamxyy * Cgxzx + l_Gamyyy * Cgxzy + l_Gamzyy * Cgxzz
                   + l_Gamxyy * Cgxxz + l_Gamyyy * Cgxyz + l_Gamzyy * Cgxzz
                   + l_Gamxyz * Cgxxy + l_Gamyyz * Cgxyy + l_Gamzyz * Cgxzy
                   + l_Gamxxy * Cgyzx + l_Gamyxy * Cgyzy + l_Gamzxy * Cgyzz)
                 + gupxz * (
                     l_Gamxxy * Cgxzz + l_Gamyxy * Cgyzz + l_Gamzxy * Cgzzz
                   + l_Gamxxz * Cgxzy + l_Gamyxz * Cgyzy + l_Gamzxz * Cgzzy
                   + l_Gamxyz * Cgxzx + l_Gamyyz * Cgxzy + l_Gamzyz * Cgxzz
                   + l_Gamxyz * Cgxxz + l_Gamyyz * Cgxyz + l_Gamzyz * Cgxzz
                   + l_Gamxzz * Cgxxy + l_Gamyzz * Cgxyy + l_Gamzzz * Cgxzy
                   + l_Gamxxy * Cgzzx + l_Gamyxy * Cgzzy + l_Gamzxy * Cgzzz)
                 + gupyy * (
                     l_Gamxyy * Cgxyz + l_Gamyyy * Cgyyz + l_Gamzyy * Cgyzz
                   + l_Gamxyz * Cgxyy + l_Gamyyz * Cgyyy + l_Gamzyz * Cgyzy
                   + l_Gamxyy * Cgyzx + l_Gamyyy * Cgyzy + l_Gamzyy * Cgyzz)
                 + gupyz * (
                     l_Gamxyy * Cgxzz + l_Gamyyy * Cgyzz + l_Gamzyy * Cgzzz
                   + l_Gamxyz * Cgxzy + l_Gamyyz * Cgyzy + l_Gamzyz * Cgzzy
                   + l_Gamxyz * Cgyzx + l_Gamyyz * Cgyzy + l_Gamzyz * Cgyzz
                   + l_Gamxyz * Cgxyz + l_Gamyyz * Cgyyz + l_Gamzyz * Cgyzz
                   + l_Gamxzz * Cgxyy + l_Gamyzz * Cgyyy + l_Gamzzz * Cgyzy
                   + l_Gamxyy * Cgzzx + l_Gamyyy * Cgzzy + l_Gamzyy * Cgzzz)
                 + gupzz * (
                     l_Gamxyz * Cgxzz + l_Gamyyz * Cgyzz + l_Gamzyz * Cgzzz
                   + l_Gamxzz * Cgxzy + l_Gamyzz * Cgyzy + l_Gamzzz * Cgzzy
                   + l_Gamxyz * Cgzzx + l_Gamyyz * Cgzzy + l_Gamzyz * Cgzzz);


    // ===================================================================
    // Chi second derivatives and contributions to Ricci
    // ===================================================================
    double chixx, chixy, chixz, chiyy, chiyz, chizz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, chi,
        &chixx, &chixy, &chixz, &chiyy, &chiyz, &chizz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);

    double D2chi_xx = chixx - l_Gamxxx * chix - l_Gamyxx * chiy - l_Gamzxx * chiz;
    double D2chi_xy = chixy - l_Gamxxy * chix - l_Gamyxy * chiy - l_Gamzxy * chiz;
    double D2chi_xz = chixz - l_Gamxxz * chix - l_Gamyxz * chiy - l_Gamzxz * chiz;
    double D2chi_yy = chiyy - l_Gamxyy * chix - l_Gamyyy * chiy - l_Gamzyy * chiz;
    double D2chi_yz = chiyz - l_Gamxyz * chix - l_Gamyyz * chiy - l_Gamzyz * chiz;
    double D2chi_zz = chizz - l_Gamxzz * chix - l_Gamyzz * chiy - l_Gamzzz * chiz;

    double f_chi_scalar = gupxx * (D2chi_xx - F3o2 * inv_chin1 * chix * chix)
                        + gupyy * (D2chi_yy - F3o2 * inv_chin1 * chiy * chiy)
                        + gupzz * (D2chi_zz - F3o2 * inv_chin1 * chiz * chiz)
         + TWO * gupxy * (D2chi_xy - F3o2 * inv_chin1 * chix * chiy)
         + TWO * gupxz * (D2chi_xz - F3o2 * inv_chin1 * chix * chiz)
         + TWO * gupyz * (D2chi_yz - F3o2 * inv_chin1 * chiy * chiz);

    l_Rxx += (D2chi_xx - chix * chix * inv_chin1 * HALF + l_gxx * f_chi_scalar) * inv_chin1 * HALF;
    l_Ryy += (D2chi_yy - chiy * chiy * inv_chin1 * HALF + l_gyy * f_chi_scalar) * inv_chin1 * HALF;
    l_Rzz += (D2chi_zz - chiz * chiz * inv_chin1 * HALF + l_gzz * f_chi_scalar) * inv_chin1 * HALF;
    l_Rxy += (D2chi_xy - chix * chiy * inv_chin1 * HALF + l_gxy * f_chi_scalar) * inv_chin1 * HALF;
    l_Rxz += (D2chi_xz - chix * chiz * inv_chin1 * HALF + l_gxz * f_chi_scalar) * inv_chin1 * HALF;
    l_Ryz += (D2chi_yz - chiy * chiz * inv_chin1 * HALF + l_gyz * f_chi_scalar) * inv_chin1 * HALF;

    // ===================================================================
    // Physical Christoffel symbols
    // ===================================================================
    double gup_chix = gupxx * chix + gupxy * chiy + gupxz * chiz;
    double gup_chiy = gupxy * chix + gupyy * chiy + gupyz * chiz;
    double gup_chiz = gupxz * chix + gupyz * chiy + gupzz * chiz;
    double gx_phy = gup_chix * inv_chin1;
    double gy_phy = gup_chiy * inv_chin1;
    double gz_phy = gup_chiz * inv_chin1;

    double phy_Gamxxx = l_Gamxxx - (HALF * TWO * chix * inv_chin1 - HALF * l_gxx * gx_phy);
    double phy_Gamyxx = l_Gamyxx - (-HALF * l_gxx * gy_phy);
    double phy_Gamzxx = l_Gamzxx - (-HALF * l_gxx * gz_phy);
    double phy_Gamxyy = l_Gamxyy - (-HALF * l_gyy * gx_phy);
    double phy_Gamyyy = l_Gamyyy - (HALF * TWO * chiy * inv_chin1 - HALF * l_gyy * gy_phy);
    double phy_Gamzyy = l_Gamzyy - (-HALF * l_gyy * gz_phy);
    double phy_Gamxzz = l_Gamxzz - (-HALF * l_gzz * gx_phy);
    double phy_Gamyzz = l_Gamyzz - (-HALF * l_gzz * gy_phy);
    double phy_Gamzzz = l_Gamzzz - (HALF * TWO * chiz * inv_chin1 - HALF * l_gzz * gz_phy);
    double phy_Gamxxy = l_Gamxxy - (HALF * chiy * inv_chin1 - HALF * l_gxy * gx_phy);
    double phy_Gamyxy = l_Gamyxy - (HALF * chix * inv_chin1 - HALF * l_gxy * gy_phy);
    double phy_Gamzxy = l_Gamzxy - (-HALF * l_gxy * gz_phy);
    double phy_Gamxxz = l_Gamxxz - (HALF * chiz * inv_chin1 - HALF * l_gxz * gx_phy);
    double phy_Gamyxz = l_Gamyxz - (-HALF * l_gxz * gy_phy);
    double phy_Gamzxz = l_Gamzxz - (HALF * chix * inv_chin1 - HALF * l_gxz * gz_phy);
    double phy_Gamxyz = l_Gamxyz - (-HALF * l_gyz * gx_phy);
    double phy_Gamyyz = l_Gamyyz - (HALF * chiz * inv_chin1 - HALF * l_gyz * gy_phy);
    double phy_Gamzyz = l_Gamzyz - (HALF * chiy * inv_chin1 - HALF * l_gyz * gz_phy);

    // ===================================================================
    // Second derivatives of Lap -> D^i D_i Lap
    // ===================================================================
    double Lapxx, Lapxy, Lapxz, Lapyy, Lapyz, Lapzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, Lap,
        &Lapxx, &Lapxy, &Lapxz, &Lapyy, &Lapyz, &Lapzz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);

    double D2Lap_xx = Lapxx - l_Gamxxx * Lapx - l_Gamyxx * Lapy - l_Gamzxx * Lapz;
    double D2Lap_yy = Lapyy - l_Gamxyy * Lapx - l_Gamyyy * Lapy - l_Gamzyy * Lapz;
    double D2Lap_zz = Lapzz - l_Gamxzz * Lapx - l_Gamyzz * Lapy - l_Gamzzz * Lapz;
    double D2Lap_xy = Lapxy - l_Gamxxy * Lapx - l_Gamyxy * Lapy - l_Gamzxy * Lapz;
    double D2Lap_xz = Lapxz - l_Gamxxz * Lapx - l_Gamyxz * Lapy - l_Gamzxz * Lapz;
    double D2Lap_yz = Lapyz - l_Gamxyz * Lapx - l_Gamyyz * Lapy - l_Gamzyz * Lapz;

    double chi_dot_Lap = gx_phy * Lapx + gy_phy * Lapy + gz_phy * Lapz;

    double fxx_Lap_phy = D2Lap_xx + chix * Lapx * inv_chin1 - HALF * l_gxx * chi_dot_Lap;
    double fyy_Lap_phy = D2Lap_yy + chiy * Lapy * inv_chin1 - HALF * l_gyy * chi_dot_Lap;
    double fzz_Lap_phy = D2Lap_zz + chiz * Lapz * inv_chin1 - HALF * l_gzz * chi_dot_Lap;
    double fxy_Lap_phy = D2Lap_xy + HALF * (chix * Lapy + chiy * Lapx) * inv_chin1 - HALF * l_gxy * chi_dot_Lap;
    double fxz_Lap_phy = D2Lap_xz + HALF * (chix * Lapz + chiz * Lapx) * inv_chin1 - HALF * l_gxz * chi_dot_Lap;
    double fyz_Lap_phy = D2Lap_yz + HALF * (chiy * Lapz + chiz * Lapy) * inv_chin1 - HALF * l_gyz * chi_dot_Lap;

    double conf_trK_rhs = gupxx * D2Lap_xx + gupyy * D2Lap_yy + gupzz * D2Lap_zz
                        + TWO * (gupxy * D2Lap_xy + gupxz * D2Lap_xz + gupyz * D2Lap_yz);
    double trK_rhs_conf = conf_trK_rhs - HALF * chi_dot_Lap;

    // ===================================================================
    // A_ij RHS: trace-free part of (alpn1 * (R_ij - 8pi S_ij) - D_i D_j Lap)
    // ===================================================================
    double fxx_P = l_alpn1 * (l_Rxx - EIGHT * PIO * l_Sxx) - fxx_Lap_phy;
    double fyy_P = l_alpn1 * (l_Ryy - EIGHT * PIO * l_Syy) - fyy_Lap_phy;
    double fzz_P = l_alpn1 * (l_Rzz - EIGHT * PIO * l_Szz) - fzz_Lap_phy;
    double fxy_P = l_alpn1 * (l_Rxy - EIGHT * PIO * l_Sxy) - fxy_Lap_phy;
    double fxz_P = l_alpn1 * (l_Rxz - EIGHT * PIO * l_Sxz) - fxz_Lap_phy;
    double fyz_P = l_alpn1 * (l_Ryz - EIGHT * PIO * l_Syz) - fyz_Lap_phy;

    double f_trace = F1o3 * (gupxx * fxx_P + gupyy * fyy_P + gupzz * fzz_P
                   + TWO * (gupxy * fxy_P + gupxz * fxz_P + gupyz * fyz_P));

    // A_il A^l_j
    double A2_xx = gupxx * l_Axx * l_Axx + gupyy * l_Axy * l_Axy + gupzz * l_Axz * l_Axz
                 + TWO * (gupxy * l_Axx * l_Axy + gupxz * l_Axx * l_Axz + gupyz * l_Axy * l_Axz);
    double A2_yy = gupxx * l_Axy * l_Axy + gupyy * l_Ayy * l_Ayy + gupzz * l_Ayz * l_Ayz
                 + TWO * (gupxy * l_Axy * l_Ayy + gupxz * l_Axy * l_Ayz + gupyz * l_Ayy * l_Ayz);
    double A2_zz = gupxx * l_Axz * l_Axz + gupyy * l_Ayz * l_Ayz + gupzz * l_Azz * l_Azz
                 + TWO * (gupxy * l_Axz * l_Ayz + gupxz * l_Axz * l_Azz + gupyz * l_Ayz * l_Azz);
    double A2_xy = gupxx * l_Axx * l_Axy + gupyy * l_Axy * l_Ayy + gupzz * l_Axz * l_Ayz
                 + gupxy * (l_Axx * l_Ayy + l_Axy * l_Axy)
                 + gupxz * (l_Axx * l_Ayz + l_Axz * l_Axy)
                 + gupyz * (l_Axy * l_Ayz + l_Axz * l_Ayy);
    double A2_xz = gupxx * l_Axx * l_Axz + gupyy * l_Axy * l_Ayz + gupzz * l_Axz * l_Azz
                 + gupxy * (l_Axx * l_Ayz + l_Axy * l_Axz)
                 + gupxz * (l_Axx * l_Azz + l_Axz * l_Axz)
                 + gupyz * (l_Axy * l_Azz + l_Axz * l_Ayz);
    double A2_yz = gupxx * l_Axy * l_Axz + gupyy * l_Ayy * l_Ayz + gupzz * l_Ayz * l_Azz
                 + gupxy * (l_Axy * l_Ayz + l_Ayy * l_Axz)
                 + gupxz * (l_Axy * l_Azz + l_Ayz * l_Axz)
                 + gupyz * (l_Ayy * l_Azz + l_Ayz * l_Ayz);

    Axx_rhs[idx] = l_chin1 * (fxx_P - l_gxx * f_trace) + l_alpn1 * (l_trK * l_Axx - TWO * A2_xx)
                 + TWO * (l_Axx * betaxx + l_Axy * betayx + l_Axz * betazx)
                 - F2o3 * l_Axx * div_beta;
    Ayy_rhs[idx] = l_chin1 * (fyy_P - l_gyy * f_trace) + l_alpn1 * (l_trK * l_Ayy - TWO * A2_yy)
                 + TWO * (l_Axy * betaxy + l_Ayy * betayy + l_Ayz * betazy)
                 - F2o3 * l_Ayy * div_beta;
    Azz_rhs[idx] = l_chin1 * (fzz_P - l_gzz * f_trace) + l_alpn1 * (l_trK * l_Azz - TWO * A2_zz)
                 + TWO * (l_Axz * betaxz + l_Ayz * betayz + l_Azz * betazz)
                 - F2o3 * l_Azz * div_beta;
    Axy_rhs[idx] = l_chin1 * (fxy_P - l_gxy * f_trace) + l_alpn1 * (l_trK * l_Axy - TWO * A2_xy)
                 + l_Axx * betaxy + l_Axz * betazy
                 + l_Ayy * betayx + l_Ayz * betazx
                 + F1o3 * l_Axy * div_beta - l_Axy * betazz;
    Ayz_rhs[idx] = l_chin1 * (fyz_P - l_gyz * f_trace) + l_alpn1 * (l_trK * l_Ayz - TWO * A2_yz)
                 + l_Axy * betaxz + l_Ayy * betayz
                 + l_Axz * betaxy + l_Azz * betazy
                 + F1o3 * l_Ayz * div_beta - l_Ayz * betaxx;
    Axz_rhs[idx] = l_chin1 * (fxz_P - l_gxz * f_trace) + l_alpn1 * (l_trK * l_Axz - TWO * A2_xz)
                 + l_Axx * betaxz + l_Axy * betayz
                 + l_Ayz * betayx + l_Azz * betazx
                 + F1o3 * l_Axz * div_beta - l_Axz * betayy;

    // ===================================================================
    // trK_rhs
    // ===================================================================
    double trA2 = l_Axx * Aup_xx + l_Ayy * Aup_yy + l_Azz * Aup_zz
                + TWO * (l_Axy * Aup_xy + l_Axz * Aup_xz + l_Ayz * Aup_yz);
    double l_S = l_chin1 * (gupxx * l_Sxx + gupyy * l_Syy + gupzz * l_Szz
               + TWO * (gupxy * l_Sxy + gupxz * l_Sxz + gupyz * l_Syz));

    trK_rhs[idx] = -l_chin1 * trK_rhs_conf
                 + l_alpn1 * (F1o3 * l_trK * l_trK + trA2 + FOUR * PIO * (l_rho + l_S));

    // Store Ricci output
    Ricxx[idx] = l_Rxx; Ricxy[idx] = l_Rxy; Ricxz[idx] = l_Rxz;
    Ricyy[idx] = l_Ryy; Ricyz[idx] = l_Ryz; Riczz[idx] = l_Rzz;


    // ===================================================================
    // Lap_rhs
    // ===================================================================
    Lap_rhs[idx] = -TWO * l_alpn1 * l_trK;

    // ===================================================================
    // Gauge section
    // ===================================================================
    if constexpr (GAUGE_PARAM == 0)
    {
        betax_rhs[idx] = FF * l_dtSfx;
        betay_rhs[idx] = FF * l_dtSfy;
        betaz_rhs[idx] = FF * l_dtSfz;
        dtSfx_rhs[idx] = val_Gamx_rhs - ETA * l_dtSfx;
        dtSfy_rhs[idx] = val_Gamy_rhs - ETA * l_dtSfy;
        dtSfz_rhs[idx] = val_Gamz_rhs - ETA * l_dtSfz;
    }
    else if constexpr (GAUGE_PARAM == 1)
    {
        betax_rhs[idx] = l_Gamx - ETA * l_betax;
        betay_rhs[idx] = l_Gamy - ETA * l_betay;
        betaz_rhs[idx] = l_Gamz - ETA * l_betaz;
        dtSfx_rhs[idx] = 0.0;
        dtSfy_rhs[idx] = 0.0;
        dtSfz_rhs[idx] = 0.0;
    }
    else if constexpr (GAUGE_PARAM == 2)
    {
        betax_rhs[idx] = FF * l_dtSfx;
        betay_rhs[idx] = FF * l_dtSfy;
        betaz_rhs[idx] = FF * l_dtSfz;
        double reta_val = gupxx * chix * chix + gupyy * chiy * chiy + gupzz * chiz * chiz
                        + TWO * (gupxy * chix * chiy + gupxz * chix * chiz + gupyz * chiy * chiz);
        double sq = sqrt(l_chin1);
        double reta = 1.31 / 2.0 * sqrt(reta_val / l_chin1) / ((1.0 - sq) * (1.0 - sq));
        dtSfx_rhs[idx] = val_Gamx_rhs - reta * l_dtSfx;
        dtSfy_rhs[idx] = val_Gamy_rhs - reta * l_dtSfy;
        dtSfz_rhs[idx] = val_Gamz_rhs - reta * l_dtSfz;
    }
    else if constexpr (GAUGE_PARAM == 3)
    {
        betax_rhs[idx] = FF * l_dtSfx;
        betay_rhs[idx] = FF * l_dtSfy;
        betaz_rhs[idx] = FF * l_dtSfz;
        double reta_val = gupxx * chix * chix + gupyy * chiy * chiy + gupzz * chiz * chiz
                        + TWO * (gupxy * chix * chiy + gupxz * chix * chiz + gupyz * chiy * chiz);
        double reta = 1.31 / 2.0 * sqrt(reta_val / l_chin1) / ((1.0 - l_chin1) * (1.0 - l_chin1));
        dtSfx_rhs[idx] = val_Gamx_rhs - reta * l_dtSfx;
        dtSfy_rhs[idx] = val_Gamy_rhs - reta * l_dtSfy;
        dtSfz_rhs[idx] = val_Gamz_rhs - reta * l_dtSfz;
    }
    else if constexpr (GAUGE_PARAM == 4)
    {
        double reta_val = gupxx * chix * chix + gupyy * chiy * chiy + gupzz * chiz * chiz
                        + TWO * (gupxy * chix * chiy + gupxz * chix * chiz + gupyz * chiy * chiz);
        double sq = sqrt(l_chin1);
        double reta = 1.31 / 2.0 * sqrt(reta_val / l_chin1) / ((1.0 - sq) * (1.0 - sq));
        betax_rhs[idx] = FF * l_Gamx - reta * l_betax;
        betay_rhs[idx] = FF * l_Gamy - reta * l_betay;
        betaz_rhs[idx] = FF * l_Gamz - reta * l_betaz;
        dtSfx_rhs[idx] = 0.0;
        dtSfy_rhs[idx] = 0.0;
        dtSfz_rhs[idx] = 0.0;
    }
    else if constexpr (GAUGE_PARAM == 5)
    {
        double reta_val = gupxx * chix * chix + gupyy * chiy * chiy + gupzz * chiz * chiz
                        + TWO * (gupxy * chix * chiy + gupxz * chix * chiz + gupyz * chiy * chiz);
        double reta = 1.31 / 2.0 * sqrt(reta_val / l_chin1) / ((1.0 - l_chin1) * (1.0 - l_chin1));
        betax_rhs[idx] = FF * l_Gamx - reta * l_betax;
        betay_rhs[idx] = FF * l_Gamy - reta * l_betay;
        betaz_rhs[idx] = FF * l_Gamz - reta * l_betaz;
        dtSfx_rhs[idx] = 0.0;
        dtSfy_rhs[idx] = 0.0;
        dtSfz_rhs[idx] = 0.0;
    }
    else if constexpr (GAUGE_PARAM == 6)
    {
        betax_rhs[idx] = FF * l_dtSfx;
        betay_rhs[idx] = FF * l_dtSfy;
        betaz_rhs[idx] = FF * l_dtSfz;
        double reta = ETA;
        if (BHN == 2) {
            double M = Mass1 + Mass2;
            double A = 2.0 / M;
            double w1 = 12.0, w2 = 12.0;
            double C1 = 1.0 / Mass1 - A;
            double C2 = 1.0 / Mass2 - A;
            double sep2 = (Porg1 - Porg4) * (Porg1 - Porg4)
                        + (Porg2 - Porg5) * (Porg2 - Porg5)
                        + (Porg3 - Porg6) * (Porg3 - Porg6);
            double r1 = ((Porg1 - X[idx]) * (Porg1 - X[idx])
                       + (Porg2 - Y[idx]) * (Porg2 - Y[idx])
                       + (Porg3 - Z[idx]) * (Porg3 - Z[idx])) / sep2;
            double r2 = ((Porg4 - X[idx]) * (Porg4 - X[idx])
                       + (Porg5 - Y[idx]) * (Porg5 - Y[idx])
                       + (Porg6 - Z[idx]) * (Porg6 - Z[idx])) / sep2;
            reta = A + C1 / (1.0 + w1 * r1) + C2 / (1.0 + w2 * r2);
        }
        dtSfx_rhs[idx] = val_Gamx_rhs - reta * l_dtSfx;
        dtSfy_rhs[idx] = val_Gamy_rhs - reta * l_dtSfy;
        dtSfz_rhs[idx] = val_Gamz_rhs - reta * l_dtSfz;
    }
    else if constexpr (GAUGE_PARAM == 7)
    {
        betax_rhs[idx] = FF * l_dtSfx;
        betay_rhs[idx] = FF * l_dtSfy;
        betaz_rhs[idx] = FF * l_dtSfz;
        double reta = ETA;
        if (BHN == 2) {
            double M = Mass1 + Mass2;
            double A = 2.0 / M;
            double w1 = 12.0, w2 = 12.0;
            double C1 = 1.0 / Mass1 - A;
            double C2 = 1.0 / Mass2 - A;
            double sep2 = (Porg1 - Porg4) * (Porg1 - Porg4)
                        + (Porg2 - Porg5) * (Porg2 - Porg5)
                        + (Porg3 - Porg6) * (Porg3 - Porg6);
            double r1 = ((Porg1 - X[idx]) * (Porg1 - X[idx])
                       + (Porg2 - Y[idx]) * (Porg2 - Y[idx])
                       + (Porg3 - Z[idx]) * (Porg3 - Z[idx])) / sep2;
            double r2 = ((Porg4 - X[idx]) * (Porg4 - X[idx])
                       + (Porg5 - Y[idx]) * (Porg5 - Y[idx])
                       + (Porg6 - Z[idx]) * (Porg6 - Z[idx])) / sep2;
            reta = A + C1 * exp(-w1 * r1) + C2 * exp(-w2 * r2);
        }
        dtSfx_rhs[idx] = val_Gamx_rhs - reta * l_dtSfx;
        dtSfy_rhs[idx] = val_Gamy_rhs - reta * l_dtSfy;
        dtSfz_rhs[idx] = val_Gamz_rhs - reta * l_dtSfz;
    }

    // Store Gam^i RHS
    Gamx_rhs[idx] = val_Gamx_rhs;
    Gamy_rhs[idx] = val_Gamy_rhs;
    Gamz_rhs[idx] = val_Gamz_rhs;


    // ===================================================================
    // Advection terms
    // ===================================================================
    double fadx, fady, fadz;

    // g_ij advection
    d_point_fderivs_shc_gw<GW>(ex, pex, dxx, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    gxx_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, gxy_f, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, ANTI, ANTI, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    gxy_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, gxz_f, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, ANTI, SYM, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    gxz_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, dyy, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    gyy_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, gyz_f, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, ANTI, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    gyz_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, dzz, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    gzz_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;

    // A_ij advection
    d_point_fderivs_shc_gw<GW>(ex, pex, Axx, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    Axx_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Axy, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, ANTI, ANTI, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    Axy_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Axz, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, ANTI, SYM, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    Axz_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Ayy, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    Ayy_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Ayz, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, ANTI, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    Ayz_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Azz, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    Azz_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;

    // chi and trK advection
    d_point_fderivs_shc_gw<GW>(ex, pex, chi, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    chi_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, trK, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    trK_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;

    // Gam^i advection
    d_point_fderivs_shc_gw<GW>(ex, pex, Gamx, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, ANTI, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    Gamx_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Gamy, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, ANTI, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    Gamy_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Gamz, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, SYM, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    Gamz_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;

    // Lap advection
    d_point_fderivs_shc_gw<GW>(ex, pex, Lap, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    Lap_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;

    // beta advection (all gauge values)
    d_point_fderivs_shc_gw<GW>(ex, pex, betax, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, ANTI, SYM, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    betax_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, betay, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, ANTI, SYM, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    betay_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    d_point_fderivs_shc_gw<GW>(ex, pex, betaz, &fadx, &fady, &fadz,
        crho, sigma, Rcrd, SYM, SYM, ANTI, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    betaz_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;

    // dtSf advection (GAUGE 0, 2, 3, 6, 7)
    if constexpr (GAUGE_PARAM == 0 || GAUGE_PARAM == 2 || GAUGE_PARAM == 3 || GAUGE_PARAM == 6 || GAUGE_PARAM == 7)
    {
        d_point_fderivs_shc_gw<GW>(ex, pex, dtSfx, &fadx, &fady, &fadz,
            crho, sigma, Rcrd, ANTI, SYM, SYM, symmetry, onoff, sst,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
        dtSfx_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
        d_point_fderivs_shc_gw<GW>(ex, pex, dtSfy, &fadx, &fady, &fadz,
            crho, sigma, Rcrd, SYM, ANTI, SYM, symmetry, onoff, sst,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
        dtSfy_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
        d_point_fderivs_shc_gw<GW>(ex, pex, dtSfz, &fadx, &fady, &fadz,
            crho, sigma, Rcrd, SYM, SYM, ANTI, symmetry, onoff, sst,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
        dtSfz_rhs[idx] += l_betax * fadx + l_betay * fady + l_betaz * fadz;
    }


    // ===================================================================
    // KO dissipation
    // ===================================================================
    if (eps > 0.0)
    {
        // SSS = (SYM, SYM, SYM): chi, trK, gxx, gyy, gzz, Axx, Ayy, Azz, Lap
        chi_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, chi, crho, sigma, Rcrd,
            SYM, SYM, SYM, symmetry, eps, sst, i, j, k);
        trK_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, trK, crho, sigma, Rcrd,
            SYM, SYM, SYM, symmetry, eps, sst, i, j, k);
        gxx_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, dxx, crho, sigma, Rcrd,
            SYM, SYM, SYM, symmetry, eps, sst, i, j, k);
        gyy_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, dyy, crho, sigma, Rcrd,
            SYM, SYM, SYM, symmetry, eps, sst, i, j, k);
        gzz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, dzz, crho, sigma, Rcrd,
            SYM, SYM, SYM, symmetry, eps, sst, i, j, k);
        Axx_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Axx, crho, sigma, Rcrd,
            SYM, SYM, SYM, symmetry, eps, sst, i, j, k);
        Ayy_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Ayy, crho, sigma, Rcrd,
            SYM, SYM, SYM, symmetry, eps, sst, i, j, k);
        Azz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Azz, crho, sigma, Rcrd,
            SYM, SYM, SYM, symmetry, eps, sst, i, j, k);
        Lap_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Lap, crho, sigma, Rcrd,
            SYM, SYM, SYM, symmetry, eps, sst, i, j, k);

        // AAS = (ANTI, ANTI, SYM): gxy, Axy
        gxy_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, gxy_f, crho, sigma, Rcrd,
            ANTI, ANTI, SYM, symmetry, eps, sst, i, j, k);
        Axy_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Axy, crho, sigma, Rcrd,
            ANTI, ANTI, SYM, symmetry, eps, sst, i, j, k);

        // ASA = (ANTI, SYM, ANTI): gxz, Axz
        gxz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, gxz_f, crho, sigma, Rcrd,
            ANTI, SYM, ANTI, symmetry, eps, sst, i, j, k);
        Axz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Axz, crho, sigma, Rcrd,
            ANTI, SYM, ANTI, symmetry, eps, sst, i, j, k);

        // SAA = (SYM, ANTI, ANTI): gyz, Ayz
        gyz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, gyz_f, crho, sigma, Rcrd,
            SYM, ANTI, ANTI, symmetry, eps, sst, i, j, k);
        Ayz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Ayz, crho, sigma, Rcrd,
            SYM, ANTI, ANTI, symmetry, eps, sst, i, j, k);

        // ASS = (ANTI, SYM, SYM): Gamx, betax, dtSfx
        Gamx_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Gamx, crho, sigma, Rcrd,
            ANTI, SYM, SYM, symmetry, eps, sst, i, j, k);
        betax_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, betax, crho, sigma, Rcrd,
            ANTI, SYM, SYM, symmetry, eps, sst, i, j, k);

        // SAS = (SYM, ANTI, SYM): Gamy, betay, dtSfy
        Gamy_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Gamy, crho, sigma, Rcrd,
            SYM, ANTI, SYM, symmetry, eps, sst, i, j, k);
        betay_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, betay, crho, sigma, Rcrd,
            SYM, ANTI, SYM, symmetry, eps, sst, i, j, k);

        // SSA = (SYM, SYM, ANTI): Gamz, betaz, dtSfz
        Gamz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, Gamz, crho, sigma, Rcrd,
            SYM, SYM, ANTI, symmetry, eps, sst, i, j, k);
        betaz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, betaz, crho, sigma, Rcrd,
            SYM, SYM, ANTI, symmetry, eps, sst, i, j, k);

        // dtSf KO dissipation (GAUGE 0, 2, 3, 6, 7)
        if constexpr (GAUGE_PARAM == 0 || GAUGE_PARAM == 2 || GAUGE_PARAM == 3 || GAUGE_PARAM == 6 || GAUGE_PARAM == 7)
        {
            dtSfx_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, dtSfx, crho, sigma, Rcrd,
                ANTI, SYM, SYM, symmetry, eps, sst, i, j, k);
            dtSfy_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, dtSfy, crho, sigma, Rcrd,
                SYM, ANTI, SYM, symmetry, eps, sst, i, j, k);
            dtSfz_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, dtSfz, crho, sigma, Rcrd,
                SYM, SYM, ANTI, symmetry, eps, sst, i, j, k);
        }
    }


    // ===================================================================
    // Constraints (only when co == 0)
    // ===================================================================
    if (co == 0)
    {
        // Hamiltonian constraint
        // ham_Res = chi * R + 2/3 * K^2 - A_ij * A^ij - 16 * pi * rho
        // where R is the physical Ricci scalar
        double conf_R_scalar = gupxx * l_Rxx + gupyy * l_Ryy + gupzz * l_Rzz
                             + TWO * (gupxy * l_Rxy + gupxz * l_Rxz + gupyz * l_Ryz);
        ham_Res[idx] = l_chin1 * conf_R_scalar + F2o3 * l_trK * l_trK - trA2 - F16 * PIO * l_rho;

        // Momentum constraints: mov_i = gup^kj * (D_k A_ij - ch^{-1} d_k chi * A_ij) - 2/3 d_i K - 8 pi S_i
        // First, compute D_k A_ij - ch^{-1} d_k chi * A_ij for all components
        double dA_xxx, dA_xxy, dA_xxz, dA_xyx, dA_xyy, dA_xyz, dA_xzx, dA_xzy, dA_xzz;
        double dA_yyx, dA_yyy, dA_yyz, dA_yzx, dA_yzy, dA_yzz, dA_zzx, dA_zzy, dA_zzz;

        // Axx first derivatives
        d_point_fderivs_shc_gw<GW>(ex, pex, Axx, &dA_xxx, &dA_xxy, &dA_xxz,
            crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
        // Axy first derivatives
        d_point_fderivs_shc_gw<GW>(ex, pex, Axy, &dA_xyx, &dA_xyy, &dA_xyz,
            crho, sigma, Rcrd, ANTI, ANTI, SYM, symmetry, onoff, sst,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
        // Axz first derivatives
        d_point_fderivs_shc_gw<GW>(ex, pex, Axz, &dA_xzx, &dA_xzy, &dA_xzz,
            crho, sigma, Rcrd, ANTI, SYM, ANTI, symmetry, onoff, sst,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
        // Ayy first derivatives
        d_point_fderivs_shc_gw<GW>(ex, pex, Ayy, &dA_yyx, &dA_yyy, &dA_yyz,
            crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
        // Ayz first derivatives
        d_point_fderivs_shc_gw<GW>(ex, pex, Ayz, &dA_yzx, &dA_yzy, &dA_yzz,
            crho, sigma, Rcrd, SYM, ANTI, ANTI, symmetry, onoff, sst,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
        // Azz first derivatives
        d_point_fderivs_shc_gw<GW>(ex, pex, Azz, &dA_zzx, &dA_zzy, &dA_zzz,
            crho, sigma, Rcrd, SYM, SYM, SYM, symmetry, onoff, sst,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

        // Covariant derivative: D_k A_ij - ch^{-1} d_k chi * A_ij
        double gxxx_c = dA_xxx - (l_Gamxxx * l_Axx + l_Gamyxx * l_Axy + l_Gamzxx * l_Axz
                       + l_Gamxxx * l_Axx + l_Gamyxx * l_Axy + l_Gamzxx * l_Axz) - chix * l_Axx * inv_chin1;
        double gxyx_c = dA_xyx - (l_Gamxxy * l_Axx + l_Gamyxy * l_Axy + l_Gamzxy * l_Axz
                       + l_Gamxxx * l_Axy + l_Gamyxx * l_Ayy + l_Gamzxx * l_Ayz) - chix * l_Axy * inv_chin1;
        double gxzx_c = dA_xzx - (l_Gamxxz * l_Axx + l_Gamyxz * l_Axy + l_Gamzxz * l_Axz
                       + l_Gamxxx * l_Axz + l_Gamyxx * l_Ayz + l_Gamzxx * l_Azz) - chix * l_Axz * inv_chin1;
        double gyyx_c = dA_yyx - (l_Gamxxy * l_Axy + l_Gamyxy * l_Ayy + l_Gamzxy * l_Ayz
                       + l_Gamxxy * l_Axy + l_Gamyxy * l_Ayy + l_Gamzxy * l_Ayz) - chix * l_Ayy * inv_chin1;
        double gyzx_c = dA_yzx - (l_Gamxxz * l_Axy + l_Gamyxz * l_Ayy + l_Gamzxz * l_Ayz
                       + l_Gamxxy * l_Axz + l_Gamyxy * l_Ayz + l_Gamzxy * l_Azz) - chix * l_Ayz * inv_chin1;
        double gzzx_c = dA_zzx - (l_Gamxxz * l_Axz + l_Gamyxz * l_Ayz + l_Gamzxz * l_Azz
                       + l_Gamxxz * l_Axz + l_Gamyxz * l_Ayz + l_Gamzxz * l_Azz) - chix * l_Azz * inv_chin1;

        double gxxy_c = dA_xxy - (l_Gamxxy * l_Axx + l_Gamyxy * l_Axy + l_Gamzxy * l_Axz
                       + l_Gamxxy * l_Axx + l_Gamyxy * l_Axy + l_Gamzxy * l_Axz) - chiy * l_Axx * inv_chin1;
        double gxyy_c = dA_xyy - (l_Gamxyy * l_Axx + l_Gamyyy * l_Axy + l_Gamzyy * l_Axz
                       + l_Gamxxy * l_Axy + l_Gamyxy * l_Ayy + l_Gamzxy * l_Ayz) - chiy * l_Axy * inv_chin1;
        double gxzy_c = dA_xzy - (l_Gamxyz * l_Axx + l_Gamyyz * l_Axy + l_Gamzyz * l_Axz
                       + l_Gamxxy * l_Axz + l_Gamyxy * l_Ayz + l_Gamzxy * l_Azz) - chiy * l_Axz * inv_chin1;
        double gyyy_c = dA_yyy - (l_Gamxyy * l_Axy + l_Gamyyy * l_Ayy + l_Gamzyy * l_Ayz
                       + l_Gamxyy * l_Axy + l_Gamyyy * l_Ayy + l_Gamzyy * l_Ayz) - chiy * l_Ayy * inv_chin1;
        double gyzy_c = dA_yzy - (l_Gamxyz * l_Axy + l_Gamyyz * l_Ayy + l_Gamzyz * l_Ayz
                       + l_Gamxyy * l_Axz + l_Gamyyy * l_Ayz + l_Gamzyy * l_Azz) - chiy * l_Ayz * inv_chin1;
        double gzzy_c = dA_zzy - (l_Gamxyz * l_Axz + l_Gamyyz * l_Ayz + l_Gamzyz * l_Azz
                       + l_Gamxyz * l_Axz + l_Gamyyz * l_Ayz + l_Gamzyz * l_Azz) - chiy * l_Azz * inv_chin1;

        double gxxz_c = dA_xxz - (l_Gamxxz * l_Axx + l_Gamyxz * l_Axy + l_Gamzxz * l_Axz
                       + l_Gamxxz * l_Axx + l_Gamyxz * l_Axy + l_Gamzxz * l_Axz) - chiz * l_Axx * inv_chin1;
        double gxyz_c = dA_xyz - (l_Gamxyz * l_Axx + l_Gamyyz * l_Axy + l_Gamzyz * l_Axz
                       + l_Gamxxz * l_Axy + l_Gamyxz * l_Ayy + l_Gamzxz * l_Ayz) - chiz * l_Axy * inv_chin1;
        double gxzz_c = dA_xzz - (l_Gamxzz * l_Axx + l_Gamyzz * l_Axy + l_Gamzzz * l_Axz
                       + l_Gamxxz * l_Axz + l_Gamyxz * l_Ayz + l_Gamzxz * l_Azz) - chiz * l_Axz * inv_chin1;
        double gyyz_c = dA_yyz - (l_Gamxyz * l_Axy + l_Gamyyz * l_Ayy + l_Gamzyz * l_Ayz
                       + l_Gamxyz * l_Axy + l_Gamyyz * l_Ayy + l_Gamzyz * l_Ayz) - chiz * l_Ayy * inv_chin1;
        double gyzz_c = dA_yzz - (l_Gamxzz * l_Axy + l_Gamyzz * l_Ayy + l_Gamzzz * l_Ayz
                       + l_Gamxyz * l_Axz + l_Gamyyz * l_Ayz + l_Gamzyz * l_Azz) - chiz * l_Ayz * inv_chin1;
        double gzzz_c = dA_zzz - (l_Gamxzz * l_Axz + l_Gamyzz * l_Ayz + l_Gamzzz * l_Azz
                       + l_Gamxzz * l_Axz + l_Gamyzz * l_Ayz + l_Gamzzz * l_Azz) - chiz * l_Azz * inv_chin1;

        double movx = gupxx * gxxx_c + gupyy * gxyy_c + gupzz * gxzz_c
                    + gupxy * gxyx_c + gupxz * gxzx_c + gupyz * gxzy_c
                    + gupxy * gxxy_c + gupxz * gxxz_c + gupyz * gxyz_c;
        double movy = gupxx * gxyx_c + gupyy * gyyy_c + gupzz * gyzz_c
                    + gupxy * gyyx_c + gupxz * gyzx_c + gupyz * gyzy_c
                    + gupxy * gxyy_c + gupxz * gxyz_c + gupyz * gyyz_c;
        double movz = gupxx * gxzx_c + gupyy * gyzy_c + gupzz * gzzz_c
                    + gupxy * gyzx_c + gupxz * gzzx_c + gupyz * gzzy_c
                    + gupxy * gxzy_c + gupxz * gxzz_c + gupyz * gyzz_c;

        movx_Res[idx] = movx - F2o3 * Kx - F8 * PIO * l_Sx;
        movy_Res[idx] = movy - F2o3 * Ky - F8 * PIO * l_Sy;
        movz_Res[idx] = movz - F2o3 * Kz - F8 * PIO * l_Sz;
    }
}


void gpu_compute_rhs_bssn_ss_launch(
    cudaStream_t &stream, int* ex, int* pex, double T,
    double* crho, double* sigma, double* R,
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
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co, int sst
)
{
    const int DIMX = 8, DIMY = 8, DIMZ = 4;
    dim3 block(DIMX, DIMY, DIMZ);
    dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );

    int ex0 = ex[0], ex1 = ex[1], ex2 = ex[2];
    int pex0 = pex[0], pex1 = pex[1], pex2 = pex[2];

    double Porg1 = PBH[0], Porg2 = PBH[1], Porg3 = PBH[2];
    double Porg4 = PBH[3], Porg5 = PBH[4], Porg6 = PBH[5];
    double Mass1 = Mass[0], Mass2 = Mass[1];
    int BHN_val = BHN;

#define LAUNCH_BSSN_RHS_SS(G) \
    bssn_rhs_ss_kernel<ghost_width, G><<<grid, block, 0, stream>>>( \
        ex0, ex1, ex2, pex0, pex1, pex2, \
        symmetry, lev, eps, co, sst, \
        BHN_val, Porg1, Porg2, Porg3, Porg4, Porg5, Porg6, Mass1, Mass2, \
        crho, sigma, R, X, Y, Z, \
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, \
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, \
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, \
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, \
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz, \
        Axx, Axy, Axz, Ayy, Ayz, Azz, \
        Gamx, Gamy, Gamz, Lap, betax, betay, betaz, \
        dtSfx, dtSfy, dtSfz, rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, \
        chi_rhs, trK_rhs, gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs, \
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs, \
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs, betax_rhs, betay_rhs, betaz_rhs, \
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs, \
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, \
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, \
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, \
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, \
        ham_Res, movx_Res, movy_Res, movz_Res, \
        Gmx_Res, Gmy_Res, Gmz_Res)

#if GAUGE == 0
    LAUNCH_BSSN_RHS_SS(0);
#elif GAUGE == 1
    LAUNCH_BSSN_RHS_SS(1);
#elif GAUGE == 2
    LAUNCH_BSSN_RHS_SS(2);
#elif GAUGE == 3
    LAUNCH_BSSN_RHS_SS(3);
#elif GAUGE == 4
    LAUNCH_BSSN_RHS_SS(4);
#elif GAUGE == 5
    LAUNCH_BSSN_RHS_SS(5);
#elif GAUGE == 6
    LAUNCH_BSSN_RHS_SS(6);
#elif GAUGE == 7
    LAUNCH_BSSN_RHS_SS(7);
#endif

#undef LAUNCH_BSSN_RHS_SS
}
