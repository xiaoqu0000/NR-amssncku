#include "macrodef.h"
#include "bssn_rhs.h"
#include "scalar_rhs.h"
#include "bssnEScalar_rhs_gpu.cuh"

#include "rhs_smem_common.cuh"
#include "diff_new_gpu.cuh"
#include "diff_new_sh_gpu.cuh"
#include "kodiss_gpu.cuh"
#include "kodiss_sh_gpu.cuh"
#include "lopsidediff_gpu.cuh"
#include "point_diff_new_sh_gpu.cuh"
#include "ricci_gamma.h"
#include "bssn_constraint.h"

#include <cuda_runtime.h>

extern "C" {
#ifdef fortran1
void set_escalar_parameter(double&, double&, double&, double&, double&);
#endif
#ifdef fortran2
void SET_ESCALAR_PARAMETER(double&, double&, double&, double&, double&);
#endif
#ifdef fortran3
void set_escalar_parameter_(double&, double&, double&, double&, double&);
#endif
}

static inline void load_escalar_parameters(double& a2, double& phi0, double& r0, double& sigma0, double& l2)
{
#ifdef fortran1
    set_escalar_parameter(a2, phi0, r0, sigma0, l2);
#endif
#ifdef fortran2
    SET_ESCALAR_PARAMETER(a2, phi0, r0, sigma0, l2);
#endif
#ifdef fortran3
    set_escalar_parameter_(a2, phi0, r0, sigma0, l2);
#endif
}

template<int GW>
__global__ void scalar_rhs_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ Sphi, const double* __restrict__ Spi,
    double* __restrict__ Sphi_rhs, double* __restrict__ Spi_rhs,
    int symmetry, int lev, double eps)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    double fxx, fxy, fxz, fyy, fyz, fzz;
    d_fdderivs_point_gw<GW>(ex, pex, Sphi, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz,
                            X, Y, Z, 1.0, 1.0, 1.0, symmetry, 0, i, j, k);

    double rhs_phi = Spi[idx];
    double rhs_pi = fxx + fyy + fzz;
    if (eps > 0.0) {
        rhs_phi += d_kodis_point<GW>(ex, pex, Sphi, X, Y, Z, 1.0, 1.0, 1.0, symmetry, eps, i, j, k);
        rhs_pi += d_kodis_point<GW>(ex, pex, Spi, X, Y, Z, 1.0, 1.0, 1.0, symmetry, eps, i, j, k);
    }
    Sphi_rhs[idx] = rhs_phi;
    Spi_rhs[idx] = rhs_pi;
}

template<int GW>
__global__ void scalar_rhs_ss_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ crho, const double* __restrict__ sigma, const double* __restrict__ R,
    const double* __restrict__ drhodx, const double* __restrict__ drhody, const double* __restrict__ drhodz,
    const double* __restrict__ dsigmadx, const double* __restrict__ dsigmady, const double* __restrict__ dsigmadz,
    const double* __restrict__ dRdx, const double* __restrict__ dRdy, const double* __restrict__ dRdz,
    const double* __restrict__ drhodxx, const double* __restrict__ drhodxy, const double* __restrict__ drhodxz,
    const double* __restrict__ drhodyy, const double* __restrict__ drhodyz, const double* __restrict__ drhodzz,
    const double* __restrict__ dsigmadxx, const double* __restrict__ dsigmadxy, const double* __restrict__ dsigmadxz,
    const double* __restrict__ dsigmadyy, const double* __restrict__ dsigmadyz, const double* __restrict__ dsigmadzz,
    const double* __restrict__ dRdxx, const double* __restrict__ dRdxy, const double* __restrict__ dRdxz,
    const double* __restrict__ dRdyy, const double* __restrict__ dRdyz, const double* __restrict__ dRdzz,
    const double* __restrict__ Sphi, const double* __restrict__ Spi,
    double* __restrict__ Sphi_rhs, double* __restrict__ Spi_rhs,
    int symmetry, int lev, double eps, int sst)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    double fxx, fxy, fxz, fyy, fyz, fzz;
    d_point_fdderivs_shc_gw<GW>(ex, pex, Sphi, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz,
                                crho, sigma, R, 1.0, 1.0, 1.0, symmetry, 0, sst,
                                drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
                                drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
                                dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
                                dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k);

    double rhs_phi = Spi[idx];
    double rhs_pi = fxx + fyy + fzz;
    if (eps > 0.0) {
        rhs_phi += d_kodis_sh_point<GW>(ex, pex, Sphi, crho, sigma, R, 1.0, 1.0, 1.0, symmetry, eps, sst, i, j, k);
        rhs_pi += d_kodis_sh_point<GW>(ex, pex, Spi, crho, sigma, R, 1.0, 1.0, 1.0, symmetry, eps, sst, i, j, k);
    }
    Sphi_rhs[idx] = rhs_phi;
    Spi_rhs[idx] = rhs_pi;
}

template<int GW>
__global__ void escalar_matter_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    const double* __restrict__ Lap,
    const double* __restrict__ Sphi, const double* __restrict__ Spi,
    double* __restrict__ Sphi_rhs, double* __restrict__ Spi_rhs,
    double* __restrict__ rho, double* __restrict__ Sx, double* __restrict__ Sy, double* __restrict__ Sz,
    double* __restrict__ Sxx, double* __restrict__ Sxy, double* __restrict__ Sxz,
    double* __restrict__ Syy, double* __restrict__ Syz, double* __restrict__ Szz,
    double a2, int symmetry, int lev)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    double chix, chiy, chiz;
    double lapx, lapy, lapz;
    double kx, ky, kz;
    double fxx, fxy, fxz, fyy, fyz, fzz;

    d_fderivs_point_gw<GW>(ex, pex, chi, &chix, &chiy, &chiz, X, Y, Z, 1.0, 1.0, 1.0, symmetry, 0, i, j, k);
    d_fderivs_point_gw<GW>(ex, pex, Lap, &lapx, &lapy, &lapz, X, Y, Z, 1.0, 1.0, 1.0, symmetry, 0, i, j, k);
    d_fderivs_point_gw<GW>(ex, pex, Sphi, &kx, &ky, &kz, X, Y, Z, 1.0, 1.0, 1.0, symmetry, 0, i, j, k);
    d_fdderivs_point_gw<GW>(ex, pex, Sphi, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz,
                            X, Y, Z, 1.0, 1.0, 1.0, symmetry, 0, i, j, k);

    const double gxx = dxx[idx] + 1.0;
    const double gyy = dyy[idx] + 1.0;
    const double gzz = dzz[idx] + 1.0;
    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    d_inverse_metric(gxx, gxy[idx], gxz[idx], gyy, gyz[idx], gzz, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);

    const double alpn1 = Lap[idx] + 1.0;
    const double chin1 = chi[idx] + 1.0;
    const double inv_chin1 = 1.0 / chin1;

    const double grad_term = gupxx * fxx + gupyy * fyy + gupzz * fzz
                           + 2.0 * (gupxy * fxy + gupxz * fxz + gupyz * fyz);
    const double gamma_grad =
        (Gamx[idx] + 0.5 * inv_chin1 * (gupxx * chix + gupxy * chiy + gupxz * chiz)) * kx +
        (Gamy[idx] + 0.5 * inv_chin1 * (gupxy * chix + gupyy * chiy + gupyz * chiz)) * ky +
        (Gamz[idx] + 0.5 * inv_chin1 * (gupxz * chix + gupyz * chiy + gupzz * chiz)) * kz;
    const double lapse_grad =
        gupxx * lapx * kx + gupxy * lapx * ky + gupxz * lapx * kz +
        gupxy * lapy * kx + gupyy * lapy * ky + gupyz * lapy * kz +
        gupxz * lapz * kx + gupyz * lapz * ky + gupzz * lapz * kz;

    double V, dVdSphi;
    d_frpotential(Sphi[idx], a2, V, dVdSphi);

    Sphi_rhs[idx] = alpn1 * Spi[idx];
    Spi_rhs[idx] = (alpn1 * (grad_term - gamma_grad) + lapse_grad) * chin1
                 + alpn1 * (trK[idx] * Spi[idx] - dVdSphi);

    const double grad_energy = chin1 * (0.5 * (gupxx * kx * kx + gupyy * ky * ky + gupzz * kz * kz)
                             + gupxy * kx * ky + gupxz * kx * kz + gupyz * ky * kz);
    const double spi = Spi[idx];
    const double local_rho = grad_energy + 0.5 * spi * spi + V;
    const double press_part = (local_rho - spi * spi) / chin1;

    rho[idx] = local_rho;
    Sx[idx] = -spi * kx;
    Sy[idx] = -spi * ky;
    Sz[idx] = -spi * kz;
    Sxx[idx] = kx * kx - press_part * gxx;
    Sxy[idx] = kx * ky - press_part * gxy[idx];
    Sxz[idx] = kx * kz - press_part * gxz[idx];
    Syy[idx] = ky * ky - press_part * gyy;
    Syz[idx] = ky * kz - press_part * gyz[idx];
    Szz[idx] = kz * kz - press_part * gzz;
}

template<int GW>
__global__ void escalar_scalar_advect_ko_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ betax, const double* __restrict__ betay, const double* __restrict__ betaz,
    const double* __restrict__ Sphi, const double* __restrict__ Spi,
    double* __restrict__ Sphi_rhs, double* __restrict__ Spi_rhs,
    int symmetry, double eps)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    const int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    Sphi_rhs[idx] += d_lopsided_point<GW>(ex, pex, Sphi, betax, betay, betaz, X, Y, Z, symmetry, 1.0, 1.0, 1.0, i, j, k);
    Spi_rhs[idx] += d_lopsided_point<GW>(ex, pex, Spi, betax, betay, betaz, X, Y, Z, symmetry, 1.0, 1.0, 1.0, i, j, k);
    if (eps > 0.0) {
        Sphi_rhs[idx] += d_kodis_point<GW>(ex, pex, Sphi, X, Y, Z, 1.0, 1.0, 1.0, symmetry, eps, i, j, k);
        Spi_rhs[idx] += d_kodis_point<GW>(ex, pex, Spi, X, Y, Z, 1.0, 1.0, 1.0, symmetry, eps, i, j, k);
    }
}

template<int GW>
__global__ void escalar_matter_ss_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ crho, const double* __restrict__ sigma, const double* __restrict__ R,
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
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    const double* __restrict__ Lap,
    const double* __restrict__ Sphi, const double* __restrict__ Spi,
    double* __restrict__ Sphi_rhs, double* __restrict__ Spi_rhs,
    double* __restrict__ rho, double* __restrict__ Sx, double* __restrict__ Sy, double* __restrict__ Sz,
    double* __restrict__ Sxx, double* __restrict__ Sxy, double* __restrict__ Sxz,
    double* __restrict__ Syy, double* __restrict__ Syz, double* __restrict__ Szz,
    double a2, int symmetry, int lev, int sst)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    double chix, chiy, chiz, lapx, lapy, lapz, kx, ky, kz;
    double fxx, fxy, fxz, fyy, fyz, fzz;

#define SHC_FDER(F,S1,S2,S3,FX,FY,FZ) \
    d_point_fderivs_shc_gw<GW>(ex, pex, F, FX, FY, FZ, crho, sigma, R, S1, S2, S3, symmetry, 0, sst, \
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k)
#define SHC_FDDER(F,S1,S2,S3,FXX,FXY,FXZ,FYY,FYZ,FZZ) \
    d_point_fdderivs_shc_gw<GW>(ex, pex, F, FXX, FXY, FXZ, FYY, FYZ, FZZ, crho, sigma, R, S1, S2, S3, symmetry, 0, sst, \
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, \
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, \
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, \
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k)

    SHC_FDER(chi, 1.0, 1.0, 1.0, &chix, &chiy, &chiz);
    SHC_FDER(Lap, 1.0, 1.0, 1.0, &lapx, &lapy, &lapz);
    SHC_FDER(Sphi, 1.0, 1.0, 1.0, &kx, &ky, &kz);
    SHC_FDDER(Sphi, 1.0, 1.0, 1.0, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz);

    const double gxx = dxx[idx] + 1.0;
    const double gyy = dyy[idx] + 1.0;
    const double gzz = dzz[idx] + 1.0;
    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    d_inverse_metric(gxx, gxy[idx], gxz[idx], gyy, gyz[idx], gzz, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);

    const double alpn1 = Lap[idx] + 1.0;
    const double chin1 = chi[idx] + 1.0;
    const double inv_chin1 = 1.0 / chin1;

    const double grad_term = gupxx * fxx + gupyy * fyy + gupzz * fzz
                           + 2.0 * (gupxy * fxy + gupxz * fxz + gupyz * fyz);
    const double gamma_grad =
        (Gamx[idx] + 0.5 * inv_chin1 * (gupxx * chix + gupxy * chiy + gupxz * chiz)) * kx +
        (Gamy[idx] + 0.5 * inv_chin1 * (gupxy * chix + gupyy * chiy + gupyz * chiz)) * ky +
        (Gamz[idx] + 0.5 * inv_chin1 * (gupxz * chix + gupyz * chiy + gupzz * chiz)) * kz;
    const double lapse_grad =
        gupxx * lapx * kx + gupxy * lapx * ky + gupxz * lapx * kz +
        gupxy * lapy * kx + gupyy * lapy * ky + gupyz * lapy * kz +
        gupxz * lapz * kx + gupyz * lapz * ky + gupzz * lapz * kz;

    double V, dVdSphi;
    d_frpotential(Sphi[idx], a2, V, dVdSphi);

    Sphi_rhs[idx] = alpn1 * Spi[idx];
    Spi_rhs[idx] = (alpn1 * (grad_term - gamma_grad) + lapse_grad) * chin1
                 + alpn1 * (trK[idx] * Spi[idx] - dVdSphi);

    const double spi = Spi[idx];
    const double local_rho = chin1 * (0.5 * (gupxx * kx * kx + gupyy * ky * ky + gupzz * kz * kz)
                           + gupxy * kx * ky + gupxz * kx * kz + gupyz * ky * kz)
                           + 0.5 * spi * spi + V;
    const double press_part = (local_rho - spi * spi) / chin1;

    rho[idx] = local_rho;
    Sx[idx] = -spi * kx; Sy[idx] = -spi * ky; Sz[idx] = -spi * kz;
    Sxx[idx] = kx * kx - press_part * gxx;
    Sxy[idx] = kx * ky - press_part * gxy[idx];
    Sxz[idx] = kx * kz - press_part * gxz[idx];
    Syy[idx] = ky * ky - press_part * gyy;
    Syz[idx] = ky * kz - press_part * gyz[idx];
    Szz[idx] = kz * kz - press_part * gzz;

#undef SHC_FDER
#undef SHC_FDDER
}

template<int GW>
__global__ void bssn_escalar_ss_eval_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ crho, const double* __restrict__ sigma, const double* __restrict__ Rcoord,
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
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Axx, const double* __restrict__ Axy, const double* __restrict__ Axz,
    const double* __restrict__ Ayy, const double* __restrict__ Ayz, const double* __restrict__ Azz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    const double* __restrict__ Lap,
    const double* __restrict__ betax, const double* __restrict__ betay, const double* __restrict__ betaz,
    const double* __restrict__ dtSfx, const double* __restrict__ dtSfy, const double* __restrict__ dtSfz,
    const double* __restrict__ rho, const double* __restrict__ Sx, const double* __restrict__ Sy, const double* __restrict__ Sz,
    const double* __restrict__ Sxx, const double* __restrict__ Sxy, const double* __restrict__ Sxz,
    const double* __restrict__ Syy, const double* __restrict__ Syz, const double* __restrict__ Szz,
    const double* __restrict__ Gamxxx, const double* __restrict__ Gamxxy, const double* __restrict__ Gamxxz,
    const double* __restrict__ Gamxyy, const double* __restrict__ Gamxyz, const double* __restrict__ Gamxzz,
    const double* __restrict__ Gamyxx, const double* __restrict__ Gamyxy, const double* __restrict__ Gamyxz,
    const double* __restrict__ Gamyyy, const double* __restrict__ Gamyyz, const double* __restrict__ Gamyzz,
    const double* __restrict__ Gamzxx, const double* __restrict__ Gamzxy, const double* __restrict__ Gamzxz,
    const double* __restrict__ Gamzyy, const double* __restrict__ Gamzyz, const double* __restrict__ Gamzzz,
    const double* __restrict__ Rxx, const double* __restrict__ Rxy, const double* __restrict__ Rxz,
    const double* __restrict__ Ryy, const double* __restrict__ Ryz, const double* __restrict__ Rzz,
    double* __restrict__ chi_rhs, double* __restrict__ trK_rhs,
    double* __restrict__ gxx_rhs, double* __restrict__ gxy_rhs, double* __restrict__ gxz_rhs,
    double* __restrict__ gyy_rhs, double* __restrict__ gyz_rhs, double* __restrict__ gzz_rhs,
    double* __restrict__ Axx_rhs, double* __restrict__ Axy_rhs, double* __restrict__ Axz_rhs,
    double* __restrict__ Ayy_rhs, double* __restrict__ Ayz_rhs, double* __restrict__ Azz_rhs,
    double* __restrict__ Gamx_rhs, double* __restrict__ Gamy_rhs, double* __restrict__ Gamz_rhs,
    double* __restrict__ Lap_rhs,
    double* __restrict__ betax_rhs, double* __restrict__ betay_rhs, double* __restrict__ betaz_rhs,
    double* __restrict__ dtSfx_rhs, double* __restrict__ dtSfy_rhs, double* __restrict__ dtSfz_rhs,
    int symmetry, int lev, int sst)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    const int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

#define SHC_FDER(F,S1,S2,S3,FX,FY,FZ) \
    d_point_fderivs_shc_gw<GW>(ex, pex, F, FX, FY, FZ, crho, sigma, Rcoord, S1, S2, S3, symmetry, 0, sst, \
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k)
#define SHC_FDDER(F,S1,S2,S3,FXX,FXY,FXZ,FYY,FYZ,FZZ) \
    d_point_fdderivs_shc_gw<GW>(ex, pex, F, FXX, FXY, FXZ, FYY, FYZ, FZZ, crho, sigma, Rcoord, S1, S2, S3, symmetry, 0, sst, \
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, \
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz, \
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz, \
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, i, j, k)

    double betaxx, betaxy, betaxz, betayx, betayy, betayz, betazx, betazy, betazz;
    double bx_xx, bx_xy, bx_xz, bx_yy, bx_yz, bx_zz;
    double by_xx, by_xy, by_xz, by_yy, by_yz, by_zz;
    double bz_xx, bz_xy, bz_xz, bz_yy, bz_yz, bz_zz;
    SHC_FDDER(betax, -1.0, 1.0, 1.0, &bx_xx, &bx_xy, &bx_xz, &bx_yy, &bx_yz, &bx_zz);
    SHC_FDDER(betay, 1.0, -1.0, 1.0, &by_xx, &by_xy, &by_xz, &by_yy, &by_yz, &by_zz);
    SHC_FDDER(betaz, 1.0, 1.0, -1.0, &bz_xx, &bz_xy, &bz_xz, &bz_yy, &bz_yz, &bz_zz);
    SHC_FDER(betax, -1.0, 1.0, 1.0, &betaxx, &betaxy, &betaxz);
    SHC_FDER(betay, 1.0, -1.0, 1.0, &betayx, &betayy, &betayz);
    SHC_FDER(betaz, 1.0, 1.0, -1.0, &betazx, &betazy, &betazz);

    double chix, chiy, chiz, Lapx, Lapy, Lapz, Kx, Ky, Kz;
    double Lapxx, Lapxy, Lapxz, Lapyy, Lapyz, Lapzz;
    SHC_FDER(chi, 1.0, 1.0, 1.0, &chix, &chiy, &chiz);
    SHC_FDDER(Lap, 1.0, 1.0, 1.0, &Lapxx, &Lapxy, &Lapxz, &Lapyy, &Lapyz, &Lapzz);
    SHC_FDER(Lap, 1.0, 1.0, 1.0, &Lapx, &Lapy, &Lapz);
    SHC_FDER(trK, 1.0, 1.0, 1.0, &Kx, &Ky, &Kz);

    double gxxx, gxxy, gxxz, gxyx, gxyy, gxyz_m, gxzx, gxzy, gxzz_m;
    double gyyx, gyyy, gyyz, gyzx, gyzy, gyzz_m, gzzx, gzzy, gzzz_m;
    SHC_FDER(dxx, 1.0, 1.0, 1.0, &gxxx, &gxxy, &gxxz);
    SHC_FDER(gxy, -1.0, -1.0, 1.0, &gxyx, &gxyy, &gxyz_m);
    SHC_FDER(gxz, -1.0, 1.0, -1.0, &gxzx, &gxzy, &gxzz_m);
    SHC_FDER(dyy, 1.0, 1.0, 1.0, &gyyx, &gyyy, &gyyz);
    SHC_FDER(gyz, 1.0, -1.0, -1.0, &gyzx, &gyzy, &gyzz_m);
    SHC_FDER(dzz, 1.0, 1.0, 1.0, &gzzx, &gzzy, &gzzz_m);

    const double l_gxx = dxx[idx] + 1.0, l_gxy = gxy[idx], l_gxz = gxz[idx];
    const double l_gyy = dyy[idx] + 1.0, l_gyz = gyz[idx], l_gzz = dzz[idx] + 1.0;
    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    d_inverse_metric(l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);

    const double alpn1 = Lap[idx] + 1.0;
    const double chin1 = chi[idx] + 1.0;
    const double inv_chin1 = 1.0 / chin1;
    const double div_beta = betaxx + betayy + betazz;
    const double val_trK = trK[idx];

    chi_rhs[idx] = ESC_F2o3 * chin1 * (alpn1 * val_trK - div_beta);

    const double l_Axx = Axx[idx], l_Axy = Axy[idx], l_Axz = Axz[idx];
    const double l_Ayy = Ayy[idx], l_Ayz = Ayz[idx], l_Azz = Azz[idx];
    gxx_rhs[idx] = -2.0 * alpn1 * l_Axx - ESC_F2o3 * l_gxx * div_beta + 2.0 * (l_gxx * betaxx + l_gxy * betayx + l_gxz * betazx);
    gyy_rhs[idx] = -2.0 * alpn1 * l_Ayy - ESC_F2o3 * l_gyy * div_beta + 2.0 * (l_gxy * betaxy + l_gyy * betayy + l_gyz * betazy);
    gzz_rhs[idx] = -2.0 * alpn1 * l_Azz - ESC_F2o3 * l_gzz * div_beta + 2.0 * (l_gxz * betaxz + l_gyz * betayz + l_gzz * betazz);
    gxy_rhs[idx] = -2.0 * alpn1 * l_Axy + ESC_F1o3 * l_gxy * div_beta + l_gxx * betaxy + l_gxz * betazy + l_gyy * betayx + l_gyz * betazx - l_gxy * betazz;
    gyz_rhs[idx] = -2.0 * alpn1 * l_Ayz + ESC_F1o3 * l_gyz * div_beta + l_gxy * betaxz + l_gyy * betayz + l_gxz * betaxy + l_gzz * betazy - l_gyz * betaxx;
    gxz_rhs[idx] = -2.0 * alpn1 * l_Axz + ESC_F1o3 * l_gxz * div_beta + l_gxx * betaxz + l_gxy * betayz + l_gyz * betayx + l_gzz * betazx - l_gxz * betayy;

    const double cGxxx = 0.5 * (gupxx * gxxx + gupxy * (2.0 * gxyx - gxxy) + gupxz * (2.0 * gxzx - gxxz));
    const double cGyxx = 0.5 * (gupxy * gxxx + gupyy * (2.0 * gxyx - gxxy) + gupyz * (2.0 * gxzx - gxxz));
    const double cGzxx = 0.5 * (gupxz * gxxx + gupyz * (2.0 * gxyx - gxxy) + gupzz * (2.0 * gxzx - gxxz));
    const double cGxyy = 0.5 * (gupxx * (2.0 * gxyy - gyyx) + gupxy * gyyy + gupxz * (2.0 * gyzy - gyyz));
    const double cGyyy = 0.5 * (gupxy * (2.0 * gxyy - gyyx) + gupyy * gyyy + gupyz * (2.0 * gyzy - gyyz));
    const double cGzyy = 0.5 * (gupxz * (2.0 * gxyy - gyyx) + gupyz * gyyy + gupzz * (2.0 * gyzy - gyyz));
    const double cGxzz = 0.5 * (gupxx * (2.0 * gxzz_m - gzzx) + gupxy * (2.0 * gyzz_m - gzzy) + gupxz * gzzz_m);
    const double cGyzz = 0.5 * (gupxy * (2.0 * gxzz_m - gzzx) + gupyy * (2.0 * gyzz_m - gzzy) + gupyz * gzzz_m);
    const double cGzzz = 0.5 * (gupxz * (2.0 * gxzz_m - gzzx) + gupyz * (2.0 * gyzz_m - gzzy) + gupzz * gzzz_m);
    const double cGxxy = 0.5 * (gupxx * gxxy + gupxy * gyyx + gupxz * (gxzy + gyzx - gxyz_m));
    const double cGyxy = 0.5 * (gupxy * gxxy + gupyy * gyyx + gupyz * (gxzy + gyzx - gxyz_m));
    const double cGzxy = 0.5 * (gupxz * gxxy + gupyz * gyyx + gupzz * (gxzy + gyzx - gxyz_m));
    const double cGxxz = 0.5 * (gupxx * gxxz + gupxy * (gxyz_m + gyzx - gxzy) + gupxz * gzzx);
    const double cGyxz = 0.5 * (gupxy * gxxz + gupyy * (gxyz_m + gyzx - gxzy) + gupyz * gzzx);
    const double cGzxz = 0.5 * (gupxz * gxxz + gupyz * (gxyz_m + gyzx - gxzy) + gupzz * gzzx);
    const double cGxyz = 0.5 * (gupxx * (gxyz_m + gxzy - gyzx) + gupxy * gyyz + gupxz * gzzy);
    const double cGyyz = 0.5 * (gupxy * (gxyz_m + gxzy - gyzx) + gupyy * gyyz + gupyz * gzzy);
    const double cGzyz = 0.5 * (gupxz * (gxyz_m + gxzy - gyzx) + gupyz * gyyz + gupzz * gzzy);

    const double pGxxx = Gamxxx[idx], pGxxy = Gamxxy[idx], pGxxz = Gamxxz[idx];
    const double pGxyy = Gamxyy[idx], pGxyz = Gamxyz[idx], pGxzz = Gamxzz[idx];
    const double pGyxx = Gamyxx[idx], pGyxy = Gamyxy[idx], pGyxz = Gamyxz[idx];
    const double pGyyy = Gamyyy[idx], pGyyz = Gamyyz[idx], pGyzz = Gamyzz[idx];
    const double pGzxx = Gamzxx[idx], pGzxy = Gamzxy[idx], pGzxz = Gamzxz[idx];
    const double pGzyy = Gamzyy[idx], pGzyz = Gamzyz[idx], pGzzz = Gamzzz[idx];

    const double D2a_xx = Lapxx - pGxxx * Lapx - pGyxx * Lapy - pGzxx * Lapz;
    const double D2a_xy = Lapxy - pGxxy * Lapx - pGyxy * Lapy - pGzxy * Lapz;
    const double D2a_xz = Lapxz - pGxxz * Lapx - pGyxz * Lapy - pGzxz * Lapz;
    const double D2a_yy = Lapyy - pGxyy * Lapx - pGyyy * Lapy - pGzyy * Lapz;
    const double D2a_yz = Lapyz - pGxyz * Lapx - pGyyz * Lapy - pGzyz * Lapz;
    const double D2a_zz = Lapzz - pGxzz * Lapx - pGyzz * Lapy - pGzzz * Lapz;
    const double trD2a = gupxx * D2a_xx + gupyy * D2a_yy + gupzz * D2a_zz + 2.0 * (gupxy * D2a_xy + gupxz * D2a_xz + gupyz * D2a_yz);

    const double Au_xx = gupxx * l_Axx + gupxy * l_Axy + gupxz * l_Axz;
    const double Au_xy = gupxx * l_Axy + gupxy * l_Ayy + gupxz * l_Ayz;
    const double Au_xz = gupxx * l_Axz + gupxy * l_Ayz + gupxz * l_Azz;
    const double Au_yx = gupxy * l_Axx + gupyy * l_Axy + gupyz * l_Axz;
    const double Au_yy = gupxy * l_Axy + gupyy * l_Ayy + gupyz * l_Ayz;
    const double Au_yz = gupxy * l_Axz + gupyy * l_Ayz + gupyz * l_Azz;
    const double Au_zx = gupxz * l_Axx + gupyz * l_Axy + gupzz * l_Axz;
    const double Au_zy = gupxz * l_Axy + gupyz * l_Ayy + gupzz * l_Ayz;
    const double Au_zz = gupxz * l_Axz + gupyz * l_Ayz + gupzz * l_Azz;

    const double A2_xx = l_Axx * Au_xx + l_Axy * Au_yx + l_Axz * Au_zx;
    const double A2_xy = l_Axx * Au_xy + l_Axy * Au_yy + l_Axz * Au_zy;
    const double A2_xz = l_Axx * Au_xz + l_Axy * Au_yz + l_Axz * Au_zz;
    const double A2_yy = l_Axy * Au_xy + l_Ayy * Au_yy + l_Ayz * Au_zy;
    const double A2_yz = l_Axy * Au_xz + l_Ayy * Au_yz + l_Ayz * Au_zz;
    const double A2_zz = l_Axz * Au_xz + l_Ayz * Au_yz + l_Azz * Au_zz;
    const double trA2 = gupxx * A2_xx + gupyy * A2_yy + gupzz * A2_zz
        + 2.0 * (gupxy * A2_xy + gupxz * A2_xz + gupyz * A2_yz);

    const double s_trace = chin1 * (gupxx * Sxx[idx] + gupyy * Syy[idx] + gupzz * Szz[idx] + 2.0 * (gupxy * Sxy[idx] + gupxz * Sxz[idx] + gupyz * Syz[idx]));
    const double src_xx = alpn1 * (Rxx[idx] - ESC_EIGHT * M_PI * Sxx[idx]) - D2a_xx;
    const double src_xy = alpn1 * (Rxy[idx] - ESC_EIGHT * M_PI * Sxy[idx]) - D2a_xy;
    const double src_xz = alpn1 * (Rxz[idx] - ESC_EIGHT * M_PI * Sxz[idx]) - D2a_xz;
    const double src_yy = alpn1 * (Ryy[idx] - ESC_EIGHT * M_PI * Syy[idx]) - D2a_yy;
    const double src_yz = alpn1 * (Ryz[idx] - ESC_EIGHT * M_PI * Syz[idx]) - D2a_yz;
    const double src_zz = alpn1 * (Rzz[idx] - ESC_EIGHT * M_PI * Szz[idx]) - D2a_zz;
    const double tf_trace = ESC_F1o3 * (gupxx * src_xx + gupyy * src_yy + gupzz * src_zz + 2.0 * (gupxy * src_xy + gupxz * src_xz + gupyz * src_yz));

    Axx_rhs[idx] = chin1 * (src_xx - l_gxx * tf_trace) + alpn1 * (val_trK * l_Axx - 2.0 * A2_xx) + 2.0 * (l_Axx * betaxx + l_Axy * betayx + l_Axz * betazx) - ESC_F2o3 * l_Axx * div_beta;
    Axy_rhs[idx] = chin1 * (src_xy - l_gxy * tf_trace) + alpn1 * (val_trK * l_Axy - 2.0 * A2_xy) + l_Axx * betaxy + l_Axz * betazy + l_Ayy * betayx + l_Ayz * betazx + ESC_F1o3 * l_Axy * div_beta - l_Axy * betazz;
    Axz_rhs[idx] = chin1 * (src_xz - l_gxz * tf_trace) + alpn1 * (val_trK * l_Axz - 2.0 * A2_xz) + l_Axx * betaxz + l_Axy * betayz + l_Ayz * betayx + l_Azz * betazx + ESC_F1o3 * l_Axz * div_beta - l_Axz * betayy;
    Ayy_rhs[idx] = chin1 * (src_yy - l_gyy * tf_trace) + alpn1 * (val_trK * l_Ayy - 2.0 * A2_yy) + 2.0 * (l_Axy * betaxy + l_Ayy * betayy + l_Ayz * betazy) - ESC_F2o3 * l_Ayy * div_beta;
    Ayz_rhs[idx] = chin1 * (src_yz - l_gyz * tf_trace) + alpn1 * (val_trK * l_Ayz - 2.0 * A2_yz) + l_Axy * betaxz + l_Ayy * betayz + l_Axz * betaxy + l_Azz * betazy + ESC_F1o3 * l_Ayz * div_beta - l_Ayz * betaxx;
    Azz_rhs[idx] = chin1 * (src_zz - l_gzz * tf_trace) + alpn1 * (val_trK * l_Azz - 2.0 * A2_zz) + 2.0 * (l_Axz * betaxz + l_Ayz * betayz + l_Azz * betazz) - ESC_F2o3 * l_Azz * div_beta;

    trK_rhs[idx] = -chin1 * trD2a + alpn1 * (ESC_F1o3 * val_trK * val_trK + trA2 + ESC_FOUR_PI * (rho[idx] + s_trace));
    Lap_rhs[idx] = -2.0 * alpn1 * val_trK;

    const double Aup_xx = gupxx * gupxx * l_Axx + gupxy * gupxy * l_Ayy + gupxz * gupxz * l_Azz
        + 2.0 * (gupxx * gupxy * l_Axy + gupxx * gupxz * l_Axz + gupxy * gupxz * l_Ayz);
    const double Aup_yy = gupxy * gupxy * l_Axx + gupyy * gupyy * l_Ayy + gupyz * gupyz * l_Azz
        + 2.0 * (gupxy * gupyy * l_Axy + gupxy * gupyz * l_Axz + gupyy * gupyz * l_Ayz);
    const double Aup_zz = gupxz * gupxz * l_Axx + gupyz * gupyz * l_Ayy + gupzz * gupzz * l_Azz
        + 2.0 * (gupxz * gupyz * l_Axy + gupxz * gupzz * l_Axz + gupyz * gupzz * l_Ayz);
    const double Aup_xy = gupxx * gupxy * l_Axx + gupxy * gupyy * l_Ayy + gupxz * gupyz * l_Azz
        + (gupxx * gupyy + gupxy * gupxy) * l_Axy
        + (gupxx * gupyz + gupxz * gupxy) * l_Axz
        + (gupxy * gupyz + gupxz * gupyy) * l_Ayz;
    const double Aup_xz = gupxx * gupxz * l_Axx + gupxy * gupyz * l_Ayy + gupxz * gupzz * l_Azz
        + (gupxx * gupyz + gupxy * gupxz) * l_Axy
        + (gupxx * gupzz + gupxz * gupxz) * l_Axz
        + (gupxy * gupzz + gupxz * gupyz) * l_Ayz;
    const double Aup_yz = gupxy * gupxz * l_Axx + gupyy * gupyz * l_Ayy + gupyz * gupzz * l_Azz
        + (gupxy * gupyz + gupyy * gupxz) * l_Axy
        + (gupxy * gupzz + gupyz * gupxz) * l_Axz
        + (gupyy * gupzz + gupyz * gupyz) * l_Ayz;
    const double Gamxa = gupxx * cGxxx + gupyy * cGxyy + gupzz * cGxzz + 2.0 * (gupxy * cGxxy + gupxz * cGxxz + gupyz * cGxyz);
    const double Gamya = gupxx * cGyxx + gupyy * cGyyy + gupzz * cGyzz + 2.0 * (gupxy * cGyxy + gupxz * cGyxz + gupyz * cGyyz);
    const double Gamza = gupxx * cGzxx + gupyy * cGzyy + gupzz * cGzzz + 2.0 * (gupxy * cGzxy + gupxz * cGzxz + gupyz * cGzyz);
    const double fxx_beta = bx_xx + by_xy + bz_xz;
    const double fxy_beta = bx_xy + by_yy + bz_yz;
    const double fxz_beta = bx_xz + by_yz + bz_zz;
    const double partial_Gamx = ESC_F1o3 * (gupxx * fxx_beta + gupxy * fxy_beta + gupxz * fxz_beta)
        + gupxx * bx_xx + gupyy * bx_yy + gupzz * bx_zz
        + 2.0 * (gupxy * bx_xy + gupxz * bx_xz + gupyz * bx_yz);
    const double partial_Gamy = ESC_F1o3 * (gupxy * fxx_beta + gupyy * fxy_beta + gupyz * fxz_beta)
        + gupxx * by_xx + gupyy * by_yy + gupzz * by_zz
        + 2.0 * (gupxy * by_xy + gupxz * by_xz + gupyz * by_yz);
    const double partial_Gamz = ESC_F1o3 * (gupxz * fxx_beta + gupyz * fxy_beta + gupzz * fxz_beta)
        + gupxx * bz_xx + gupyy * bz_yy + gupzz * bz_zz
        + 2.0 * (gupxy * bz_xy + gupxz * bz_xz + gupyz * bz_yz);

    double gxrhs = -2.0 * (Lapx * Aup_xx + Lapy * Aup_xy + Lapz * Aup_xz)
        + 2.0 * alpn1 * (-ESC_F3o2 * inv_chin1 * (chix * Aup_xx + chiy * Aup_xy + chiz * Aup_xz)
        - gupxx * (ESC_F2o3 * Kx + ESC_EIGHT * M_PI * Sx[idx])
        - gupxy * (ESC_F2o3 * Ky + ESC_EIGHT * M_PI * Sy[idx])
        - gupxz * (ESC_F2o3 * Kz + ESC_EIGHT * M_PI * Sz[idx])
        + cGxxx * Aup_xx + cGxyy * Aup_yy + cGxzz * Aup_zz + 2.0 * (cGxxy * Aup_xy + cGxxz * Aup_xz + cGxyz * Aup_yz));
    double gyrhs = -2.0 * (Lapx * Aup_xy + Lapy * Aup_yy + Lapz * Aup_yz)
        + 2.0 * alpn1 * (-ESC_F3o2 * inv_chin1 * (chix * Aup_xy + chiy * Aup_yy + chiz * Aup_yz)
        - gupxy * (ESC_F2o3 * Kx + ESC_EIGHT * M_PI * Sx[idx])
        - gupyy * (ESC_F2o3 * Ky + ESC_EIGHT * M_PI * Sy[idx])
        - gupyz * (ESC_F2o3 * Kz + ESC_EIGHT * M_PI * Sz[idx])
        + cGyxx * Aup_xx + cGyyy * Aup_yy + cGyzz * Aup_zz + 2.0 * (cGyxy * Aup_xy + cGyxz * Aup_xz + cGyyz * Aup_yz));
    double gzrhs = -2.0 * (Lapx * Aup_xz + Lapy * Aup_yz + Lapz * Aup_zz)
        + 2.0 * alpn1 * (-ESC_F3o2 * inv_chin1 * (chix * Aup_xz + chiy * Aup_yz + chiz * Aup_zz)
        - gupxz * (ESC_F2o3 * Kx + ESC_EIGHT * M_PI * Sx[idx])
        - gupyz * (ESC_F2o3 * Ky + ESC_EIGHT * M_PI * Sy[idx])
        - gupzz * (ESC_F2o3 * Kz + ESC_EIGHT * M_PI * Sz[idx])
        + cGzxx * Aup_xx + cGzyy * Aup_yy + cGzzz * Aup_zz + 2.0 * (cGzxy * Aup_xy + cGzxz * Aup_xz + cGzyz * Aup_yz));
    gxrhs += ESC_F2o3 * Gamxa * div_beta - (Gamxa * betaxx + Gamya * betaxy + Gamza * betaxz)
          + partial_Gamx;
    gyrhs += ESC_F2o3 * Gamya * div_beta - (Gamxa * betayx + Gamya * betayy + Gamza * betayz)
          + partial_Gamy;
    gzrhs += ESC_F2o3 * Gamza * div_beta - (Gamxa * betazx + Gamya * betazy + Gamza * betazz)
          + partial_Gamz;
    Gamx_rhs[idx] = gxrhs; Gamy_rhs[idx] = gyrhs; Gamz_rhs[idx] = gzrhs;

    betax_rhs[idx] = ESC_FF * dtSfx[idx];
    betay_rhs[idx] = ESC_FF * dtSfy[idx];
    betaz_rhs[idx] = ESC_FF * dtSfz[idx];
    dtSfx_rhs[idx] = gxrhs - ESC_ETA * dtSfx[idx];
    dtSfy_rhs[idx] = gyrhs - ESC_ETA * dtSfy[idx];
    dtSfz_rhs[idx] = gzrhs - ESC_ETA * dtSfz[idx];

#undef SHC_FDER
#undef SHC_FDDER
}

template<int GW>
__global__ void bssn_escalar_ss_advect_ko_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ crho, const double* __restrict__ sigma, const double* __restrict__ Rcoord,
    const double* __restrict__ drhodx, const double* __restrict__ drhody, const double* __restrict__ drhodz,
    const double* __restrict__ dsigmadx, const double* __restrict__ dsigmady, const double* __restrict__ dsigmadz,
    const double* __restrict__ dRdx, const double* __restrict__ dRdy, const double* __restrict__ dRdz,
    const double* __restrict__ betax, const double* __restrict__ betay, const double* __restrict__ betaz,
    const double* __restrict__ f, double* __restrict__ f_rhs,
    double sym1, double sym2, double sym3, int symmetry, int lev, double eps, int sst)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    const int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    double fx, fy, fz;
    d_point_fderivs_shc_gw<GW>(ex, pex, f, &fx, &fy, &fz, crho, sigma, Rcoord,
                               sym1, sym2, sym3, symmetry, 0, sst,
                               drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);
    f_rhs[idx] += betax[idx] * fx + betay[idx] * fy + betaz[idx] * fz;
    if (eps > 0.0) {
        f_rhs[idx] += d_kodis_sh_point<GW>(ex, pex, f, crho, sigma, Rcoord, sym1, sym2, sym3, symmetry, eps, sst, i, j, k);
    }
}

static inline void scalar_grid(int* ex, dim3& grid, dim3& block)
{
    block = dim3(8, 8, 4);
    grid = dim3((ex[0] + block.x - 1) / block.x,
                (ex[1] + block.y - 1) / block.y,
                (ex[2] + block.z - 1) / block.z);
}

void gpu_compute_rhs_scalar_launch(cudaStream_t& stream, int* ex, int* pex, double,
                                   double* X, double* Y, double* Z,
                                   double* Sphi, double* Spi,
                                   double* Sphi_rhs, double* Spi_rhs,
                                   int symmetry, int lev, double eps)
{
    dim3 grid, block;
    scalar_grid(ex, grid, block);
    scalar_rhs_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        X, Y, Z, Sphi, Spi, Sphi_rhs, Spi_rhs, symmetry, lev, eps);
}

void gpu_compute_rhs_scalar_ss_launch(cudaStream_t& stream, int* ex, int* pex, double,
                                      double* crho, double* sigma, double* R,
                                      double*, double*, double*,
                                      double* drhodx, double* drhody, double* drhodz,
                                      double* dsigmadx, double* dsigmady, double* dsigmadz,
                                      double* dRdx, double* dRdy, double* dRdz,
                                      double* drhodxx, double* drhodxy, double* drhodxz,
                                      double* drhodyy, double* drhodyz, double* drhodzz,
                                      double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
                                      double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
                                      double* dRdxx, double* dRdxy, double* dRdxz,
                                      double* dRdyy, double* dRdyz, double* dRdzz,
                                      double* Sphi, double* Spi,
                                      double* Sphi_rhs, double* Spi_rhs,
                                      int symmetry, int lev, double eps, int sst)
{
    dim3 grid, block;
    scalar_grid(ex, grid, block);
    scalar_rhs_ss_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        crho, sigma, R,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        Sphi, Spi, Sphi_rhs, Spi_rhs, symmetry, lev, eps, sst);
}

void gpu_compute_rhs_bssn_escalar_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* Sphi, double* Spi,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* Sphi_rhs, double* Spi_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co)
{
    dim3 grid, block;
    scalar_grid(ex, grid, block);
    double a2, phi0, r0, sigma0, l2;
    load_escalar_parameters(a2, phi0, r0, sigma0, l2);
    escalar_matter_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz, Lap,
        Sphi, Spi, Sphi_rhs, Spi_rhs, rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        a2, symmetry, lev);

    gpu_compute_rhs_bssn_launch(
        stream, ex, pex, T, X, Y, Z, chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, Gamx, Gamy, Gamz, Lap, betax, betay, betaz,
        dtSfx, dtSfy, dtSfz, chi_rhs, trK_rhs, gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs, Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs, dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz, Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, ham_Res, movx_Res, movy_Res, movz_Res,
        Gmx_Res, Gmy_Res, Gmz_Res, symmetry, lev, eps, co);

    escalar_scalar_advect_ko_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], X, Y, Z,
        betax, betay, betaz, Sphi, Spi, Sphi_rhs, Spi_rhs, symmetry, eps);
}

void gpu_compute_rhs_bssn_escalar_ss_launch(
    cudaStream_t &stream, int* ex, int* pex, double,
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
    double* Sphi, double* Spi,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* Sphi_rhs, double* Spi_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int sst, int co)
{
    dim3 grid, block;
    scalar_grid(ex, grid, block);
    double a2, phi0, r0, sigma0, l2;
    load_escalar_parameters(a2, phi0, r0, sigma0, l2);

    escalar_matter_ss_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        crho, sigma, R,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz, Lap,
        Sphi, Spi, Sphi_rhs, Spi_rhs, rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        a2, symmetry, lev, sst);

    gpu_ricci_gamma_ss_launch(
        stream, ex, pex, crho, sigma, R, chi,
        dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, symmetry, sst);

    bssn_escalar_ss_eval_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        crho, sigma, R,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, Gamx, Gamy, Gamz, Lap,
        betax, betay, betaz, dtSfx, dtSfy, dtSfz,
        rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        chi_rhs, trK_rhs, gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs, dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        symmetry, lev, sst);

#define LAUNCH_SH_ADV(F, RHS, SX, SY, SZ) \
    bssn_escalar_ss_advect_ko_kernel<ghost_width><<<grid, block, 0, stream>>>( \
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], crho, sigma, R, \
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, \
        betax, betay, betaz, F, RHS, SX, SY, SZ, symmetry, lev, eps, sst)

    LAUNCH_SH_ADV(chi, chi_rhs, 1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(trK, trK_rhs, 1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(dxx, gxx_rhs, 1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(gxy, gxy_rhs, -1.0, -1.0, 1.0);
    LAUNCH_SH_ADV(gxz, gxz_rhs, -1.0, 1.0, -1.0);
    LAUNCH_SH_ADV(dyy, gyy_rhs, 1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(gyz, gyz_rhs, 1.0, -1.0, -1.0);
    LAUNCH_SH_ADV(dzz, gzz_rhs, 1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(Axx, Axx_rhs, 1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(Axy, Axy_rhs, -1.0, -1.0, 1.0);
    LAUNCH_SH_ADV(Axz, Axz_rhs, -1.0, 1.0, -1.0);
    LAUNCH_SH_ADV(Ayy, Ayy_rhs, 1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(Ayz, Ayz_rhs, 1.0, -1.0, -1.0);
    LAUNCH_SH_ADV(Azz, Azz_rhs, 1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(Gamx, Gamx_rhs, -1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(Gamy, Gamy_rhs, 1.0, -1.0, 1.0);
    LAUNCH_SH_ADV(Gamz, Gamz_rhs, 1.0, 1.0, -1.0);
    LAUNCH_SH_ADV(Lap, Lap_rhs, 1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(betax, betax_rhs, -1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(betay, betay_rhs, 1.0, -1.0, 1.0);
    LAUNCH_SH_ADV(betaz, betaz_rhs, 1.0, 1.0, -1.0);
    LAUNCH_SH_ADV(dtSfx, dtSfx_rhs, -1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(dtSfy, dtSfy_rhs, 1.0, -1.0, 1.0);
    LAUNCH_SH_ADV(dtSfz, dtSfz_rhs, 1.0, 1.0, -1.0);
    LAUNCH_SH_ADV(Sphi, Sphi_rhs, 1.0, 1.0, 1.0);
    LAUNCH_SH_ADV(Spi, Spi_rhs, 1.0, 1.0, 1.0);
#undef LAUNCH_SH_ADV

    if (co == 0) {
        gpu_constraint_bssn_ss_launch(
            stream, ex, pex, crho, sigma, R, X, Y, Z,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
            drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
            dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
            dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
            chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
            Axx, Axy, Axz, Ayy, Ayz, Azz, Gamx, Gamy, Gamz,
            Lap, dtSfx, dtSfy, dtSfz, rho, Sx, Sy, Sz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
            ham_Res, movx_Res, movy_Res, movz_Res, Gmx_Res, Gmy_Res, Gmz_Res,
            symmetry, lev, sst);
    }
}

// ============================================================================
// EScalar fR-constraint (Cons_fR) — full GPU pipeline
//
// Mirrors bssnEScalar_rhs.f90::compute_constraint_fr:
//   call get4ricciscalar(...)              ! 4D Ricci scalar -> RR
//   call frfprim(RR, fprim)                ! f'(R) = 1 + 2*a2*R (for CC=2/3)
//   Cons_fR = Sphi - sqrt(3/pi)/4 * log(fprim)
//
// Caller supplies a scratch device buffer of the block's padded extent for RR.
// We use Cons_fR itself as the fprim destination, so the finalize kernel reads
// it and writes back in place (pointwise — safe).
// ============================================================================

// Forward declarations of host-side launch helpers defined in other TUs.
extern "C" void gpu_get4ricciscalar_launch(
    cudaStream_t stream, int ex[3], int pex[3],
    double* chi, double* trK, double* rho,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* RR);

extern "C" void gpu_frfprim_launch(
    cudaStream_t stream, int ex[3], int pex[3],
    double* RR, double* fprim,
    double a2, double r0, double phi0, double sigma, double l2);

__global__ void finalize_cons_fr_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* __restrict__ Sphi, double* __restrict__ Cons_fR)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int idx = IDX3D_PAD(i, j, k, pex0, pex1);
    // Cons_fR currently holds fprim from gpu_frfprim_launch.
    Cons_fR[idx] = Sphi[idx] - sqrt(3.0 / M_PI) * 0.25 * log(Cons_fR[idx]);
}

extern "C" void gpu_compute_constraint_fr_launch(
    cudaStream_t stream, int ex[3], int pex[3],
    double* chi, double* trK, double* rho, double* Sphi,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Cons_fR, double* d_scratch_RR)
{
    double a2, phi0_p, r0, sigma0, l2;
    load_escalar_parameters(a2, phi0_p, r0, sigma0, l2);

    // 1. RR = 4D Ricci scalar (written into the scratch buffer)
    gpu_get4ricciscalar_launch(stream, ex, pex, chi, trK, rho,
                               dxx, gxy, gxz, dyy, gyz, dzz,
                               Axx, Axy, Axz, Ayy, Ayz, Azz,
                               Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
                               Sxx, Sxy, Sxz, Syy, Syz, Szz,
                               d_scratch_RR);

    // 2. f'(R) -> Cons_fR (initial)
    gpu_frfprim_launch(stream, ex, pex, d_scratch_RR, Cons_fR,
                       a2, r0, phi0_p, sigma0, l2);

    // 3. Finalise: Cons_fR = Sphi - sqrt(3/pi)/4 * log(fprim)
    dim3 block(8, 8, 4);
    dim3 grid((ex[0] + block.x - 1) / block.x,
              (ex[1] + block.y - 1) / block.y,
              (ex[2] + block.z - 1) / block.z);
    finalize_cons_fr_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], Sphi, Cons_fR);
}
