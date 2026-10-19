#ifndef BSSN_ESCALAR_RHS_GPU_CUH
#define BSSN_ESCALAR_RHS_GPU_CUH

#include <math.h>

#ifndef MACRODEF_FH_INCLUDED
#define MACRODEF_FH_INCLUDED
#include "macrodef.fh"
#endif

#ifndef IDX3D_PAD
#define IDX3D_PAD(i, j, k, pnx, pny) ((i) + (pnx) * ((j) + (pny) * (k)))
#endif

constexpr double ESC_F1o3 = 1.0 / 3.0;
constexpr double ESC_F2o3 = 2.0 / 3.0;
constexpr double ESC_F3o2 = 1.5;
constexpr double ESC_EIGHT = 8.0;
constexpr double ESC_FOUR_PI = 4.0 * M_PI;
constexpr double ESC_FF = 0.75;
constexpr double ESC_ETA = 2.0;

__device__ __forceinline__ void d_frpotential(double sphi, double a2, double& v, double& dvdsphi)
{
#if (EScalar_CC == 2 || EScalar_CC == 3)
    const double pi = M_PI;
    const double e4 = exp(4.0 * sqrt(pi / 3.0) * sphi);
    const double e8 = exp(-8.0 * sqrt(pi / 3.0) * sphi);
    v = e8 * (1.0 - e4) * (1.0 - e4) / (32.0 * pi * a2);
    dvdsphi = sqrt(3.0 / pi) * e8 * (-1.0 + e4) / (12.0 * a2);
#else
    (void)sphi;
    v = 0.0;
    dvdsphi = 0.0;
#endif
}

__device__ __forceinline__ void d_inverse_metric(
    double gxx, double gxy, double gxz,
    double gyy, double gyz, double gzz,
    double& gupxx, double& gupxy, double& gupxz,
    double& gupyy, double& gupyz, double& gupzz)
{
    const double det = gxx * gyy * gzz + 2.0 * gxy * gyz * gxz
                     - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz;
    const double inv_det = 1.0 / det;
    gupxx = (gyy * gzz - gyz * gyz) * inv_det;
    gupxy = (gyz * gxz - gxy * gzz) * inv_det;
    gupxz = (gxy * gyz - gyy * gxz) * inv_det;
    gupyy = (gxx * gzz - gxz * gxz) * inv_det;
    gupyz = (gxy * gxz - gxx * gyz) * inv_det;
    gupzz = (gxx * gyy - gxy * gxy) * inv_det;
}

#endif /* BSSN_ESCALAR_RHS_GPU_CUH */
