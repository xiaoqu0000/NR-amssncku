#ifndef BSSN2ADM_GPU_CUH
#define BSSN2ADM_GPU_CUH

#include <cmath>

// Convert BSSN conformal variables to ADM physical variables.
// Ported from bssn2adm.f90.
// All inputs are conformal (tilde) BSSN quantities; outputs are ADM.
__device__ __forceinline__ void d_bssn2adm(
    double chi,
    double trK,
    double gxx, double gxy, double gxz,
    double gyy, double gyz, double gzz,
    double Axx, double Axy, double Axz,
    double Ayy, double Ayz, double Azz,
    double& adm_gxx, double& adm_gxy, double& adm_gxz,
    double& adm_gyy, double& adm_gyz, double& adm_gzz,
    double& Kxx, double& Kxy, double& Kxz,
    double& Kyy, double& Kyz, double& Kzz
) {
    const double F1o3 = 1.0 / 3.0;
    const double inv_chi = 1.0 / chi;

    adm_gxx = gxx * inv_chi;
    adm_gxy = gxy * inv_chi;
    adm_gxz = gxz * inv_chi;
    adm_gyy = gyy * inv_chi;
    adm_gyz = gyz * inv_chi;
    adm_gzz = gzz * inv_chi;

    Kxx = Axx * inv_chi + F1o3 * trK * adm_gxx;
    Kxy = Axy * inv_chi + F1o3 * trK * adm_gxy;
    Kxz = Axz * inv_chi + F1o3 * trK * adm_gxz;
    Kyy = Ayy * inv_chi + F1o3 * trK * adm_gyy;
    Kyz = Ayz * inv_chi + F1o3 * trK * adm_gyz;
    Kzz = Azz * inv_chi + F1o3 * trK * adm_gzz;
}

#endif // BSSN2ADM_GPU_CUH
