#ifndef BSSN_CONSTRAINT_GPU_CUH
#define BSSN_CONSTRAINT_GPU_CUH

#include "adm_ricci_gamma_gpu.cuh"
#include "bssn2adm_gpu.cuh"
#include "diff_new_gpu.cuh"
#include "gpu_common.cuh"

#include <math_constants.h>

namespace bssn_constraint_gpu {

constexpr double ONE = 1.0;
constexpr double TWO = 2.0;
constexpr double F2O3 = 2.0 / 3.0;
constexpr double F8 = 8.0;
constexpr double F16 = 16.0;
constexpr double SYM = 1.0;
constexpr double ANTI = -1.0;

__device__ __forceinline__ size_t idx3(int i, int j, int k, const int pex[3])
{
    return static_cast<size_t>(i) + static_cast<size_t>(pex[0]) *
        (static_cast<size_t>(j) + static_cast<size_t>(pex[1]) * static_cast<size_t>(k));
}

__device__ __forceinline__ void invert_metric(
    double gxx, double gxy, double gxz, double gyy, double gyz, double gzz,
    double& gupxx, double& gupxy, double& gupxz,
    double& gupyy, double& gupyz, double& gupzz)
{
    const double det = gxx * gyy * gzz + TWO * gxy * gyz * gxz
                     - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz;
    const double inv_det = ONE / det;
    gupxx =  (gyy * gzz - gyz * gyz) * inv_det;
    gupxy = -(gxy * gzz - gyz * gxz) * inv_det;
    gupxz =  (gxy * gyz - gyy * gxz) * inv_det;
    gupyy =  (gxx * gzz - gxz * gxz) * inv_det;
    gupyz = -(gxx * gyz - gxy * gxz) * inv_det;
    gupzz =  (gxx * gyy - gxy * gxy) * inv_det;
}

__device__ __forceinline__ double fetch_lower_symmetric(
    const int ex[3], const int pex[3], const double* f,
    int i, int j, int k, double sym1, double sym2, double sym3)
{
    double factor = ONE;
    if (i < 0) { i = -i - 1; factor *= sym1; }
    if (j < 0) { j = -j - 1; factor *= sym2; }
    if (k < 0) { k = -k - 1; factor *= sym3; }
    i = (i < 0) ? 0 : ((i >= ex[0]) ? ex[0] - 1 : i);
    j = (j < 0) ? 0 : ((j >= ex[1]) ? ex[1] - 1 : j);
    k = (k < 0) ? 0 : ((k >= ex[2]) ? ex[2] - 1 : k);
    return factor * f[idx3(i, j, k, pex)];
}

__device__ __forceinline__ bool lower_symmetry_enabled(
    int axis, int symmetry, const double* X, const double* Y, const double* Z)
{
    if (axis == 2) return symmetry > 0 && fabs(Z[0]) < fabs(Z[1] - Z[0]);
    if (axis == 0) return symmetry > 1 && fabs(X[0]) < fabs(X[1] - X[0]);
    return symmetry > 1 && fabs(Y[0]) < fabs(Y[1] - Y[0]);
}

template<int Axis>
__device__ __forceinline__ double deriv_axis_safe(
    const int ex[3], const int pex[3], const double* f,
    int i, int j, int k, double h, bool lower_sym, bool upper_sym,
    double sym1, double sym2, double sym3)
{
    const int n = ex[Axis];
    const int q = (Axis == 0) ? i : ((Axis == 1) ? j : k);
    const auto val = [&](int qq) {
        int ii = i, jj = j, kk = k;
        if constexpr (Axis == 0) ii = qq;
        if constexpr (Axis == 1) jj = qq;
        if constexpr (Axis == 2) kk = qq;
        if (qq >= n && upper_sym) {
            qq = 2 * n - qq - 1;
            if constexpr (Axis == 0) { ii = qq; return sym1 * fetch_lower_symmetric(ex, pex, f, ii, jj, kk, sym1, sym2, sym3); }
            if constexpr (Axis == 1) { jj = qq; return sym2 * fetch_lower_symmetric(ex, pex, f, ii, jj, kk, sym1, sym2, sym3); }
            if constexpr (Axis == 2) { kk = qq; return sym3 * fetch_lower_symmetric(ex, pex, f, ii, jj, kk, sym1, sym2, sym3); }
        }
        return fetch_lower_symmetric(ex, pex, f, ii, jj, kk, sym1, sym2, sym3);
    };

    if (n <= 1) return 0.0;
    if (q >= 2 && q + 2 < n) {
        return (val(q - 2) - 8.0 * val(q - 1) + 8.0 * val(q + 1) - val(q + 2)) / (12.0 * h);
    }
    if (q > 0 && q + 1 < n) {
        return (val(q + 1) - val(q - 1)) / (2.0 * h);
    }
    if (q == n - 1 && upper_sym && n > 1) {
        return (val(n) - val(n - 2)) / (2.0 * h);
    }
    if (q == 0 && lower_sym && n > 1) {
        return (val(1) - val(-1)) / (2.0 * h);
    }
    if (q == 0) {
        if (n >= 3) return (-3.0 * val(0) + 4.0 * val(1) - val(2)) / (2.0 * h);
        return (val(1) - val(0)) / h;
    }
    if (n >= 3) return (3.0 * val(n - 1) - 4.0 * val(n - 2) + val(n - 3)) / (2.0 * h);
    return (val(n - 1) - val(n - 2)) / h;
}

template<int GW>
__device__ __forceinline__ void first_derivs_cart(
    const int ex[3], const int pex[3], const double* f,
    double& fx, double& fy, double& fz,
    const double* X, const double* Y, const double* Z,
    double sym1, double sym2, double sym3, int symmetry,
    int i, int j, int k)
{
    d_fderivs_point_gw<GW>(
        ex, pex, f, &fx, &fy, &fz,
        X, Y, Z, sym1, sym2, sym3, symmetry, 0, i, j, k);
}

__device__ __forceinline__ bool shell_lower_symmetry_enabled(
    int axis, int symmetry, int sst, const double* crho, const double* sigma)
{
    if (symmetry == 2) {
        if (axis == 0) return fabs(crho[0]) < fabs(crho[1] - crho[0]);
        if (axis == 1) return fabs(sigma[0]) < fabs(sigma[1] - sigma[0]);
    }
    if (symmetry == 1 && axis == 1) {
        return (sst == 2 || sst == 4) && fabs(sigma[0]) < fabs(sigma[1] - sigma[0]);
    }
    return false;
}

__device__ __forceinline__ bool shell_upper_symmetry_enabled(
    int axis, int symmetry, int sst, const int ex[3], const double* sigma)
{
    return symmetry == 1 && axis == 1 && (sst == 3 || sst == 5) &&
           fabs(sigma[ex[1] - 1]) < fabs(sigma[1] - sigma[0]);
}

template<int GW>
__device__ __forceinline__ void first_derivs_shell(
    const int ex[3], const int pex[3], const double* f,
    double& fx, double& fy, double& fz,
    const double* crho, const double* sigma, const double* R,
    const double* drhodx, const double* drhody, const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx, const double* dRdy, const double* dRdz,
    double sym1, double sym2, double sym3, int symmetry,
    int sst, int i, int j, int k)
{
    double frho, fsigma, fR;
    const double drho = crho[1] - crho[0];
    const double dsigma = sigma[1] - sigma[0];
    const double dRad = R[1] - R[0];
    frho = deriv_axis_safe<0>(ex, pex, f, i, j, k, drho,
                              shell_lower_symmetry_enabled(0, symmetry, sst, crho, sigma),
                              false, sym1, sym2, sym3);
    fsigma = deriv_axis_safe<1>(ex, pex, f, i, j, k, dsigma,
                                shell_lower_symmetry_enabled(1, symmetry, sst, crho, sigma),
                                shell_upper_symmetry_enabled(1, symmetry, sst, ex, sigma),
                                sym1, sym2, sym3);
    fR = deriv_axis_safe<2>(ex, pex, f, i, j, k, dRad, false, false, sym1, sym2, sym3);
    const size_t id = idx3(i, j, k, pex);
    fx = frho * drhodx[id] + fsigma * dsigmadx[id] + fR * dRdx[id];
    fy = frho * drhody[id] + fsigma * dsigmady[id] + fR * dRdy[id];
    fz = frho * drhodz[id] + fsigma * dsigmadz[id] + fR * dRdz[id];
}

__device__ __forceinline__ double contract_symmetric_square(
    double gupxx, double gupxy, double gupxz, double gupyy, double gupyz, double gupzz,
    double Axx, double Axy, double Axz, double Ayy, double Ayz, double Azz)
{
    return gupxx * (
               gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz +
               TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz)) +
           gupyy * (
               gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz +
               TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz)) +
           gupzz * (
               gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz +
               TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz)) +
           TWO * (
               gupxy * (
                   gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz +
                   gupxy * (Axx * Ayy + Axy * Axy) +
                   gupxz * (Axx * Ayz + Axz * Axy) +
                   gupyz * (Axy * Ayz + Axz * Ayy)) +
               gupxz * (
                   gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz +
                   gupxy * (Axx * Ayz + Axy * Axz) +
                   gupxz * (Axx * Azz + Axz * Axz) +
                   gupyz * (Axy * Azz + Axz * Ayz)) +
               gupyz * (
                   gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz +
                   gupxy * (Axy * Ayz + Ayy * Axz) +
                   gupxz * (Axy * Azz + Ayz * Axz) +
                   gupyz * (Ayy * Azz + Ayz * Ayz)));
}

__device__ __forceinline__ void gamma_constraint_residual(
    double gupxx, double gupxy, double gupxz, double gupyy, double gupyz, double gupzz,
    double gxxx, double gxyx, double gxzx, double gyyx, double gyzx, double gzzx,
    double gxxy, double gxyy, double gxzy, double gyyy, double gyzy, double gzzy,
    double gxxz, double gxyz, double gxzz, double gyyz, double gyzz, double gzzz,
    double Gmx, double Gmy, double Gmz,
    double& Gmx_res, double& Gmy_res, double& Gmz_res)
{
    Gmx_res = Gmx - (gupxx*(gupxx*gxxx+gupxy*gxyx+gupxz*gxzx)
                    +gupxy*(gupxx*gxyx+gupxy*gyyx+gupxz*gyzx)
                    +gupxz*(gupxx*gxzx+gupxy*gyzx+gupxz*gzzx)
                    +gupxx*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)
                    +gupxy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)
                    +gupxz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)
                    +gupxx*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)
                    +gupxy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)
                    +gupxz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz));
    Gmy_res = Gmy - (gupxx*(gupxy*gxxx+gupyy*gxyx+gupyz*gxzx)
                    +gupxy*(gupxy*gxyx+gupyy*gyyx+gupyz*gyzx)
                    +gupxz*(gupxy*gxzx+gupyy*gyzx+gupyz*gzzx)
                    +gupxy*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)
                    +gupyy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)
                    +gupyz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)
                    +gupxy*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)
                    +gupyy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)
                    +gupyz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz));
    Gmz_res = Gmz - (gupxx*(gupxz*gxxx+gupyz*gxyx+gupzz*gxzx)
                    +gupxy*(gupxz*gxyx+gupyz*gyyx+gupzz*gyzx)
                    +gupxz*(gupxz*gxzx+gupyz*gyzx+gupzz*gzzx)
                    +gupxy*(gupxz*gxxy+gupyz*gxyy+gupzz*gxzy)
                    +gupyy*(gupxz*gxyy+gupyz*gyyy+gupzz*gyzy)
                    +gupyz*(gupxz*gxzy+gupyz*gyzy+gupzz*gzzy)
                    +gupxz*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)
                    +gupyz*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)
                    +gupzz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz));
}

__device__ __forceinline__ void covariant_momentum(
    double gupxx, double gupxy, double gupxz, double gupyy, double gupyz, double gupzz,
    double Txx, double Txy, double Txz, double Tyy, double Tyz, double Tzz,
    double Txxx, double Txyx, double Txzx, double Tyyx, double Tyzx, double Tzzx,
    double Txxy, double Txyy, double Txzy, double Tyyy, double Tyzy, double Tzzy,
    double Txxz, double Txyz, double Txzz, double Tyyz, double Tyzz, double Tzzz,
    double Gamxxx, double Gamxxy, double Gamxxz, double Gamxyy, double Gamxyz, double Gamxzz,
    double Gamyxx, double Gamyxy, double Gamyxz, double Gamyyy, double Gamyyz, double Gamyzz,
    double Gamzxx, double Gamzxy, double Gamzxz, double Gamzyy, double Gamzyz, double Gamzzz,
    double chix, double chiy, double chiz, double inv_chi_factor,
    double& movx, double& movy, double& movz)
{
    Txxx = Txxx - (Gamxxx*Txx + Gamyxx*Txy + Gamzxx*Txz + Gamxxx*Txx + Gamyxx*Txy + Gamzxx*Txz) - chix*Txx*inv_chi_factor;
    Txyx = Txyx - (Gamxxy*Txx + Gamyxy*Txy + Gamzxy*Txz + Gamxxx*Txy + Gamyxx*Tyy + Gamzxx*Tyz) - chix*Txy*inv_chi_factor;
    Txzx = Txzx - (Gamxxz*Txx + Gamyxz*Txy + Gamzxz*Txz + Gamxxx*Txz + Gamyxx*Tyz + Gamzxx*Tzz) - chix*Txz*inv_chi_factor;
    Tyyx = Tyyx - (Gamxxy*Txy + Gamyxy*Tyy + Gamzxy*Tyz + Gamxxy*Txy + Gamyxy*Tyy + Gamzxy*Tyz) - chix*Tyy*inv_chi_factor;
    Tyzx = Tyzx - (Gamxxz*Txy + Gamyxz*Tyy + Gamzxz*Tyz + Gamxxy*Txz + Gamyxy*Tyz + Gamzxy*Tzz) - chix*Tyz*inv_chi_factor;
    Tzzx = Tzzx - (Gamxxz*Txz + Gamyxz*Tyz + Gamzxz*Tzz + Gamxxz*Txz + Gamyxz*Tyz + Gamzxz*Tzz) - chix*Tzz*inv_chi_factor;
    Txxy = Txxy - (Gamxxy*Txx + Gamyxy*Txy + Gamzxy*Txz + Gamxxy*Txx + Gamyxy*Txy + Gamzxy*Txz) - chiy*Txx*inv_chi_factor;
    Txyy = Txyy - (Gamxyy*Txx + Gamyyy*Txy + Gamzyy*Txz + Gamxxy*Txy + Gamyxy*Tyy + Gamzxy*Tyz) - chiy*Txy*inv_chi_factor;
    Txzy = Txzy - (Gamxyz*Txx + Gamyyz*Txy + Gamzyz*Txz + Gamxxy*Txz + Gamyxy*Tyz + Gamzxy*Tzz) - chiy*Txz*inv_chi_factor;
    Tyyy = Tyyy - (Gamxyy*Txy + Gamyyy*Tyy + Gamzyy*Tyz + Gamxyy*Txy + Gamyyy*Tyy + Gamzyy*Tyz) - chiy*Tyy*inv_chi_factor;
    Tyzy = Tyzy - (Gamxyz*Txy + Gamyyz*Tyy + Gamzyz*Tyz + Gamxyy*Txz + Gamyyy*Tyz + Gamzyy*Tzz) - chiy*Tyz*inv_chi_factor;
    Tzzy = Tzzy - (Gamxyz*Txz + Gamyyz*Tyz + Gamzyz*Tzz + Gamxyz*Txz + Gamyyz*Tyz + Gamzyz*Tzz) - chiy*Tzz*inv_chi_factor;
    Txxz = Txxz - (Gamxxz*Txx + Gamyxz*Txy + Gamzxz*Txz + Gamxxz*Txx + Gamyxz*Txy + Gamzxz*Txz) - chiz*Txx*inv_chi_factor;
    Txyz = Txyz - (Gamxyz*Txx + Gamyyz*Txy + Gamzyz*Txz + Gamxxz*Txy + Gamyxz*Tyy + Gamzxz*Tyz) - chiz*Txy*inv_chi_factor;
    Txzz = Txzz - (Gamxzz*Txx + Gamyzz*Txy + Gamzzz*Txz + Gamxxz*Txz + Gamyxz*Tyz + Gamzxz*Tzz) - chiz*Txz*inv_chi_factor;
    Tyyz = Tyyz - (Gamxyz*Txy + Gamyyz*Tyy + Gamzyz*Tyz + Gamxyz*Txy + Gamyyz*Tyy + Gamzyz*Tyz) - chiz*Tyy*inv_chi_factor;
    Tyzz = Tyzz - (Gamxzz*Txy + Gamyzz*Tyy + Gamzzz*Tyz + Gamxyz*Txz + Gamyyz*Tyz + Gamzyz*Tzz) - chiz*Tyz*inv_chi_factor;
    Tzzz = Tzzz - (Gamxzz*Txz + Gamyzz*Tyz + Gamzzz*Tzz + Gamxzz*Txz + Gamyzz*Tyz + Gamzzz*Tzz) - chiz*Tzz*inv_chi_factor;

    movx = gupxx*Txxx + gupyy*Txyy + gupzz*Txzz
         + gupxy*Txyx + gupxz*Txzx + gupyz*Txzy
         + gupxy*Txxy + gupxz*Txxz + gupyz*Txyz;
    movy = gupxx*Txyx + gupyy*Tyyy + gupzz*Tyzz
         + gupxy*Tyyx + gupxz*Tyzx + gupyz*Tyzy
         + gupxy*Txyy + gupxz*Txyz + gupyz*Tyyz;
    movz = gupxx*Txzx + gupyy*Tyzy + gupzz*Tzzz
         + gupxy*Tyzx + gupxz*Tzzx + gupyz*Tzzy
         + gupxy*Txzy + gupxz*Txzz + gupyz*Tyzz;
}

} // namespace bssn_constraint_gpu

#endif
