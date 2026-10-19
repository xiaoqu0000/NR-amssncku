#include "macrodef.h"
#include "bssn_constraint.h"
#include "bssn_constraint_gpu.cuh"
#include "gpu_manager.h"
#include "macrodef.fh"

#include <cuda_runtime.h>

using namespace bssn_constraint_gpu;

namespace {

template<bool Shell>
__device__ __forceinline__ void first_derivs(
    const int ex[3], const int pex[3], const double* f,
    double& fx, double& fy, double& fz,
    const double* X, const double* Y, const double* Z,
    const double* crho, const double* sigma, const double* Rcoord,
    const double* drhodx, const double* drhody, const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx, const double* dRdy, const double* dRdz,
    double sym1, double sym2, double sym3, int symmetry,
    int sst, int i, int j, int k)
{
    if constexpr (Shell) {
        first_derivs_shell<ghost_width>(ex, pex, f, fx, fy, fz, crho, sigma, Rcoord,
                                        drhodx, drhody, drhodz,
                                        dsigmadx, dsigmady, dsigmadz,
                                        dRdx, dRdy, dRdz,
                                        sym1, sym2, sym3, symmetry, sst, i, j, k);
    } else {
        first_derivs_cart<ghost_width>(ex, pex, f, fx, fy, fz, X, Y, Z,
                                       sym1, sym2, sym3, symmetry, i, j, k);
    }
}

template<bool Shell>
__global__ void constraint_bssn_kernel(
    const int ex0, const int ex1, const int ex2,
    const int pex0, const int pex1, const int pex2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ crho, const double* __restrict__ sigma, const double* __restrict__ Rcoord,
    const double* __restrict__ drhodx, const double* __restrict__ drhody, const double* __restrict__ drhodz,
    const double* __restrict__ dsigmadx, const double* __restrict__ dsigmady, const double* __restrict__ dsigmadz,
    const double* __restrict__ dRdx, const double* __restrict__ dRdy, const double* __restrict__ dRdz,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Axx, const double* __restrict__ Axy, const double* __restrict__ Axz,
    const double* __restrict__ Ayy, const double* __restrict__ Ayz, const double* __restrict__ Azz,
    const double* __restrict__ Gmx, const double* __restrict__ Gmy, const double* __restrict__ Gmz,
    const double* __restrict__ rho, const double* __restrict__ Sx, const double* __restrict__ Sy, const double* __restrict__ Sz,
    const double* __restrict__ Gamxxx_in, const double* __restrict__ Gamxxy_in, const double* __restrict__ Gamxxz_in,
    const double* __restrict__ Gamxyy_in, const double* __restrict__ Gamxyz_in, const double* __restrict__ Gamxzz_in,
    const double* __restrict__ Gamyxx_in, const double* __restrict__ Gamyxy_in, const double* __restrict__ Gamyxz_in,
    const double* __restrict__ Gamyyy_in, const double* __restrict__ Gamyyz_in, const double* __restrict__ Gamyzz_in,
    const double* __restrict__ Gamzxx_in, const double* __restrict__ Gamzxy_in, const double* __restrict__ Gamzxz_in,
    const double* __restrict__ Gamzyy_in, const double* __restrict__ Gamzyz_in, const double* __restrict__ Gamzzz_in,
    const double* __restrict__ Rxx_in, const double* __restrict__ Rxy_in, const double* __restrict__ Rxz_in,
    const double* __restrict__ Ryy_in, const double* __restrict__ Ryz_in, const double* __restrict__ Rzz_in,
    double* __restrict__ ham_Res,
    double* __restrict__ movx_Res, double* __restrict__ movy_Res, double* __restrict__ movz_Res,
    double* __restrict__ Gmx_Res, double* __restrict__ Gmy_Res, double* __restrict__ Gmz_Res,
    const int symmetry, const int sst)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int n = ex0 * ex1 * ex2;
    if (tid >= n) return;

    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);

    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const size_t id = idx3(i, j, k, pex);

    double gxxx, gxxy, gxxz;
    double gxyx, gxyy, gxyz;
    double gxzx, gxzy, gxzz;
    double gyyx, gyyy, gyyz;
    double gyzx, gyzy, gyzz;
    double gzzx, gzzy, gzzz;

    first_derivs<Shell>(ex, pex, dxx, gxxx, gxxy, gxxz, X, Y, Z, crho, sigma, Rcoord,
                        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
                        dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, gxy, gxyx, gxyy, gxyz, X, Y, Z, crho, sigma, Rcoord,
                        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
                        dRdx, dRdy, dRdz, ANTI, ANTI, SYM, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, gxz, gxzx, gxzy, gxzz, X, Y, Z, crho, sigma, Rcoord,
                        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
                        dRdx, dRdy, dRdz, ANTI, SYM, ANTI, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, dyy, gyyx, gyyy, gyyz, X, Y, Z, crho, sigma, Rcoord,
                        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
                        dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, gyz, gyzx, gyzy, gyzz, X, Y, Z, crho, sigma, Rcoord,
                        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
                        dRdx, dRdy, dRdz, SYM, ANTI, ANTI, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, dzz, gzzx, gzzy, gzzz, X, Y, Z, crho, sigma, Rcoord,
                        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
                        dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    const double l_gxx = dxx[id] + ONE;
    const double l_gxy = gxy[id];
    const double l_gxz = gxz[id];
    const double l_gyy = dyy[id] + ONE;
    const double l_gyz = gyz[id];
    const double l_gzz = dzz[id] + ONE;
    invert_metric(l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz,
                  gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);

    gamma_constraint_residual(
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        gxxx, gxyx, gxzx, gyyx, gyzx, gzzx,
        gxxy, gxyy, gxzy, gyyy, gyzy, gzzy,
        gxxz, gxyz, gxzz, gyyz, gyzz, gzzz,
        Gmx[id], Gmy[id], Gmz[id], Gmx_Res[id], Gmy_Res[id], Gmz_Res[id]);

#if (ABV == 1)
    double adm_gxx, adm_gxy, adm_gxz, adm_gyy, adm_gyz, adm_gzz;
    double Kxx, Kxy, Kxz, Kyy, Kyz, Kzz;
    d_bssn2adm(chi[id] + ONE, trK[id], l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz,
               Axx[id], Axy[id], Axz[id], Ayy[id], Ayz[id], Azz[id],
               adm_gxx, adm_gxy, adm_gxz, adm_gyy, adm_gyz, adm_gzz,
               Kxx, Kxy, Kxz, Kyy, Kyz, Kzz);
    const double cdxx = adm_gxx - ONE;
    const double cdyy = adm_gyy - ONE;
    const double cdzz = adm_gzz - ONE;
    double Rxx = Rxx_in[id], Rxy = Rxy_in[id], Rxz = Rxz_in[id];
    double Ryy = Ryy_in[id], Ryz = Ryz_in[id], Rzz = Rzz_in[id];
    double Gamxxx = Gamxxx_in[id], Gamxxy = Gamxxy_in[id], Gamxxz = Gamxxz_in[id];
    double Gamxyy = Gamxyy_in[id], Gamxyz = Gamxyz_in[id], Gamxzz = Gamxzz_in[id];
    double Gamyxx = Gamyxx_in[id], Gamyxy = Gamyxy_in[id], Gamyxz = Gamyxz_in[id];
    double Gamyyy = Gamyyy_in[id], Gamyyz = Gamyyz_in[id], Gamyzz = Gamyzz_in[id];
    double Gamzxx = Gamzxx_in[id], Gamzxy = Gamzxy_in[id], Gamzxz = Gamzxz_in[id];
    double Gamzyy = Gamzyy_in[id], Gamzyz = Gamzyz_in[id], Gamzzz = Gamzzz_in[id];
    if constexpr (!Shell) {
        d_adm_ricci_gamma_point<ghost_width>(
            ex, pex, i, j, k, id, X, Y, Z,
            dxx, gxy, gxz, dyy, gyz, dzz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, symmetry);
    }
    invert_metric(cdxx + ONE, adm_gxy, adm_gxz, cdyy + ONE, adm_gyz, cdzz + ONE,
                  gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);
    const double trK_adm = gupxx*Kxx + gupyy*Kyy + gupzz*Kzz
                         + TWO*(gupxy*Kxy + gupxz*Kxz + gupyz*Kyz);
    ham_Res[id] = gupxx*Rxx + gupyy*Ryy + gupzz*Rzz + TWO*(gupxy*Rxy + gupxz*Rxz + gupyz*Ryz)
                + trK_adm*trK_adm
                - contract_symmetric_square(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Kxx, Kxy, Kxz, Kyy, Kyz, Kzz)
                - F16 * CUDART_PI * rho[id];
    double Txxx, Txyx, Txzx, Tyyx, Tyzx, Tzzx, Txxy, Txyy, Txzy, Tyyy, Tyzy, Tzzy, Txxz, Txyz, Txzz, Tyyz, Tyzz, Tzzz;
    first_derivs<Shell>(ex, pex, Axx, Txxx, Txxy, Txxz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, Axy, Txyx, Txyy, Txyz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, ANTI, ANTI, SYM, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, Axz, Txzx, Txzy, Txzz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, ANTI, SYM, ANTI, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, Ayy, Tyyx, Tyyy, Tyyz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, Ayz, Tyzx, Tyzy, Tyzz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, ANTI, ANTI, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, Azz, Tzzx, Tzzy, Tzzz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
    double Kx, Ky, Kz;
    first_derivs<Shell>(ex, pex, trK, Kx, Ky, Kz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
    covariant_momentum(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
                       Kxx, Kxy, Kxz, Kyy, Kyz, Kzz,
                       Txxx, Txyx, Txzx, Tyyx, Tyzx, Tzzx, Txxy, Txyy, Txzy, Tyyy, Tyzy, Tzzy, Txxz, Txyz, Txzz, Tyyz, Tyzz, Tzzz,
                       Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
                       Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
                       Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
                       0.0, 0.0, 0.0, 0.0, movx_Res[id], movy_Res[id], movz_Res[id]);
    movx_Res[id] -= Kx + F8 * CUDART_PI * Sx[id];
    movy_Res[id] -= Ky + F8 * CUDART_PI * Sy[id];
    movz_Res[id] -= Kz + F8 * CUDART_PI * Sz[id];
#else
    const double chin1 = chi[id] + ONE;
    ham_Res[id] = chin1 * (gupxx*Rxx_in[id] + gupyy*Ryy_in[id] + gupzz*Rzz_in[id]
                + TWO*(gupxy*Rxy_in[id] + gupxz*Rxz_in[id] + gupyz*Ryz_in[id]))
                + F2O3 * trK[id] * trK[id]
                - contract_symmetric_square(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
                                            Axx[id], Axy[id], Axz[id], Ayy[id], Ayz[id], Azz[id])
                - F16 * CUDART_PI * rho[id];

    double chix, chiy, chiz;
    first_derivs<Shell>(ex, pex, chi, chix, chiy, chiz, X, Y, Z, crho, sigma, Rcoord,
                        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
                        dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);

    double Axxx, Axyx, Axzx, Ayyx, Ayzx, Azzx, Axxy, Axyy, Axzy, Ayyy, Ayzy, Azzy, Axxz, Axyz, Axzz, Ayyz, Ayzz, Azzz;
    first_derivs<Shell>(ex, pex, Axx, Axxx, Axxy, Axxz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, Axy, Axyx, Axyy, Axyz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, ANTI, ANTI, SYM, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, Axz, Axzx, Axzy, Axzz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, ANTI, SYM, ANTI, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, Ayy, Ayyx, Ayyy, Ayyz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, Ayz, Ayzx, Ayzy, Ayzz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, ANTI, ANTI, symmetry, sst, i, j, k);
    first_derivs<Shell>(ex, pex, Azz, Azzx, Azzy, Azzz, X, Y, Z, crho, sigma, Rcoord, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);

    covariant_momentum(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
                       Axx[id], Axy[id], Axz[id], Ayy[id], Ayz[id], Azz[id],
                       Axxx, Axyx, Axzx, Ayyx, Ayzx, Azzx, Axxy, Axyy, Axzy, Ayyy, Ayzy, Azzy, Axxz, Axyz, Axzz, Ayyz, Ayzz, Azzz,
                       Gamxxx_in[id], Gamxxy_in[id], Gamxxz_in[id], Gamxyy_in[id], Gamxyz_in[id], Gamxzz_in[id],
                       Gamyxx_in[id], Gamyxy_in[id], Gamyxz_in[id], Gamyyy_in[id], Gamyyz_in[id], Gamyzz_in[id],
                       Gamzxx_in[id], Gamzxy_in[id], Gamzxz_in[id], Gamzyy_in[id], Gamzyz_in[id], Gamzzz_in[id],
                       chix, chiy, chiz, ONE / chin1,
                       movx_Res[id], movy_Res[id], movz_Res[id]);
    double Kx, Ky, Kz;
    first_derivs<Shell>(ex, pex, trK, Kx, Ky, Kz, X, Y, Z, crho, sigma, Rcoord,
                        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
                        dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
    movx_Res[id] -= F2O3 * Kx + F8 * CUDART_PI * Sx[id];
    movy_Res[id] -= F2O3 * Ky + F8 * CUDART_PI * Sy[id];
    movz_Res[id] -= F2O3 * Kz + F8 * CUDART_PI * Sz[id];
#endif
}

#if (ABV != 1)
void launch_common(
    cudaStream_t stream, bool shell, int ex[3], int pex[3],
    double* X, double* Y, double* Z,
    double* crho, double* sigma, double* Rcoord,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gmx, double* Gmy, double* Gmz,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int sst)
{
    const int n = ex[0] * ex[1] * ex[2];
    const int block = 256;
    const int grid = (n + block - 1) / block;
    if (shell) {
        constraint_bssn_kernel<true><<<grid, block, 0, stream>>>(
            ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], X, Y, Z, crho, sigma, Rcoord,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
            chi, trK, dxx, gxy, gxz, dyy, gyz, dzz, Axx, Axy, Axz, Ayy, Ayz, Azz,
            Gmx, Gmy, Gmz, rho, Sx, Sy, Sz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
            ham_Res, movx_Res, movy_Res, movz_Res, Gmx_Res, Gmy_Res, Gmz_Res, symmetry, sst);
    } else {
        constraint_bssn_kernel<false><<<grid, block, 0, stream>>>(
            ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], X, Y, Z, crho, sigma, Rcoord,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
            chi, trK, dxx, gxy, gxz, dyy, gyz, dzz, Axx, Axy, Axz, Ayy, Ayz, Azz,
            Gmx, Gmy, Gmz, rho, Sx, Sy, Sz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
            ham_Res, movx_Res, movy_Res, movz_Res, Gmx_Res, Gmy_Res, Gmz_Res, symmetry, sst);
    }
}
#endif

} // namespace


namespace {

enum AdmComp {
    ADM_DXX = 0, ADM_GXY, ADM_GXZ, ADM_DYY, ADM_GYZ, ADM_DZZ,
    ADM_KXX, ADM_KXY, ADM_KXZ, ADM_KYY, ADM_KYZ, ADM_KZZ,
    ADM_TRK
};

struct BssnInput {
    const double* chi;
    const double* trK;
    const double* dxx;
    const double* gxy;
    const double* gxz;
    const double* dyy;
    const double* gyz;
    const double* dzz;
    const double* Axx;
    const double* Axy;
    const double* Axz;
    const double* Ayy;
    const double* Ayz;
    const double* Azz;
};

__device__ __forceinline__ void reflect_lower(
    int ex0, int ex1, int ex2, int& i, int& j, int& k,
    double& factor, double sx, double sy, double sz)
{
    if (i < 0) { i = -i - 1; factor *= sx; }
    if (j < 0) { j = -j - 1; factor *= sy; }
    if (k < 0) { k = -k - 1; factor *= sz; }
    i = (i < 0) ? 0 : ((i >= ex0) ? ex0 - 1 : i);
    j = (j < 0) ? 0 : ((j >= ex1) ? ex1 - 1 : j);
    k = (k < 0) ? 0 : ((k >= ex2) ? ex2 - 1 : k);
}

__device__ __forceinline__ double comp_sym(int comp, int axis)
{
    if (comp == ADM_GXY || comp == ADM_KXY) return (axis == 0 || axis == 1) ? ANTI : SYM;
    if (comp == ADM_GXZ || comp == ADM_KXZ) return (axis == 0 || axis == 2) ? ANTI : SYM;
    if (comp == ADM_GYZ || comp == ADM_KYZ) return (axis == 1 || axis == 2) ? ANTI : SYM;
    return SYM;
}

__device__ __forceinline__ double adm_component_at(
    const int ex[3], const int pex[3], const BssnInput& in,
    int i, int j, int k, int comp)
{
    double factor = ONE;
    reflect_lower(ex[0], ex[1], ex[2], i, j, k, factor,
                  comp_sym(comp, 0), comp_sym(comp, 1), comp_sym(comp, 2));

    const size_t id = idx3(i, j, k, pex);
    const double gxx = in.dxx[id] + ONE;
    const double gyy = in.dyy[id] + ONE;
    const double gzz = in.dzz[id] + ONE;
    double adm_gxx, adm_gxy, adm_gxz, adm_gyy, adm_gyz, adm_gzz;
    double Kxx, Kxy, Kxz, Kyy, Kyz, Kzz;
    d_bssn2adm(in.chi[id] + ONE, in.trK[id],
               gxx, in.gxy[id], in.gxz[id], gyy, in.gyz[id], gzz,
               in.Axx[id], in.Axy[id], in.Axz[id], in.Ayy[id], in.Ayz[id], in.Azz[id],
               adm_gxx, adm_gxy, adm_gxz, adm_gyy, adm_gyz, adm_gzz,
               Kxx, Kxy, Kxz, Kyy, Kyz, Kzz);

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    if (comp == ADM_TRK) {
        invert_metric(adm_gxx, adm_gxy, adm_gxz, adm_gyy, adm_gyz, adm_gzz,
                      gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);
        return factor * (gupxx*Kxx + gupyy*Kyy + gupzz*Kzz +
                         TWO*(gupxy*Kxy + gupxz*Kxz + gupyz*Kyz));
    }

    double v = 0.0;
    switch (comp) {
        case ADM_DXX: v = adm_gxx - ONE; break;
        case ADM_GXY: v = adm_gxy; break;
        case ADM_GXZ: v = adm_gxz; break;
        case ADM_DYY: v = adm_gyy - ONE; break;
        case ADM_GYZ: v = adm_gyz; break;
        case ADM_DZZ: v = adm_gzz - ONE; break;
        case ADM_KXX: v = Kxx; break;
        case ADM_KXY: v = Kxy; break;
        case ADM_KXZ: v = Kxz; break;
        case ADM_KYY: v = Kyy; break;
        case ADM_KYZ: v = Kyz; break;
        case ADM_KZZ: v = Kzz; break;
        default: v = 0.0; break;
    }
    return factor * v;
}

template<int Axis>
__device__ __forceinline__ double deriv1_comp(
    const int ex[3], const int pex[3], const BssnInput& in,
    int i, int j, int k, double h, int comp)
{
    const int n = ex[Axis];
    const int q = (Axis == 0) ? i : ((Axis == 1) ? j : k);
    // Match the CPU ADM-constraint boundary contract: derivatives on the
    // outermost high-side diagnostic layer are not evaluated with extrapolated
    // one-sided stencils.  This keeps ABV=1 boundary constraints finite and
    // numerically close to the CPU reference path.
    if (i >= ex[0] - 1 || j >= ex[1] - 1 || k >= ex[2] - 1) return 0.0;
    const auto val = [&](int qq) {
        int ii = i, jj = j, kk = k;
        if constexpr (Axis == 0) ii = qq;
        if constexpr (Axis == 1) jj = qq;
        if constexpr (Axis == 2) kk = qq;
        return adm_component_at(ex, pex, in, ii, jj, kk, comp);
    };

    if (n <= 1) return 0.0;
    if constexpr (ghost_width == 2) {
        if (q >= 1 && q + 1 < n) return (val(q + 1) - val(q - 1)) / (2.0*h);
        return 0.0;
    } else if constexpr (ghost_width == 3) {
        if (i >= 2 && i + 2 < ex[0] &&
            j >= 2 && j + 2 < ex[1] &&
            k >= 2 && k + 2 < ex[2]) {
            return (val(q - 2) - 8.0*val(q - 1) + 8.0*val(q + 1) - val(q + 2)) / (12.0*h);
        }
        if (i >= 1 && i + 1 < ex[0] &&
            j >= 1 && j + 1 < ex[1] &&
            k >= 1 && k + 1 < ex[2]) {
            return (val(q + 1) - val(q - 1)) / (2.0*h);
        }
        return 0.0;
    }
    if (q + 1 >= n) return 0.0;
    if constexpr (ghost_width >= 5) {
        if (q >= 4 && q + 4 < n) {
            return (3.0*val(q - 4) - 32.0*val(q - 3) + 168.0*val(q - 2) - 672.0*val(q - 1)
                    + 672.0*val(q + 1) - 168.0*val(q + 2) + 32.0*val(q + 3) - 3.0*val(q + 4)) / (840.0*h);
        }
    }
    if constexpr (ghost_width >= 4) {
        if (q >= 3 && q + 3 < n) {
            return (-val(q - 3) + 9.0*val(q - 2) - 45.0*val(q - 1)
                    + 45.0*val(q + 1) - 9.0*val(q + 2) + val(q + 3)) / (60.0*h);
        }
    }
    if (q >= 2 && q + 2 < n) return (val(q - 2) - 8.0*val(q - 1) + 8.0*val(q + 1) - val(q + 2)) / (12.0*h);
    if (q >= 1 && q + 1 < n) return (val(q + 1) - val(q - 1)) / (2.0*h);
    return 0.0;
}

template<int Axis>
__device__ __forceinline__ double deriv2_comp(
    const int ex[3], const int pex[3], const BssnInput& in,
    int i, int j, int k, double h, int comp)
{
    const int n = ex[Axis];
    const int q = (Axis == 0) ? i : ((Axis == 1) ? j : k);
    // See deriv1_comp: ABV=1 ADM diagnostic boundaries intentionally mirror
    // the CPU path instead of extrapolating high-side one-sided stencils.
    if (i >= ex[0] - 1 || j >= ex[1] - 1 || k >= ex[2] - 1) return 0.0;
    const auto val = [&](int qq) {
        int ii = i, jj = j, kk = k;
        if constexpr (Axis == 0) ii = qq;
        if constexpr (Axis == 1) jj = qq;
        if constexpr (Axis == 2) kk = qq;
        return adm_component_at(ex, pex, in, ii, jj, kk, comp);
    };

    if (n <= 2) return 0.0;
    if constexpr (ghost_width == 2) {
        if (q >= 1 && q + 1 < n) return (val(q - 1) - 2.0*val(q) + val(q + 1)) / (h*h);
        return 0.0;
    } else if constexpr (ghost_width == 3) {
        if (i >= 2 && i + 2 < ex[0] &&
            j >= 2 && j + 2 < ex[1] &&
            k >= 2 && k + 2 < ex[2]) {
            return (-val(q - 2) + 16.0*val(q - 1) - 30.0*val(q) + 16.0*val(q + 1) - val(q + 2)) / (12.0*h*h);
        }
        if (i >= 1 && i + 1 < ex[0] &&
            j >= 1 && j + 1 < ex[1] &&
            k >= 1 && k + 1 < ex[2]) {
            return (val(q - 1) - 2.0*val(q) + val(q + 1)) / (h*h);
        }
        return 0.0;
    }
    if (q + 1 >= n) return 0.0;
    if constexpr (ghost_width >= 5) {
        if (q >= 4 && q + 4 < n) {
            return (-9.0*val(q - 4) + 128.0*val(q - 3) - 1008.0*val(q - 2) + 8064.0*val(q - 1)
                    - 14350.0*val(q) + 8064.0*val(q + 1) - 1008.0*val(q + 2) + 128.0*val(q + 3)
                    - 9.0*val(q + 4)) / (5040.0*h*h);
        }
    }
    if constexpr (ghost_width >= 4) {
        if (q >= 3 && q + 3 < n) {
            return (2.0*val(q - 3) - 27.0*val(q - 2) + 270.0*val(q - 1) - 490.0*val(q)
                    + 270.0*val(q + 1) - 27.0*val(q + 2) + 2.0*val(q + 3)) / (180.0*h*h);
        }
    }
    if (q >= 2 && q + 2 < n) {
        return (-val(q - 2) + 16.0*val(q - 1) - 30.0*val(q) + 16.0*val(q + 1) - val(q + 2)) / (12.0*h*h);
    }
    if (q >= 1 && q + 1 < n) return (val(q - 1) - 2.0*val(q) + val(q + 1)) / (h*h);
    return 0.0;
}

template<int AxisA, int AxisB>
__device__ __forceinline__ double cross2_comp(
    const int ex[3], const int pex[3], const BssnInput& in,
    int i, int j, int k, double ha, double hb, int comp)
{
    const auto val = [&](int qa, int qb) {
        int ii = i, jj = j, kk = k;
        if constexpr (AxisA == 0) ii = qa;
        if constexpr (AxisA == 1) jj = qa;
        if constexpr (AxisA == 2) kk = qa;
        if constexpr (AxisB == 0) ii = qb;
        if constexpr (AxisB == 1) jj = qb;
        if constexpr (AxisB == 2) kk = qb;
        return adm_component_at(ex, pex, in, ii, jj, kk, comp);
    };
    const int qa = (AxisA == 0) ? i : ((AxisA == 1) ? j : k);
    const int qb = (AxisB == 0) ? i : ((AxisB == 1) ? j : k);
    const int na = ex[AxisA], nb = ex[AxisB];
    if (i >= ex[0] - 1 || j >= ex[1] - 1 || k >= ex[2] - 1) return 0.0;
    if (na <= 2 || nb <= 2) return 0.0;
    if constexpr (ghost_width >= 5) {
        if (qa >= 4 && qa + 4 < na && qb >= 4 && qb + 4 < nb) {
            double sum = 0.0;
            const int off[8] = {-4, -3, -2, -1, 1, 2, 3, 4};
            const double c[8] = {3.0, -32.0, 168.0, -672.0, 672.0, -168.0, 32.0, -3.0};
            #pragma unroll
            for (int b = 0; b < 8; ++b) {
                #pragma unroll
                for (int a = 0; a < 8; ++a) {
                    sum += c[a] * c[b] * val(qa + off[a], qb + off[b]);
                }
            }
            return sum / (705600.0*ha*hb);
        }
    }
    if constexpr (ghost_width >= 4) {
        if (qa >= 3 && qa + 3 < na && qb >= 3 && qb + 3 < nb) {
            double sum = 0.0;
            const int off[6] = {-3, -2, -1, 1, 2, 3};
            const double c[6] = {-1.0, 9.0, -45.0, 45.0, -9.0, 1.0};
            #pragma unroll
            for (int b = 0; b < 6; ++b) {
                #pragma unroll
                for (int a = 0; a < 6; ++a) {
                    sum += c[a] * c[b] * val(qa + off[a], qb + off[b]);
                }
            }
            return sum / (3600.0*ha*hb);
        }
    }
    if constexpr (ghost_width == 3) {
        if (i >= 2 && i + 2 < ex[0] &&
            j >= 2 && j + 2 < ex[1] &&
            k >= 2 && k + 2 < ex[2]) {
            double sum = 0.0;
            const int off[4] = {-2, -1, 1, 2};
            const double c[4] = {1.0, -8.0, 8.0, -1.0};
            #pragma unroll
            for (int b = 0; b < 4; ++b) {
                #pragma unroll
                for (int a = 0; a < 4; ++a) {
                    sum += c[a] * c[b] * val(qa + off[a], qb + off[b]);
                }
            }
            return sum / (144.0*ha*hb);
        }
        if (!(i >= 1 && i + 1 < ex[0] &&
              j >= 1 && j + 1 < ex[1] &&
              k >= 1 && k + 1 < ex[2])) return 0.0;
    }
    if (qa <= 0 || qa + 1 >= na || qb <= 0 || qb + 1 >= nb) return 0.0;
    return (val(qa + 1, qb + 1) - val(qa + 1, qb - 1) -
            val(qa - 1, qb + 1) + val(qa - 1, qb - 1)) / (4.0*ha*hb);
}

__device__ __forceinline__ void adm_metric_derivs_cart(
    const int ex[3], const int pex[3], const BssnInput& in,
    const double* X, const double* Y, const double* Z,
    int i, int j, int k, int comp,
    double& fx, double& fy, double& fz,
    double& fxx, double& fxy, double& fxz, double& fyy, double& fyz, double& fzz)
{
    const double hx = X[1] - X[0], hy = Y[1] - Y[0], hz = Z[1] - Z[0];
    fx = deriv1_comp<0>(ex, pex, in, i, j, k, hx, comp);
    fy = deriv1_comp<1>(ex, pex, in, i, j, k, hy, comp);
    fz = deriv1_comp<2>(ex, pex, in, i, j, k, hz, comp);
    fxx = deriv2_comp<0>(ex, pex, in, i, j, k, hx, comp);
    fyy = deriv2_comp<1>(ex, pex, in, i, j, k, hy, comp);
    fzz = deriv2_comp<2>(ex, pex, in, i, j, k, hz, comp);
    fxy = cross2_comp<0, 1>(ex, pex, in, i, j, k, hx, hy, comp);
    fxz = cross2_comp<0, 2>(ex, pex, in, i, j, k, hx, hz, comp);
    fyz = cross2_comp<1, 2>(ex, pex, in, i, j, k, hy, hz, comp);
}

__device__ __forceinline__ void adm_metric_derivs_shell(
    const int ex[3], const int pex[3], const BssnInput& in,
    const double* crho, const double* sigma, const double* R,
    const double* drhodx, const double* drhody, const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx, const double* dRdy, const double* dRdz,
    const double* drhodxx, const double* drhodxy, const double* drhodxz,
    const double* drhodyy, const double* drhodyz, const double* drhodzz,
    const double* dsigmadxx, const double* dsigmadxy, const double* dsigmadxz,
    const double* dsigmadyy, const double* dsigmadyz, const double* dsigmadzz,
    const double* dRdxx, const double* dRdxy, const double* dRdxz,
    const double* dRdyy, const double* dRdyz, const double* dRdzz,
    int i, int j, int k, int comp,
    double& fx, double& fy, double& fz,
    double& fxx, double& fxy, double& fxz, double& fyy, double& fyz, double& fzz)
{
    const double hr = crho[1] - crho[0], hs = sigma[1] - sigma[0], hR = R[1] - R[0];
    const double fr = deriv1_comp<0>(ex, pex, in, i, j, k, hr, comp);
    const double fs = deriv1_comp<1>(ex, pex, in, i, j, k, hs, comp);
    const double fR = deriv1_comp<2>(ex, pex, in, i, j, k, hR, comp);
    const double frr = deriv2_comp<0>(ex, pex, in, i, j, k, hr, comp);
    const double fss = deriv2_comp<1>(ex, pex, in, i, j, k, hs, comp);
    const double fRR = deriv2_comp<2>(ex, pex, in, i, j, k, hR, comp);
    const double frs = cross2_comp<0, 1>(ex, pex, in, i, j, k, hr, hs, comp);
    const double frR = cross2_comp<0, 2>(ex, pex, in, i, j, k, hr, hR, comp);
    const double fsR = cross2_comp<1, 2>(ex, pex, in, i, j, k, hs, hR, comp);
    const size_t id = idx3(i, j, k, pex);

    const double rx = drhodx[id], ry = drhody[id], rz = drhodz[id];
    const double sx = dsigmadx[id], sy = dsigmady[id], sz = dsigmadz[id];
    const double Rx = dRdx[id], Ry = dRdy[id], Rz = dRdz[id];

    fx = fr*rx + fs*sx + fR*Rx;
    fy = fr*ry + fs*sy + fR*Ry;
    fz = fr*rz + fs*sz + fR*Rz;

    fxx = frr*rx*rx + fss*sx*sx + fRR*Rx*Rx + 2.0*(frs*rx*sx + frR*rx*Rx + fsR*sx*Rx)
        + fr*drhodxx[id] + fs*dsigmadxx[id] + fR*dRdxx[id];
    fxy = frr*rx*ry + fss*sx*sy + fRR*Rx*Ry + frs*(rx*sy + ry*sx) + frR*(rx*Ry + ry*Rx) + fsR*(sx*Ry + sy*Rx)
        + fr*drhodxy[id] + fs*dsigmadxy[id] + fR*dRdxy[id];
    fxz = frr*rx*rz + fss*sx*sz + fRR*Rx*Rz + frs*(rx*sz + rz*sx) + frR*(rx*Rz + rz*Rx) + fsR*(sx*Rz + sz*Rx)
        + fr*drhodxz[id] + fs*dsigmadxz[id] + fR*dRdxz[id];
    fyy = frr*ry*ry + fss*sy*sy + fRR*Ry*Ry + 2.0*(frs*ry*sy + frR*ry*Ry + fsR*sy*Ry)
        + fr*drhodyy[id] + fs*dsigmadyy[id] + fR*dRdyy[id];
    fyz = frr*ry*rz + fss*sy*sz + fRR*Ry*Rz + frs*(ry*sz + rz*sy) + frR*(ry*Rz + rz*Ry) + fsR*(sy*Rz + sz*Ry)
        + fr*drhodyz[id] + fs*dsigmadyz[id] + fR*dRdyz[id];
    fzz = frr*rz*rz + fss*sz*sz + fRR*Rz*Rz + 2.0*(frs*rz*sz + frR*rz*Rz + fsR*sz*Rz)
        + fr*drhodzz[id] + fs*dsigmadzz[id] + fR*dRdzz[id];
}

__device__ __forceinline__ void adm_metric_first_derivs_shell(
    const int ex[3], const int pex[3], const BssnInput& in,
    const double* crho, const double* sigma, const double* R,
    const double* drhodx, const double* drhody, const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx, const double* dRdy, const double* dRdz,
    int i, int j, int k, int comp,
    double& fx, double& fy, double& fz)
{
    const double hr = crho[1] - crho[0], hs = sigma[1] - sigma[0], hR = R[1] - R[0];
    const double fr = deriv1_comp<0>(ex, pex, in, i, j, k, hr, comp);
    const double fs = deriv1_comp<1>(ex, pex, in, i, j, k, hs, comp);
    const double fR = deriv1_comp<2>(ex, pex, in, i, j, k, hR, comp);
    const size_t id = idx3(i, j, k, pex);

    fx = fr*drhodx[id] + fs*dsigmadx[id] + fR*dRdx[id];
    fy = fr*drhody[id] + fs*dsigmady[id] + fR*dRdy[id];
    fz = fr*drhodz[id] + fs*dsigmadz[id] + fR*dRdz[id];
}

template<bool Shell>
__device__ __forceinline__ void adm_constraint_point(
    const int ex[3], const int pex[3], const BssnInput& in,
    const double* X, const double* Y, const double* Z,
    const double* crho, const double* sigma, const double* R,
    const double* drhodx, const double* drhody, const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx, const double* dRdy, const double* dRdz,
    const double* drhodxx, const double* drhodxy, const double* drhodxz,
    const double* drhodyy, const double* drhodyz, const double* drhodzz,
    const double* dsigmadxx, const double* dsigmadxy, const double* dsigmadxz,
    const double* dsigmadyy, const double* dsigmadyz, const double* dsigmadzz,
    const double* dRdxx, const double* dRdxy, const double* dRdxz,
    const double* dRdyy, const double* dRdyz, const double* dRdzz,
    const double* rho, const double* Sx, const double* Sy, const double* Sz,
    int i, int j, int k,
    double& ham, double& mx, double& my, double& mz)
{
    const size_t id = idx3(i, j, k, pex);
    const double dxx = adm_component_at(ex, pex, in, i, j, k, ADM_DXX);
    const double gxy = adm_component_at(ex, pex, in, i, j, k, ADM_GXY);
    const double gxz = adm_component_at(ex, pex, in, i, j, k, ADM_GXZ);
    const double dyy = adm_component_at(ex, pex, in, i, j, k, ADM_DYY);
    const double gyz = adm_component_at(ex, pex, in, i, j, k, ADM_GYZ);
    const double dzz = adm_component_at(ex, pex, in, i, j, k, ADM_DZZ);
    const double Kxx = adm_component_at(ex, pex, in, i, j, k, ADM_KXX);
    const double Kxy = adm_component_at(ex, pex, in, i, j, k, ADM_KXY);
    const double Kxz = adm_component_at(ex, pex, in, i, j, k, ADM_KXZ);
    const double Kyy = adm_component_at(ex, pex, in, i, j, k, ADM_KYY);
    const double Kyz = adm_component_at(ex, pex, in, i, j, k, ADM_KYZ);
    const double Kzz = adm_component_at(ex, pex, in, i, j, k, ADM_KZZ);

    double md[6][9];
    for (int c = 0; c < 6; ++c) {
        if constexpr (Shell) {
            adm_metric_derivs_shell(ex, pex, in, crho, sigma, R,
                drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
                drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
                dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
                dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
                i, j, k, c, md[c][0], md[c][1], md[c][2], md[c][3], md[c][4], md[c][5], md[c][6], md[c][7], md[c][8]);
        } else {
            adm_metric_derivs_cart(ex, pex, in, X, Y, Z, i, j, k, c,
                md[c][0], md[c][1], md[c][2], md[c][3], md[c][4], md[c][5], md[c][6], md[c][7], md[c][8]);
        }
    }

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    invert_metric(dxx + ONE, gxy, gxz, dyy + ONE, gyz, dzz + ONE,
                  gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);

    double aGxxx, aGxxy, aGxxz, aGxyy, aGxyz, aGxzz;
    double aGyxx, aGyxy, aGyxz, aGyyy, aGyyz, aGyzz;
    double aGzxx, aGzxy, aGzxz, aGzyy, aGzyz, aGzzz;
    d_kind1_connection(
        md[0][0], md[1][0], md[2][0], md[3][0], md[4][0], md[5][0],
        md[0][1], md[1][1], md[2][1], md[3][1], md[4][1], md[5][1],
        md[0][2], md[1][2], md[2][2], md[3][2], md[4][2], md[5][2],
        aGxxx, aGxxy, aGxxz, aGxyy, aGxyz, aGxzz,
        aGyxx, aGyxy, aGyxz, aGyyy, aGyyz, aGyzz,
        aGzxx, aGzxy, aGzxz, aGzyy, aGzyz, aGzzz);

    double Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz;
    double Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz;
    double Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz;
    d_kind2_connection(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        aGxxx, aGxxy, aGxxz, aGxyy, aGxyz, aGxzz,
        aGyxx, aGyxy, aGyxz, aGyyy, aGyyz, aGyzz,
        aGzxx, aGzxy, aGzxz, aGzyy, aGzyz, aGzzz,
        Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz,
        Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz,
        Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz);

    double Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz;
    d_adm_riemann(
        md[0][3], md[0][4], md[0][5], md[0][6], md[0][7], md[0][8],
        md[1][3], md[1][4], md[1][5], md[1][6], md[1][7], md[1][8],
        md[2][3], md[2][4], md[2][5], md[2][6], md[2][7], md[2][8],
        md[3][3], md[3][4], md[3][5], md[3][6], md[3][7], md[3][8],
        md[4][3], md[4][4], md[4][5], md[4][6], md[4][7], md[4][8],
        md[5][3], md[5][4], md[5][5], md[5][6], md[5][7], md[5][8],
        Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz,
        Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz,
        Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz,
        aGxxx, aGxxy, aGxxz, aGxyy, aGxyz, aGxzz,
        aGyxx, aGyxy, aGyxz, aGyyy, aGyyz, aGyzz,
        aGzxx, aGzxy, aGzxz, aGzyy, aGzyz, aGzzz,
        Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz);

    double Rxx, Rxy, Rxz, Ryy, Ryz, Rzz;
    d_adm_ricci(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
                Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz,
                Rxx, Rxy, Rxz, Ryy, Ryz, Rzz);

    const double trK = gupxx*Kxx + gupyy*Kyy + gupzz*Kzz + TWO*(gupxy*Kxy + gupxz*Kxz + gupyz*Kyz);
    ham = gupxx*Rxx + gupyy*Ryy + gupzz*Rzz + TWO*(gupxy*Rxy + gupxz*Rxz + gupyz*Ryz)
        + trK*trK - contract_symmetric_square(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, Kxx, Kxy, Kxz, Kyy, Kyz, Kzz)
        - F16*CUDART_PI*rho[id];

    double kdx[6][3];
    for (int c = 0; c < 6; ++c) {
        const int kc = ADM_KXX + c;
        if constexpr (Shell) {
            adm_metric_first_derivs_shell(ex, pex, in, crho, sigma, R,
                drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
                i, j, k, kc, kdx[c][0], kdx[c][1], kdx[c][2]);
        } else {
            kdx[c][0] = deriv1_comp<0>(ex, pex, in, i, j, k, X[1] - X[0], kc);
            kdx[c][1] = deriv1_comp<1>(ex, pex, in, i, j, k, Y[1] - Y[0], kc);
            kdx[c][2] = deriv1_comp<2>(ex, pex, in, i, j, k, Z[1] - Z[0], kc);
        }
    }
    double trKx, trKy, trKz;
    if constexpr (Shell) {
        adm_metric_first_derivs_shell(ex, pex, in, crho, sigma, R,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
            i, j, k, ADM_TRK, trKx, trKy, trKz);
    } else {
        trKx = deriv1_comp<0>(ex, pex, in, i, j, k, X[1] - X[0], ADM_TRK);
        trKy = deriv1_comp<1>(ex, pex, in, i, j, k, Y[1] - Y[0], ADM_TRK);
        trKz = deriv1_comp<2>(ex, pex, in, i, j, k, Z[1] - Z[0], ADM_TRK);
    }

    covariant_momentum(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        Kxx, Kxy, Kxz, Kyy, Kyz, Kzz,
        kdx[0][0], kdx[1][0], kdx[2][0], kdx[3][0], kdx[4][0], kdx[5][0],
        kdx[0][1], kdx[1][1], kdx[2][1], kdx[3][1], kdx[4][1], kdx[5][1],
        kdx[0][2], kdx[1][2], kdx[2][2], kdx[3][2], kdx[4][2], kdx[5][2],
        Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz,
        Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz,
        Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz,
        0.0, 0.0, 0.0, 0.0, mx, my, mz);
    mx -= trKx + F8*CUDART_PI*Sx[id];
    my -= trKy + F8*CUDART_PI*Sy[id];
    mz -= trKz + F8*CUDART_PI*Sz[id];
}

template<bool Shell>
__global__ void constraint_bssn_abv1_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* X, const double* Y, const double* Z,
    const double* crho, const double* sigma, const double* R,
    const double* drhodx, const double* drhody, const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx, const double* dRdy, const double* dRdz,
    const double* drhodxx, const double* drhodxy, const double* drhodxz,
    const double* drhodyy, const double* drhodyz, const double* drhodzz,
    const double* dsigmadxx, const double* dsigmadxy, const double* dsigmadxz,
    const double* dsigmadyy, const double* dsigmadyz, const double* dsigmadzz,
    const double* dRdxx, const double* dRdxy, const double* dRdxz,
    const double* dRdyy, const double* dRdyz, const double* dRdzz,
    BssnInput in,
    const double* Gmx, const double* Gmy, const double* Gmz,
    const double* rho, const double* Sx, const double* Sy, const double* Sz,
    double* ham, double* mx, double* my, double* mz,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int sst)
{
    const int n = ex0 * ex1 * ex2;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const size_t id = idx3(i, j, k, pex);

    double gxxx, gxxy, gxxz, gxyx, gxyy, gxyz, gxzx, gxzy, gxzz;
    double gyyx, gyyy, gyyz, gyzx, gyzy, gyzz, gzzx, gzzy, gzzz;
    if constexpr (Shell) {
        first_derivs_shell<ghost_width>(ex, pex, in.dxx, gxxx, gxxy, gxxz, crho, sigma, R,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
        first_derivs_shell<ghost_width>(ex, pex, in.gxy, gxyx, gxyy, gxyz, crho, sigma, R,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, ANTI, ANTI, SYM, symmetry, sst, i, j, k);
        first_derivs_shell<ghost_width>(ex, pex, in.gxz, gxzx, gxzy, gxzz, crho, sigma, R,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, ANTI, SYM, ANTI, symmetry, sst, i, j, k);
        first_derivs_shell<ghost_width>(ex, pex, in.dyy, gyyx, gyyy, gyyz, crho, sigma, R,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
        first_derivs_shell<ghost_width>(ex, pex, in.gyz, gyzx, gyzy, gyzz, crho, sigma, R,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, ANTI, ANTI, symmetry, sst, i, j, k);
        first_derivs_shell<ghost_width>(ex, pex, in.dzz, gzzx, gzzy, gzzz, crho, sigma, R,
            drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, SYM, SYM, SYM, symmetry, sst, i, j, k);
    } else {
        first_derivs_cart<ghost_width>(ex, pex, in.dxx, gxxx, gxxy, gxxz, X, Y, Z, SYM, SYM, SYM, symmetry, i, j, k);
        first_derivs_cart<ghost_width>(ex, pex, in.gxy, gxyx, gxyy, gxyz, X, Y, Z, ANTI, ANTI, SYM, symmetry, i, j, k);
        first_derivs_cart<ghost_width>(ex, pex, in.gxz, gxzx, gxzy, gxzz, X, Y, Z, ANTI, SYM, ANTI, symmetry, i, j, k);
        first_derivs_cart<ghost_width>(ex, pex, in.dyy, gyyx, gyyy, gyyz, X, Y, Z, SYM, SYM, SYM, symmetry, i, j, k);
        first_derivs_cart<ghost_width>(ex, pex, in.gyz, gyzx, gyzy, gyzz, X, Y, Z, SYM, ANTI, ANTI, symmetry, i, j, k);
        first_derivs_cart<ghost_width>(ex, pex, in.dzz, gzzx, gzzy, gzzz, X, Y, Z, SYM, SYM, SYM, symmetry, i, j, k);
    }

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    invert_metric(in.dxx[id] + ONE, in.gxy[id], in.gxz[id], in.dyy[id] + ONE, in.gyz[id], in.dzz[id] + ONE,
                  gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);
    gamma_constraint_residual(
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        gxxx, gxyx, gxzx, gyyx, gyzx, gzzx,
        gxxy, gxyy, gxzy, gyyy, gyzy, gzzy,
        gxxz, gxyz, gxzz, gyyz, gyzz, gzzz,
        Gmx ? Gmx[id] : 0.0, Gmy ? Gmy[id] : 0.0, Gmz ? Gmz[id] : 0.0,
        Gmx_Res[id], Gmy_Res[id], Gmz_Res[id]);

    adm_constraint_point<Shell>(ex, pex, in, X, Y, Z, crho, sigma, R,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        rho, Sx, Sy, Sz, i, j, k, ham[id], mx[id], my[id], mz[id]);
}

} // namespace

__global__ void bssn_abv1_adm_arrays_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    BssnInput in,
    double* adm_dxx, double* adm_gxy, double* adm_gxz,
    double* adm_dyy, double* adm_gyz, double* adm_dzz,
    double* Kxx, double* Kxy, double* Kxz,
    double* Kyy, double* Kyz, double* Kzz);

__global__ void bssn_abv1_hamiltonian_from_adm_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* X, const double* Y, const double* Z,
    const double* adm_dxx, const double* adm_gxy, const double* adm_gxz,
    const double* adm_dyy, const double* adm_gyz, const double* adm_dzz,
    const double* Kxx, const double* Kxy, const double* Kxz,
    const double* Kyy, const double* Kyz, const double* Kzz,
    const double* rho, double* ham, int symmetry);

extern "C" void gpu_constraint_bssn_abv1_launch(
    cudaStream_t stream, int ex[3], int pex[3],
    double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gmx, double* Gmy, double* Gmz,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* ham, double* mx, double* my, double* mz,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry)
{
    const int n = ex[0] * ex[1] * ex[2];
    const int block = 128;
    const int grid = (n + block - 1) / block;
    const size_t num_elements = static_cast<size_t>(pex[0]) * pex[1] * pex[2];
    GPUManager& gpu_manager = GPUManager::getInstance();
    BssnInput in{chi, trK, dxx, gxy, gxz, dyy, gyz, dzz, Axx, Axy, Axz, Ayy, Ayz, Azz};
    constraint_bssn_abv1_kernel<false><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        X, Y, Z, X, Y, Z,
        nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
        in, Gmx, Gmy, Gmz, rho, Sx, Sy, Sz,
        ham, mx, my, mz, Gmx_Res, Gmy_Res, Gmz_Res, symmetry, 0);

    constexpr int workspace_fields = 12;
    double* workspace = gpu_manager.allocate_device_memory<double>(workspace_fields * num_elements, stream);
    double* adm_dxx = workspace + 0 * num_elements;
    double* adm_gxy = workspace + 1 * num_elements;
    double* adm_gxz = workspace + 2 * num_elements;
    double* adm_dyy = workspace + 3 * num_elements;
    double* adm_gyz = workspace + 4 * num_elements;
    double* adm_dzz = workspace + 5 * num_elements;
    double* Kxx_tmp = workspace + 6 * num_elements;
    double* Kxy_tmp = workspace + 7 * num_elements;
    double* Kxz_tmp = workspace + 8 * num_elements;
    double* Kyy_tmp = workspace + 9 * num_elements;
    double* Kyz_tmp = workspace + 10 * num_elements;
    double* Kzz_tmp = workspace + 11 * num_elements;
    bssn_abv1_adm_arrays_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], in,
        adm_dxx, adm_gxy, adm_gxz, adm_dyy, adm_gyz, adm_dzz,
        Kxx_tmp, Kxy_tmp, Kxz_tmp, Kyy_tmp, Kyz_tmp, Kzz_tmp);
    bssn_abv1_hamiltonian_from_adm_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        X, Y, Z, adm_dxx, adm_gxy, adm_gxz, adm_dyy, adm_gyz, adm_dzz,
        Kxx_tmp, Kxy_tmp, Kxz_tmp, Kyy_tmp, Kyz_tmp, Kzz_tmp, rho, ham, symmetry);
    gpu_manager.free_device_memory(workspace, workspace_fields * num_elements, stream);
}

extern "C" void gpu_constraint_bssn_abv1_ss_launch(
    cudaStream_t stream, int ex[3], int pex[3],
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
    double* Gmx, double* Gmy, double* Gmz,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* ham, double* mx, double* my, double* mz,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int sst)
{
    const int n = ex[0] * ex[1] * ex[2];
    const int block = 128;
    const int grid = (n + block - 1) / block;
    BssnInput in{chi, trK, dxx, gxy, gxz, dyy, gyz, dzz, Axx, Axy, Axz, Ayy, Ayz, Azz};
    constraint_bssn_abv1_kernel<true><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        crho, sigma, R, crho, sigma, R,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        in, Gmx, Gmy, Gmz, rho, Sx, Sy, Sz,
        ham, mx, my, mz, Gmx_Res, Gmy_Res, Gmz_Res, symmetry, sst);
}

__global__ void bssn_abv1_adm_arrays_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    BssnInput in,
    double* adm_dxx, double* adm_gxy, double* adm_gxz,
    double* adm_dyy, double* adm_gyz, double* adm_dzz,
    double* Kxx, double* Kxy, double* Kxz,
    double* Kyy, double* Kyz, double* Kzz)
{
    const int n = ex0 * ex1 * ex2;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    const int pex[3] = {pex0, pex1, 0};
    const size_t id = idx3(i, j, k, pex);

    double agxx, agxy, agxz, agyy, agyz, agzz;
    double kxx, kxy, kxz, kyy, kyz, kzz;
    d_bssn2adm(in.chi[id] + ONE, in.trK[id],
               in.dxx[id] + ONE, in.gxy[id], in.gxz[id],
               in.dyy[id] + ONE, in.gyz[id], in.dzz[id] + ONE,
               in.Axx[id], in.Axy[id], in.Axz[id],
               in.Ayy[id], in.Ayz[id], in.Azz[id],
               agxx, agxy, agxz, agyy, agyz, agzz,
               kxx, kxy, kxz, kyy, kyz, kzz);
    adm_dxx[id] = agxx - ONE;
    adm_gxy[id] = agxy;
    adm_gxz[id] = agxz;
    adm_dyy[id] = agyy - ONE;
    adm_gyz[id] = agyz;
    adm_dzz[id] = agzz - ONE;
    Kxx[id] = kxx;
    Kxy[id] = kxy;
    Kxz[id] = kxz;
    Kyy[id] = kyy;
    Kyz[id] = kyz;
    Kzz[id] = kzz;
}

__global__ void bssn_abv1_hamiltonian_from_adm_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* X, const double* Y, const double* Z,
    const double* adm_dxx, const double* adm_gxy, const double* adm_gxz,
    const double* adm_dyy, const double* adm_gyz, const double* adm_dzz,
    const double* Kxx, const double* Kxy, const double* Kxz,
    const double* Kyy, const double* Kyz, const double* Kzz,
    const double* rho, double* ham, int symmetry)
{
    const int n = ex0 * ex1 * ex2;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const size_t id = idx3(i, j, k, pex);

    double Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz;
    double Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz;
    double Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz;
    double Rxx, Rxy, Rxz, Ryy, Ryz, Rzz;
    d_adm_ricci_gamma_point<ghost_width>(
        ex, pex, i, j, k, id, X, Y, Z,
        adm_dxx, adm_gxy, adm_gxz, adm_dyy, adm_gyz, adm_dzz,
        Gxxx, Gxxy, Gxxz, Gxyy, Gxyz, Gxzz,
        Gyxx, Gyxy, Gyxz, Gyyy, Gyyz, Gyzz,
        Gzxx, Gzxy, Gzxz, Gzyy, Gzyz, Gzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, symmetry);

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    invert_metric(adm_dxx[id] + ONE, adm_gxy[id], adm_gxz[id],
                  adm_dyy[id] + ONE, adm_gyz[id], adm_dzz[id] + ONE,
                  gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);
    const double trK = gupxx*Kxx[id] + gupyy*Kyy[id] + gupzz*Kzz[id]
                     + TWO*(gupxy*Kxy[id] + gupxz*Kxz[id] + gupyz*Kyz[id]);
    ham[id] = gupxx*Rxx + gupyy*Ryy + gupzz*Rzz + TWO*(gupxy*Rxy + gupxz*Rxz + gupyz*Ryz)
            + trK*trK
            - contract_symmetric_square(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
                                        Kxx[id], Kxy[id], Kxz[id], Kyy[id], Kyz[id], Kzz[id])
            - F16*CUDART_PI*rho[id];
}
void gpu_constraint_bssn_launch(
    cudaStream_t stream, int ex[3], int pex[3], double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gmx, double* Gmy, double* Gmz,
    double* Lap, double* Sfx, double* Sfy, double* Sfz,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry)
{
    (void)Lap; (void)Sfx; (void)Sfy; (void)Sfz;
#if (ABV == 1)
    (void)Gamxxx; (void)Gamxxy; (void)Gamxxz; (void)Gamxyy; (void)Gamxyz; (void)Gamxzz;
    (void)Gamyxx; (void)Gamyxy; (void)Gamyxz; (void)Gamyyy; (void)Gamyyz; (void)Gamyzz;
    (void)Gamzxx; (void)Gamzxy; (void)Gamzxz; (void)Gamzyy; (void)Gamzyz; (void)Gamzzz;
    (void)Rxx; (void)Rxy; (void)Rxz; (void)Ryy; (void)Ryz; (void)Rzz;
    gpu_constraint_bssn_abv1_launch(stream, ex, pex, X, Y, Z,
                                    chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
                                    Axx, Axy, Axz, Ayy, Ayz, Azz,
                                    Gmx, Gmy, Gmz, rho, Sx, Sy, Sz,
                                    ham_Res, movx_Res, movy_Res, movz_Res,
                                    Gmx_Res, Gmy_Res, Gmz_Res, symmetry);
#else
    launch_common(stream, false, ex, pex, X, Y, Z, X, Y, Z,
                  nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                  chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
                  Axx, Axy, Axz, Ayy, Ayz, Azz, Gmx, Gmy, Gmz,
                  rho, Sx, Sy, Sz,
                  Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
                  Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
                  Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
                  Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
                  ham_Res, movx_Res, movy_Res, movz_Res, Gmx_Res, Gmy_Res, Gmz_Res,
                  symmetry, 0);
#endif
}

void gpu_constraint_bssn_ss_launch(
    cudaStream_t stream, int ex[3], int pex[3],
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
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gmx, double* Gmy, double* Gmz,
    double* Lap, double* Sfx, double* Sfy, double* Sfz,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, int sst)
{
#if (ABV == 1)
    (void)X; (void)Y; (void)Z; (void)Lap; (void)Sfx; (void)Sfy; (void)Sfz;
    (void)Gamxxx; (void)Gamxxy; (void)Gamxxz; (void)Gamxyy; (void)Gamxyz; (void)Gamxzz;
    (void)Gamyxx; (void)Gamyxy; (void)Gamyxz; (void)Gamyyy; (void)Gamyyz; (void)Gamyzz;
    (void)Gamzxx; (void)Gamzxy; (void)Gamzxz; (void)Gamzyy; (void)Gamzyz; (void)Gamzzz;
    (void)Rxx; (void)Rxy; (void)Rxz; (void)Ryy; (void)Ryz; (void)Rzz; (void)lev;
    gpu_constraint_bssn_abv1_ss_launch(stream, ex, pex, crho, sigma, Rcoord,
                                       drhodx, drhody, drhodz,
                                       dsigmadx, dsigmady, dsigmadz,
                                       dRdx, dRdy, dRdz,
                                       drhodxx, drhodxy, drhodxz,
                                       drhodyy, drhodyz, drhodzz,
                                       dsigmadxx, dsigmadxy, dsigmadxz,
                                       dsigmadyy, dsigmadyz, dsigmadzz,
                                       dRdxx, dRdxy, dRdxz,
                                       dRdyy, dRdyz, dRdzz,
                                       chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
                                       Axx, Axy, Axz, Ayy, Ayz, Azz,
                                       Gmx, Gmy, Gmz, rho, Sx, Sy, Sz,
                                       ham_Res, movx_Res, movy_Res, movz_Res,
                                       Gmx_Res, Gmy_Res, Gmz_Res,
                                       symmetry, sst);
#else
    (void)X; (void)Y; (void)Z; (void)Lap; (void)Sfx; (void)Sfy; (void)Sfz;
    (void)drhodxx; (void)drhodxy; (void)drhodxz; (void)drhodyy; (void)drhodyz; (void)drhodzz;
    (void)dsigmadxx; (void)dsigmadxy; (void)dsigmadxz; (void)dsigmadyy; (void)dsigmadyz; (void)dsigmadzz;
    (void)dRdxx; (void)dRdxy; (void)dRdxz; (void)dRdyy; (void)dRdyz; (void)dRdzz;
    (void)lev;
    launch_common(stream, true, ex, pex, crho, sigma, Rcoord, crho, sigma, Rcoord,
                  drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
                  chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
                  Axx, Axy, Axz, Ayy, Ayz, Azz, Gmx, Gmy, Gmz,
                  rho, Sx, Sy, Sz,
                  Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
                  Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
                  Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
                  Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
                  ham_Res, movx_Res, movy_Res, movz_Res, Gmx_Res, Gmy_Res, Gmz_Res,
                  symmetry, sst);
#endif
}
