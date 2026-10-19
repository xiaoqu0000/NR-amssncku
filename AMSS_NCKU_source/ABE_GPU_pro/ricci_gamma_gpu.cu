#include "macrodef.h"
#include <cuda_runtime.h>
#include <math.h>

#include "ricci_gamma.h"
#include "gpu_manager.h"
#include "gpu_common.cuh"
#include "diff_new_gpu.cuh"
#include "diff_new_sh_gpu.cuh"

// ============================================================================
// Per-point device function: Ricci tensor + Christoffel symbols (Cartesian)
// ============================================================================
__device__ __forceinline__ void d_ricci_gamma_point(
    const int ex[3], const int pex[3],
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    int symmetry,
    int i, int j, int k
) {
    const double ONE = 1.0, TWO = 2.0, HALF = 0.5, F3o2 = 1.5;
    const double SYM = 1.0, ANTI = -1.0;

    // ========================================================================
    // Phase 1: Build metric and compute 1st derivatives
    // ========================================================================

    // g_ij = d_ij + delta_ij  (at this point)
    int idx = (k * pex[1] + j) * pex[0] + i;
    double gxx_v = dxx[idx] + ONE;
    double gyy_v = dyy[idx] + ONE;
    double gzz_v = dzz[idx] + ONE;
    double gxy_v = gxy[idx];
    double gxz_v = gxz[idx];
    double gyz_v = gyz[idx];
    double chi_v = chi[idx];
    double chin1 = fmax(chi_v + ONE, 1e-10);

    // 1st derivatives of chi
    double chix, chiy, chiz;
    d_fderivs_point(ex, pex, chi, &chix, &chiy, &chiz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);

    // 1st derivatives of metric
    double gxxx, gxxy, gxxz;
    d_fderivs_point(ex, pex, dxx, &gxxx, &gxxy, &gxxz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);

    double gxyx, gxyy, gxyz_d;
    d_fderivs_point(ex, pex, gxy, &gxyx, &gxyy, &gxyz_d, X, Y, Z, ANTI, ANTI, SYM, symmetry, 0, i, j, k);

    double gxzx, gxzy, gxzz_d;
    d_fderivs_point(ex, pex, gxz, &gxzx, &gxzy, &gxzz_d, X, Y, Z, ANTI, SYM, ANTI, symmetry, 0, i, j, k);

    double gyyx, gyyy, gyyz;
    d_fderivs_point(ex, pex, dyy, &gyyx, &gyyy, &gyyz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);

    double gyzx, gyzy, gyzz_d;
    d_fderivs_point(ex, pex, gyz, &gyzx, &gyzy, &gyzz_d, X, Y, Z, SYM, ANTI, ANTI, symmetry, 0, i, j, k);

    double gzzx, gzzy, gzzz_d;
    d_fderivs_point(ex, pex, dzz, &gzzx, &gzzy, &gzzz_d, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);

    // Inverse metric
    double detg = gxx_v * gyy_v * gzz_v + gxy_v * gyz_v * gxz_v + gxz_v * gxy_v * gyz_v
                - gxz_v * gyy_v * gxz_v - gxy_v * gxy_v * gzz_v - gxx_v * gyz_v * gyz_v;
    double gupxx = (gyy_v * gzz_v - gyz_v * gyz_v) / detg;
    double gupxy = -(gxy_v * gzz_v - gyz_v * gxz_v) / detg;
    double gupxz = (gxy_v * gyz_v - gyy_v * gxz_v) / detg;
    double gupyy = (gxx_v * gzz_v - gxz_v * gxz_v) / detg;
    double gupyz = -(gxx_v * gyz_v - gxy_v * gxz_v) / detg;
    double gupzz = (gxx_v * gyy_v - gxy_v * gxy_v) / detg;

    // ========================================================================
    // Phase 2: Second kind Christoffel symbols (18 components)
    // ========================================================================

    double l_Gamxxx = HALF*( gupxx*gxxx + gupxy*(TWO*gxyx - gxxy ) + gupxz*(TWO*gxzx - gxxz ));
    double l_Gamyxx = HALF*( gupxy*gxxx + gupyy*(TWO*gxyx - gxxy ) + gupyz*(TWO*gxzx - gxxz ));
    double l_Gamzxx = HALF*( gupxz*gxxx + gupyz*(TWO*gxyx - gxxy ) + gupzz*(TWO*gxzx - gxxz ));

    double l_Gamxyy = HALF*( gupxx*(TWO*gxyy - gyyx ) + gupxy*gyyy + gupxz*(TWO*gyzy - gyyz ));
    double l_Gamyyy = HALF*( gupxy*(TWO*gxyy - gyyx ) + gupyy*gyyy + gupyz*(TWO*gyzy - gyyz ));
    double l_Gamzyy = HALF*( gupxz*(TWO*gxyy - gyyx ) + gupyz*gyyy + gupzz*(TWO*gyzy - gyyz ));

    double l_Gamxzz = HALF*( gupxx*(TWO*gxzz_d - gzzx ) + gupxy*(TWO*gyzz_d - gzzy ) + gupxz*gzzz_d);
    double l_Gamyzz = HALF*( gupxy*(TWO*gxzz_d - gzzx ) + gupyy*(TWO*gyzz_d - gzzy ) + gupyz*gzzz_d);
    double l_Gamzzz = HALF*( gupxz*(TWO*gxzz_d - gzzx ) + gupyz*(TWO*gyzz_d - gzzy ) + gupzz*gzzz_d);

    double l_Gamxxy = HALF*( gupxx*gxxy + gupxy*gyyx + gupxz*( gxzy + gyzx - gxyz_d ) );
    double l_Gamyxy = HALF*( gupxy*gxxy + gupyy*gyyx + gupyz*( gxzy + gyzx - gxyz_d ) );
    double l_Gamzxy = HALF*( gupxz*gxxy + gupyz*gyyx + gupzz*( gxzy + gyzx - gxyz_d ) );

    double l_Gamxxz = HALF*( gupxx*gxxz + gupxy*( gxyz_d + gyzx - gxzy ) + gupxz*gzzx );
    double l_Gamyxz = HALF*( gupxy*gxxz + gupyy*( gxyz_d + gyzx - gxzy ) + gupyz*gzzx );
    double l_Gamzxz = HALF*( gupxz*gxxz + gupyz*( gxyz_d + gyzx - gxzy ) + gupzz*gzzx );

    double l_Gamxyz = HALF*( gupxx*( gxyz_d + gxzy - gyzx ) + gupxy*gyyz + gupxz*gzzy );
    double l_Gamyyz = HALF*( gupxy*( gxyz_d + gxzy - gyzx ) + gupyy*gyyz + gupyz*gzzy );
    double l_Gamzyz = HALF*( gupxz*( gxyz_d + gxzy - gyzx ) + gupyz*gyyz + gupzz*gzzy );

    // ========================================================================
    // Phase 3: Contracted Christoffel + derivatives of Gam
    // ========================================================================

    double Gamxa =       gupxx * l_Gamxxx + gupyy * l_Gamxyy + gupzz * l_Gamxzz +
                  TWO*( gupxy * l_Gamxxy + gupxz * l_Gamxxz + gupyz * l_Gamxyz );
    double Gamya =       gupxx * l_Gamyxx + gupyy * l_Gamyyy + gupzz * l_Gamyzz +
                  TWO*( gupxy * l_Gamyxy + gupxz * l_Gamyxz + gupyz * l_Gamyyz );
    double Gamza =       gupxx * l_Gamzxx + gupyy * l_Gamzyy + gupzz * l_Gamzzz +
                  TWO*( gupxy * l_Gamzxy + gupxz * l_Gamzxz + gupyz * l_Gamzyz );

    double Gamxx, Gamxy, Gamxz;
    d_fderivs_point(ex, pex, Gamx, &Gamxx, &Gamxy, &Gamxz, X, Y, Z, ANTI, SYM, SYM, symmetry, 0, i, j, k);
    double Gamyx, Gamyy, Gamyz;
    d_fderivs_point(ex, pex, Gamy, &Gamyx, &Gamyy, &Gamyz, X, Y, Z, SYM, ANTI, SYM, symmetry, 0, i, j, k);
    double Gamzx, Gamzy, Gamzz;
    d_fderivs_point(ex, pex, Gamz, &Gamzx, &Gamzy, &Gamzz, X, Y, Z, SYM, SYM, ANTI, symmetry, 0, i, j, k);

    // First kind connections g_ij,k
    double gxxx_c = gxx_v * l_Gamxxx + gxy_v * l_Gamyxx + gxz_v * l_Gamzxx;
    double gxyx_c = gxx_v * l_Gamxxy + gxy_v * l_Gamyxy + gxz_v * l_Gamzxy;
    double gxzx_c = gxx_v * l_Gamxxz + gxy_v * l_Gamyxz + gxz_v * l_Gamzxz;
    double gyyx_c = gxx_v * l_Gamxyy + gxy_v * l_Gamyyy + gxz_v * l_Gamzyy;
    double gyzx_c = gxx_v * l_Gamxyz + gxy_v * l_Gamyyz + gxz_v * l_Gamzyz;
    double gzzx_c = gxx_v * l_Gamxzz + gxy_v * l_Gamyzz + gxz_v * l_Gamzzz;

    double gxxy_c = gxy_v * l_Gamxxx + gyy_v * l_Gamyxx + gyz_v * l_Gamzxx;
    double gxyy_c = gxy_v * l_Gamxxy + gyy_v * l_Gamyxy + gyz_v * l_Gamzxy;
    double gxzy_c = gxy_v * l_Gamxxz + gyy_v * l_Gamyxz + gyz_v * l_Gamzxz;
    double gyyy_c = gxy_v * l_Gamxyy + gyy_v * l_Gamyyy + gyz_v * l_Gamzyy;
    double gyzy_c = gxy_v * l_Gamxyz + gyy_v * l_Gamyyz + gyz_v * l_Gamzyz;
    double gzzy_c = gxy_v * l_Gamxzz + gyy_v * l_Gamyzz + gyz_v * l_Gamzzz;

    double gxxz_c = gxz_v * l_Gamxxx + gyz_v * l_Gamyxx + gzz_v * l_Gamzxx;
    double gxyz_c = gxz_v * l_Gamxxy + gyz_v * l_Gamyxy + gzz_v * l_Gamzxy;
    double gxzz_c = gxz_v * l_Gamxxz + gyz_v * l_Gamyxz + gzz_v * l_Gamzxz;
    double gyyz_c = gxz_v * l_Gamxyy + gyz_v * l_Gamyyy + gzz_v * l_Gamzyy;
    double gyzz_c = gxz_v * l_Gamxyz + gyz_v * l_Gamyyz + gzz_v * l_Gamzyz;
    double gzzz_c = gxz_v * l_Gamxzz + gyz_v * l_Gamyzz + gzz_v * l_Gamzzz;

    // ========================================================================
    // Phase 4: Ricci tensor — second derivative part
    // ========================================================================

    double fxx, fxy, fxz, fyy, fyz, fzz;

    // Rxx part from dxx 2nd derivatives
    d_fdderivs_point(ex, pex, dxx, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
    double l_Rxx = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    // Ryy from dyy 2nd derivatives
    d_fdderivs_point(ex, pex, dyy, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
    double l_Ryy = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    // Rzz from dzz 2nd derivatives
    d_fdderivs_point(ex, pex, dzz, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
    double l_Rzz = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    // Rxy from gxy 2nd derivatives
    d_fdderivs_point(ex, pex, gxy, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, ANTI, ANTI, SYM, symmetry, 0, i, j, k);
    double l_Rxy = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    // Rxz from gxz 2nd derivatives
    d_fdderivs_point(ex, pex, gxz, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, ANTI, SYM, ANTI, symmetry, 0, i, j, k);
    double l_Rxz = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    // Ryz from gyz 2nd derivatives
    d_fdderivs_point(ex, pex, gyz, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, SYM, ANTI, ANTI, symmetry, 0, i, j, k);
    double l_Ryz = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    // ========================================================================
    // Phase 4 continued: Add Christoffel products to Ricci
    // ========================================================================

    // Rxx
    l_Rxx = -HALF * l_Rxx +
            gxx_v * Gamxx + gxy_v * Gamyx + gxz_v * Gamzx +
            Gamxa * gxxx_c + Gamya * gxyx_c + Gamza * gxzx_c +
            gupxx * (
                TWO*(l_Gamxxx*gxxx_c + l_Gamyxx*gxyx_c + l_Gamzxx*gxzx_c) +
                     l_Gamxxx*gxxx_c + l_Gamyxx*gxxy_c + l_Gamzxx*gxxz_c) +
            gupxy * (
                TWO*(l_Gamxxx*gxyx_c + l_Gamyxx*gyyx_c + l_Gamzxx*gyzx_c +
                     l_Gamxxy*gxxx_c + l_Gamyxy*gxyx_c + l_Gamzxy*gxzx_c) +
                     l_Gamxxy*gxxx_c + l_Gamyxy*gxxy_c + l_Gamzxy*gxxz_c +
                     l_Gamxxx*gxyx_c + l_Gamyxx*gxyy_c + l_Gamzxx*gxyz_c) +
            gupxz * (
                TWO*(l_Gamxxx*gxzx_c + l_Gamyxx*gyzx_c + l_Gamzxx*gzzx_c +
                     l_Gamxxz*gxxx_c + l_Gamyxz*gxyx_c + l_Gamzxz*gxzx_c) +
                     l_Gamxxz*gxxx_c + l_Gamyxz*gxxy_c + l_Gamzxz*gxxz_c +
                     l_Gamxxx*gxzx_c + l_Gamyxx*gxzy_c + l_Gamzxx*gxzz_c) +
            gupyy * (
                TWO*(l_Gamxxy*gxyx_c + l_Gamyxy*gyyx_c + l_Gamzxy*gyzx_c) +
                     l_Gamxxy*gxyx_c + l_Gamyxy*gxyy_c + l_Gamzxy*gxyz_c) +
            gupyz * (
                TWO*(l_Gamxxy*gxzx_c + l_Gamyxy*gyzx_c + l_Gamzxy*gzzx_c +
                     l_Gamxxz*gxyx_c + l_Gamyxz*gyyx_c + l_Gamzxz*gyzx_c) +
                     l_Gamxxz*gxyx_c + l_Gamyxz*gxyy_c + l_Gamzxz*gxyz_c +
                     l_Gamxxy*gxzx_c + l_Gamyxy*gxzy_c + l_Gamzxy*gxzz_c) +
            gupzz * (
                TWO*(l_Gamxxz*gxzx_c + l_Gamyxz*gyzx_c + l_Gamzxz*gzzx_c) +
                     l_Gamxxz*gxzx_c + l_Gamyxz*gxzy_c + l_Gamzxz*gxzz_c);

    // Ryy
    l_Ryy = -HALF * l_Ryy +
            gxy_v * Gamxy + gyy_v * Gamyy + gyz_v * Gamzy +
            Gamxa * gxyy_c + Gamya * gyyy_c + Gamza * gyzy_c +
            gupxx * (
                TWO*(l_Gamxxy*gxxy_c + l_Gamyxy*gxyy_c + l_Gamzxy*gxzy_c) +
                     l_Gamxxy*gxyx_c + l_Gamyxy*gxyy_c + l_Gamzxy*gxyz_c) +
            gupxy * (
                TWO*(l_Gamxxy*gxyy_c + l_Gamyxy*gyyy_c + l_Gamzxy*gyzy_c +
                     l_Gamxyy*gxxy_c + l_Gamyyy*gxyy_c + l_Gamzyy*gxzy_c) +
                     l_Gamxyy*gxyx_c + l_Gamyyy*gxyy_c + l_Gamzyy*gxyz_c +
                     l_Gamxxy*gyyx_c + l_Gamyxy*gyyy_c + l_Gamzxy*gyyz_c) +
            gupxz * (
                TWO*(l_Gamxxy*gxzy_c + l_Gamyxy*gyzy_c + l_Gamzxy*gzzy_c +
                     l_Gamxyz*gxxy_c + l_Gamyyz*gxyy_c + l_Gamzyz*gxzy_c) +
                     l_Gamxyz*gxyx_c + l_Gamyyz*gxyy_c + l_Gamzyz*gxyz_c +
                     l_Gamxxy*gyzx_c + l_Gamyxy*gyzy_c + l_Gamzxy*gyzz_c) +
            gupyy * (
                TWO*(l_Gamxyy*gxyy_c + l_Gamyyy*gyyy_c + l_Gamzyy*gyzy_c) +
                     l_Gamxyy*gyyx_c + l_Gamyyy*gyyy_c + l_Gamzyy*gyyz_c) +
            gupyz * (
                TWO*(l_Gamxyy*gxzy_c + l_Gamyyy*gyzy_c + l_Gamzyy*gzzy_c +
                     l_Gamxyz*gxyy_c + l_Gamyyz*gyyy_c + l_Gamzyz*gyzy_c) +
                     l_Gamxyz*gyyx_c + l_Gamyyz*gyyy_c + l_Gamzyz*gyyz_c +
                     l_Gamxyy*gyzx_c + l_Gamyyy*gyzy_c + l_Gamzyy*gyzz_c) +
            gupzz * (
                TWO*(l_Gamxyz*gxzy_c + l_Gamyyz*gyzy_c + l_Gamzyz*gzzy_c) +
                     l_Gamxyz*gyzx_c + l_Gamyyz*gyzy_c + l_Gamzyz*gyzz_c);

    // Rzz
    l_Rzz = -HALF * l_Rzz +
            gxz_v * Gamxz + gyz_v * Gamyz + gzz_v * Gamzz +
            Gamxa * gxzz_c + Gamya * gyzz_c + Gamza * gzzz_c +
            gupxx * (
                TWO*(l_Gamxxz*gxxz_c + l_Gamyxz*gxyz_c + l_Gamzxz*gxzz_c) +
                     l_Gamxxz*gxzx_c + l_Gamyxz*gxzy_c + l_Gamzxz*gxzz_c) +
            gupxy * (
                TWO*(l_Gamxxz*gxyz_c + l_Gamyxz*gyyz_c + l_Gamzxz*gyzz_c +
                     l_Gamxyz*gxxz_c + l_Gamyyz*gxyz_c + l_Gamzyz*gxzz_c) +
                     l_Gamxyz*gxzx_c + l_Gamyyz*gxzy_c + l_Gamzyz*gxzz_c +
                     l_Gamxxz*gyzx_c + l_Gamyxz*gyzy_c + l_Gamzxz*gyzz_c) +
            gupxz * (
                TWO*(l_Gamxxz*gxzz_c + l_Gamyxz*gyzz_c + l_Gamzxz*gzzz_c +
                     l_Gamxzz*gxxz_c + l_Gamyzz*gxyz_c + l_Gamzzz*gxzz_c) +
                     l_Gamxzz*gxzx_c + l_Gamyzz*gxzy_c + l_Gamzzz*gxzz_c +
                     l_Gamxxz*gzzx_c + l_Gamyxz*gzzy_c + l_Gamzxz*gzzz_c) +
            gupyy * (
                TWO*(l_Gamxyz*gxyz_c + l_Gamyyz*gyyz_c + l_Gamzyz*gyzz_c) +
                     l_Gamxyz*gyzx_c + l_Gamyyz*gyzy_c + l_Gamzyz*gyzz_c) +
            gupyz * (
                TWO*(l_Gamxyz*gxzz_c + l_Gamyyz*gyzz_c + l_Gamzyz*gzzz_c +
                     l_Gamxzz*gxyz_c + l_Gamyzz*gyyz_c + l_Gamzzz*gyzz_c) +
                     l_Gamxzz*gyzx_c + l_Gamyzz*gyzy_c + l_Gamzzz*gyzz_c +
                     l_Gamxyz*gzzx_c + l_Gamyyz*gzzy_c + l_Gamzyz*gzzz_c) +
            gupzz * (
                TWO*(l_Gamxzz*gxzz_c + l_Gamyzz*gyzz_c + l_Gamzzz*gzzz_c) +
                     l_Gamxzz*gzzx_c + l_Gamyzz*gzzy_c + l_Gamzzz*gzzz_c);

    // Rxy
    l_Rxy = HALF*( -l_Rxy +
            gxx_v * Gamxy + gxy_v * Gamyy + gxz_v * Gamzy +
            gxy_v * Gamxx + gyy_v * Gamyx + gyz_v * Gamzx +
            Gamxa * gxyx_c + Gamya * gyyx_c + Gamza * gyzx_c +
            Gamxa * gxxy_c + Gamya * gxyy_c + Gamza * gxzy_c) +
            gupxx * (
                l_Gamxxx*gxxy_c + l_Gamyxx*gxyy_c + l_Gamzxx*gxzy_c +
                l_Gamxxy*gxxx_c + l_Gamyxy*gxyx_c + l_Gamzxy*gxzx_c +
                l_Gamxxx*gxyx_c + l_Gamyxx*gxyy_c + l_Gamzxx*gxyz_c) +
            gupxy * (
                l_Gamxxx*gxyy_c + l_Gamyxx*gyyy_c + l_Gamzxx*gyzy_c +
                l_Gamxxy*gxyx_c + l_Gamyxy*gyyx_c + l_Gamzxy*gyzx_c +
                l_Gamxxy*gxyx_c + l_Gamyxy*gxyy_c + l_Gamzxy*gxyz_c +
                l_Gamxxy*gxxy_c + l_Gamyxy*gxyy_c + l_Gamzxy*gxzy_c +
                l_Gamxyy*gxxx_c + l_Gamyyy*gxyx_c + l_Gamzyy*gxzx_c +
                l_Gamxxx*gyyx_c + l_Gamyxx*gyyy_c + l_Gamzxx*gyyz_c) +
            gupxz * (
                l_Gamxxx*gxzy_c + l_Gamyxx*gyzy_c + l_Gamzxx*gzzy_c +
                l_Gamxxy*gxzx_c + l_Gamyxy*gyzx_c + l_Gamzxy*gzzx_c +
                l_Gamxxz*gxyx_c + l_Gamyxz*gxyy_c + l_Gamzxz*gxyz_c +
                l_Gamxxz*gxxy_c + l_Gamyxz*gxyy_c + l_Gamzxz*gxzy_c +
                l_Gamxyz*gxxx_c + l_Gamyyz*gxyx_c + l_Gamzyz*gxzx_c +
                l_Gamxxx*gyzx_c + l_Gamyxx*gyzy_c + l_Gamzxx*gyzz_c) +
            gupyy * (
                l_Gamxxy*gxyy_c + l_Gamyxy*gyyy_c + l_Gamzxy*gyzy_c +
                l_Gamxyy*gxyx_c + l_Gamyyy*gyyx_c + l_Gamzyy*gyzx_c +
                l_Gamxxy*gyyx_c + l_Gamyxy*gyyy_c + l_Gamzxy*gyyz_c) +
            gupyz * (
                l_Gamxxy*gxzy_c + l_Gamyxy*gyzy_c + l_Gamzxy*gzzy_c +
                l_Gamxyy*gxzx_c + l_Gamyyy*gyzx_c + l_Gamzyy*gzzx_c +
                l_Gamxxz*gyyx_c + l_Gamyxz*gyyy_c + l_Gamzxz*gyyz_c +
                l_Gamxxz*gxyy_c + l_Gamyxz*gyyy_c + l_Gamzxz*gyzy_c +
                l_Gamxyz*gxyx_c + l_Gamyyz*gyyx_c + l_Gamzyz*gyzx_c +
                l_Gamxxy*gyzx_c + l_Gamyxy*gyzy_c + l_Gamzxy*gyzz_c) +
            gupzz * (
                l_Gamxxz*gxzy_c + l_Gamyxz*gyzy_c + l_Gamzxz*gzzy_c +
                l_Gamxyz*gxzx_c + l_Gamyyz*gyzx_c + l_Gamzyz*gzzx_c +
                l_Gamxxz*gyzx_c + l_Gamyxz*gyzy_c + l_Gamzxz*gyzz_c);

    // Rxz
    l_Rxz = HALF*( -l_Rxz +
            gxx_v * Gamxz + gxy_v * Gamyz + gxz_v * Gamzz +
            gxz_v * Gamxx + gyz_v * Gamyx + gzz_v * Gamzx +
            Gamxa * gxzx_c + Gamya * gyzx_c + Gamza * gzzx_c +
            Gamxa * gxxz_c + Gamya * gxyz_c + Gamza * gxzz_c) +
            gupxx * (
                l_Gamxxx*gxxz_c + l_Gamyxx*gxyz_c + l_Gamzxx*gxzz_c +
                l_Gamxxz*gxxx_c + l_Gamyxz*gxyx_c + l_Gamzxz*gxzx_c +
                l_Gamxxx*gxzx_c + l_Gamyxx*gxzy_c + l_Gamzxx*gxzz_c) +
            gupxy * (
                l_Gamxxx*gxyz_c + l_Gamyxx*gyyz_c + l_Gamzxx*gyzz_c +
                l_Gamxxz*gxyx_c + l_Gamyxz*gyyx_c + l_Gamzxz*gyzx_c +
                l_Gamxxy*gxzx_c + l_Gamyxy*gxzy_c + l_Gamzxy*gxzz_c +
                l_Gamxxy*gxxz_c + l_Gamyxy*gxyz_c + l_Gamzxy*gxzz_c +
                l_Gamxyz*gxxx_c + l_Gamyyz*gxyx_c + l_Gamzyz*gxzx_c +
                l_Gamxxx*gyzx_c + l_Gamyxx*gyzy_c + l_Gamzxx*gyzz_c) +
            gupxz * (
                l_Gamxxx*gxzz_c + l_Gamyxx*gyzz_c + l_Gamzxx*gzzz_c +
                l_Gamxxz*gxzx_c + l_Gamyxz*gyzx_c + l_Gamzxz*gzzx_c +
                l_Gamxxz*gxzx_c + l_Gamyxz*gxzy_c + l_Gamzxz*gxzz_c +
                l_Gamxxz*gxxz_c + l_Gamyxz*gxyz_c + l_Gamzxz*gxzz_c +
                l_Gamxzz*gxxx_c + l_Gamyzz*gxyx_c + l_Gamzzz*gxzx_c +
                l_Gamxxx*gzzx_c + l_Gamyxx*gzzy_c + l_Gamzxx*gzzz_c) +
            gupyy * (
                l_Gamxxy*gxyz_c + l_Gamyxy*gyyz_c + l_Gamzxy*gyzz_c +
                l_Gamxyz*gxyx_c + l_Gamyyz*gyyx_c + l_Gamzyz*gyzx_c +
                l_Gamxxy*gyzx_c + l_Gamyxy*gyzy_c + l_Gamzxy*gyzz_c) +
            gupyz * (
                l_Gamxxy*gxzz_c + l_Gamyxy*gyzz_c + l_Gamzxy*gzzz_c +
                l_Gamxyz*gxzx_c + l_Gamyyz*gyzx_c + l_Gamzyz*gzzx_c +
                l_Gamxxz*gyzx_c + l_Gamyxz*gyzy_c + l_Gamzxz*gyzz_c +
                l_Gamxxz*gxyz_c + l_Gamyxz*gyyz_c + l_Gamzxz*gyzz_c +
                l_Gamxzz*gxyx_c + l_Gamyzz*gyyx_c + l_Gamzzz*gyzx_c +
                l_Gamxxy*gzzx_c + l_Gamyxy*gzzy_c + l_Gamzxy*gzzz_c) +
            gupzz * (
                l_Gamxxz*gxzz_c + l_Gamyxz*gyzz_c + l_Gamzxz*gzzz_c +
                l_Gamxzz*gxzx_c + l_Gamyzz*gyzx_c + l_Gamzzz*gzzx_c +
                l_Gamxxz*gzzx_c + l_Gamyxz*gzzy_c + l_Gamzxz*gzzz_c);

    // Ryz
    l_Ryz = HALF*( -l_Ryz +
            gxy_v * Gamxz + gyy_v * Gamyz + gyz_v * Gamzz +
            gxz_v * Gamxy + gyz_v * Gamyy + gzz_v * Gamzy +
            Gamxa * gxzy_c + Gamya * gyzy_c + Gamza * gzzy_c +
            Gamxa * gxyz_c + Gamya * gyyz_c + Gamza * gyzz_c) +
            gupxx * (
                l_Gamxxy*gxxz_c + l_Gamyxy*gxyz_c + l_Gamzxy*gxzz_c +
                l_Gamxxz*gxxy_c + l_Gamyxz*gxyy_c + l_Gamzxz*gxzy_c +
                l_Gamxxy*gxzx_c + l_Gamyxy*gxzy_c + l_Gamzxy*gxzz_c) +
            gupxy * (
                l_Gamxxy*gxyz_c + l_Gamyxy*gyyz_c + l_Gamzxy*gyzz_c +
                l_Gamxxz*gxyy_c + l_Gamyxz*gyyy_c + l_Gamzxz*gyzy_c +
                l_Gamxyy*gxzx_c + l_Gamyyy*gxzy_c + l_Gamzyy*gxzz_c +
                l_Gamxyy*gxxz_c + l_Gamyyy*gxyz_c + l_Gamzyy*gxzz_c +
                l_Gamxyz*gxxy_c + l_Gamyyz*gxyy_c + l_Gamzyz*gxzy_c +
                l_Gamxxy*gyzx_c + l_Gamyxy*gyzy_c + l_Gamzxy*gyzz_c) +
            gupxz * (
                l_Gamxxy*gxzz_c + l_Gamyxy*gyzz_c + l_Gamzxy*gzzz_c +
                l_Gamxxz*gxzy_c + l_Gamyxz*gyzy_c + l_Gamzxz*gzzy_c +
                l_Gamxyz*gxzx_c + l_Gamyyz*gxzy_c + l_Gamzyz*gxzz_c +
                l_Gamxyz*gxxz_c + l_Gamyyz*gxyz_c + l_Gamzyz*gxzz_c +
                l_Gamxzz*gxxy_c + l_Gamyzz*gxyy_c + l_Gamzzz*gxzy_c +
                l_Gamxxy*gzzx_c + l_Gamyxy*gzzy_c + l_Gamzxy*gzzz_c) +
            gupyy * (
                l_Gamxyy*gxyz_c + l_Gamyyy*gyyz_c + l_Gamzyy*gyzz_c +
                l_Gamxyz*gxyy_c + l_Gamyyz*gyyy_c + l_Gamzyz*gyzy_c +
                l_Gamxyy*gyzx_c + l_Gamyyy*gyzy_c + l_Gamzyy*gyzz_c) +
            gupyz * (
                l_Gamxyy*gxzz_c + l_Gamyyy*gyzz_c + l_Gamzyy*gzzz_c +
                l_Gamxyz*gxzy_c + l_Gamyyz*gyzy_c + l_Gamzyz*gzzy_c +
                l_Gamxyz*gyzx_c + l_Gamyyz*gyzy_c + l_Gamzyz*gyzz_c +
                l_Gamxyz*gxyz_c + l_Gamyyz*gyyy_c + l_Gamzyz*gyzz_c +
                l_Gamxzz*gxyy_c + l_Gamyzz*gyyy_c + l_Gamzzz*gyzy_c +
                l_Gamxyy*gzzx_c + l_Gamyyy*gzzy_c + l_Gamzyy*gzzz_c) +
            gupzz * (
                l_Gamxyz*gxzz_c + l_Gamyyz*gyzz_c + l_Gamzyz*gzzz_c +
                l_Gamxzz*gxzy_c + l_Gamyzz*gyzy_c + l_Gamzzz*gzzy_c +
                l_Gamxyz*gzzx_c + l_Gamyyz*gzzy_c + l_Gamzyz*gzzz_c);

    // ========================================================================
    // Phase 5: Chi contribution to Ricci
    // ========================================================================

    d_fdderivs_point(ex, pex, chi, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);

    fxx = fxx - l_Gamxxx * chix - l_Gamyxx * chiy - l_Gamzxx * chiz;
    fxy = fxy - l_Gamxxy * chix - l_Gamyxy * chiy - l_Gamzxy * chiz;
    fxz = fxz - l_Gamxxz * chix - l_Gamyxz * chiy - l_Gamzxz * chiz;
    fyy = fyy - l_Gamxyy * chix - l_Gamyyy * chiy - l_Gamzyy * chiz;
    fyz = fyz - l_Gamxyz * chix - l_Gamyyz * chiy - l_Gamzyz * chiz;
    fzz = fzz - l_Gamxzz * chix - l_Gamyzz * chiy - l_Gamzzz * chiz;

    double f = gupxx*(fxx - F3o2/chin1*chix*chix) +
               gupyy*(fyy - F3o2/chin1*chiy*chiy) +
               gupzz*(fzz - F3o2/chin1*chiz*chiz) +
         TWO*(gupxy*(fxy - F3o2/chin1*chix*chiy) +
              gupxz*(fxz - F3o2/chin1*chix*chiz) +
              gupyz*(fyz - F3o2/chin1*chiy*chiz));

    l_Rxx = l_Rxx + (fxx - chix*chix/chin1/TWO + gxx_v*f)/chin1/TWO;
    l_Ryy = l_Ryy + (fyy - chiy*chiy/chin1/TWO + gyy_v*f)/chin1/TWO;
    l_Rzz = l_Rzz + (fzz - chiz*chiz/chin1/TWO + gzz_v*f)/chin1/TWO;
    l_Rxy = l_Rxy + (fxy - chix*chiy/chin1/TWO + gxy_v*f)/chin1/TWO;
    l_Rxz = l_Rxz + (fxz - chix*chiz/chin1/TWO + gxz_v*f)/chin1/TWO;
    l_Ryz = l_Ryz + (fyz - chiy*chiz/chin1/TWO + gyz_v*f)/chin1/TWO;

    // ========================================================================
    // Phase 6: Physical Christoffel symbols
    // ========================================================================

    double gxxx_ch = (gupxx*chix + gupxy*chiy + gupxz*chiz)/chin1;
    double gxxy_ch = (gupxy*chix + gupyy*chiy + gupyz*chiz)/chin1;
    double gxxz_ch = (gupxz*chix + gupyz*chiy + gupzz*chiz)/chin1;

    *Gamxxx = l_Gamxxx - ((chix+chix)/chin1 - gxx_v*gxxx_ch)*HALF;
    *Gamyxx = l_Gamyxx - (                    - gxx_v*gxxy_ch)*HALF;
    *Gamzxx = l_Gamzxx - (                    - gxx_v*gxxz_ch)*HALF;
    *Gamxyy = l_Gamxyy - (                    - gyy_v*gxxx_ch)*HALF;
    *Gamyyy = l_Gamyyy - ((chiy+chiy)/chin1 - gyy_v*gxxy_ch)*HALF;
    *Gamzyy = l_Gamzyy - (                    - gyy_v*gxxz_ch)*HALF;
    *Gamxzz = l_Gamxzz - (                    - gzz_v*gxxx_ch)*HALF;
    *Gamyzz = l_Gamyzz - (                    - gzz_v*gxxy_ch)*HALF;
    *Gamzzz = l_Gamzzz - ((chiz+chiz)/chin1 - gzz_v*gxxz_ch)*HALF;
    *Gamxxy = l_Gamxxy - ( chiy       /chin1 - gxy_v*gxxx_ch)*HALF;
    *Gamyxy = l_Gamyxy - (        chix/chin1 - gxy_v*gxxy_ch)*HALF;
    *Gamzxy = l_Gamzxy - (                    - gxy_v*gxxz_ch)*HALF;
    *Gamxxz = l_Gamxxz - ( chiz       /chin1 - gxz_v*gxxx_ch)*HALF;
    *Gamyxz = l_Gamyxz - (                    - gxz_v*gxxy_ch)*HALF;
    *Gamzxz = l_Gamzxz - (        chix/chin1 - gxz_v*gxxz_ch)*HALF;
    *Gamxyz = l_Gamxyz - (                    - gyz_v*gxxx_ch)*HALF;
    *Gamyyz = l_Gamyyz - ( chiz       /chin1 - gyz_v*gxxy_ch)*HALF;
    *Gamzyz = l_Gamzyz - (        chiy/chin1 - gyz_v*gxxz_ch)*HALF;

    *Rxx = l_Rxx;
    *Ryy = l_Ryy;
    *Rzz = l_Rzz;
    *Rxy = l_Rxy;
    *Rxz = l_Rxz;
    *Ryz = l_Ryz;
}

// ============================================================================
// Global kernel: Ricci tensor + Christoffel symbols
// ============================================================================
__global__ void ricci_gamma_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    double* __restrict__ Gamxxx, double* __restrict__ Gamxxy, double* __restrict__ Gamxxz,
    double* __restrict__ Gamxyy, double* __restrict__ Gamxyz, double* __restrict__ Gamxzz,
    double* __restrict__ Gamyxx, double* __restrict__ Gamyxy, double* __restrict__ Gamyxz,
    double* __restrict__ Gamyyy, double* __restrict__ Gamyyz, double* __restrict__ Gamyzz,
    double* __restrict__ Gamzxx, double* __restrict__ Gamzxy, double* __restrict__ Gamzxz,
    double* __restrict__ Gamzyy, double* __restrict__ Gamzyz, double* __restrict__ Gamzzz,
    double* __restrict__ Rxx, double* __restrict__ Rxy, double* __restrict__ Rxz,
    double* __restrict__ Ryy, double* __restrict__ Ryz, double* __restrict__ Rzz,
    int symmetry
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total = ex0 * ex1 * ex2;
    if (tid >= total) return;

    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);

    int idx = (k * pex1 + j) * pex0 + i;

    int ex[3]  = {ex0, ex1, ex2};
    int pex[3] = {pex0, pex1, pex2};

    double l_Gamxxx, l_Gamxxy, l_Gamxxz;
    double l_Gamxyy, l_Gamxyz, l_Gamxzz;
    double l_Gamyxx, l_Gamyxy, l_Gamyxz;
    double l_Gamyyy, l_Gamyyz, l_Gamyzz;
    double l_Gamzxx, l_Gamzxy, l_Gamzxz;
    double l_Gamzyy, l_Gamzyz, l_Gamzzz;
    double l_Rxx, l_Rxy, l_Rxz, l_Ryy, l_Ryz, l_Rzz;

    d_ricci_gamma_point(
        ex, pex,
        X, Y, Z,
        chi, dxx, gxy, gxz, dyy, gyz, dzz,
        Gamx, Gamy, Gamz,
        &l_Gamxxx, &l_Gamxxy, &l_Gamxxz,
        &l_Gamxyy, &l_Gamxyz, &l_Gamxzz,
        &l_Gamyxx, &l_Gamyxy, &l_Gamyxz,
        &l_Gamyyy, &l_Gamyyz, &l_Gamyzz,
        &l_Gamzxx, &l_Gamzxy, &l_Gamzxz,
        &l_Gamzyy, &l_Gamzyz, &l_Gamzzz,
        &l_Rxx, &l_Rxy, &l_Rxz,
        &l_Ryy, &l_Ryz, &l_Rzz,
        symmetry, i, j, k
    );

    Gamxxx[idx] = l_Gamxxx; Gamxxy[idx] = l_Gamxxy; Gamxxz[idx] = l_Gamxxz;
    Gamxyy[idx] = l_Gamxyy; Gamxyz[idx] = l_Gamxyz; Gamxzz[idx] = l_Gamxzz;
    Gamyxx[idx] = l_Gamyxx; Gamyxy[idx] = l_Gamyxy; Gamyxz[idx] = l_Gamyxz;
    Gamyyy[idx] = l_Gamyyy; Gamyyz[idx] = l_Gamyyz; Gamyzz[idx] = l_Gamyzz;
    Gamzxx[idx] = l_Gamzxx; Gamzxy[idx] = l_Gamzxy; Gamzxz[idx] = l_Gamzxz;
    Gamzyy[idx] = l_Gamzyy; Gamzyz[idx] = l_Gamzyz; Gamzzz[idx] = l_Gamzzz;

    Rxx[idx] = l_Rxx; Rxy[idx] = l_Rxy; Rxz[idx] = l_Rxz;
    Ryy[idx] = l_Ryy; Ryz[idx] = l_Ryz; Rzz[idx] = l_Rzz;
}

// ============================================================================
// Host launch wrapper
// ============================================================================
void gpu_ricci_gamma_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Gamx, const double* d_Gamy, const double* d_Gamz,
    double* d_Gamxxx, double* d_Gamxxy, double* d_Gamxxz,
    double* d_Gamxyy, double* d_Gamxyz, double* d_Gamxzz,
    double* d_Gamyxx, double* d_Gamyxy, double* d_Gamyxz,
    double* d_Gamyyy, double* d_Gamyyz, double* d_Gamyzz,
    double* d_Gamzxx, double* d_Gamzxy, double* d_Gamzxz,
    double* d_Gamzyy, double* d_Gamzyz, double* d_Gamzzz,
    double* d_Rxx, double* d_Rxy, double* d_Rxz,
    double* d_Ryy, double* d_Ryz, double* d_Rzz,
    int symmetry
) {
    int total = ex[0] * ex[1] * ex[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);

    ricci_gamma_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2],
        pex[0], pex[1], pex[2],
        d_X, d_Y, d_Z,
        d_chi,
        d_dxx, d_gxy, d_gxz,
        d_dyy, d_gyz, d_dzz,
        d_Gamx, d_Gamy, d_Gamz,
        d_Gamxxx, d_Gamxxy, d_Gamxxz,
        d_Gamxyy, d_Gamxyz, d_Gamxzz,
        d_Gamyxx, d_Gamyxy, d_Gamyxz,
        d_Gamyyy, d_Gamyyz, d_Gamyzz,
        d_Gamzxx, d_Gamzxy, d_Gamzxz,
        d_Gamzyy, d_Gamzyz, d_Gamzzz,
        d_Rxx, d_Rxy, d_Rxz,
        d_Ryy, d_Ryz, d_Rzz,
        symmetry
    );
}

// ============================================================================
// Shell coordinate variant: d_ricci_gamma_ss_point
// Identical computation to Cartesian, but uses d_fderivs_sh_point / d_fdderivs_sh_point
// which handle shell coordinate boundary conditions (upper+lower reflection)
// ============================================================================
__device__ __forceinline__ void d_ricci_gamma_ss_point(
    const int ex[3], const int pex[3],
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    int symmetry, int sst,
    int i, int j, int k
) {
    const double ONE = 1.0, TWO = 2.0, HALF = 0.5, F3o2 = 1.5;
    const double SYM = 1.0, ANTI = -1.0;

    int idx = (k * pex[1] + j) * pex[0] + i;
    double gxx_v = dxx[idx] + ONE;
    double gyy_v = dyy[idx] + ONE;
    double gzz_v = dzz[idx] + ONE;
    double gxy_v = gxy[idx];
    double gxz_v = gxz[idx];
    double gyz_v = gyz[idx];
    double chi_v = chi[idx];
    double chin1 = fmax(chi_v + ONE, 1e-10);

    // 1st derivatives (shell)
    double chix, chiy, chiz;
    d_fderivs_sh_point(ex, pex, chi, &chix, &chiy, &chiz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);

    double gxxx, gxxy, gxxz;
    d_fderivs_sh_point(ex, pex, dxx, &gxxx, &gxxy, &gxxz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);

    double gxyx, gxyy, gxyz_d;
    d_fderivs_sh_point(ex, pex, gxy, &gxyx, &gxyy, &gxyz_d, X, Y, Z, ANTI, ANTI, SYM, symmetry, 0, sst, i, j, k);

    double gxzx, gxzy, gxzz_d;
    d_fderivs_sh_point(ex, pex, gxz, &gxzx, &gxzy, &gxzz_d, X, Y, Z, ANTI, SYM, ANTI, symmetry, 0, sst, i, j, k);

    double gyyx, gyyy, gyyz;
    d_fderivs_sh_point(ex, pex, dyy, &gyyx, &gyyy, &gyyz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);

    double gyzx, gyzy, gyzz_d;
    d_fderivs_sh_point(ex, pex, gyz, &gyzx, &gyzy, &gyzz_d, X, Y, Z, SYM, ANTI, ANTI, symmetry, 0, sst, i, j, k);

    double gzzx, gzzy, gzzz_d;
    d_fderivs_sh_point(ex, pex, dzz, &gzzx, &gzzy, &gzzz_d, X, Y, Z, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);

    // Inverse metric
    double detg = gxx_v * gyy_v * gzz_v + gxy_v * gyz_v * gxz_v + gxz_v * gxy_v * gyz_v
                - gxz_v * gyy_v * gxz_v - gxy_v * gxy_v * gzz_v - gxx_v * gyz_v * gyz_v;
    double gupxx = (gyy_v * gzz_v - gyz_v * gyz_v) / detg;
    double gupxy = -(gxy_v * gzz_v - gyz_v * gxz_v) / detg;
    double gupxz = (gxy_v * gyz_v - gyy_v * gxz_v) / detg;
    double gupyy = (gxx_v * gzz_v - gxz_v * gxz_v) / detg;
    double gupyz = -(gxx_v * gyz_v - gxy_v * gxz_v) / detg;
    double gupzz = (gxx_v * gyy_v - gxy_v * gxy_v) / detg;

    // Second kind Christoffel symbols (18 components) — same formulas as Cartesian
    double l_Gamxxx = HALF*( gupxx*gxxx + gupxy*(TWO*gxyx - gxxy ) + gupxz*(TWO*gxzx - gxxz ));
    double l_Gamyxx = HALF*( gupxy*gxxx + gupyy*(TWO*gxyx - gxxy ) + gupyz*(TWO*gxzx - gxxz ));
    double l_Gamzxx = HALF*( gupxz*gxxx + gupyz*(TWO*gxyx - gxxy ) + gupzz*(TWO*gxzx - gxxz ));

    double l_Gamxyy = HALF*( gupxx*(TWO*gxyy - gyyx ) + gupxy*gyyy + gupxz*(TWO*gyzy - gyyz ));
    double l_Gamyyy = HALF*( gupxy*(TWO*gxyy - gyyx ) + gupyy*gyyy + gupyz*(TWO*gyzy - gyyz ));
    double l_Gamzyy = HALF*( gupxz*(TWO*gxyy - gyyx ) + gupyz*gyyy + gupzz*(TWO*gyzy - gyyz ));

    double l_Gamxzz = HALF*( gupxx*(TWO*gxzz_d - gzzx ) + gupxy*(TWO*gyzz_d - gzzy ) + gupxz*gzzz_d);
    double l_Gamyzz = HALF*( gupxy*(TWO*gxzz_d - gzzx ) + gupyy*(TWO*gyzz_d - gzzy ) + gupyz*gzzz_d);
    double l_Gamzzz = HALF*( gupxz*(TWO*gxzz_d - gzzx ) + gupyz*(TWO*gyzz_d - gzzy ) + gupzz*gzzz_d);

    double l_Gamxxy = HALF*( gupxx*gxxy + gupxy*gyyx + gupxz*( gxzy + gyzx - gxyz_d ) );
    double l_Gamyxy = HALF*( gupxy*gxxy + gupyy*gyyx + gupyz*( gxzy + gyzx - gxyz_d ) );
    double l_Gamzxy = HALF*( gupxz*gxxy + gupyz*gyyx + gupzz*( gxzy + gyzx - gxyz_d ) );

    double l_Gamxxz = HALF*( gupxx*gxxz + gupxy*( gxyz_d + gyzx - gxzy ) + gupxz*gzzx );
    double l_Gamyxz = HALF*( gupxy*gxxz + gupyy*( gxyz_d + gyzx - gxzy ) + gupyz*gzzx );
    double l_Gamzxz = HALF*( gupxz*gxxz + gupyz*( gxyz_d + gyzx - gxzy ) + gupzz*gzzx );

    double l_Gamxyz = HALF*( gupxx*( gxyz_d + gxzy - gyzx ) + gupxy*gyyz + gupxz*gzzy );
    double l_Gamyyz = HALF*( gupxy*( gxyz_d + gxzy - gyzx ) + gupyy*gyyz + gupyz*gzzy );
    double l_Gamzyz = HALF*( gupxz*( gxyz_d + gxzy - gyzx ) + gupyz*gyyz + gupzz*gzzy );

    // Contracted Christoffel + Gam derivatives (shell)
    double Gamxa =       gupxx * l_Gamxxx + gupyy * l_Gamxyy + gupzz * l_Gamxzz +
                  TWO*( gupxy * l_Gamxxy + gupxz * l_Gamxxz + gupyz * l_Gamxyz );
    double Gamya =       gupxx * l_Gamyxx + gupyy * l_Gamyyy + gupzz * l_Gamyzz +
                  TWO*( gupxy * l_Gamyxy + gupxz * l_Gamyxz + gupyz * l_Gamyyz );
    double Gamza =       gupxx * l_Gamzxx + gupyy * l_Gamzyy + gupzz * l_Gamzzz +
                  TWO*( gupxy * l_Gamzxy + gupxz * l_Gamzxz + gupyz * l_Gamzyz );

    double Gamxx, Gamxy, Gamxz;
    d_fderivs_sh_point(ex, pex, Gamx, &Gamxx, &Gamxy, &Gamxz, X, Y, Z, ANTI, SYM, SYM, symmetry, 0, sst, i, j, k);
    double Gamyx, Gamyy, Gamyz;
    d_fderivs_sh_point(ex, pex, Gamy, &Gamyx, &Gamyy, &Gamyz, X, Y, Z, SYM, ANTI, SYM, symmetry, 0, sst, i, j, k);
    double Gamzx, Gamzy, Gamzz;
    d_fderivs_sh_point(ex, pex, Gamz, &Gamzx, &Gamzy, &Gamzz, X, Y, Z, SYM, SYM, ANTI, symmetry, 0, sst, i, j, k);

    // First kind connections
    double gxxx_c = gxx_v * l_Gamxxx + gxy_v * l_Gamyxx + gxz_v * l_Gamzxx;
    double gxyx_c = gxx_v * l_Gamxxy + gxy_v * l_Gamyxy + gxz_v * l_Gamzxy;
    double gxzx_c = gxx_v * l_Gamxxz + gxy_v * l_Gamyxz + gxz_v * l_Gamzxz;
    double gyyx_c = gxx_v * l_Gamxyy + gxy_v * l_Gamyyy + gxz_v * l_Gamzyy;
    double gyzx_c = gxx_v * l_Gamxyz + gxy_v * l_Gamyyz + gxz_v * l_Gamzyz;
    double gzzx_c = gxx_v * l_Gamxzz + gxy_v * l_Gamyzz + gxz_v * l_Gamzzz;

    double gxxy_c = gxy_v * l_Gamxxx + gyy_v * l_Gamyxx + gyz_v * l_Gamzxx;
    double gxyy_c = gxy_v * l_Gamxxy + gyy_v * l_Gamyxy + gyz_v * l_Gamzxy;
    double gxzy_c = gxy_v * l_Gamxxz + gyy_v * l_Gamyxz + gyz_v * l_Gamzxz;
    double gyyy_c = gxy_v * l_Gamxyy + gyy_v * l_Gamyyy + gyz_v * l_Gamzyy;
    double gyzy_c = gxy_v * l_Gamxyz + gyy_v * l_Gamyyz + gyz_v * l_Gamzyz;
    double gzzy_c = gxy_v * l_Gamxzz + gyy_v * l_Gamyzz + gyz_v * l_Gamzzz;

    double gxxz_c = gxz_v * l_Gamxxx + gyz_v * l_Gamyxx + gzz_v * l_Gamzxx;
    double gxyz_c = gxz_v * l_Gamxxy + gyz_v * l_Gamyxy + gzz_v * l_Gamzxy;
    double gxzz_c = gxz_v * l_Gamxxz + gyz_v * l_Gamyxz + gzz_v * l_Gamzxz;
    double gyyz_c = gxz_v * l_Gamxyy + gyz_v * l_Gamyyy + gzz_v * l_Gamzyy;
    double gyzz_c = gxz_v * l_Gamxyz + gyz_v * l_Gamyyz + gzz_v * l_Gamzyz;
    double gzzz_c = gxz_v * l_Gamxzz + gyz_v * l_Gamyzz + gzz_v * l_Gamzzz;

    // Ricci tensor — second derivatives (shell)
    double fxx, fxy, fxz, fyy, fyz, fzz;

    d_fdderivs_sh_point(ex, pex, dxx, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);
    double l_Rxx = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    d_fdderivs_sh_point(ex, pex, dyy, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);
    double l_Ryy = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    d_fdderivs_sh_point(ex, pex, dzz, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);
    double l_Rzz = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    d_fdderivs_sh_point(ex, pex, gxy, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, ANTI, ANTI, SYM, symmetry, 0, sst, i, j, k);
    double l_Rxy = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    d_fdderivs_sh_point(ex, pex, gxz, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, ANTI, SYM, ANTI, symmetry, 0, sst, i, j, k);
    double l_Rxz = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    d_fdderivs_sh_point(ex, pex, gyz, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, SYM, ANTI, ANTI, symmetry, 0, sst, i, j, k);
    double l_Ryz = gupxx*fxx + gupyy*fyy + gupzz*fzz + TWO*(gupxy*fxy + gupxz*fxz + gupyz*fyz);

    // Rxx Christoffel products (same as Cartesian)
    l_Rxx = -HALF * l_Rxx + gxx_v * Gamxx + gxy_v * Gamyx + gxz_v * Gamzx + Gamxa * gxxx_c + Gamya * gxyx_c + Gamza * gxzx_c +
        gupxx*(TWO*(l_Gamxxx*gxxx_c+l_Gamyxx*gxyx_c+l_Gamzxx*gxzx_c)+l_Gamxxx*gxxx_c+l_Gamyxx*gxxy_c+l_Gamzxx*gxxz_c) +
        gupxy*(TWO*(l_Gamxxx*gxyx_c+l_Gamyxx*gyyx_c+l_Gamzxx*gyzx_c+l_Gamxxy*gxxx_c+l_Gamyxy*gxyx_c+l_Gamzxy*gxzx_c)+l_Gamxxy*gxxx_c+l_Gamyxy*gxxy_c+l_Gamzxy*gxxz_c+l_Gamxxx*gxyx_c+l_Gamyxx*gxyy_c+l_Gamzxx*gxyz_c) +
        gupxz*(TWO*(l_Gamxxx*gxzx_c+l_Gamyxx*gyzx_c+l_Gamzxx*gzzx_c+l_Gamxxz*gxxx_c+l_Gamyxz*gxyx_c+l_Gamzxz*gxzx_c)+l_Gamxxz*gxxx_c+l_Gamyxz*gxxy_c+l_Gamzxz*gxxz_c+l_Gamxxx*gxzx_c+l_Gamyxx*gxzy_c+l_Gamzxx*gxzz_c) +
        gupyy*(TWO*(l_Gamxxy*gxyx_c+l_Gamyxy*gyyx_c+l_Gamzxy*gyzx_c)+l_Gamxxy*gxyx_c+l_Gamyxy*gxyy_c+l_Gamzxy*gxyz_c) +
        gupyz*(TWO*(l_Gamxxy*gxzx_c+l_Gamyxy*gyzx_c+l_Gamzxy*gzzx_c+l_Gamxxz*gxyx_c+l_Gamyxz*gyyx_c+l_Gamzxz*gyzx_c)+l_Gamxxz*gxyx_c+l_Gamyxz*gxyy_c+l_Gamzxz*gxyz_c+l_Gamxxy*gxzx_c+l_Gamyxy*gxzy_c+l_Gamzxy*gxzz_c) +
        gupzz*(TWO*(l_Gamxxz*gxzx_c+l_Gamyxz*gyzx_c+l_Gamzxz*gzzx_c)+l_Gamxxz*gxzx_c+l_Gamyxz*gxzy_c+l_Gamzxz*gxzz_c);

    // Ryy (same as Cartesian)
    l_Ryy = -HALF * l_Ryy + gxy_v * Gamxy + gyy_v * Gamyy + gyz_v * Gamzy + Gamxa * gxyy_c + Gamya * gyyy_c + Gamza * gyzy_c +
        gupxx*(TWO*(l_Gamxxy*gxxy_c+l_Gamyxy*gxyy_c+l_Gamzxy*gxzy_c)+l_Gamxxy*gxyx_c+l_Gamyxy*gxyy_c+l_Gamzxy*gxyz_c) +
        gupxy*(TWO*(l_Gamxxy*gxyy_c+l_Gamyxy*gyyy_c+l_Gamzxy*gyzy_c+l_Gamxyy*gxxy_c+l_Gamyyy*gxyy_c+l_Gamzyy*gxzy_c)+l_Gamxyy*gxyx_c+l_Gamyyy*gxyy_c+l_Gamzyy*gxyz_c+l_Gamxxy*gyyx_c+l_Gamyxy*gyyy_c+l_Gamzxy*gyyz_c) +
        gupxz*(TWO*(l_Gamxxy*gxzy_c+l_Gamyxy*gyzy_c+l_Gamzxy*gzzy_c+l_Gamxyz*gxxy_c+l_Gamyyz*gxyy_c+l_Gamzyz*gxzy_c)+l_Gamxyz*gxyx_c+l_Gamyyz*gxyy_c+l_Gamzyz*gxyz_c+l_Gamxxy*gyzx_c+l_Gamyxy*gyzy_c+l_Gamzxy*gyzz_c) +
        gupyy*(TWO*(l_Gamxyy*gxyy_c+l_Gamyyy*gyyy_c+l_Gamzyy*gyzy_c)+l_Gamxyy*gyyx_c+l_Gamyyy*gyyy_c+l_Gamzyy*gyyz_c) +
        gupyz*(TWO*(l_Gamxyy*gxzy_c+l_Gamyyy*gyzy_c+l_Gamzyy*gzzy_c+l_Gamxyz*gxyy_c+l_Gamyyz*gyyy_c+l_Gamzyz*gyzy_c)+l_Gamxyz*gyyx_c+l_Gamyyz*gyyy_c+l_Gamzyz*gyyz_c+l_Gamxyy*gyzx_c+l_Gamyyy*gyzy_c+l_Gamzyy*gyzz_c) +
        gupzz*(TWO*(l_Gamxyz*gxzy_c+l_Gamyyz*gyzy_c+l_Gamzyz*gzzy_c)+l_Gamxyz*gyzx_c+l_Gamyyz*gyzy_c+l_Gamzyz*gyzz_c);

    // Rzz (same as Cartesian)
    l_Rzz = -HALF * l_Rzz + gxz_v * Gamxz + gyz_v * Gamyz + gzz_v * Gamzz + Gamxa * gxzz_c + Gamya * gyzz_c + Gamza * gzzz_c +
        gupxx*(TWO*(l_Gamxxz*gxxz_c+l_Gamyxz*gxyz_c+l_Gamzxz*gxzz_c)+l_Gamxxz*gxzx_c+l_Gamyxz*gxzy_c+l_Gamzxz*gxzz_c) +
        gupxy*(TWO*(l_Gamxxz*gxyz_c+l_Gamyxz*gyyz_c+l_Gamzxz*gyzz_c+l_Gamxyz*gxxz_c+l_Gamyyz*gxyz_c+l_Gamzyz*gxzz_c)+l_Gamxyz*gxzx_c+l_Gamyyz*gxzy_c+l_Gamzyz*gxzz_c+l_Gamxxz*gyzx_c+l_Gamyxz*gyzy_c+l_Gamzxz*gyzz_c) +
        gupxz*(TWO*(l_Gamxxz*gxzz_c+l_Gamyxz*gyzz_c+l_Gamzxz*gzzz_c+l_Gamxzz*gxxz_c+l_Gamyzz*gxyz_c+l_Gamzzz*gxzz_c)+l_Gamxzz*gxzx_c+l_Gamyzz*gxzy_c+l_Gamzzz*gxzz_c+l_Gamxxz*gzzx_c+l_Gamyxz*gzzy_c+l_Gamzxz*gzzz_c) +
        gupyy*(TWO*(l_Gamxyz*gxyz_c+l_Gamyyz*gyyz_c+l_Gamzyz*gyzz_c)+l_Gamxyz*gyzx_c+l_Gamyyz*gyzy_c+l_Gamzyz*gyzz_c) +
        gupyz*(TWO*(l_Gamxyz*gxzz_c+l_Gamyyz*gyzz_c+l_Gamzyz*gzzz_c+l_Gamxzz*gxyz_c+l_Gamyzz*gyyz_c+l_Gamzzz*gyzz_c)+l_Gamxzz*gyzx_c+l_Gamyzz*gyzy_c+l_Gamzzz*gyzz_c+l_Gamxyz*gzzx_c+l_Gamyyz*gzzy_c+l_Gamzyz*gzzz_c) +
        gupzz*(TWO*(l_Gamxzz*gxzz_c+l_Gamyzz*gyzz_c+l_Gamzzz*gzzz_c)+l_Gamxzz*gzzx_c+l_Gamyzz*gzzy_c+l_Gamzzz*gzzz_c);

    // Rxy (same as Cartesian)
    l_Rxy = HALF*(-l_Rxy + gxx_v*Gamxy+gxy_v*Gamyy+gxz_v*Gamzy+gxy_v*Gamxx+gyy_v*Gamyx+gyz_v*Gamzx+Gamxa*gxyx_c+Gamya*gyyx_c+Gamza*gyzx_c+Gamxa*gxxy_c+Gamya*gxyy_c+Gamza*gxzy_c) +
        gupxx*(l_Gamxxx*gxxy_c+l_Gamyxx*gxyy_c+l_Gamzxx*gxzy_c+l_Gamxxy*gxxx_c+l_Gamyxy*gxyx_c+l_Gamzxy*gxzx_c+l_Gamxxx*gxyx_c+l_Gamyxx*gxyy_c+l_Gamzxx*gxyz_c) +
        gupxy*(l_Gamxxx*gxyy_c+l_Gamyxx*gyyy_c+l_Gamzxx*gyzy_c+l_Gamxxy*gxyx_c+l_Gamyxy*gyyx_c+l_Gamzxy*gyzx_c+l_Gamxxy*gxyx_c+l_Gamyxy*gxyy_c+l_Gamzxy*gxyz_c+l_Gamxxy*gxxy_c+l_Gamyxy*gxyy_c+l_Gamzxy*gxzy_c+l_Gamxyy*gxxx_c+l_Gamyyy*gxyx_c+l_Gamzyy*gxzx_c+l_Gamxxx*gyyx_c+l_Gamyxx*gyyy_c+l_Gamzxx*gyyz_c) +
        gupxz*(l_Gamxxx*gxzy_c+l_Gamyxx*gyzy_c+l_Gamzxx*gzzy_c+l_Gamxxy*gxzx_c+l_Gamyxy*gyzx_c+l_Gamzxy*gzzx_c+l_Gamxxz*gxyx_c+l_Gamyxz*gxyy_c+l_Gamzxz*gxyz_c+l_Gamxxz*gxxy_c+l_Gamyxz*gxyy_c+l_Gamzxz*gxzy_c+l_Gamxyz*gxxx_c+l_Gamyyz*gxyx_c+l_Gamzyz*gxzx_c+l_Gamxxx*gyzx_c+l_Gamyxx*gyzy_c+l_Gamzxx*gyzz_c) +
        gupyy*(l_Gamxxy*gxyy_c+l_Gamyxy*gyyy_c+l_Gamzxy*gyzy_c+l_Gamxyy*gxyx_c+l_Gamyyy*gyyx_c+l_Gamzyy*gyzx_c+l_Gamxxy*gyyx_c+l_Gamyxy*gyyy_c+l_Gamzxy*gyyz_c) +
        gupyz*(l_Gamxxy*gxzy_c+l_Gamyxy*gyzy_c+l_Gamzxy*gzzy_c+l_Gamxyy*gxzx_c+l_Gamyyy*gyzx_c+l_Gamzyy*gzzx_c+l_Gamxxz*gyyx_c+l_Gamyxz*gyyy_c+l_Gamzxz*gyyz_c+l_Gamxxz*gxyy_c+l_Gamyxz*gyyy_c+l_Gamzxz*gyzy_c+l_Gamxyz*gxyx_c+l_Gamyyz*gyyx_c+l_Gamzyz*gyzx_c+l_Gamxxy*gyzx_c+l_Gamyxy*gyzy_c+l_Gamzxy*gyzz_c) +
        gupzz*(l_Gamxxz*gxzy_c+l_Gamyxz*gyzy_c+l_Gamzxz*gzzy_c+l_Gamxyz*gxzx_c+l_Gamyyz*gyzx_c+l_Gamzyz*gzzx_c+l_Gamxxz*gyzx_c+l_Gamyxz*gyzy_c+l_Gamzxz*gyzz_c);

    // Rxz (same as Cartesian)
    l_Rxz = HALF*(-l_Rxz + gxx_v*Gamxz+gxy_v*Gamyz+gxz_v*Gamzz+gxz_v*Gamxx+gyz_v*Gamyx+gzz_v*Gamzx+Gamxa*gxzx_c+Gamya*gyzx_c+Gamza*gzzx_c+Gamxa*gxxz_c+Gamya*gxyz_c+Gamza*gxzz_c) +
        gupxx*(l_Gamxxx*gxxz_c+l_Gamyxx*gxyz_c+l_Gamzxx*gxzz_c+l_Gamxxz*gxxx_c+l_Gamyxz*gxyx_c+l_Gamzxz*gxzx_c+l_Gamxxx*gxzx_c+l_Gamyxx*gxzy_c+l_Gamzxx*gxzz_c) +
        gupxy*(l_Gamxxx*gxyz_c+l_Gamyxx*gyyz_c+l_Gamzxx*gyzz_c+l_Gamxxz*gxyx_c+l_Gamyxz*gyyx_c+l_Gamzxz*gyzx_c+l_Gamxxy*gxzx_c+l_Gamyxy*gxzy_c+l_Gamzxy*gxzz_c+l_Gamxxy*gxxz_c+l_Gamyxy*gxyz_c+l_Gamzxy*gxzz_c+l_Gamxyz*gxxx_c+l_Gamyyz*gxyx_c+l_Gamzyz*gxzx_c+l_Gamxxx*gyzx_c+l_Gamyxx*gyzy_c+l_Gamzxx*gyzz_c) +
        gupxz*(l_Gamxxx*gxzz_c+l_Gamyxx*gyzz_c+l_Gamzxx*gzzz_c+l_Gamxxz*gxzx_c+l_Gamyxz*gyzx_c+l_Gamzxz*gzzx_c+l_Gamxxz*gxzx_c+l_Gamyxz*gxzy_c+l_Gamzxz*gxzz_c+l_Gamxxz*gxxz_c+l_Gamyxz*gxyz_c+l_Gamzxz*gxzz_c+l_Gamxzz*gxxx_c+l_Gamyzz*gxyx_c+l_Gamzzz*gxzx_c+l_Gamxxx*gzzx_c+l_Gamyxx*gzzy_c+l_Gamzxx*gzzz_c) +
        gupyy*(l_Gamxxy*gxyz_c+l_Gamyxy*gyyz_c+l_Gamzxy*gyzz_c+l_Gamxyz*gxyx_c+l_Gamyyz*gyyx_c+l_Gamzyz*gyzx_c+l_Gamxxy*gyzx_c+l_Gamyxy*gyzy_c+l_Gamzxy*gyzz_c) +
        gupyz*(l_Gamxxy*gxzz_c+l_Gamyxy*gyzz_c+l_Gamzxy*gzzz_c+l_Gamxyz*gxzx_c+l_Gamyyz*gyzx_c+l_Gamzyz*gzzx_c+l_Gamxxz*gyzx_c+l_Gamyxz*gyzy_c+l_Gamzxz*gyzz_c+l_Gamxxz*gxyz_c+l_Gamyxz*gyyz_c+l_Gamzxz*gyzz_c+l_Gamxzz*gxyx_c+l_Gamyzz*gyyx_c+l_Gamzzz*gyzx_c+l_Gamxxy*gzzx_c+l_Gamyxy*gzzy_c+l_Gamzxy*gzzz_c) +
        gupzz*(l_Gamxxz*gxzz_c+l_Gamyxz*gyzz_c+l_Gamzxz*gzzz_c+l_Gamxzz*gxzx_c+l_Gamyzz*gyzx_c+l_Gamzzz*gzzx_c+l_Gamxxz*gzzx_c+l_Gamyxz*gzzy_c+l_Gamzxz*gzzz_c);

    // Ryz (same as Cartesian)
    l_Ryz = HALF*(-l_Ryz + gxy_v*Gamxz+gyy_v*Gamyz+gyz_v*Gamzz+gxz_v*Gamxy+gyz_v*Gamyy+gzz_v*Gamzy+Gamxa*gxzy_c+Gamya*gyzy_c+Gamza*gzzy_c+Gamxa*gxyz_c+Gamya*gyyz_c+Gamza*gyzz_c) +
        gupxx*(l_Gamxxy*gxxz_c+l_Gamyxy*gxyz_c+l_Gamzxy*gxzz_c+l_Gamxxz*gxxy_c+l_Gamyxz*gxyy_c+l_Gamzxz*gxzy_c+l_Gamxxy*gxzx_c+l_Gamyxy*gxzy_c+l_Gamzxy*gxzz_c) +
        gupxy*(l_Gamxxy*gxyz_c+l_Gamyxy*gyyz_c+l_Gamzxy*gyzz_c+l_Gamxxz*gxyy_c+l_Gamyxz*gyyy_c+l_Gamzxz*gyzy_c+l_Gamxyy*gxzx_c+l_Gamyyy*gxzy_c+l_Gamzyy*gxzz_c+l_Gamxyy*gxxz_c+l_Gamyyy*gxyz_c+l_Gamzyy*gxzz_c+l_Gamxyz*gxxy_c+l_Gamyyz*gxyy_c+l_Gamzyz*gxzy_c+l_Gamxxy*gyzx_c+l_Gamyxy*gyzy_c+l_Gamzxy*gyzz_c) +
        gupxz*(l_Gamxxy*gxzz_c+l_Gamyxy*gyzz_c+l_Gamzxy*gzzz_c+l_Gamxxz*gxzy_c+l_Gamyxz*gyzy_c+l_Gamzxz*gzzy_c+l_Gamxyz*gxzx_c+l_Gamyyz*gxzy_c+l_Gamzyz*gxzz_c+l_Gamxyz*gxxz_c+l_Gamyyz*gxyz_c+l_Gamzyz*gxzz_c+l_Gamxzz*gxxy_c+l_Gamyzz*gxyy_c+l_Gamzzz*gxzy_c+l_Gamxxy*gzzx_c+l_Gamyxy*gzzy_c+l_Gamzxy*gzzz_c) +
        gupyy*(l_Gamxyy*gxyz_c+l_Gamyyy*gyyz_c+l_Gamzyy*gyzz_c+l_Gamxyz*gxyy_c+l_Gamyyz*gyyy_c+l_Gamzyz*gyzy_c+l_Gamxyy*gyzx_c+l_Gamyyy*gyzy_c+l_Gamzyy*gyzz_c) +
        gupyz*(l_Gamxyy*gxzz_c+l_Gamyyy*gyzz_c+l_Gamzyy*gzzz_c+l_Gamxyz*gxzy_c+l_Gamyyz*gyzy_c+l_Gamzyz*gzzy_c+l_Gamxyz*gyzx_c+l_Gamyyz*gyzy_c+l_Gamzyz*gyzz_c+l_Gamxyz*gxyz_c+l_Gamyyz*gyyy_c+l_Gamzyz*gyzz_c+l_Gamxzz*gxyy_c+l_Gamyzz*gyyy_c+l_Gamzzz*gyzy_c+l_Gamxyy*gzzx_c+l_Gamyyy*gzzy_c+l_Gamzyy*gzzz_c) +
        gupzz*(l_Gamxyz*gxzz_c+l_Gamyyz*gyzz_c+l_Gamzyz*gzzz_c+l_Gamxzz*gxzy_c+l_Gamyzz*gyzy_c+l_Gamzzz*gzzy_c+l_Gamxyz*gzzx_c+l_Gamyyz*gzzy_c+l_Gamzyz*gzzz_c);

    // Chi contribution (shell 2nd derivative)
    d_fdderivs_sh_point(ex, pex, chi, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, sst, i, j, k);
    fxx = fxx - l_Gamxxx*chix - l_Gamyxx*chiy - l_Gamzxx*chiz;
    fxy = fxy - l_Gamxxy*chix - l_Gamyxy*chiy - l_Gamzxy*chiz;
    fxz = fxz - l_Gamxxz*chix - l_Gamyxz*chiy - l_Gamzxz*chiz;
    fyy = fyy - l_Gamxyy*chix - l_Gamyyy*chiy - l_Gamzyy*chiz;
    fyz = fyz - l_Gamxyz*chix - l_Gamyyz*chiy - l_Gamzyz*chiz;
    fzz = fzz - l_Gamxzz*chix - l_Gamyzz*chiy - l_Gamzzz*chiz;

    double f = gupxx*(fxx-F3o2/chin1*chix*chix) + gupyy*(fyy-F3o2/chin1*chiy*chiy) + gupzz*(fzz-F3o2/chin1*chiz*chiz) +
        TWO*(gupxy*(fxy-F3o2/chin1*chix*chiy) + gupxz*(fxz-F3o2/chin1*chix*chiz) + gupyz*(fyz-F3o2/chin1*chiy*chiz));

    l_Rxx += (fxx-chix*chix/chin1/TWO+gxx_v*f)/chin1/TWO;
    l_Ryy += (fyy-chiy*chiy/chin1/TWO+gyy_v*f)/chin1/TWO;
    l_Rzz += (fzz-chiz*chiz/chin1/TWO+gzz_v*f)/chin1/TWO;
    l_Rxy += (fxy-chix*chiy/chin1/TWO+gxy_v*f)/chin1/TWO;
    l_Rxz += (fxz-chix*chiz/chin1/TWO+gxz_v*f)/chin1/TWO;
    l_Ryz += (fyz-chiy*chiz/chin1/TWO+gyz_v*f)/chin1/TWO;

    // Physical Christoffel symbols
    double gxxx_ch = (gupxx*chix+gupxy*chiy+gupxz*chiz)/chin1;
    double gxxy_ch = (gupxy*chix+gupyy*chiy+gupyz*chiz)/chin1;
    double gxxz_ch = (gupxz*chix+gupyz*chiy+gupzz*chiz)/chin1;

    *Gamxxx = l_Gamxxx - ((chix+chix)/chin1 - gxx_v*gxxx_ch)*HALF;
    *Gamyxx = l_Gamyxx - (                     - gxx_v*gxxy_ch)*HALF;
    *Gamzxx = l_Gamzxx - (                     - gxx_v*gxxz_ch)*HALF;
    *Gamxyy = l_Gamxyy - (                     - gyy_v*gxxx_ch)*HALF;
    *Gamyyy = l_Gamyyy - ((chiy+chiy)/chin1 - gyy_v*gxxy_ch)*HALF;
    *Gamzyy = l_Gamzyy - (                     - gyy_v*gxxz_ch)*HALF;
    *Gamxzz = l_Gamxzz - (                     - gzz_v*gxxx_ch)*HALF;
    *Gamyzz = l_Gamyzz - (                     - gzz_v*gxxy_ch)*HALF;
    *Gamzzz = l_Gamzzz - ((chiz+chiz)/chin1 - gzz_v*gxxz_ch)*HALF;
    *Gamxxy = l_Gamxxy - ( chiy       /chin1 - gxy_v*gxxx_ch)*HALF;
    *Gamyxy = l_Gamyxy - (        chix/chin1 - gxy_v*gxxy_ch)*HALF;
    *Gamzxy = l_Gamzxy - (                     - gxy_v*gxxz_ch)*HALF;
    *Gamxxz = l_Gamxxz - ( chiz       /chin1 - gxz_v*gxxx_ch)*HALF;
    *Gamyxz = l_Gamyxz - (                     - gxz_v*gxxy_ch)*HALF;
    *Gamzxz = l_Gamzxz - (        chix/chin1 - gxz_v*gxxz_ch)*HALF;
    *Gamxyz = l_Gamxyz - (                     - gyz_v*gxxx_ch)*HALF;
    *Gamyyz = l_Gamyyz - ( chiz       /chin1 - gyz_v*gxxy_ch)*HALF;
    *Gamzyz = l_Gamzyz - (        chiy/chin1 - gyz_v*gxxz_ch)*HALF;

    *Rxx = l_Rxx; *Ryy = l_Ryy; *Rzz = l_Rzz;
    *Rxy = l_Rxy; *Rxz = l_Rxz; *Ryz = l_Ryz;
}

// ============================================================================
// Global kernel: ricci_gamma_ss (shell coordinate)
// ============================================================================
__global__ void ricci_gamma_ss_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    double* __restrict__ Gamxxx, double* __restrict__ Gamxxy, double* __restrict__ Gamxxz,
    double* __restrict__ Gamxyy, double* __restrict__ Gamxyz, double* __restrict__ Gamxzz,
    double* __restrict__ Gamyxx, double* __restrict__ Gamyxy, double* __restrict__ Gamyxz,
    double* __restrict__ Gamyyy, double* __restrict__ Gamyyz, double* __restrict__ Gamyzz,
    double* __restrict__ Gamzxx, double* __restrict__ Gamzxy, double* __restrict__ Gamzxz,
    double* __restrict__ Gamzyy, double* __restrict__ Gamzyz, double* __restrict__ Gamzzz,
    double* __restrict__ Rxx, double* __restrict__ Rxy, double* __restrict__ Rxz,
    double* __restrict__ Ryy, double* __restrict__ Ryz, double* __restrict__ Rzz,
    int symmetry, int sst
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total = ex0 * ex1 * ex2;
    if (tid >= total) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = (k * pex1 + j) * pex0 + i;

    int ex[3]  = {ex0, ex1, ex2};
    int pex[3] = {pex0, pex1, pex2};

    double l_Gamxxx, l_Gamxxy, l_Gamxxz, l_Gamxyy, l_Gamxyz, l_Gamxzz;
    double l_Gamyxx, l_Gamyxy, l_Gamyxz, l_Gamyyy, l_Gamyyz, l_Gamyzz;
    double l_Gamzxx, l_Gamzxy, l_Gamzxz, l_Gamzyy, l_Gamzyz, l_Gamzzz;
    double l_Rxx, l_Rxy, l_Rxz, l_Ryy, l_Ryz, l_Rzz;

    d_ricci_gamma_ss_point(ex, pex, X, Y, Z, chi, dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz,
        &l_Gamxxx, &l_Gamxxy, &l_Gamxxz, &l_Gamxyy, &l_Gamxyz, &l_Gamxzz,
        &l_Gamyxx, &l_Gamyxy, &l_Gamyxz, &l_Gamyyy, &l_Gamyyz, &l_Gamyzz,
        &l_Gamzxx, &l_Gamzxy, &l_Gamzxz, &l_Gamzyy, &l_Gamzyz, &l_Gamzzz,
        &l_Rxx, &l_Rxy, &l_Rxz, &l_Ryy, &l_Ryz, &l_Rzz,
        symmetry, sst, i, j, k);

    Gamxxx[idx]=l_Gamxxx; Gamxxy[idx]=l_Gamxxy; Gamxxz[idx]=l_Gamxxz;
    Gamxyy[idx]=l_Gamxyy; Gamxyz[idx]=l_Gamxyz; Gamxzz[idx]=l_Gamxzz;
    Gamyxx[idx]=l_Gamyxx; Gamyxy[idx]=l_Gamyxy; Gamyxz[idx]=l_Gamyxz;
    Gamyyy[idx]=l_Gamyyy; Gamyyz[idx]=l_Gamyyz; Gamyzz[idx]=l_Gamyzz;
    Gamzxx[idx]=l_Gamzxx; Gamzxy[idx]=l_Gamzxy; Gamzxz[idx]=l_Gamzxz;
    Gamzyy[idx]=l_Gamzyy; Gamzyz[idx]=l_Gamzyz; Gamzzz[idx]=l_Gamzzz;
    Rxx[idx]=l_Rxx; Rxy[idx]=l_Rxy; Rxz[idx]=l_Rxz;
    Ryy[idx]=l_Ryy; Ryz[idx]=l_Ryz; Rzz[idx]=l_Rzz;
}

// ============================================================================
// Host launch wrapper: ricci_gamma_ss
// ============================================================================
void gpu_ricci_gamma_ss_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Gamx, const double* d_Gamy, const double* d_Gamz,
    double* d_Gamxxx, double* d_Gamxxy, double* d_Gamxxz,
    double* d_Gamxyy, double* d_Gamxyz, double* d_Gamxzz,
    double* d_Gamyxx, double* d_Gamyxy, double* d_Gamyxz,
    double* d_Gamyyy, double* d_Gamyyz, double* d_Gamyzz,
    double* d_Gamzxx, double* d_Gamzxy, double* d_Gamzxz,
    double* d_Gamzyy, double* d_Gamzyz, double* d_Gamzzz,
    double* d_Rxx, double* d_Rxy, double* d_Rxz,
    double* d_Ryy, double* d_Ryz, double* d_Rzz,
    int symmetry, int sst
) {
    int total = ex[0] * ex[1] * ex[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);

    ricci_gamma_ss_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2],
        pex[0], pex[1], pex[2],
        d_X, d_Y, d_Z,
        d_chi,
        d_dxx, d_gxy, d_gxz,
        d_dyy, d_gyz, d_dzz,
        d_Gamx, d_Gamy, d_Gamz,
        d_Gamxxx, d_Gamxxy, d_Gamxxz,
        d_Gamxyy, d_Gamxyz, d_Gamxzz,
        d_Gamyxx, d_Gamyxy, d_Gamyxz,
        d_Gamyyy, d_Gamyyz, d_Gamyzz,
        d_Gamzxx, d_Gamzxy, d_Gamzxz,
        d_Gamzyy, d_Gamzyz, d_Gamzzz,
        d_Rxx, d_Rxy, d_Rxz,
        d_Ryy, d_Ryz, d_Rzz,
        symmetry, sst
    );
}
