#include "macrodef.h"
#include "empart.h"
#include "diff_new_gpu.cuh"
#include "point_diff_new_sh_gpu.cuh"
#include "lopsidediff_gpu.cuh"
#include "kodiss_gpu.cuh"
#include "kodiss_sh_gpu.cuh"
#include "rhs_smem_common.cuh"
#include "gpu_manager.h"
#include "gpu_common.cuh"
#include "gpu_info.h"
#include "macrodef.fh"
#include <cuda_runtime.h>
#include <math_constants.h>
#include <cmath>

#ifndef IDX3D_PAD
#define IDX3D_PAD(i, j, k, pnx, pny) ((i) + (pnx) * ((j) + (pny) * (k)))
#endif

namespace {
    constexpr double SYM  =  1.0;
    constexpr double ANTI = -1.0;
    constexpr double ONE  =  1.0;
    constexpr double TWO  =  2.0;
    constexpr double FOUR =  4.0;
    constexpr double F3o2 =  1.5;
    constexpr double KAPPA = 1.0;
    constexpr int ONOFF = 0;
}

// ============================================================================
// Cartesian kernel: compute_rhs_empart
// SMEM-tiled first-derivative load + per-point physics assembly.
// Replaces the previous per-point d_fderivs_point_gw model — every derivative
// is now reused from shared memory tiles loaded by the whole block in 7
// batches (cf. bssn_ricci_kernel). Ghost zones must be pre-mirrored by
// gpu_apply_symmetry_bd_empart_launch (see bssnEM_gpu_class::Step_GPU) so
// that the smem load can read them directly.
// ============================================================================
template<int GW>
__launch_bounds__(256, 1)
__global__ void empart_rhs_kernel(
    const int ex0, const int ex1, const int ex2,
    const int pex0, const int pex1, const int pex2,
    const int symmetry, const int lev,
    const double eps,
    const double* __restrict__ X,
    const double* __restrict__ Y,
    const double* __restrict__ Z,
    const double* __restrict__ chi,
    const double* __restrict__ Jx,
    const double* __restrict__ Jy,
    const double* __restrict__ Jz,
    const double* __restrict__ qchar,
    const double* __restrict__ dxx,
    const double* __restrict__ dxy,
    const double* __restrict__ dxz,
    const double* __restrict__ dyy,
    const double* __restrict__ dyz,
    const double* __restrict__ dzz,
    const double* __restrict__ Lap,
    const double* __restrict__ betax,
    const double* __restrict__ betay,
    const double* __restrict__ betaz,
    const double* __restrict__ trK,
    const double* __restrict__ Ex,
    const double* __restrict__ Ey,
    const double* __restrict__ Ez,
    const double* __restrict__ Bx,
    const double* __restrict__ By,
    const double* __restrict__ Bz,
    const double* __restrict__ Kpsi,
    const double* __restrict__ Kphi,
    double* __restrict__ Ex_rhs,
    double* __restrict__ Ey_rhs,
    double* __restrict__ Ez_rhs,
    double* __restrict__ Bx_rhs,
    double* __restrict__ By_rhs,
    double* __restrict__ Bz_rhs,
    double* __restrict__ Kpsi_rhs,
    double* __restrict__ Kphi_rhs,
    double* __restrict__ rho,
    double* __restrict__ Sx,
    double* __restrict__ Sy,
    double* __restrict__ Sz,
    double* __restrict__ Sxx,
    double* __restrict__ Sxy,
    double* __restrict__ Sxz,
    double* __restrict__ Syy,
    double* __restrict__ Syz,
    double* __restrict__ Szz
) {
    // ---- Thread index ----
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;

    // ---- ex/pex as arrays for derivative function calls ----
    const int ex[3]  = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    // All threads in the block must participate in the cooperative SMEM
    // loads below, so we *don't* early-return out-of-bounds threads yet —
    // we only mark them invalid and skip the physics writes at the end.
    const bool is_valid = (i < ex0 && j < ex1 && k < ex2);
    const int idx = is_valid ? IDX3D_PAD(i, j, k, pex0, pex1) : 0;

    // ---- Derivative bookkeeping (mirrors BSSN ricci kernel) ----
    const int imax = ex0 - 1;
    const int jmax = ex1 - 1;
    const int kmax = ex2 - 1;
    int imin = 0, jmin = 0, kmin = 0;
    constexpr int ORD = GW - 1;
    const double dX_d = X[1] - X[0];
    const double dY_d = Y[1] - Y[0];
    const double dZ_d = Z[1] - Z[0];
    if (symmetry > 0 && fabs(Z[0]) < dZ_d) kmin = -ORD;
    if (symmetry > 1 && fabs(X[0]) < dX_d) imin = -ORD;
    if (symmetry > 1 && fabs(Y[0]) < dY_d) jmin = -ORD;

    const int block_i = blockIdx.x * blockDim.x;
    const int block_j = blockIdx.y * blockDim.y;
    const int block_k = blockIdx.z * blockDim.z;

    // ---- Constants ----
    constexpr double PI = CUDART_PI;
    constexpr double FOUR_PI = 4.0 * CUDART_PI;
    constexpr double EIT = 8.0;

    // ---- Load local field values (safe at idx=0 if invalid; result unused) ----
    const double val_chi  = chi[idx];
    const double val_Jx   = Jx[idx];
    const double val_Jy   = Jy[idx];
    const double val_Jz   = Jz[idx];
    const double val_qchar = qchar[idx];

    const double val_dxx = dxx[idx];
    const double val_dxy = dxy[idx];
    const double val_dxz = dxz[idx];
    const double val_dyy = dyy[idx];
    const double val_dyz = dyz[idx];
    const double val_dzz = dzz[idx];

    const double val_Lap   = Lap[idx];
    const double val_trK   = trK[idx];

    const double val_Ex = Ex[idx];
    const double val_Ey = Ey[idx];
    const double val_Ez = Ez[idx];
    const double val_Bx = Bx[idx];
    const double val_By = By[idx];
    const double val_Bz = Bz[idx];
    const double val_Kpsi = Kpsi[idx];
    const double val_Kphi = Kphi[idx];

    // ---- alpn1 and chin1 (point-wise) ----
    const double alpn1 = val_Lap + ONE;
    // Clamp chin1 to prevent sqrt(negative) -> NaN and 1/0 -> inf when
    // chi is near the puncture (chitiny safeguard matching CPU lowerboundset).
    // On CPU, f_compute_rhs_bssn calls lowerboundset(chi, chitiny) before
    // the RHS evaluation; the GPU EM kernel runs independently and must
    // apply the same safeguard.  chitiny = -1 + delta, so chin1 >= delta > 0.
    // Match Fortran: chin1 = chi + ONE, no clamping (CPU applies lowerboundset before RHS)
    const double chin1 = val_chi + ONE;
    const double chi3o2 = chin1 * sqrt(chin1);
    const double inv_chin1 = 1.0 / chin1;

    // ---- Conformal metric components ----
    double gxx = val_dxx + ONE;
    double gyy = val_dyy + ONE;
    double gzz = val_dzz + ONE;
    double gxy = val_dxy;
    double gxz = val_dxz;
    double gyz = val_dyz;

    // ======================================================================
    // Stage 3: First derivatives of the 19 fields — collaborative SMEM tile.
    // Layout: 3 dynamic-smem tiles, reused across 7 load batches.
    // ======================================================================
    extern __shared__ double smem_raw[];
    double (*smem3)[TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X] =
        reinterpret_cast<double(*)[TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X]>(smem_raw);

    // ---- Batch 1: Lap, chi, betax ----
    load_3_fields_to_smem<GW>(smem3, Lap, chi, betax, pex, block_i, block_j, block_k);
    __syncthreads();
    double Lapx, Lapy, Lapz;
    compute_first_derivs_smem<GW>(smem3[0], &Lapx, &Lapy, &Lapz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double chix, chiy, chiz;
    compute_first_derivs_smem<GW>(smem3[1], &chix, &chiy, &chiz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double betaxx, betaxy, betaxz;
    compute_first_derivs_smem<GW>(smem3[2], &betaxx, &betaxy, &betaxz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    __syncthreads();

    // ---- Batch 2: betay, betaz, dxx ----
    load_3_fields_to_smem<GW>(smem3, betay, betaz, dxx, pex, block_i, block_j, block_k);
    __syncthreads();
    double betayx, betayy, betayz;
    compute_first_derivs_smem<GW>(smem3[0], &betayx, &betayy, &betayz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double betazx, betazy, betazz;
    compute_first_derivs_smem<GW>(smem3[1], &betazx, &betazy, &betazz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double gxxx, gxxy, gxxz;
    compute_first_derivs_smem<GW>(smem3[2], &gxxx, &gxxy, &gxxz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    __syncthreads();

    // ---- Batch 3: dxy, dxz, dyy ----
    load_3_fields_to_smem<GW>(smem3, dxy, dxz, dyy, pex, block_i, block_j, block_k);
    __syncthreads();
    double gxyx, gxyy, gxyz;
    compute_first_derivs_smem<GW>(smem3[0], &gxyx, &gxyy, &gxyz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double gxzx, gxzy, gxzz;
    compute_first_derivs_smem<GW>(smem3[1], &gxzx, &gxzy, &gxzz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double gyyx, gyyy, gyyz;
    compute_first_derivs_smem<GW>(smem3[2], &gyyx, &gyyy, &gyyz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    __syncthreads();

    // ---- Batch 4: dyz, dzz, Kpsi ----
    load_3_fields_to_smem<GW>(smem3, dyz, dzz, Kpsi, pex, block_i, block_j, block_k);
    __syncthreads();
    double gyzx, gyzy, gyzz;
    compute_first_derivs_smem<GW>(smem3[0], &gyzx, &gyzy, &gyzz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double gzzx, gzzy, gzzz;
    compute_first_derivs_smem<GW>(smem3[1], &gzzx, &gzzy, &gzzz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double Kpsix, Kpsiy, Kpsiz;
    compute_first_derivs_smem<GW>(smem3[2], &Kpsix, &Kpsiy, &Kpsiz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    __syncthreads();

    // ---- Batch 5: Kphi, Ex, Ey ----
    load_3_fields_to_smem<GW>(smem3, Kphi, Ex, Ey, pex, block_i, block_j, block_k);
    __syncthreads();
    double Kphix, Kphiy, Kphiz;
    compute_first_derivs_smem<GW>(smem3[0], &Kphix, &Kphiy, &Kphiz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double Exx, Exy, Exz;
    compute_first_derivs_smem<GW>(smem3[1], &Exx, &Exy, &Exz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double Eyx, Eyy, Eyz;
    compute_first_derivs_smem<GW>(smem3[2], &Eyx, &Eyy, &Eyz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    __syncthreads();

    // ---- Batch 6: Ez, Bx, By ----
    load_3_fields_to_smem<GW>(smem3, Ez, Bx, By, pex, block_i, block_j, block_k);
    __syncthreads();
    double Ezx, Ezy, Ezz;
    compute_first_derivs_smem<GW>(smem3[0], &Ezx, &Ezy, &Ezz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double Bxx, Bxy, Bxz;
    compute_first_derivs_smem<GW>(smem3[1], &Bxx, &Bxy, &Bxz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    double Byx, Byy, Byz;
    compute_first_derivs_smem<GW>(smem3[2], &Byx, &Byy, &Byz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    __syncthreads();

    // ---- Batch 7: Bz alone ----
    load_field_to_smem<GW>(smem3[0], Bz, ex, pex, block_i, block_j, block_k, ANTI, ANTI, SYM);
    __syncthreads();
    double Bzx, Bzy, Bzz;
    compute_first_derivs_smem<GW>(smem3[0], &Bzx, &Bzy, &Bzz, dX_d, dY_d, dZ_d, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    // No __syncthreads() needed — smem is no longer touched.

    // From this point on no thread re-reads the tile; bail out the OOB ones.
    if (!is_valid) return;

    // ======================================================================
    // Stage 4: Physical metric reconstruction
    // ======================================================================

    // Physical metric: gij_phys = gij_conf / chin1
    gxx /= chin1;
    gxy /= chin1;
    gxz /= chin1;
    gyy /= chin1;
    gyz /= chin1;
    gzz /= chin1;

    // Physical metric derivatives: gij,k = (gij,k_conf - gij_phys * chi,k) / chin1
    gxxx = (gxxx - gxx * chix) * inv_chin1;
    gxxy = (gxxy - gxx * chiy) * inv_chin1;
    gxxz = (gxxz - gxx * chiz) * inv_chin1;

    gxyx = (gxyx - gxy * chix) * inv_chin1;
    gxyy = (gxyy - gxy * chiy) * inv_chin1;
    gxyz = (gxyz - gxy * chiz) * inv_chin1;

    gxzx = (gxzx - gxz * chix) * inv_chin1;
    gxzy = (gxzy - gxz * chiy) * inv_chin1;
    gxzz = (gxzz - gxz * chiz) * inv_chin1;

    gyyx = (gyyx - gyy * chix) * inv_chin1;
    gyyy = (gyyy - gyy * chiy) * inv_chin1;
    gyyz = (gyyz - gyy * chiz) * inv_chin1;

    gyzx = (gyzx - gyz * chix) * inv_chin1;
    gyzy = (gyzy - gyz * chiy) * inv_chin1;
    gyzz = (gyzz - gyz * chiz) * inv_chin1;

    gzzx = (gzzx - gzz * chix) * inv_chin1;
    gzzy = (gzzy - gzz * chiy) * inv_chin1;
    gzzz = (gzzz - gzz * chiz) * inv_chin1;

    // Inverse physical metric via Cramer's rule
    const double det = gxx * gyy * gzz
        + gxy * gyz * gxz
        + gxz * gxy * gyz
        - gxz * gyy * gxz
        - gxy * gxy * gzz
        - gxx * gyz * gyz;
    // Match Fortran: direct inverse, no fabs/fmax clamping
    const double inv_det = 1.0 / det;

    const double gupxx = (gyy * gzz - gyz * gyz) * inv_det;
    const double gupxy = -(gxy * gzz - gyz * gxz) * inv_det;
    const double gupxz = (gxy * gyz - gyy * gxz) * inv_det;
    const double gupyy = (gxx * gzz - gxz * gxz) * inv_det;
    const double gupyz = -(gxx * gyz - gxy * gxz) * inv_det;
    const double gupzz = (gxx * gyy - gxy * gxy) * inv_det;

    // ======================================================================
    // Stage 5: RHS source computation
    // ======================================================================

    double val_Ex_rhs, val_Ey_rhs, val_Ez_rhs;
    double val_Bx_rhs, val_By_rhs, val_Bz_rhs;
    double val_Kpsi_rhs, val_Kphi_rhs;

    // --- Ex_rhs ---
    val_Ex_rhs = alpn1 * val_trK * val_Ex
        - (val_Ex * betaxx + val_Ey * betaxy + val_Ez * betaxz)
        - FOUR_PI * alpn1 * val_Jx
        - alpn1 * (gupxx * Kpsix + gupxy * Kpsiy + gupxz * Kpsiz)
        + chi3o2 * (
            ((gxz * val_Bx + gyz * val_By + gzz * val_Bz) * Lapy
                + alpn1 * (gxz * Bxy + gyz * Byy + gzz * Bzy)
                + alpn1 * (val_Bx * gxzy + val_By * gyzy + val_Bz * gzzy))
          - ((gxy * val_Bx + gyy * val_By + gyz * val_Bz) * Lapz
                + alpn1 * (gxy * Bxz + gyy * Byz + gyz * Bzz)
                + alpn1 * (val_Bx * gxyz + val_By * gyyz + val_Bz * gyzz))
        );

    // --- Ey_rhs ---
    val_Ey_rhs = alpn1 * val_trK * val_Ey
        - (val_Ex * betayx + val_Ey * betayy + val_Ez * betayz)
        - FOUR_PI * alpn1 * val_Jy
        - alpn1 * (gupxy * Kpsix + gupyy * Kpsiy + gupyz * Kpsiz)
        + chi3o2 * (
            ((gxx * val_Bx + gxy * val_By + gxz * val_Bz) * Lapz
                + alpn1 * (gxx * Bxz + gxy * Byz + gxz * Bzz)
                + alpn1 * (val_Bx * gxxz + val_By * gxyz + val_Bz * gxzz))
          - ((gxz * val_Bx + gyz * val_By + gzz * val_Bz) * Lapx
                + alpn1 * (gxz * Bxx + gyz * Byx + gzz * Bzx)
                + alpn1 * (val_Bx * gxzx + val_By * gyzx + val_Bz * gzzx))
        );

    // --- Ez_rhs ---
    val_Ez_rhs = alpn1 * val_trK * val_Ez
        - (val_Ex * betazx + val_Ey * betazy + val_Ez * betazz)
        - FOUR_PI * alpn1 * val_Jz
        - alpn1 * (gupxz * Kpsix + gupyz * Kpsiy + gupzz * Kpsiz)
        + chi3o2 * (
            ((gxy * val_Bx + gyy * val_By + gyz * val_Bz) * Lapx
                + alpn1 * (gxy * Bxx + gyy * Byx + gyz * Bzx)
                + alpn1 * (val_Bx * gxyx + val_By * gyyx + val_Bz * gyzx))
          - ((gxx * val_Bx + gxy * val_By + gxz * val_Bz) * Lapy
                + alpn1 * (gxx * Bxy + gxy * Byy + gxz * Bzy)
                + alpn1 * (val_Bx * gxxy + val_By * gxyy + val_Bz * gxzy))
        );

    // --- Bx_rhs ---
    val_Bx_rhs = alpn1 * val_trK * val_Bx
        - (val_Bx * betaxx + val_By * betaxy + val_Bz * betaxz)
        - alpn1 * (gupxx * Kphix + gupxy * Kphiy + gupxz * Kphiz)
        - chi3o2 * (
            ((gxz * val_Ex + gyz * val_Ey + gzz * val_Ez) * Lapy
                + alpn1 * (gxz * Exy + gyz * Eyy + gzz * Ezy)
                + alpn1 * (val_Ex * gxzy + val_Ey * gyzy + val_Ez * gzzy))
          - ((gxy * val_Ex + gyy * val_Ey + gyz * val_Ez) * Lapz
                + alpn1 * (gxy * Exz + gyy * Eyz + gyz * Ezz)
                + alpn1 * (val_Ex * gxyz + val_Ey * gyyz + val_Ez * gyzz))
        );

    // --- By_rhs ---
    val_By_rhs = alpn1 * val_trK * val_By
        - (val_Bx * betayx + val_By * betayy + val_Bz * betayz)
        - alpn1 * (gupxy * Kphix + gupyy * Kphiy + gupyz * Kphiz)
        - chi3o2 * (
            ((gxx * val_Ex + gxy * val_Ey + gxz * val_Ez) * Lapz
                + alpn1 * (gxx * Exz + gxy * Eyz + gxz * Ezz)
                + alpn1 * (val_Ex * gxxz + val_Ey * gxyz + val_Ez * gxzz))
          - ((gxz * val_Ex + gyz * val_Ey + gzz * val_Ez) * Lapx
                + alpn1 * (gxz * Exx + gyz * Eyx + gzz * Ezx)
                + alpn1 * (val_Ex * gxzx + val_Ey * gyzx + val_Ez * gzzx))
        );

    // --- Bz_rhs ---
    val_Bz_rhs = alpn1 * val_trK * val_Bz
        - (val_Bx * betazx + val_By * betazy + val_Bz * betazz)
        - alpn1 * (gupxz * Kphix + gupyz * Kphiy + gupzz * Kphiz)
        - chi3o2 * (
            ((gxy * val_Ex + gyy * val_Ey + gyz * val_Ez) * Lapx
                + alpn1 * (gxy * Exx + gyy * Eyx + gyz * Ezx)
                + alpn1 * (val_Ex * gxyx + val_Ey * gyyx + val_Ez * gyzx))
          - ((gxx * val_Ex + gxy * val_Ey + gxz * val_Ez) * Lapy
                + alpn1 * (gxx * Exy + gxy * Eyy + gxz * Ezy)
                + alpn1 * (val_Ex * gxxy + val_Ey * gxyy + val_Ez * gxzy))
        );

    // --- Kpsi_rhs ---
    val_Kpsi_rhs = FOUR * PI * alpn1 * val_qchar
        - alpn1 * KAPPA * val_Kpsi
        - alpn1 * (Exx + Eyy + Ezz
            - F3o2 * inv_chin1 * (chix * val_Ex + chiy * val_Ey + chiz * val_Ez));

    // --- Kphi_rhs ---
    val_Kphi_rhs = -alpn1 * KAPPA * val_Kphi
        - alpn1 * (Bxx + Byy + Bzz
            - F3o2 * inv_chin1 * (chix * val_Bx + chiy * val_By + chiz * val_Bz));

    // Stage 6 (lopsided advection + Kreiss-Oliger dissipation) has been
    // split out into empart_advection_dissipation_kernel — see below. The
    // current kernel only emits the "source" RHS for each field, the
    // separate kernel finishes the sum with the upwind/KO contribution.

    // ======================================================================
    // Stage 7: Stress-energy tensor
    // ======================================================================

    // --- rho ---
    const double val_rho = (
        gxx * (val_Ex * val_Ex + val_Bx * val_Bx)
        + gyy * (val_Ey * val_Ey + val_By * val_By)
        + gzz * (val_Ez * val_Ez + val_Bz * val_Bz)
        + TWO * (
            gxy * (val_Ex * val_Ey + val_Bx * val_By)
            + gxz * (val_Ex * val_Ez + val_Bx * val_Bz)
            + gyz * (val_Ey * val_Ez + val_By * val_Bz)
        )
    ) / (EIT * PI);

    // --- Sx, Sy, Sz (Poynting vector) ---
    const double val_Sx = (val_Ey * val_Bz - val_Ez * val_By) / (FOUR * PI * chi3o2);
    const double val_Sy = (val_Ez * val_Bx - val_Ex * val_Bz) / (FOUR * PI * chi3o2);
    const double val_Sz = (val_Ex * val_By - val_Ey * val_Bx) / (FOUR * PI * chi3o2);

    // --- Lowered-index electric field (physical metric) ---
    const double lEx = gxx * val_Ex + gxy * val_Ey + gxz * val_Ez;
    const double lEy = gxy * val_Ex + gyy * val_Ey + gyz * val_Ez;
    const double lEz = gxz * val_Ex + gyz * val_Ey + gzz * val_Ez;

    // --- Lowered-index magnetic field (physical metric) ---
    const double lBx = gxx * val_Bx + gxy * val_By + gxz * val_Bz;
    const double lBy = gxy * val_Bx + gyy * val_By + gyz * val_Bz;
    const double lBz = gxz * val_Bx + gyz * val_By + gzz * val_Bz;

    // --- Stress tensor S_ij ---
    const double val_Sxx = val_rho * gxx - (lEx * lEx + lBx * lBx) / (FOUR * PI);
    const double val_Sxy = val_rho * gxy - (lEx * lEy + lBx * lBy) / (FOUR * PI);
    const double val_Sxz = val_rho * gxz - (lEx * lEz + lBx * lBz) / (FOUR * PI);
    const double val_Syy = val_rho * gyy - (lEy * lEy + lBy * lBy) / (FOUR * PI);
    const double val_Syz = val_rho * gyz - (lEy * lEz + lBy * lBz) / (FOUR * PI);
    const double val_Szz = val_rho * gzz - (lEz * lEz + lBz * lBz) / (FOUR * PI);

    // ======================================================================
    // Write outputs
    // ======================================================================
    Ex_rhs[idx]   = val_Ex_rhs;
    Ey_rhs[idx]   = val_Ey_rhs;
    Ez_rhs[idx]   = val_Ez_rhs;
    Bx_rhs[idx]   = val_Bx_rhs;
    By_rhs[idx]   = val_By_rhs;
    Bz_rhs[idx]   = val_Bz_rhs;
    Kpsi_rhs[idx] = val_Kpsi_rhs;
    Kphi_rhs[idx] = val_Kphi_rhs;

    rho[idx] = val_rho;
    Sx[idx]  = val_Sx;
    Sy[idx]  = val_Sy;
    Sz[idx]  = val_Sz;
    Sxx[idx] = val_Sxx;
    Sxy[idx] = val_Sxy;
    Sxz[idx] = val_Sxz;
    Syy[idx] = val_Syy;
    Syz[idx] = val_Syz;
    Szz[idx] = val_Szz;

    (void)lev;
}

// ============================================================================
// Cartesian advection + Kreiss-Oliger dissipation kernel for EM fields
// 8 fields (Kpsi, Kphi, Ex, Ey, Ez, Bx, By, Bz) over a fused 1D weight
// stencil (lopsided upwind + KO). Mirrors bssn_advection_dissipation_kernel
// in spirit but for the EM-only field set and without GAUGE skips.
// Ghost zones for Kpsi..Bz must have been pre-mirrored (see Step_GPU).
// ============================================================================
#include "advection_dissipation_smem_gpu.cuh"

template<int GW>
__launch_bounds__(256, 2)
__global__ void empart_advection_dissipation_kernel(
    const int ex0, const int ex1, const int ex2,
    const int pex0, const int pex1, const int pex2,
    const int symmetry, const int lev, const double eps,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ betax, const double* __restrict__ betay, const double* __restrict__ betaz,
    const double* __restrict__ Kpsi, const double* __restrict__ Kphi,
    const double* __restrict__ Ex,   const double* __restrict__ Ey,   const double* __restrict__ Ez,
    const double* __restrict__ Bx,   const double* __restrict__ By,   const double* __restrict__ Bz,
    double* __restrict__ Kpsi_rhs, double* __restrict__ Kphi_rhs,
    double* __restrict__ Ex_rhs,   double* __restrict__ Ey_rhs,   double* __restrict__ Ez_rhs,
    double* __restrict__ Bx_rhs,   double* __restrict__ By_rhs,   double* __restrict__ Bz_rhs
) {
    constexpr int MAX_DIST = adv_diss::max_dist<GW>();
    constexpr int WSIZE = 2 * MAX_DIST + 1;
    constexpr int WCTR = MAX_DIST;

    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;

    const int imax = ex0 - 1;
    const int jmax = ex1 - 1;
    const int kmax = ex2 - 1;
    if (i > imax || j > jmax || k > kmax) return;

    const int idx = i + pex0 * (j + pex1 * k);

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    // Note: imin = -GW (not -ORD = -(GW-1)) — the upwind+KO stencil reaches
    // one cell further than the centered stencil, matching the BSSN base.
    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > 0 && fabs(Z[0]) < dZ) kmin = -GW;
    if (symmetry > 1 && fabs(X[0]) < dX) imin = -GW;
    if (symmetry > 1 && fabs(Y[0]) < dY) jmin = -GW;

    const double vx = __ldg(&betax[idx]);
    const double vy = __ldg(&betay[idx]);
    const double vz = __ldg(&betaz[idx]);

    double Wx[WSIZE] = {0};
    double Wy[WSIZE] = {0};
    double Wz[WSIZE] = {0};

    // Per-axis lopsided / centered upwind stencil weights (from helper).
    adv_diss::assemble_advect_weights_1d<GW>(Wx, vx, dX, i, imin, imax);
    adv_diss::assemble_advect_weights_1d<GW>(Wy, vy, dY, j, jmin, jmax);
    adv_diss::assemble_advect_weights_1d<GW>(Wz, vz, dZ, k, kmin, kmax);

    // Fuse Kreiss-Oliger dissipation into the same weight array — same gating
    // as d_kodis_point (stencil must fit entirely inside [imin, imax]).
    const bool do_dissipation = (eps > 0.0
        && adv_diss::ko_fits<GW>(i, imin, imax)
        && adv_diss::ko_fits<GW>(j, jmin, jmax)
        && adv_diss::ko_fits<GW>(k, kmin, kmax));
    if (do_dissipation) {
        constexpr double cof  = adv_diss::ko_cof<GW>();
        constexpr double sgn  = adv_diss::ko_sign<GW>();
        const double edx = sgn * eps / (cof * dX);
        const double edy = sgn * eps / (cof * dY);
        const double edz = sgn * eps / (cof * dZ);
        adv_diss::assemble_ko_weights_1d<GW>(Wx, edx);
        adv_diss::assemble_ko_weights_1d<GW>(Wy, edy);
        adv_diss::assemble_ko_weights_1d<GW>(Wz, edz);
    }

    // 8-field batched __ldg sweep (Kpsi, Kphi, Ex, Ey, Ez, Bx, By, Bz).
    const int dx_s = 1;
    const int dy_s = pex0;
    const int dz_s = pex0 * pex1;
    const double center_weight = Wx[WCTR] + Wy[WCTR] + Wz[WCTR];

    const double *f;
    double *f_rhs;

    #define GET_EMPART_STEP_POINTERS(step_idx) \
        f = nullptr; f_rhs = nullptr; \
        switch (step_idx) { \
            case 0: f = Kpsi; f_rhs = Kpsi_rhs; break; \
            case 1: f = Kphi; f_rhs = Kphi_rhs; break; \
            case 2: f = Ex;   f_rhs = Ex_rhs;   break; \
            case 3: f = Ey;   f_rhs = Ey_rhs;   break; \
            case 4: f = Ez;   f_rhs = Ez_rhs;   break; \
            case 5: f = Bx;   f_rhs = Bx_rhs;   break; \
            case 6: f = By;   f_rhs = By_rhs;   break; \
            case 7: f = Bz;   f_rhs = Bz_rhs;   break; \
        }

    #pragma unroll 1
    for (int step = 0; step < 8; ++step) {
        GET_EMPART_STEP_POINTERS(step);
        double center_val = __ldg(&f[idx]);
        double res = center_weight * center_val;

        #pragma unroll
        for (int d = 1; d <= MAX_DIST; ++d) {
            res += Wx[WCTR-d] * __ldg(&f[idx - d * dx_s]) + Wx[WCTR+d] * __ldg(&f[idx + d * dx_s]);
            res += Wy[WCTR-d] * __ldg(&f[idx - d * dy_s]) + Wy[WCTR+d] * __ldg(&f[idx + d * dy_s]);
            res += Wz[WCTR-d] * __ldg(&f[idx - d * dz_s]) + Wz[WCTR+d] * __ldg(&f[idx + d * dz_s]);
        }

        f_rhs[idx] += res;
    }

    #undef GET_EMPART_STEP_POINTERS
    (void)lev;
}

// ============================================================================
// Shell kernel: compute_rhs_empart_ss per-point
// ============================================================================
template<int GW>
__launch_bounds__(256, 2)
__global__ void empart_rhs_ss_kernel(
    const int ex0, const int ex1, const int ex2,
    const int pex0, const int pex1, const int pex2,
    const int symmetry, const int lev, const int sst,
    const double eps,
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
    const double* __restrict__ Jx,
    const double* __restrict__ Jy,
    const double* __restrict__ Jz,
    const double* __restrict__ qchar,
    const double* __restrict__ dxx,
    const double* __restrict__ dxy,
    const double* __restrict__ dxz,
    const double* __restrict__ dyy,
    const double* __restrict__ dyz,
    const double* __restrict__ dzz,
    const double* __restrict__ Lap,
    const double* __restrict__ betax,
    const double* __restrict__ betay,
    const double* __restrict__ betaz,
    const double* __restrict__ trK,
    const double* __restrict__ Ex,
    const double* __restrict__ Ey,
    const double* __restrict__ Ez,
    const double* __restrict__ Bx,
    const double* __restrict__ By,
    const double* __restrict__ Bz,
    const double* __restrict__ Kpsi,
    const double* __restrict__ Kphi,
    double* __restrict__ Ex_rhs,
    double* __restrict__ Ey_rhs,
    double* __restrict__ Ez_rhs,
    double* __restrict__ Bx_rhs,
    double* __restrict__ By_rhs,
    double* __restrict__ Bz_rhs,
    double* __restrict__ Kpsi_rhs,
    double* __restrict__ Kphi_rhs,
    double* __restrict__ rho,
    double* __restrict__ Sx,
    double* __restrict__ Sy,
    double* __restrict__ Sz,
    double* __restrict__ Sxx,
    double* __restrict__ Sxy,
    double* __restrict__ Sxz,
    double* __restrict__ Syy,
    double* __restrict__ Syz,
    double* __restrict__ Szz
) {
    // ---- Thread index ----
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;

    // ---- ex/pex as arrays for derivative function calls ----
    const int ex[3]  = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    // Match BSSN kernel bounds: compute on all domain points [0, ex0-1]
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int idx = IDX3D_PAD(i, j, k, pex0, pex1);

    // ---- Constants ----
    constexpr double PI = CUDART_PI;
    constexpr double FOUR_PI = 4.0 * CUDART_PI;
    constexpr double EIT = 8.0;

    // ---- Load local field values ----
    const double val_chi   = chi[idx];
    const double val_Jx    = Jx[idx];
    const double val_Jy    = Jy[idx];
    const double val_Jz    = Jz[idx];
    const double val_qchar = qchar[idx];

    const double val_dxx = dxx[idx];
    const double val_dxy = dxy[idx];
    const double val_dxz = dxz[idx];
    const double val_dyy = dyy[idx];
    const double val_dyz = dyz[idx];
    const double val_dzz = dzz[idx];

    const double val_Lap   = Lap[idx];
    const double val_betax = betax[idx];
    const double val_betay = betay[idx];
    const double val_betaz = betaz[idx];
    const double val_trK   = trK[idx];

    const double val_Ex = Ex[idx];
    const double val_Ey = Ey[idx];
    const double val_Ez = Ez[idx];
    const double val_Bx = Bx[idx];
    const double val_By = By[idx];
    const double val_Bz = Bz[idx];
    const double val_Kpsi = Kpsi[idx];
    const double val_Kphi = Kphi[idx];

    // ---- alpn1 and chin1 (point-wise) ----
    const double alpn1 = val_Lap + ONE;
    // Clamp chin1 to prevent sqrt(negative) -> NaN and 1/0 -> inf when
    // chi is near the puncture (chitiny safeguard matching CPU lowerboundset).
    // Match Fortran: chin1 = chi + ONE, no clamping (CPU applies lowerboundset before RHS)
    const double chin1 = val_chi + ONE;
    const double chi3o2 = chin1 * sqrt(chin1);
    const double inv_chin1 = 1.0 / chin1;

    // ---- Conformal metric components ----
    double gxx = val_dxx + ONE;
    double gyy = val_dyy + ONE;
    double gzz = val_dzz + ONE;
    double gxy = val_dxy;
    double gxz = val_dxz;
    double gyz = val_dyz;

    // ======================================================================
    // Stage 3: Shell coordinate-transformed first derivatives of 19 fields
    // ======================================================================

    // --- Lapse ---
    double Lapx, Lapy, Lapz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Lap, &Lapx, &Lapy, &Lapz,
        crho, sigma, R, SYM, SYM, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    // --- betax ---
    double betaxx, betaxy, betaxz;
    d_point_fderivs_shc_gw<GW>(ex, pex, betax, &betaxx, &betaxy, &betaxz,
        crho, sigma, R, ANTI, SYM, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    // --- betay ---
    double betayx, betayy, betayz;
    d_point_fderivs_shc_gw<GW>(ex, pex, betay, &betayx, &betayy, &betayz,
        crho, sigma, R, SYM, ANTI, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    // --- betaz ---
    double betazx, betazy, betazz;
    d_point_fderivs_shc_gw<GW>(ex, pex, betaz, &betazx, &betazy, &betazz,
        crho, sigma, R, SYM, SYM, ANTI, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    // --- Conformal factor ---
    double chix, chiy, chiz;
    d_point_fderivs_shc_gw<GW>(ex, pex, chi, &chix, &chiy, &chiz,
        crho, sigma, R, SYM, SYM, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    // --- Conformal metric perturbation derivatives ---
    double gxxx, gxxy, gxxz;
    d_point_fderivs_shc_gw<GW>(ex, pex, dxx, &gxxx, &gxxy, &gxxz,
        crho, sigma, R, SYM, SYM, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double gxyx, gxyy, gxyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, dxy, &gxyx, &gxyy, &gxyz,
        crho, sigma, R, ANTI, ANTI, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double gxzx, gxzy, gxzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, dxz, &gxzx, &gxzy, &gxzz,
        crho, sigma, R, ANTI, SYM, ANTI, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double gyyx, gyyy, gyyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, dyy, &gyyx, &gyyy, &gyyz,
        crho, sigma, R, SYM, SYM, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double gyzx, gyzy, gyzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, dyz, &gyzx, &gyzy, &gyzz,
        crho, sigma, R, SYM, ANTI, ANTI, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double gzzx, gzzy, gzzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, dzz, &gzzx, &gzzy, &gzzz,
        crho, sigma, R, SYM, SYM, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    // --- Constraint damping scalars ---
    double Kpsix, Kpsiy, Kpsiz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Kpsi, &Kpsix, &Kpsiy, &Kpsiz,
        crho, sigma, R, SYM, SYM, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double Kphix, Kphiy, Kphiz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Kphi, &Kphix, &Kphiy, &Kphiz,
        crho, sigma, R, SYM, SYM, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    // --- Electric field ---
    double Exx, Exy, Exz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Ex, &Exx, &Exy, &Exz,
        crho, sigma, R, ANTI, SYM, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double Eyx, Eyy, Eyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Ey, &Eyx, &Eyy, &Eyz,
        crho, sigma, R, SYM, ANTI, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double Ezx, Ezy, Ezz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Ez, &Ezx, &Ezy, &Ezz,
        crho, sigma, R, SYM, SYM, ANTI, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    // --- Magnetic field (axial vector symmetry: SYM,ANTI,ANTI etc.) ---
    double Bxx, Bxy, Bxz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Bx, &Bxx, &Bxy, &Bxz,
        crho, sigma, R, SYM, ANTI, ANTI, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double Byx, Byy, Byz;
    d_point_fderivs_shc_gw<GW>(ex, pex, By, &Byx, &Byy, &Byz,
        crho, sigma, R, ANTI, SYM, ANTI, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    double Bzx, Bzy, Bzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Bz, &Bzx, &Bzy, &Bzz,
        crho, sigma, R, ANTI, ANTI, SYM, symmetry, ONOFF, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);

    // ======================================================================
    // Stage 4: Physical metric reconstruction (same as Cartesian)
    // ======================================================================

    gxx /= chin1;
    gxy /= chin1;
    gxz /= chin1;
    gyy /= chin1;
    gyz /= chin1;
    gzz /= chin1;

    gxxx = (gxxx - gxx * chix) * inv_chin1;
    gxxy = (gxxy - gxx * chiy) * inv_chin1;
    gxxz = (gxxz - gxx * chiz) * inv_chin1;

    gxyx = (gxyx - gxy * chix) * inv_chin1;
    gxyy = (gxyy - gxy * chiy) * inv_chin1;
    gxyz = (gxyz - gxy * chiz) * inv_chin1;

    gxzx = (gxzx - gxz * chix) * inv_chin1;
    gxzy = (gxzy - gxz * chiy) * inv_chin1;
    gxzz = (gxzz - gxz * chiz) * inv_chin1;

    gyyx = (gyyx - gyy * chix) * inv_chin1;
    gyyy = (gyyy - gyy * chiy) * inv_chin1;
    gyyz = (gyyz - gyy * chiz) * inv_chin1;

    gyzx = (gyzx - gyz * chix) * inv_chin1;
    gyzy = (gyzy - gyz * chiy) * inv_chin1;
    gyzz = (gyzz - gyz * chiz) * inv_chin1;

    gzzx = (gzzx - gzz * chix) * inv_chin1;
    gzzy = (gzzy - gzz * chiy) * inv_chin1;
    gzzz = (gzzz - gzz * chiz) * inv_chin1;

    // Inverse physical metric
    const double det = gxx * gyy * gzz
        + gxy * gyz * gxz
        + gxz * gxy * gyz
        - gxz * gyy * gxz
        - gxy * gxy * gzz
        - gxx * gyz * gyz;
    // Match Fortran: direct inverse, no fabs/fmax clamping
    const double inv_det = 1.0 / det;

    const double gupxx = (gyy * gzz - gyz * gyz) * inv_det;
    const double gupxy = -(gxy * gzz - gyz * gxz) * inv_det;
    const double gupxz = (gxy * gyz - gyy * gxz) * inv_det;
    const double gupyy = (gxx * gzz - gxz * gxz) * inv_det;
    const double gupyz = -(gxx * gyz - gxy * gxz) * inv_det;
    const double gupzz = (gxx * gyy - gxy * gxy) * inv_det;

    // ======================================================================
    // Stage 5: RHS source computation
    // ======================================================================

    double val_Ex_rhs, val_Ey_rhs, val_Ez_rhs;
    double val_Bx_rhs, val_By_rhs, val_Bz_rhs;
    double val_Kpsi_rhs, val_Kphi_rhs;

    // --- Ex_rhs ---
    val_Ex_rhs = alpn1 * val_trK * val_Ex
        - (val_Ex * betaxx + val_Ey * betaxy + val_Ez * betaxz)
        - FOUR_PI * alpn1 * val_Jx
        - alpn1 * (gupxx * Kpsix + gupxy * Kpsiy + gupxz * Kpsiz)
        + chi3o2 * (
            ((gxz * val_Bx + gyz * val_By + gzz * val_Bz) * Lapy
                + alpn1 * (gxz * Bxy + gyz * Byy + gzz * Bzy)
                + alpn1 * (val_Bx * gxzy + val_By * gyzy + val_Bz * gzzy))
          - ((gxy * val_Bx + gyy * val_By + gyz * val_Bz) * Lapz
                + alpn1 * (gxy * Bxz + gyy * Byz + gyz * Bzz)
                + alpn1 * (val_Bx * gxyz + val_By * gyyz + val_Bz * gyzz))
        );

    // --- Ey_rhs ---
    val_Ey_rhs = alpn1 * val_trK * val_Ey
        - (val_Ex * betayx + val_Ey * betayy + val_Ez * betayz)
        - FOUR_PI * alpn1 * val_Jy
        - alpn1 * (gupxy * Kpsix + gupyy * Kpsiy + gupyz * Kpsiz)
        + chi3o2 * (
            ((gxx * val_Bx + gxy * val_By + gxz * val_Bz) * Lapz
                + alpn1 * (gxx * Bxz + gxy * Byz + gxz * Bzz)
                + alpn1 * (val_Bx * gxxz + val_By * gxyz + val_Bz * gxzz))
          - ((gxz * val_Bx + gyz * val_By + gzz * val_Bz) * Lapx
                + alpn1 * (gxz * Bxx + gyz * Byx + gzz * Bzx)
                + alpn1 * (val_Bx * gxzx + val_By * gyzx + val_Bz * gzzx))
        );

    // --- Ez_rhs ---
    val_Ez_rhs = alpn1 * val_trK * val_Ez
        - (val_Ex * betazx + val_Ey * betazy + val_Ez * betazz)
        - FOUR_PI * alpn1 * val_Jz
        - alpn1 * (gupxz * Kpsix + gupyz * Kpsiy + gupzz * Kpsiz)
        + chi3o2 * (
            ((gxy * val_Bx + gyy * val_By + gyz * val_Bz) * Lapx
                + alpn1 * (gxy * Bxx + gyy * Byx + gyz * Bzx)
                + alpn1 * (val_Bx * gxyx + val_By * gyyx + val_Bz * gyzx))
          - ((gxx * val_Bx + gxy * val_By + gxz * val_Bz) * Lapy
                + alpn1 * (gxx * Bxy + gxy * Byy + gxz * Bzy)
                + alpn1 * (val_Bx * gxxy + val_By * gxyy + val_Bz * gxzy))
        );

    // --- Bx_rhs ---
    val_Bx_rhs = alpn1 * val_trK * val_Bx
        - (val_Bx * betaxx + val_By * betaxy + val_Bz * betaxz)
        - alpn1 * (gupxx * Kphix + gupxy * Kphiy + gupxz * Kphiz)
        - chi3o2 * (
            ((gxz * val_Ex + gyz * val_Ey + gzz * val_Ez) * Lapy
                + alpn1 * (gxz * Exy + gyz * Eyy + gzz * Ezy)
                + alpn1 * (val_Ex * gxzy + val_Ey * gyzy + val_Ez * gzzy))
          - ((gxy * val_Ex + gyy * val_Ey + gyz * val_Ez) * Lapz
                + alpn1 * (gxy * Exz + gyy * Eyz + gyz * Ezz)
                + alpn1 * (val_Ex * gxyz + val_Ey * gyyz + val_Ez * gyzz))
        );

    // --- By_rhs ---
    val_By_rhs = alpn1 * val_trK * val_By
        - (val_Bx * betayx + val_By * betayy + val_Bz * betayz)
        - alpn1 * (gupxy * Kphix + gupyy * Kphiy + gupyz * Kphiz)
        - chi3o2 * (
            ((gxx * val_Ex + gxy * val_Ey + gxz * val_Ez) * Lapz
                + alpn1 * (gxx * Exz + gxy * Eyz + gxz * Ezz)
                + alpn1 * (val_Ex * gxxz + val_Ey * gxyz + val_Ez * gxzz))
          - ((gxz * val_Ex + gyz * val_Ey + gzz * val_Ez) * Lapx
                + alpn1 * (gxz * Exx + gyz * Eyx + gzz * Ezx)
                + alpn1 * (val_Ex * gxzx + val_Ey * gyzx + val_Ez * gzzx))
        );

    // --- Bz_rhs ---
    val_Bz_rhs = alpn1 * val_trK * val_Bz
        - (val_Bx * betazx + val_By * betazy + val_Bz * betazz)
        - alpn1 * (gupxz * Kphix + gupyz * Kphiy + gupzz * Kphiz)
        - chi3o2 * (
            ((gxy * val_Ex + gyy * val_Ey + gyz * val_Ez) * Lapx
                + alpn1 * (gxy * Exx + gyy * Eyx + gyz * Ezx)
                + alpn1 * (val_Ex * gxyx + val_Ey * gyyx + val_Ez * gyzx))
          - ((gxx * val_Ex + gxy * val_Ey + gxz * val_Ez) * Lapy
                + alpn1 * (gxx * Exy + gxy * Eyy + gxz * Ezy)
                + alpn1 * (val_Ex * gxxy + val_Ey * gxyy + val_Ez * gxzy))
        );

    // --- Kpsi_rhs ---
    val_Kpsi_rhs = FOUR * PI * alpn1 * val_qchar
        - alpn1 * KAPPA * val_Kpsi
        - alpn1 * (Exx + Eyy + Ezz
            - F3o2 * inv_chin1 * (chix * val_Ex + chiy * val_Ey + chiz * val_Ez));

    // --- Kphi_rhs ---
    val_Kphi_rhs = -alpn1 * KAPPA * val_Kphi
        - alpn1 * (Bxx + Byy + Bzz
            - F3o2 * inv_chin1 * (chix * val_Bx + chiy * val_By + chiz * val_Bz));

    // ======================================================================
    // Stage 6: Advection (direct) + Kreiss-Oliger dissipation (shell)
    // ======================================================================

    // Advection: f_rhs += betax*fx + betay*fy + betaz*fz
    val_Kpsi_rhs += val_betax * Kpsix + val_betay * Kpsiy + val_betaz * Kpsiz;
    val_Kphi_rhs += val_betax * Kphix + val_betay * Kphiy + val_betaz * Kphiz;

    val_Ex_rhs += val_betax * Exx + val_betay * Exy + val_betaz * Exz;
    val_Ey_rhs += val_betax * Eyx + val_betay * Eyy + val_betaz * Eyz;
    val_Ez_rhs += val_betax * Ezx + val_betay * Ezy + val_betaz * Ezz;

    val_Bx_rhs += val_betax * Bxx + val_betay * Bxy + val_betaz * Bxz;
    val_By_rhs += val_betax * Byx + val_betay * Byy + val_betaz * Byz;
    val_Bz_rhs += val_betax * Bzx + val_betay * Bzy + val_betaz * Bzz;

    // Kreiss-Oliger dissipation
    if (eps > 0.0) {
        val_Kpsi_rhs += d_kodis_sh_point<GW>(ex, pex, Kpsi,
            crho, sigma, R, SYM, SYM, SYM, symmetry, eps, sst, i, j, k);

        val_Kphi_rhs += d_kodis_sh_point<GW>(ex, pex, Kphi,
            crho, sigma, R, SYM, SYM, SYM, symmetry, eps, sst, i, j, k);

        val_Ex_rhs += d_kodis_sh_point<GW>(ex, pex, Ex,
            crho, sigma, R, ANTI, SYM, SYM, symmetry, eps, sst, i, j, k);

        val_Ey_rhs += d_kodis_sh_point<GW>(ex, pex, Ey,
            crho, sigma, R, SYM, ANTI, SYM, symmetry, eps, sst, i, j, k);

        val_Ez_rhs += d_kodis_sh_point<GW>(ex, pex, Ez,
            crho, sigma, R, SYM, SYM, ANTI, symmetry, eps, sst, i, j, k);

        val_Bx_rhs += d_kodis_sh_point<GW>(ex, pex, Bx,
            crho, sigma, R, SYM, ANTI, ANTI, symmetry, eps, sst, i, j, k);

        val_By_rhs += d_kodis_sh_point<GW>(ex, pex, By,
            crho, sigma, R, ANTI, SYM, ANTI, symmetry, eps, sst, i, j, k);

        val_Bz_rhs += d_kodis_sh_point<GW>(ex, pex, Bz,
            crho, sigma, R, ANTI, ANTI, SYM, symmetry, eps, sst, i, j, k);
    }

    // ======================================================================
    // Stage 7: Stress-energy tensor
    // ======================================================================

    // --- rho ---
    const double val_rho = (
        gxx * (val_Ex * val_Ex + val_Bx * val_Bx)
        + gyy * (val_Ey * val_Ey + val_By * val_By)
        + gzz * (val_Ez * val_Ez + val_Bz * val_Bz)
        + TWO * (
            gxy * (val_Ex * val_Ey + val_Bx * val_By)
            + gxz * (val_Ex * val_Ez + val_Bx * val_Bz)
            + gyz * (val_Ey * val_Ez + val_By * val_Bz)
        )
    ) / (EIT * PI);

    // --- Sx, Sy, Sz (Poynting vector) ---
    const double val_Sx = (val_Ey * val_Bz - val_Ez * val_By) / (FOUR * PI * chi3o2);
    const double val_Sy = (val_Ez * val_Bx - val_Ex * val_Bz) / (FOUR * PI * chi3o2);
    const double val_Sz = (val_Ex * val_By - val_Ey * val_Bx) / (FOUR * PI * chi3o2);

    // --- Lowered-index electric field (physical metric) ---
    const double lEx = gxx * val_Ex + gxy * val_Ey + gxz * val_Ez;
    const double lEy = gxy * val_Ex + gyy * val_Ey + gyz * val_Ez;
    const double lEz = gxz * val_Ex + gyz * val_Ey + gzz * val_Ez;

    // --- Lowered-index magnetic field (physical metric) ---
    const double lBx = gxx * val_Bx + gxy * val_By + gxz * val_Bz;
    const double lBy = gxy * val_Bx + gyy * val_By + gyz * val_Bz;
    const double lBz = gxz * val_Bx + gyz * val_By + gzz * val_Bz;

    // --- Stress tensor S_ij ---
    const double val_Sxx = val_rho * gxx - (lEx * lEx + lBx * lBx) / (FOUR * PI);
    const double val_Sxy = val_rho * gxy - (lEx * lEy + lBx * lBy) / (FOUR * PI);
    const double val_Sxz = val_rho * gxz - (lEx * lEz + lBx * lBz) / (FOUR * PI);
    const double val_Syy = val_rho * gyy - (lEy * lEy + lBy * lBy) / (FOUR * PI);
    const double val_Syz = val_rho * gyz - (lEy * lEz + lBy * lBz) / (FOUR * PI);
    const double val_Szz = val_rho * gzz - (lEz * lEz + lBz * lBz) / (FOUR * PI);

    // ======================================================================
    // Write outputs
    // ======================================================================
    Ex_rhs[idx]   = val_Ex_rhs;
    Ey_rhs[idx]   = val_Ey_rhs;
    Ez_rhs[idx]   = val_Ez_rhs;
    Bx_rhs[idx]   = val_Bx_rhs;
    By_rhs[idx]   = val_By_rhs;
    Bz_rhs[idx]   = val_Bz_rhs;
    Kpsi_rhs[idx] = val_Kpsi_rhs;
    Kphi_rhs[idx] = val_Kphi_rhs;

    rho[idx] = val_rho;
    Sx[idx]  = val_Sx;
    Sy[idx]  = val_Sy;
    Sz[idx]  = val_Sz;
    Sxx[idx] = val_Sxx;
    Sxy[idx] = val_Sxy;
    Sxz[idx] = val_Sxz;
    Syy[idx] = val_Syy;
    Syz[idx] = val_Syz;
    Szz[idx] = val_Szz;

    (void)lev;
}

// ============================================================================
// Host launch: Cartesian variant
// ============================================================================
void gpu_compute_rhs_empart_launch(
    cudaStream_t &stream,
    int* ex, int* pex,
    double* X, double* Y, double* Z,
    double* chi, double* Jx, double* Jy, double* Jz, double* qchar,
    double* dxx, double* dxy, double* dxz, double* dyy, double* dyz, double* dzz,
    double* Lap, double* betax, double* betay, double* betaz, double* trK,
    double* Ex, double* Ey, double* Ez, double* Bx, double* By, double* Bz,
    double* Kpsi, double* Kphi,
    double* Ex_rhs, double* Ey_rhs, double* Ez_rhs,
    double* Bx_rhs, double* By_rhs, double* Bz_rhs,
    double* Kpsi_rhs, double* Kphi_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz,
    double* Syy, double* Syz, double* Szz,
    int symmetry, int lev, double eps
) {
    const dim3 block(DIMX, DIMY, DIMZ);
    const dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );

    // 3 SMEM tiles of TileParams<GW>::SIZE doubles each — exceeds the 48KB
    // static smem limit for GW>=4, so we must opt in via cudaFuncSetAttribute.
    constexpr int smem_size = 3 * TileParams<ghost_width>::SIZE * sizeof(double);
    static bool smem_attr_set = false;
    if (!smem_attr_set) {
        cudaFuncSetAttribute(empart_rhs_kernel<ghost_width>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        smem_attr_set = true;
    }

    empart_rhs_kernel<ghost_width><<<grid, block, smem_size, stream>>>(
        ex[0], ex[1], ex[2],
        pex[0], pex[1], pex[2],
        symmetry, lev, eps,
        X, Y, Z,
        chi, Jx, Jy, Jz, qchar,
        dxx, dxy, dxz, dyy, dyz, dzz,
        Lap, betax, betay, betaz, trK,
        Ex, Ey, Ez, Bx, By, Bz, Kpsi, Kphi,
        Ex_rhs, Ey_rhs, Ez_rhs,
        Bx_rhs, By_rhs, Bz_rhs,
        Kpsi_rhs, Kphi_rhs,
        rho, Sx, Sy, Sz,
        Sxx, Sxy, Sxz,
        Syy, Syz, Szz
    );

    // Second stage: fused lopsided advection + KO dissipation on the same
    // stream. Block shape (32, 8, 1) keeps the x-warp coalesced for the
    // ldg-heavy stencil sweep — different from the (8, 8, 4) SMEM-source
    // RHS kernel above. Same stream → serial after RHS.
    const dim3 block_adv(32, 8, 1);
    const dim3 grid_adv(
        (ex[0] + block_adv.x - 1) / block_adv.x,
        (ex[1] + block_adv.y - 1) / block_adv.y,
        (ex[2] + block_adv.z - 1) / block_adv.z
    );
    empart_advection_dissipation_kernel<ghost_width><<<grid_adv, block_adv, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        symmetry, lev, eps,
        X, Y, Z,
        betax, betay, betaz,
        Kpsi, Kphi,
        Ex, Ey, Ez, Bx, By, Bz,
        Kpsi_rhs, Kphi_rhs,
        Ex_rhs, Ey_rhs, Ez_rhs,
        Bx_rhs, By_rhs, Bz_rhs
    );

    CUDA_CHECK(cudaGetLastError());
}

// ============================================================================
// Host launch: Shell variant
// ============================================================================
void gpu_compute_rhs_empart_ss_launch(
    cudaStream_t &stream,
    int* ex, int* pex,
    double* crho, double* sigma, double* R,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* chi, double* Jx, double* Jy, double* Jz, double* qchar,
    double* dxx, double* dxy, double* dxz, double* dyy, double* dyz, double* dzz,
    double* Lap, double* betax, double* betay, double* betaz, double* trK,
    double* Ex, double* Ey, double* Ez, double* Bx, double* By, double* Bz,
    double* Kpsi, double* Kphi,
    double* Ex_rhs, double* Ey_rhs, double* Ez_rhs,
    double* Bx_rhs, double* By_rhs, double* Bz_rhs,
    double* Kpsi_rhs, double* Kphi_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz,
    double* Syy, double* Syz, double* Szz,
    int symmetry, int lev, double eps, int sst
) {
    const dim3 block(8, 8, 4);
    const dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );

    empart_rhs_ss_kernel<ghost_width><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2],
        pex[0], pex[1], pex[2],
        symmetry, lev, sst, eps,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        chi, Jx, Jy, Jz, qchar,
        dxx, dxy, dxz, dyy, dyz, dzz,
        Lap, betax, betay, betaz, trK,
        Ex, Ey, Ez, Bx, By, Bz, Kpsi, Kphi,
        Ex_rhs, Ey_rhs, Ez_rhs,
        Bx_rhs, By_rhs, Bz_rhs,
        Kpsi_rhs, Kphi_rhs,
        rho, Sx, Sy, Sz,
        Sxx, Sxy, Sxz,
        Syy, Syz, Szz
    );

    CUDA_CHECK(cudaGetLastError());
}

// ============================================================================
// Fused symmetry boundary condition kernel for EM fields
// ============================================================================
__global__ void empart_apply_symmetry_bd_kernel(
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    double* __restrict__ Ex,   double* __restrict__ Ey,   double* __restrict__ Ez,
    double* __restrict__ Bx,   double* __restrict__ By,   double* __restrict__ Bz,
    double* __restrict__ Kpsi, double* __restrict__ Kphi
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= pextc0 * pextc1 * pextc2) return;
    int i_pad, j_pad, k_pad;
    linear_to_ijk(tid, pextc0, pextc1, i_pad, j_pad, k_pad);

    // i_pad/j_pad/k_pad walk the full padded buffer; (i,j,k) are the logical
    // coordinates the d_fgfs[..] pointers are aligned to (Block.C arranges
    // d_fgfs[i] = d_fgfs_base + offset so the pointer sits at logical (0,0,0)
    // and negative-index reads land in the ghost zone).
    int i = i_pad - GPU_PAD_X_LEFT;
    int j = j_pad - ghost_width;
    int k = k_pad - ghost_width;

    if (i < -ghost_width) return; // skip the x-axis alignment padding
    if (i >= 0 && j >= 0 && k >= 0) return;
    if (i >= extc0 || j >= extc1 || k >= extc2) return;
    if (i + extc0 <= 0 || j + extc1 <= 0 || k + extc2 <= 0) return;

    // out_idx / src_idx are in LOGICAL coordinates (matching the
    // apply_symmetry_bd_fused_kernel in fmisc_gpu.cu, line 1511) so that
    // Ex[out_idx] addresses the correct cell relative to d_fgfs[i]'s
    // logical-(0,0,0) origin.
    int out_idx = (k * pextc1 + j) * pextc0 + i;

    bool i_neg = (i < 0);
    bool j_neg = (j < 0);
    bool k_neg = (k < 0);
#ifdef Vertex
    int ri = i_neg ? (-i) : i;
    int rj = j_neg ? (-j) : j;
    int rk = k_neg ? (-k) : k;
#else
    int ri = i_neg ? (-i - 1) : i;
    int rj = j_neg ? (-j - 1) : j;
    int rk = k_neg ? (-k - 1) : k;
#endif
    int src_idx = (rk * pextc1 + rj) * pextc0 + ri;

    double fi = i_neg ? -1.0 : 1.0;
    double fj = j_neg ? -1.0 : 1.0;
    double fk = k_neg ? -1.0 : 1.0;

    // Ex: ANTI,SYM,SYM; Ey: SYM,ANTI,SYM; Ez: SYM,SYM,ANTI
    // Bx: SYM,ANTI,ANTI; By: ANTI,SYM,ANTI; Bz: ANTI,ANTI,SYM
    // Kpsi: SYM,SYM,SYM; Kphi: SYM,SYM,SYM
    Ex[out_idx]   = Ex[src_idx]   * fi;
    Ey[out_idx]   = Ey[src_idx]   * fj;
    Ez[out_idx]   = Ez[src_idx]   * fk;
    Bx[out_idx]   = Bx[src_idx]   * fj * fk;
    By[out_idx]   = By[src_idx]   * fi * fk;
    Bz[out_idx]   = Bz[src_idx]   * fi * fj;
    Kpsi[out_idx] = Kpsi[src_idx];
    Kphi[out_idx] = Kphi[src_idx];
}

void gpu_apply_symmetry_bd_empart_launch(
    cudaStream_t stream,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    double* Ex,  double* Ey,  double* Ez,
    double* Bx,  double* By,  double* Bz,
    double* Kpsi, double* Kphi
) {
    int total = pextc0 * pextc1 * pextc2;
    dim3 block(256);
    dim3 grid((total + 255) / 256);

    empart_apply_symmetry_bd_kernel<<<grid, block, 0, stream>>>(
        extc0, extc1, extc2, pextc0, pextc1, pextc2,
        Ex, Ey, Ez, Bx, By, Bz, Kpsi, Kphi
    );

    CUDA_CHECK(cudaGetLastError());
}
