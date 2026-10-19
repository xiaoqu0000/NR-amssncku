#ifndef ADM_RICCI_GAMMA_GPU_CUH
#define ADM_RICCI_GAMMA_GPU_CUH

#include "diff_new_gpu.cuh"
#include "diff_new_sh_gpu.cuh"

// ============================================================
// Device helper functions ported from getnp4old.f90
// ============================================================

// Kind 1 (Christoffel symbol of the first kind, lowered)
// ass_Gam_ijk = 0.5 * (g_ij,k + g_ki,j - g_jk,i)
__device__ __forceinline__ void d_kind1_connection(
    double gxxx, double gxyx, double gxzx, double gyyx, double gyzx, double gzzx,
    double gxxy, double gxyy, double gxzy, double gyyy, double gyzy, double gzzy,
    double gxxz, double gxyz, double gxzz, double gyyz, double gyzz, double gzzz,
    double& ass_Gamxxx, double& ass_Gamxxy, double& ass_Gamxxz,
    double& ass_Gamxyy, double& ass_Gamxyz, double& ass_Gamxzz,
    double& ass_Gamyxx, double& ass_Gamyxy, double& ass_Gamyxz,
    double& ass_Gamyyy, double& ass_Gamyyz, double& ass_Gamyzz,
    double& ass_Gamzxx, double& ass_Gamzxy, double& ass_Gamzxz,
    double& ass_Gamzyy, double& ass_Gamzyz, double& ass_Gamzzz
) {
    const double HLF = 0.5;

    ass_Gamxxx = HLF * (gxxx);
    ass_Gamyxx = HLF * (gxyx + gxyx - gxxy);
    ass_Gamzxx = HLF * (gxzx + gxzx - gxxz);
    ass_Gamxyy = HLF * (gxyy + gxyy - gyyx);
    ass_Gamyyy = HLF * (gyyy);
    ass_Gamzyy = HLF * (gyzy + gyzy - gyyz);
    ass_Gamxzz = HLF * (gxzz + gxzz - gzzx);
    ass_Gamyzz = HLF * (gyzz + gyzz - gzzy);
    ass_Gamzzz = HLF * (gzzz);
    ass_Gamxxy = HLF * (gxxy + gxyx - gxyx);
    ass_Gamyxy = HLF * (gxyy + gyyx - gxyy);
    ass_Gamzxy = HLF * (gxzy + gyzx - gxyz);
    ass_Gamxxz = HLF * (gxxz + gxzx - gxzx);
    ass_Gamyxz = HLF * (gxyz + gyzx - gxzy);
    ass_Gamzxz = HLF * (gxzz + gzzx - gxzz);
    ass_Gamxyz = HLF * (gxyz + gxzy - gyzx);
    ass_Gamyyz = HLF * (gyyz + gyzy - gyzy);
    ass_Gamzyz = HLF * (gyzz + gzzy - gyzz);
}

// Kind 2 (Christoffel symbol of the second kind)
// Gam^i_jk = gup^il * ass_Gam_ljk
__device__ __forceinline__ void d_kind2_connection(
    double gupxx, double gupxy, double gupxz, double gupyy, double gupyz, double gupzz,
    double ass_Gamxxx, double ass_Gamxxy, double ass_Gamxxz,
    double ass_Gamxyy, double ass_Gamxyz, double ass_Gamxzz,
    double ass_Gamyxx, double ass_Gamyxy, double ass_Gamyxz,
    double ass_Gamyyy, double ass_Gamyyz, double ass_Gamyzz,
    double ass_Gamzxx, double ass_Gamzxy, double ass_Gamzxz,
    double ass_Gamzyy, double ass_Gamzyz, double ass_Gamzzz,
    double& Gamxxx, double& Gamxxy, double& Gamxxz,
    double& Gamxyy, double& Gamxyz, double& Gamxzz,
    double& Gamyxx, double& Gamyxy, double& Gamyxz,
    double& Gamyyy, double& Gamyyz, double& Gamyzz,
    double& Gamzxx, double& Gamzxy, double& Gamzxz,
    double& Gamzyy, double& Gamzyz, double& Gamzzz
) {
    Gamxxx = gupxx * ass_Gamxxx + gupxy * ass_Gamyxx + gupxz * ass_Gamzxx;
    Gamxxy = gupxx * ass_Gamxxy + gupxy * ass_Gamyxy + gupxz * ass_Gamzxy;
    Gamxxz = gupxx * ass_Gamxxz + gupxy * ass_Gamyxz + gupxz * ass_Gamzxz;
    Gamxyy = gupxx * ass_Gamxyy + gupxy * ass_Gamyyy + gupxz * ass_Gamzyy;
    Gamxyz = gupxx * ass_Gamxyz + gupxy * ass_Gamyyz + gupxz * ass_Gamzyz;
    Gamxzz = gupxx * ass_Gamxzz + gupxy * ass_Gamyzz + gupxz * ass_Gamzzz;

    Gamyxx = gupxy * ass_Gamxxx + gupyy * ass_Gamyxx + gupyz * ass_Gamzxx;
    Gamyxy = gupxy * ass_Gamxxy + gupyy * ass_Gamyxy + gupyz * ass_Gamzxy;
    Gamyxz = gupxy * ass_Gamxxz + gupyy * ass_Gamyxz + gupyz * ass_Gamzxz;
    Gamyyy = gupxy * ass_Gamxyy + gupyy * ass_Gamyyy + gupyz * ass_Gamzyy;
    Gamyyz = gupxy * ass_Gamxyz + gupyy * ass_Gamyyz + gupyz * ass_Gamzyz;
    Gamyzz = gupxy * ass_Gamxzz + gupyy * ass_Gamyzz + gupyz * ass_Gamzzz;

    Gamzxx = gupxz * ass_Gamxxx + gupyz * ass_Gamyxx + gupzz * ass_Gamzxx;
    Gamzxy = gupxz * ass_Gamxxy + gupyz * ass_Gamyxy + gupzz * ass_Gamzxy;
    Gamzxz = gupxz * ass_Gamxxz + gupyz * ass_Gamyxz + gupzz * ass_Gamzxz;
    Gamzyy = gupxz * ass_Gamxyy + gupyz * ass_Gamyyy + gupzz * ass_Gamzyy;
    Gamzyz = gupxz * ass_Gamxyz + gupyz * ass_Gamyyz + gupzz * ass_Gamzyz;
    Gamzzz = gupxz * ass_Gamxzz + gupyz * ass_Gamyzz + gupzz * ass_Gamzzz;
}

// ADM Riemann tensor (independent components for 3D)
__device__ __forceinline__ void d_adm_riemann(
    double gxxxx, double gxxxy, double gxxxz, double gxxyy, double gxxyz, double gxxzz,
    double gxyxx, double gxyxy, double gxyxz, double gxyyy, double gxyyz, double gxyzz,
    double gxzxx, double gxzxy, double gxzxz, double gxzyy, double gxzyz, double gxzzz,
    double gyyxx, double gyyxy, double gyyxz, double gyyyy, double gyyyz, double gyyzz,
    double gyzxx, double gyzxy, double gyzxz, double gyzyy, double gyzyz, double gyzzz,
    double gzzxx, double gzzxy, double gzzxz, double gzzyy, double gzzyz, double gzzzz,
    double Gamxxx, double Gamxxy, double Gamxxz,
    double Gamxyy, double Gamxyz, double Gamxzz,
    double Gamyxx, double Gamyxy, double Gamyxz,
    double Gamyyy, double Gamyyz, double Gamyzz,
    double Gamzxx, double Gamzxy, double Gamzxz,
    double Gamzyy, double Gamzyz, double Gamzzz,
    double ass_Gamxxx, double ass_Gamxxy, double ass_Gamxxz,
    double ass_Gamxyy, double ass_Gamxyz, double ass_Gamxzz,
    double ass_Gamyxx, double ass_Gamyxy, double ass_Gamyxz,
    double ass_Gamyyy, double ass_Gamyyz, double ass_Gamyzz,
    double ass_Gamzxx, double ass_Gamzxy, double ass_Gamzxz,
    double ass_Gamzyy, double ass_Gamzyz, double ass_Gamzzz,
    double& Rxyxy, double& Rxyxz, double& Rxyyz,
    double& Rxzxz, double& Rxzyz, double& Ryzyz
) {
    const double HLF = 0.5;

    Rxyxy = HLF * (gxyxy + gxyxy - gxxyy - gyyxx)
          + (ass_Gamxxy * Gamxxy + ass_Gamyxy * Gamyxy + ass_Gamzxy * Gamzxy)
          - (ass_Gamxyy * Gamxxx + ass_Gamyyy * Gamyxx + ass_Gamzyy * Gamzxx);

    Rxyxz = HLF * (gxzxy + gxyxz - gxxyz - gyzxx)
          + (ass_Gamxxy * Gamxxz + ass_Gamyxy * Gamyxz + ass_Gamzxy * Gamzxz)
          - (ass_Gamxyz * Gamxxx + ass_Gamyyz * Gamyxx + ass_Gamzyz * Gamzxx);

    Rxyyz = HLF * (gxzyy + gyyxz - gxyyz - gyzxy)
          + (ass_Gamxyy * Gamxxz + ass_Gamyyy * Gamyxz + ass_Gamzyy * Gamzxz)
          - (ass_Gamxyz * Gamxxy + ass_Gamyyz * Gamyxy + ass_Gamzyz * Gamzxy);

    Rxzxz = HLF * (gxzxz + gxzxz - gxxzz - gzzxx)
          + (ass_Gamxxz * Gamxxz + ass_Gamyxz * Gamyxz + ass_Gamzxz * Gamzxz)
          - (ass_Gamxzz * Gamxxx + ass_Gamyzz * Gamyxx + ass_Gamzzz * Gamzxx);

    Rxzyz = HLF * (gxzyz + gyzxz - gxyzz - gzzxy)
          + (ass_Gamxyz * Gamxxz + ass_Gamyyz * Gamyxz + ass_Gamzyz * Gamzxz)
          - (ass_Gamxzz * Gamxxy + ass_Gamyzz * Gamyxy + ass_Gamzzz * Gamzxy);

    Ryzyz = HLF * (gyzyz + gyzyz - gyyzz - gzzyy)
          + (ass_Gamxyz * Gamxyz + ass_Gamyyz * Gamyyz + ass_Gamzyz * Gamzyz)
          - (ass_Gamxzz * Gamxyy + ass_Gamyzz * Gamyyy + ass_Gamzzz * Gamzyy);
}

// ADM Ricci tensor: R_ij = gup^kl * R_ikjl
__device__ __forceinline__ void d_adm_ricci(
    double gupxx, double gupxy, double gupxz,
    double gupyy, double gupyz, double gupzz,
    double Rxyxy, double Rxyxz, double Rxyyz,
    double Rxzxz, double Rxzyz, double Ryzyz,
    double& Rxx, double& Rxy, double& Rxz,
    double& Ryy, double& Ryz, double& Rzz
) {
    Rxx = gupyy * Rxyxy + gupyz * Rxyxz + gupyz * Rxyxz + gupzz * Rxzxz;
    Rxy = -gupxy * Rxyxy + gupyz * Rxyyz - gupxz * Rxyxz + gupzz * Rxzyz;
    Rxz = -gupxy * Rxyxz - gupyy * Rxyyz - gupxz * Rxzxz - gupyz * Rxzyz;
    Ryy = gupxx * Rxyxy - gupxz * Rxyyz - gupxz * Rxyyz + gupzz * Ryzyz;
    Ryz = gupxx * Rxyxz + gupxy * Rxyyz - gupxz * Rxzyz - gupyz * Ryzyz;
    Rzz = gupxx * Rxzxz + gupxy * Rxzyz + gupxy * Rxzyz + gupyy * Ryzyz;
}

// ============================================================
// Inline device function: compute Christoffel symbols and Ricci
// tensor from metric arrays at a single grid point.
// The metric arrays use the convention: dxx = gxx - 1 (tilde metric).
// Works with any metric arrays (BSSN physical or ADM).
// ============================================================
template<int GW>
__device__ __forceinline__ void d_adm_ricci_gamma_point(
    const int ex[3], const int pex[3],
    int i, int j, int k, size_t idx,
    const double* X, const double* Y, const double* Z,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    double& Gamxxx, double& Gamxxy, double& Gamxxz,
    double& Gamxyy, double& Gamxyz, double& Gamxzz,
    double& Gamyxx, double& Gamyxy, double& Gamyxz,
    double& Gamyyy, double& Gamyyz, double& Gamyzz,
    double& Gamzxx, double& Gamzxy, double& Gamzxz,
    double& Gamzyy, double& Gamzyz, double& Gamzzz,
    double& Rxx, double& Rxy, double& Rxz,
    double& Ryy, double& Ryz, double& Rzz,
    int symmetry
) {
    const double SYM = 1.0;
    const double ANTI = -1.0;
    const double ONE = 1.0;

    // First derivatives of metric (using d_fderivs_point_gw with padded arrays)
    double gxxx, gxxy, gxxz;
    d_fderivs_point_gw<GW>(ex, pex, dxx, &gxxx, &gxxy, &gxxz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, i, j, k);

    double gxyx, gxyy, gxyz;
    d_fderivs_point_gw<GW>(ex, pex, gxy, &gxyx, &gxyy, &gxyz, X, Y, Z, ANTI, ANTI, SYM,  symmetry, 0, i, j, k);

    double gxzx, gxzy, gxzz;
    d_fderivs_point_gw<GW>(ex, pex, gxz, &gxzx, &gxzy, &gxzz, X, Y, Z, ANTI, SYM,  ANTI, symmetry, 0, i, j, k);

    double gyyx, gyyy, gyyz;
    d_fderivs_point_gw<GW>(ex, pex, dyy, &gyyx, &gyyy, &gyyz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, i, j, k);

    double gyzx, gyzy, gyzz;
    d_fderivs_point_gw<GW>(ex, pex, gyz, &gyzx, &gyzy, &gyzz, X, Y, Z, SYM,  ANTI, ANTI, symmetry, 0, i, j, k);

    double gzzx, gzzy, gzzz;
    d_fderivs_point_gw<GW>(ex, pex, dzz, &gzzx, &gzzy, &gzzz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, i, j, k);

    // Kind 1 connection (Christoffel of first kind)
    double ass_Gamxxx, ass_Gamxxy, ass_Gamxxz;
    double ass_Gamxyy, ass_Gamxyz, ass_Gamxzz;
    double ass_Gamyxx, ass_Gamyxy, ass_Gamyxz;
    double ass_Gamyyy, ass_Gamyyz, ass_Gamyzz;
    double ass_Gamzxx, ass_Gamzxy, ass_Gamzxz;
    double ass_Gamzyy, ass_Gamzyz, ass_Gamzzz;

    d_kind1_connection(
        gxxx, gxyx, gxzx, gyyx, gyzx, gzzx,
        gxxy, gxyy, gxzy, gyyy, gyzy, gzzy,
        gxxz, gxyz, gxzz, gyyz, gyzz, gzzz,
        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz,
        ass_Gamxyy, ass_Gamxyz, ass_Gamxzz,
        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz,
        ass_Gamyyy, ass_Gamyyz, ass_Gamyzz,
        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz,
        ass_Gamzyy, ass_Gamzyz, ass_Gamzzz
    );

    // Load physical metric at this point: g_phys = dxx + 1
    double l_gxx = dxx[idx] + ONE;
    double l_gyy = dyy[idx] + ONE;
    double l_gzz = dzz[idx] + ONE;
    double l_gxy = gxy[idx];
    double l_gxz = gxz[idx];
    double l_gyz = gyz[idx];

    double detG = l_gxx * l_gyy * l_gzz + l_gxy * l_gyz * l_gxz + l_gxz * l_gxy * l_gyz
                - l_gxz * l_gyy * l_gxz - l_gxy * l_gxy * l_gzz - l_gxx * l_gyz * l_gyz;

    double gupxx =  (l_gyy * l_gzz - l_gyz * l_gyz) / detG;
    double gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) / detG;
    double gupxz =  (l_gxy * l_gyz - l_gyy * l_gxz) / detG;
    double gupyy =  (l_gxx * l_gzz - l_gxz * l_gxz) / detG;
    double gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) / detG;
    double gupzz =  (l_gxx * l_gyy - l_gxy * l_gxy) / detG;

    // Kind 2 connection (Christoffel of second kind)
    d_kind2_connection(
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz,
        ass_Gamxyy, ass_Gamxyz, ass_Gamxzz,
        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz,
        ass_Gamyyy, ass_Gamyyz, ass_Gamyzz,
        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz,
        ass_Gamzyy, ass_Gamzyz, ass_Gamzzz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz
    );

    // Second derivatives of metric
    double gxxxx, gxxxy, gxxxz, gxxyy, gxxyz, gxxzz;
    d_fdderivs_point_gw<GW>(ex, pex, dxx, &gxxxx, &gxxxy, &gxxxz, &gxxyy, &gxxyz, &gxxzz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, i, j, k);

    double gxyxx, gxyxy, gxyxz, gxyyy, gxyyz, gxyzz;
    d_fdderivs_point_gw<GW>(ex, pex, gxy, &gxyxx, &gxyxy, &gxyxz, &gxyyy, &gxyyz, &gxyzz, X, Y, Z, ANTI, ANTI, SYM,  symmetry, 0, i, j, k);

    double gxzxx, gxzxy, gxzxz, gxzyy, gxzyz, gxzzz;
    d_fdderivs_point_gw<GW>(ex, pex, gxz, &gxzxx, &gxzxy, &gxzxz, &gxzyy, &gxzyz, &gxzzz, X, Y, Z, ANTI, SYM,  ANTI, symmetry, 0, i, j, k);

    double gyyxx, gyyxy, gyyxz, gyyyy, gyyyz, gyyzz;
    d_fdderivs_point_gw<GW>(ex, pex, dyy, &gyyxx, &gyyxy, &gyyxz, &gyyyy, &gyyyz, &gyyzz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, i, j, k);

    double gyzxx, gyzxy, gyzxz, gyzyy, gyzyz, gyzzz;
    d_fdderivs_point_gw<GW>(ex, pex, gyz, &gyzxx, &gyzxy, &gyzxz, &gyzyy, &gyzyz, &gyzzz, X, Y, Z, SYM,  ANTI, ANTI, symmetry, 0, i, j, k);

    double gzzxx, gzzxy, gzzxz, gzzyy, gzzyz, gzzzz;
    d_fdderivs_point_gw<GW>(ex, pex, dzz, &gzzxx, &gzzxy, &gzzxz, &gzzyy, &gzzyz, &gzzzz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, i, j, k);

    // Riemann tensor
    double Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz;
    d_adm_riemann(
        gxxxx, gxxxy, gxxxz, gxxyy, gxxyz, gxxzz,
        gxyxx, gxyxy, gxyxz, gxyyy, gxyyz, gxyzz,
        gxzxx, gxzxy, gxzxz, gxzyy, gxzyz, gxzzz,
        gyyxx, gyyxy, gyyxz, gyyyy, gyyyz, gyyzz,
        gyzxx, gyzxy, gyzxz, gyzyy, gyzyz, gyzzz,
        gzzxx, gzzxy, gzzxz, gzzyy, gzzyz, gzzzz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz,
        ass_Gamxyy, ass_Gamxyz, ass_Gamxzz,
        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz,
        ass_Gamyyy, ass_Gamyyz, ass_Gamyzz,
        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz,
        ass_Gamzyy, ass_Gamzyz, ass_Gamzzz,
        Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz
    );

    // Ricci tensor
    d_adm_ricci(
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz
    );
}

// ============================================================
// Shell variant: compute Christoffel and Ricci at one point
// using shell-coordinate derivatives (crho, sigma, R).
// ============================================================
template<int GW>
__device__ __forceinline__ void d_adm_ricci_gamma_point_ss(
    const int ex[3], const int pex[3],
    int i, int j, int k, size_t idx,
    const double* crho, const double* sigma, const double* Rcoord,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    double& Gamxxx, double& Gamxxy, double& Gamxxz,
    double& Gamxyy, double& Gamxyz, double& Gamxzz,
    double& Gamyxx, double& Gamyxy, double& Gamyxz,
    double& Gamyyy, double& Gamyyz, double& Gamyzz,
    double& Gamzxx, double& Gamzxy, double& Gamzxz,
    double& Gamzyy, double& Gamzyz, double& Gamzzz,
    double& Rxx, double& Rxy, double& Rxz,
    double& Ryy, double& Ryz, double& Rzz,
    int symmetry, int sst
) {
    const double SYM = 1.0;
    const double ANTI = -1.0;
    const double ONE = 1.0;

    // First derivatives using shell stencil
    double gxxx, gxxy, gxxz;
    d_fderivs_sh_point_gw<GW>(ex, pex, dxx, &gxxx, &gxxy, &gxxz, crho, sigma, Rcoord,
                               SYM, SYM, SYM, symmetry, 0, sst, i, j, k);

    double gxyx, gxyy, gxyz;
    d_fderivs_sh_point_gw<GW>(ex, pex, gxy, &gxyx, &gxyy, &gxyz, crho, sigma, Rcoord,
                               ANTI, ANTI, SYM, symmetry, 0, sst, i, j, k);

    double gxzx, gxzy, gxzz;
    d_fderivs_sh_point_gw<GW>(ex, pex, gxz, &gxzx, &gxzy, &gxzz, crho, sigma, Rcoord,
                               ANTI, SYM, ANTI, symmetry, 0, sst, i, j, k);

    double gyyx, gyyy, gyyz;
    d_fderivs_sh_point_gw<GW>(ex, pex, dyy, &gyyx, &gyyy, &gyyz, crho, sigma, Rcoord,
                               SYM, SYM, SYM, symmetry, 0, sst, i, j, k);

    double gyzx, gyzy, gyzz;
    d_fderivs_sh_point_gw<GW>(ex, pex, gyz, &gyzx, &gyzy, &gyzz, crho, sigma, Rcoord,
                               SYM, ANTI, ANTI, symmetry, 0, sst, i, j, k);

    double gzzx, gzzy, gzzz;
    d_fderivs_sh_point_gw<GW>(ex, pex, dzz, &gzzx, &gzzy, &gzzz, crho, sigma, Rcoord,
                               SYM, SYM, SYM, symmetry, 0, sst, i, j, k);

    // Kind 1 connection (same algebra)
    double ass_Gamxxx, ass_Gamxxy, ass_Gamxxz;
    double ass_Gamxyy, ass_Gamxyz, ass_Gamxzz;
    double ass_Gamyxx, ass_Gamyxy, ass_Gamyxz;
    double ass_Gamyyy, ass_Gamyyz, ass_Gamyzz;
    double ass_Gamzxx, ass_Gamzxy, ass_Gamzxz;
    double ass_Gamzyy, ass_Gamzyz, ass_Gamzzz;

    d_kind1_connection(
        gxxx, gxyx, gxzx, gyyx, gyzx, gzzx,
        gxxy, gxyy, gxzy, gyyy, gyzy, gzzy,
        gxxz, gxyz, gxzz, gyyz, gyzz, gzzz,
        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz,
        ass_Gamxyy, ass_Gamxyz, ass_Gamxzz,
        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz,
        ass_Gamyyy, ass_Gamyyz, ass_Gamyzz,
        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz,
        ass_Gamzyy, ass_Gamzyz, ass_Gamzzz
    );

    // Metric and inverse
    double l_gxx = dxx[idx] + ONE;
    double l_gyy = dyy[idx] + ONE;
    double l_gzz = dzz[idx] + ONE;
    double l_gxy = gxy[idx];
    double l_gxz = gxz[idx];
    double l_gyz = gyz[idx];

    double detG = l_gxx * l_gyy * l_gzz + l_gxy * l_gyz * l_gxz + l_gxz * l_gxy * l_gyz
                - l_gxz * l_gyy * l_gxz - l_gxy * l_gxy * l_gzz - l_gxx * l_gyz * l_gyz;

    double gupxx =  (l_gyy * l_gzz - l_gyz * l_gyz) / detG;
    double gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) / detG;
    double gupxz =  (l_gxy * l_gyz - l_gyy * l_gxz) / detG;
    double gupyy =  (l_gxx * l_gzz - l_gxz * l_gxz) / detG;
    double gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) / detG;
    double gupzz =  (l_gxx * l_gyy - l_gxy * l_gxy) / detG;

    // Kind 2 connection
    d_kind2_connection(
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz,
        ass_Gamxyy, ass_Gamxyz, ass_Gamxzz,
        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz,
        ass_Gamyyy, ass_Gamyyz, ass_Gamyzz,
        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz,
        ass_Gamzyy, ass_Gamzyz, ass_Gamzzz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz
    );

    // Second derivatives using shell stencil
    double gxxxx, gxxxy, gxxxz, gxxyy, gxxyz, gxxzz;
    d_fdderivs_sh_point_gw<GW>(ex, pex, dxx,
        &gxxxx, &gxxxy, &gxxxz, &gxxyy, &gxxyz, &gxxzz,
        crho, sigma, Rcoord,
        SYM, SYM, SYM, symmetry, 0, sst, i, j, k);

    double gxyxx, gxyxy, gxyxz, gxyyy, gxyyz, gxyzz;
    d_fdderivs_sh_point_gw<GW>(ex, pex, gxy,
        &gxyxx, &gxyxy, &gxyxz, &gxyyy, &gxyyz, &gxyzz,
        crho, sigma, Rcoord,
        ANTI, ANTI, SYM, symmetry, 0, sst, i, j, k);

    double gxzxx, gxzxy, gxzxz, gxzyy, gxzyz, gxzzz;
    d_fdderivs_sh_point_gw<GW>(ex, pex, gxz,
        &gxzxx, &gxzxy, &gxzxz, &gxzyy, &gxzyz, &gxzzz,
        crho, sigma, Rcoord,
        ANTI, SYM, ANTI, symmetry, 0, sst, i, j, k);

    double gyyxx, gyyxy, gyyxz, gyyyy, gyyyz, gyyzz;
    d_fdderivs_sh_point_gw<GW>(ex, pex, dyy,
        &gyyxx, &gyyxy, &gyyxz, &gyyyy, &gyyyz, &gyyzz,
        crho, sigma, Rcoord,
        SYM, SYM, SYM, symmetry, 0, sst, i, j, k);

    double gyzxx, gyzxy, gyzxz, gyzyy, gyzyz, gyzzz;
    d_fdderivs_sh_point_gw<GW>(ex, pex, gyz,
        &gyzxx, &gyzxy, &gyzxz, &gyzyy, &gyzyz, &gyzzz,
        crho, sigma, Rcoord,
        SYM, ANTI, ANTI, symmetry, 0, sst, i, j, k);

    double gzzxx, gzzxy, gzzxz, gzzyy, gzzyz, gzzzz;
    d_fdderivs_sh_point_gw<GW>(ex, pex, dzz,
        &gzzxx, &gzzxy, &gzzxz, &gzzyy, &gzzyz, &gzzzz,
        crho, sigma, Rcoord,
        SYM, SYM, SYM, symmetry, 0, sst, i, j, k);

    // Riemann
    double Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz;
    d_adm_riemann(
        gxxxx, gxxxy, gxxxz, gxxyy, gxxyz, gxxzz,
        gxyxx, gxyxy, gxyxz, gxyyy, gxyyz, gxyzz,
        gxzxx, gxzxy, gxzxz, gxzyy, gxzyz, gxzzz,
        gyyxx, gyyxy, gyyxz, gyyyy, gyyyz, gyyzz,
        gyzxx, gyzxy, gyzxz, gyzyy, gyzyz, gyzzz,
        gzzxx, gzzxy, gzzxz, gzzyy, gzzyz, gzzzz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz,
        ass_Gamxyy, ass_Gamxyz, ass_Gamxzz,
        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz,
        ass_Gamyyy, ass_Gamyyz, ass_Gamyzz,
        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz,
        ass_Gamzyy, ass_Gamzyz, ass_Gamzzz,
        Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz
    );

    // Ricci
    d_adm_ricci(
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz
    );
}

#endif // ADM_RICCI_GAMMA_GPU_CUH
