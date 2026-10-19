#include "getnp4.h"

#include <cmath>
#include "derivatives.h"
#include "macrodef.h"
#include "diff_new_gpu.cuh"
#include "point_diff_new_sh_gpu.cuh"
#include "bssn2adm_gpu.cuh"
#include "adm_ricci_gamma_gpu.cuh"
#include "point_diff_new_sh_gpu.cuh"

// ---------------------------------------------------------
// Pre-processing kernel for ABV==1 path:
// Convert BSSN to ADM metric and compute ADM Christoffel + Ricci.
// This kernel must be called before getnp4_kernel when ABV==1.
// Writes Christoffel symbols and Ricci tensor into Gam/R output arrays.
// ---------------------------------------------------------
// Phase 1 kernel: BSSN -> ADM metric conversion.
// Writes ADM metric in (dxx+1, gxy, ...) convention to temp output arrays.
template<int GW>
__global__ void bssn2adm_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* chi, const double* trK,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    const double* Axx, const double* Axy, const double* Axz,
    const double* Ayy, const double* Ayz, const double* Azz,
    double* temp_adm_dxx, double* temp_adm_gxy, double* temp_adm_gxz,
    double* temp_adm_dyy, double* temp_adm_gyz, double* temp_adm_dzz
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    size_t idx = i + j * pex0 + k * (size_t)pex0 * pex1;

    double l_gxx = dxx[idx] + 1.0;
    double l_gyy = dyy[idx] + 1.0;
    double l_gzz = dzz[idx] + 1.0;
    double l_gxy = gxy[idx];
    double l_gxz = gxz[idx];
    double l_gyz = gyz[idx];

    double adm_gxx, adm_gxy, adm_gxz, adm_gyy, adm_gyz, adm_gzz;
    double Kxx, Kxy, Kxz, Kyy, Kyz, Kzz;
    d_bssn2adm(chi[idx] + 1.0, trK[idx],
        l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz,
        Axx[idx], Axy[idx], Axz[idx], Ayy[idx], Ayz[idx], Azz[idx],
        adm_gxx, adm_gxy, adm_gxz, adm_gyy, adm_gyz, adm_gzz,
        Kxx, Kxy, Kxz, Kyy, Kyz, Kzz);

    // Store ADM metric in (dxx+1, gxy, ...) convention
    temp_adm_dxx[idx] = adm_gxx - 1.0;
    temp_adm_gxy[idx] = adm_gxy;
    temp_adm_gxz[idx] = adm_gxz;
    temp_adm_dyy[idx] = adm_gyy - 1.0;
    temp_adm_gyz[idx] = adm_gyz;
    temp_adm_dzz[idx] = adm_gzz - 1.0;
}

template<int GW>
__global__ void adm_ricci_from_adm_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* X, const double* Y, const double* Z,
    const double* temp_dxx, const double* temp_gxy, const double* temp_gxz,
    const double* temp_dyy, const double* temp_gyz, const double* temp_dzz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    int symmetry
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    size_t idx = i + j * pex0 + k * (size_t)pex0 * pex1;

    const int ex[3]  = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    d_adm_ricci_gamma_point<GW>(
        ex, pex, i, j, k, idx,
        X, Y, Z,
        temp_dxx, temp_gxy, temp_gxz,
        temp_dyy, temp_gyz, temp_dzz,
        Gamxxx[idx], Gamxxy[idx], Gamxxz[idx],
        Gamxyy[idx], Gamxyz[idx], Gamxzz[idx],
        Gamyxx[idx], Gamyxy[idx], Gamyxz[idx],
        Gamyyy[idx], Gamyyz[idx], Gamyzz[idx],
        Gamzxx[idx], Gamzxy[idx], Gamzxz[idx],
        Gamzyy[idx], Gamzyz[idx], Gamzzz[idx],
        Rxx[idx], Rxy[idx], Rxz[idx],
        Ryy[idx], Ryz[idx], Rzz[idx],
        symmetry
    );
}

// ---------------------------------------------------------
// Main getnp4 kernel
// ---------------------------------------------------------
template<int GW, int TETRAD>
__global__ void getnp4_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* X, const double* Y, const double* Z,
    const double* chi, const double* trK,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    const double* Axx, const double* Axy, const double* Axz,
    const double* Ayy, const double* Ayz, const double* Azz,
    const double* Gamxxx, const double* Gamxxy, const double* Gamxxz,
    const double* Gamxyy, const double* Gamxyz, const double* Gamxzz,
    const double* Gamyxx, const double* Gamyxy, const double* Gamyxz,
    const double* Gamyyy, const double* Gamyyz, const double* Gamyzz,
    const double* Gamzxx, const double* Gamzxy, const double* Gamzxz,
    const double* Gamzyy, const double* Gamzyz, const double* Gamzzz,
    const double* Rxx, const double* Rxy, const double* Rxz,
    const double* Ryy, const double* Ryz, const double* Rzz,
    double* Rpsi4, double* Ipsi4,
    int symmetry
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    size_t idx = i + j * pex0 + k * (size_t)pex0 * pex1;

    const int ex[3]  = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    const double ZEO = 0.0;
    const double ONE = 1.0;
    const double TWO = 2.0;
    const double F1o3 = 1.0 / 3.0;
    const double SYM = 1.0;
    const double ANTI = -1.0;
    const double TINYRR = 1e-14;

    double cx = X[i];
    double cy = Y[j];
    double cz = Z[k];

    // Load physical metric
    double l_gxx = dxx[idx] + ONE;
    double l_gyy = dyy[idx] + ONE;
    double l_gzz = dzz[idx] + ONE;
    double l_gxy = gxy[idx];
    double l_gxz = gxz[idx];
    double l_gyz = gyz[idx];

    double l_chi = chi[idx];
    double chipn1 = l_chi + ONE;
    double f = ONE / chipn1;

    // Inverse metric
    double detG = l_gxx * l_gyy * l_gzz + l_gxy * l_gyz * l_gxz + l_gxz * l_gxy * l_gyz -
                  l_gxz * l_gyy * l_gxz - l_gxy * l_gxy * l_gzz - l_gxx * l_gyz * l_gyz;

    double l_gupxx =  (l_gyy * l_gzz - l_gyz * l_gyz) / detG;
    double l_gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) / detG;
    double l_gupxz =  (l_gxy * l_gyz - l_gyy * l_gxz) / detG;
    double l_gupyy =  (l_gxx * l_gzz - l_gxz * l_gxz) / detG;
    double l_gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) / detG;
    double l_gupzz =  (l_gxx * l_gyy - l_gxy * l_gxy) / detG;

    // Initialize tetrad vectors U, V, W
    double vx, vy, vz, ux, uy, uz, wx, wy, wz;

    if (fabs(cx) < TINYRR && fabs(cy) < TINYRR && fabs(cz) < TINYRR) {
        vx = TINYRR; vy = TINYRR; vz = TINYRR;
    } else {
        vx = cx; vy = cy; vz = cz;
    }

    if (fabs(cx) < TINYRR && fabs(cy) < TINYRR) {
        ux = -TINYRR; uy = TINYRR; uz = ZEO;
        wx = TINYRR * cz; wy = TINYRR * cz; wz = -TWO * TINYRR * TINYRR;
    } else {
        ux = -cy; uy = cx; uz = ZEO;
        wx = cx * cz; wy = cy * cz; wz = -(cx * cx + cy * cy);
    }

    // Gram-Schmidt orthogonalization (order depends on tetradtype)
    if constexpr (TETRAD == 2) {
        double fxt = vx, fyt = vy, fzt = vz;
        vx = l_gupxx * fxt + l_gupxy * fyt + l_gupxz * fzt;
        vy = l_gupxy * fxt + l_gupyy * fyt + l_gupyz * fzt;
        vz = l_gupxz * fxt + l_gupyz * fyt + l_gupzz * fzt;

        double temp;
        temp = l_gxx*vx*vx + l_gyy*vy*vy + l_gzz*vz*vz + (l_gxy*vx*vy + l_gxz*vx*vz + l_gyz*vy*vz)*TWO;
        temp = sqrt(temp*f);
        vx /= temp; vy /= temp; vz /= temp;

        temp = l_gxx*vx*ux + l_gxy*vx*uy + l_gxz*vx*uz +
               l_gxy*vy*ux + l_gyy*vy*uy + l_gyz*vy*uz +
               l_gxz*vz*ux + l_gyz*vz*uy + l_gzz*vz*uz;
        temp *= f;
        ux -= temp*vx; uy -= temp*vy; uz -= temp*vz;

        temp = l_gxx*ux*ux + l_gyy*uy*uy + l_gzz*uz*uz + (l_gxy*ux*uy + l_gxz*ux*uz + l_gyz*uy*uz)*TWO;
        temp = sqrt(temp*f);
        ux /= temp; uy /= temp; uz /= temp;

        temp = l_gxx*vx*wx + l_gxy*vx*wy + l_gxz*vx*wz +
               l_gxy*vy*wx + l_gyy*vy*wy + l_gyz*vy*wz +
               l_gxz*vz*wx + l_gyz*vz*wy + l_gzz*vz*wz;
        temp *= f;
        wx -= temp*vx; wy -= temp*vy; wz -= temp*vz;

        temp = l_gxx*ux*wx + l_gxy*ux*wy + l_gxz*ux*wz +
               l_gxy*uy*wx + l_gyy*uy*wy + l_gyz*uy*wz +
               l_gxz*uz*wx + l_gyz*uz*wy + l_gzz*uz*wz;
        temp *= f;
        wx -= temp*ux; wy -= temp*uy; wz -= temp*uz;

        temp = l_gxx*wx*wx + l_gyy*wy*wy + l_gzz*wz*wz + (l_gxy*wx*wy + l_gxz*wx*wz + l_gyz*wy*wz)*TWO;
        temp = sqrt(temp*f);
        wx /= temp; wy /= temp; wz /= temp;

    } else if constexpr (TETRAD == 0) {
        double temp;
        temp = l_gxx*vx*vx + l_gyy*vy*vy + l_gzz*vz*vz + (l_gxy*vx*vy + l_gxz*vx*vz + l_gyz*vy*vz)*TWO;
        temp = sqrt(temp*f);
        vx /= temp; vy /= temp; vz /= temp;

        temp = l_gxx*vx*ux + l_gxy*vx*uy + l_gxz*vx*uz +
               l_gxy*vy*ux + l_gyy*vy*uy + l_gyz*vy*uz +
               l_gxz*vz*ux + l_gyz*vz*uy + l_gzz*vz*uz;
        temp *= f;
        ux -= temp*vx; uy -= temp*vy; uz -= temp*vz;

        temp = l_gxx*ux*ux + l_gyy*uy*uy + l_gzz*uz*uz + (l_gxy*ux*uy + l_gxz*ux*uz + l_gyz*uy*uz)*TWO;
        temp = sqrt(temp*f);
        ux /= temp; uy /= temp; uz /= temp;

        temp = l_gxx*vx*wx + l_gxy*vx*wy + l_gxz*vx*wz +
               l_gxy*vy*wx + l_gyy*vy*wy + l_gyz*vy*wz +
               l_gxz*vz*wx + l_gyz*vz*wy + l_gzz*vz*wz;
        temp *= f;
        wx -= temp*vx; wy -= temp*vy; wz -= temp*vz;

        temp = l_gxx*ux*wx + l_gxy*ux*wy + l_gxz*ux*wz +
               l_gxy*uy*wx + l_gyy*uy*wy + l_gyz*uy*wz +
               l_gxz*uz*wx + l_gyz*uz*wy + l_gzz*uz*wz;
        temp *= f;
        wx -= temp*ux; wy -= temp*uy; wz -= temp*uz;

        temp = l_gxx*wx*wx + l_gyy*wy*wy + l_gzz*wz*wz + (l_gxy*wx*wy + l_gxz*wx*wz + l_gyz*wy*wz)*TWO;
        temp = sqrt(temp*f);
        wx /= temp; wy /= temp; wz /= temp;

    } else if constexpr (TETRAD == 1) {
        double temp;

        temp = l_gxx*wx*wx + l_gyy*wy*wy + l_gzz*wz*wz + (l_gxy*wx*wy + l_gxz*wx*wz + l_gyz*wy*wz)*TWO;
        temp = sqrt(temp*f);
        wx /= temp; wy /= temp; wz /= temp;

        temp = l_gxx*wx*ux + l_gxy*wx*uy + l_gxz*wx*uz +
               l_gxy*wy*ux + l_gyy*wy*uy + l_gyz*wy*uz +
               l_gxz*wz*ux + l_gyz*wz*uy + l_gzz*wz*uz;
        temp *= f;
        ux -= temp*wx; uy -= temp*wy; uz -= temp*wz;

        temp = l_gxx*ux*ux + l_gyy*uy*uy + l_gzz*uz*uz + (l_gxy*ux*uy + l_gxz*ux*uz + l_gyz*uy*uz)*TWO;
        temp = sqrt(temp*f);
        ux /= temp; uy /= temp; uz /= temp;

        temp = l_gxx*vx*wx + l_gxy*vx*wy + l_gxz*vx*wz +
               l_gxy*vy*wx + l_gyy*vy*wy + l_gyz*vy*wz +
               l_gxz*vz*wx + l_gyz*vz*wy + l_gzz*vz*wz;
        temp *= f;
        vx -= temp*wx; vy -= temp*wy; vz -= temp*wz;

        temp = l_gxx*ux*vx + l_gxy*ux*vy + l_gxz*ux*vz +
               l_gxy*uy*vx + l_gyy*uy*vy + l_gyz*uy*vz +
               l_gxz*uz*vx + l_gyz*uz*vy + l_gzz*uz*vz;
        temp *= f;
        vx -= temp*ux; vy -= temp*uy; vz -= temp*uz;

        temp = l_gxx*vx*vx + l_gyy*vy*vy + l_gzz*vz*vz + (l_gxy*vx*vy + l_gxz*vx*vz + l_gyz*vy*vz)*TWO;
        temp = sqrt(temp*f);
        vx /= temp; vy /= temp; vz /= temp;
    }

    // Compute derivatives of A_ij, chi, trK
    double Axxx, Axxy, Axxz;
    d_fderivs_point_gw<GW>(ex, pex, Axx, &Axxx, &Axxy, &Axxz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, i, j, k);

    double Axyx, Axyy, Axyz;
    d_fderivs_point_gw<GW>(ex, pex, Axy, &Axyx, &Axyy, &Axyz, X, Y, Z, ANTI, ANTI, SYM,  symmetry, 0, i, j, k);

    double Axzx, Axzy, Axzz;
    d_fderivs_point_gw<GW>(ex, pex, Axz, &Axzx, &Axzy, &Axzz, X, Y, Z, ANTI, SYM,  ANTI, symmetry, 0, i, j, k);

    double Ayyx, Ayyy, Ayyz;
    d_fderivs_point_gw<GW>(ex, pex, Ayy, &Ayyx, &Ayyy, &Ayyz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, i, j, k);

    double Ayzx, Ayzy, Ayzz;
    d_fderivs_point_gw<GW>(ex, pex, Ayz, &Ayzx, &Ayzy, &Ayzz, X, Y, Z, SYM,  ANTI, ANTI, symmetry, 0, i, j, k);

    double Azzx, Azzy, Azzz;
    d_fderivs_point_gw<GW>(ex, pex, Azz, &Azzx, &Azzy, &Azzz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, i, j, k);

    double chix, chiy, chiz;
    d_fderivs_point_gw<GW>(ex, pex, chi, &chix, &chiy, &chiz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);

    double fx_trk, fy_trk, fz_trk;
    d_fderivs_point_gw<GW>(ex, pex, trK, &fx_trk, &fy_trk, &fz_trk, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);

    // Load traceless extrinsic curvature at current point
    double l_Axx = Axx[idx];
    double l_Axy = Axy[idx];
    double l_Axz = Axz[idx];
    double l_Ayy = Ayy[idx];
    double l_Ayz = Ayz[idx];
    double l_Azz = Azz[idx];

    // Compute D_k K_ij (covariant derivative corrections using Christoffel symbols)
    Axxx -= (Gamxxx[idx]*l_Axx + Gamyxx[idx]*l_Axy + Gamzxx[idx]*l_Axz)*TWO + chix/chipn1*l_Axx - F1o3*l_gxx*fx_trk;
    Axxy -= (Gamxxy[idx]*l_Axx + Gamyxy[idx]*l_Axy + Gamzxy[idx]*l_Axz)*TWO + chiy/chipn1*l_Axx - F1o3*l_gxx*fy_trk;
    Axxz -= (Gamxxz[idx]*l_Axx + Gamyxz[idx]*l_Axy + Gamzxz[idx]*l_Axz)*TWO + chiz/chipn1*l_Axx - F1o3*l_gxx*fz_trk;

    Ayyx -= (Gamxxy[idx]*l_Axy + Gamyxy[idx]*l_Ayy + Gamzxy[idx]*l_Ayz)*TWO + chix/chipn1*l_Ayy - F1o3*l_gyy*fx_trk;
    Ayyy -= (Gamxyy[idx]*l_Axy + Gamyyy[idx]*l_Ayy + Gamzyy[idx]*l_Ayz)*TWO + chiy/chipn1*l_Ayy - F1o3*l_gyy*fy_trk;
    Ayyz -= (Gamxyz[idx]*l_Axy + Gamyyz[idx]*l_Ayy + Gamzyz[idx]*l_Ayz)*TWO + chiz/chipn1*l_Ayy - F1o3*l_gyy*fz_trk;

    Azzx -= (Gamxxz[idx]*l_Axz + Gamyxz[idx]*l_Ayz + Gamzxz[idx]*l_Azz)*TWO + chix/chipn1*l_Azz - F1o3*l_gzz*fx_trk;
    Azzy -= (Gamxyz[idx]*l_Axz + Gamyyz[idx]*l_Ayz + Gamzyz[idx]*l_Azz)*TWO + chiy/chipn1*l_Azz - F1o3*l_gzz*fy_trk;
    Azzz -= (Gamxzz[idx]*l_Axz + Gamyzz[idx]*l_Ayz + Gamzzz[idx]*l_Azz)*TWO + chiz/chipn1*l_Azz - F1o3*l_gzz*fz_trk;

    Axyx -= (Gamxxy[idx]*l_Axx + Gamyxy[idx]*l_Axy + Gamzxy[idx]*l_Axz + Gamxxx[idx]*l_Axy + Gamyxx[idx]*l_Ayy + Gamzxx[idx]*l_Ayz) + chix/chipn1*l_Axy - F1o3*l_gxy*fx_trk;
    Axyy -= (Gamxyy[idx]*l_Axx + Gamyyy[idx]*l_Axy + Gamzyy[idx]*l_Axz + Gamxxy[idx]*l_Axy + Gamyxy[idx]*l_Ayy + Gamzxy[idx]*l_Ayz) + chiy/chipn1*l_Axy - F1o3*l_gxy*fy_trk;
    Axyz -= (Gamxyz[idx]*l_Axx + Gamyyz[idx]*l_Axy + Gamzyz[idx]*l_Axz + Gamxxz[idx]*l_Axy + Gamyxz[idx]*l_Ayy + Gamzxz[idx]*l_Ayz) + chiz/chipn1*l_Axy - F1o3*l_gxy*fz_trk;

    Axzx -= (Gamxxz[idx]*l_Axx + Gamyxz[idx]*l_Axy + Gamzxz[idx]*l_Axz + Gamxxx[idx]*l_Axz + Gamyxx[idx]*l_Ayz + Gamzxx[idx]*l_Azz) + chix/chipn1*l_Axz - F1o3*l_gxz*fx_trk;
    Axzy -= (Gamxyz[idx]*l_Axx + Gamyyz[idx]*l_Axy + Gamzyz[idx]*l_Axz + Gamxxy[idx]*l_Axz + Gamyxy[idx]*l_Ayz + Gamzxy[idx]*l_Azz) + chiy/chipn1*l_Axz - F1o3*l_gxz*fy_trk;
    Axzz -= (Gamxzz[idx]*l_Axx + Gamyzz[idx]*l_Axy + Gamzzz[idx]*l_Axz + Gamxxz[idx]*l_Axz + Gamyxz[idx]*l_Ayz + Gamzxz[idx]*l_Azz) + chiz/chipn1*l_Axz - F1o3*l_gxz*fz_trk;

    Ayzx -= (Gamxxz[idx]*l_Axy + Gamyxz[idx]*l_Ayy + Gamzxz[idx]*l_Ayz + Gamxxy[idx]*l_Axz + Gamyxy[idx]*l_Ayz + Gamzxy[idx]*l_Azz) + chix/chipn1*l_Ayz - F1o3*l_gyz*fx_trk;
    Ayzy -= (Gamxyz[idx]*l_Axy + Gamyyz[idx]*l_Ayy + Gamzyz[idx]*l_Ayz + Gamxyy[idx]*l_Axz + Gamyyy[idx]*l_Ayz + Gamzyy[idx]*l_Azz) + chiy/chipn1*l_Ayz - F1o3*l_gyz*fy_trk;
    Ayzz -= (Gamxzz[idx]*l_Axy + Gamyzz[idx]*l_Ayy + Gamzzz[idx]*l_Ayz + Gamxyz[idx]*l_Axz + Gamyyz[idx]*l_Ayz + Gamzyz[idx]*l_Azz) + chiz/chipn1*l_Ayz - F1o3*l_gyz*fz_trk;

    // Symmetrize B_ij = v^k (D_k K_ij - D_j K_ik)
    double Bxx = (vy*(Axxy - Axyx) + vz*(Axxz - Axzx)) * f;
    double Byy = (vx*(Ayyx - Axyy) + vz*(Ayyz - Ayzy)) * f;
    double Bzz = (vx*(Azzx - Axzz) + vy*(Azzy - Ayzz)) * f;
    double Bxy = (vx*(Axyx - (Axxy+Axyx)/TWO) + vy*(Axyy-Ayyx)/TWO + vz*(Axyz - (Axzy+Ayzx)/TWO)) * f;
    double Bxz = (vx*(Axzx - (Axxz+Axzx)/TWO) + vy*(Axzy - (Axyz+Ayzx)/TWO) + vz*(Axzz-Azzx)/TWO) * f;
    double Byz = (vx*(Ayzx - (Axyz+Axzy)/TWO) + vy*(Ayzy - (Ayyz+Ayzy)/TWO) + vz*(Ayzz-Azzy)/TWO) * f;

    double l_trK = trK[idx];

    // Physical K_ij = A_ij + (1/3)*trK*g_ij
    double Kxx_c = l_Axx + F1o3 * l_trK * l_gxx;
    double Kxy_c = l_Axy + F1o3 * l_trK * l_gxy;
    double Kxz_c = l_Axz + F1o3 * l_trK * l_gxz;
    double Kyy_c = l_Ayy + F1o3 * l_trK * l_gyy;
    double Kyz_c = l_Ayz + F1o3 * l_trK * l_gyz;
    double Kzz_c = l_Azz + F1o3 * l_trK * l_gzz;

    // Compute E_ij = K_ik K^k_j
    double Exx = l_gupxx*Kxx_c*Kxx_c + l_gupyy*Kxy_c*Kxy_c + l_gupzz*Kxz_c*Kxz_c +
                 TWO*(l_gupxy*Kxx_c*Kxy_c + l_gupxz*Kxx_c*Kxz_c + l_gupyz*Kxy_c*Kxz_c);
    double Eyy = l_gupxx*Kxy_c*Kxy_c + l_gupyy*Kyy_c*Kyy_c + l_gupzz*Kyz_c*Kyz_c +
                 TWO*(l_gupxy*Kxy_c*Kyy_c + l_gupxz*Kxy_c*Kyz_c + l_gupyz*Kyy_c*Kyz_c);
    double Ezz = l_gupxx*Kxz_c*Kxz_c + l_gupyy*Kyz_c*Kyz_c + l_gupzz*Kzz_c*Kzz_c +
                 TWO*(l_gupxy*Kxz_c*Kyz_c + l_gupxz*Kxz_c*Kzz_c + l_gupyz*Kyz_c*Kzz_c);

    double Exy = l_gupxx*Kxx_c*Kxy_c + l_gupyy*Kxy_c*Kyy_c + l_gupzz*Kxz_c*Kyz_c +
                 l_gupxy*(Kxx_c*Kyy_c + Kxy_c*Kxy_c) +
                 l_gupxz*(Kxx_c*Kyz_c + Kxz_c*Kxy_c) +
                 l_gupyz*(Kxy_c*Kyz_c + Kxz_c*Kyy_c);
    double Exz = l_gupxx*Kxx_c*Kxz_c + l_gupyy*Kxy_c*Kyz_c + l_gupzz*Kxz_c*Kzz_c +
                 l_gupxy*(Kxx_c*Kyz_c + Kxy_c*Kxz_c) +
                 l_gupxz*(Kxx_c*Kzz_c + Kxz_c*Kxz_c) +
                 l_gupyz*(Kxy_c*Kzz_c + Kxz_c*Kyz_c);
    double Eyz = l_gupxx*Kxy_c*Kxz_c + l_gupyy*Kyy_c*Kyz_c + l_gupzz*Kyz_c*Kzz_c +
                 l_gupxy*(Kxy_c*Kyz_c + Kyy_c*Kxz_c) +
                 l_gupxz*(Kxy_c*Kzz_c + Kyz_c*Kxz_c) +
                 l_gupyz*(Kyy_c*Kzz_c + Kyz_c*Kyz_c);

    Exx = Rxx[idx] - (Exx - Kxx_c*l_trK)*f - Bxx;
    Exy = Rxy[idx] - (Exy - Kxy_c*l_trK)*f - Bxy;
    Exz = Rxz[idx] - (Exz - Kxz_c*l_trK)*f - Bxz;
    Eyy = Ryy[idx] - (Eyy - Kyy_c*l_trK)*f - Byy;
    Eyz = Ryz[idx] - (Eyz - Kyz_c*l_trK)*f - Byz;
    Ezz = Rzz[idx] - (Ezz - Kzz_c*l_trK)*f - Bzz;

    // Projection tensors: uuww^ij = u^i u^j - w^i w^j, uw^ij = u^i w^j + w^i u^j
    double uuwwxx = ux * ux - wx * wx;
    double uuwwxy = ux * uy - wx * wy;
    double uuwwxz = ux * uz - wx * wz;
    double uuwwyy = uy * uy - wy * wy;
    double uuwwyz = uy * uz - wy * wz;
    double uuwwzz = uz * uz - wz * wz;

    double uwxx = ux * wx + wx * ux;
    double uwxy = ux * wy + wx * uy;
    double uwxz = ux * wz + wx * uz;
    double uwyy = uy * wy + wy * uy;
    double uwyz = uy * wz + wy * uz;
    double uwzz = uz * wz + wz * uz;

    // Combine Psi4
    double l_Rpsi4 = Exx*uuwwxx + Eyy*uuwwyy + Ezz*uuwwzz + (Exy*uuwwxy + Exz*uuwwxz + Eyz*uuwwyz)*TWO;
    double l_Ipsi4 = Exx*uwxx + Eyy*uwyy + Ezz*uwzz + (Exy*uwxy + Exz*uwxz + Eyz*uwyz)*TWO;

    Rpsi4[idx] = -l_Rpsi4 / TWO;
    Ipsi4[idx] = -l_Ipsi4 / TWO;
}

// ---------------------------------------------------------
// Template launch function
// ---------------------------------------------------------
template<int GW, int TETRAD>
void gpu_getnp4_launch_impl(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    const double* d_Gamxxx, const double* d_Gamxxy, const double* d_Gamxxz,
    const double* d_Gamxyy, const double* d_Gamxyz, const double* d_Gamxzz,
    const double* d_Gamyxx, const double* d_Gamyxy, const double* d_Gamyxz,
    const double* d_Gamyyy, const double* d_Gamyyz, const double* d_Gamyzz,
    const double* d_Gamzxx, const double* d_Gamzxy, const double* d_Gamzxz,
    const double* d_Gamzyy, const double* d_Gamzyz, const double* d_Gamzzz,
    const double* d_Rxx, const double* d_Rxy, const double* d_Rxz,
    const double* d_Ryy, const double* d_Ryz, const double* d_Rzz,
    double* d_Rpsi4, double* d_Ipsi4,
    int symmetry
) {
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    getnp4_kernel<GW, TETRAD><<<gridDim, blockDim, 0, stream>>>(
        ex[0], ex[1], ex[2],
        pex[0], pex[1], pex[2],
        d_X, d_Y, d_Z,
        d_chi, d_trK,
        d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
        d_Axx, d_Axy, d_Axz, d_Ayy, d_Ayz, d_Azz,
        d_Gamxxx, d_Gamxxy, d_Gamxxz,
        d_Gamxyy, d_Gamxyz, d_Gamxzz,
        d_Gamyxx, d_Gamyxy, d_Gamyxz,
        d_Gamyyy, d_Gamyyz, d_Gamyzz,
        d_Gamzxx, d_Gamzxy, d_Gamzxz,
        d_Gamzyy, d_Gamzyz, d_Gamzzz,
        d_Rxx, d_Rxy, d_Rxz, d_Ryy, d_Ryz, d_Rzz,
        d_Rpsi4, d_Ipsi4,
        symmetry
    );
}

// Backward-compatible wrapper: default to GW=3, tetradtype=2
void gpu_getnp4_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    const double* d_Gamxxx, const double* d_Gamxxy, const double* d_Gamxxz,
    const double* d_Gamxyy, const double* d_Gamxyz, const double* d_Gamxzz,
    const double* d_Gamyxx, const double* d_Gamyxy, const double* d_Gamyxz,
    const double* d_Gamyyy, const double* d_Gamyyz, const double* d_Gamyzz,
    const double* d_Gamzxx, const double* d_Gamzxy, const double* d_Gamzxz,
    const double* d_Gamzyy, const double* d_Gamzyz, const double* d_Gamzzz,
    const double* d_Rxx, const double* d_Rxy, const double* d_Rxz,
    const double* d_Ryy, const double* d_Ryz, const double* d_Rzz,
    double* d_Rpsi4, double* d_Ipsi4,
    int symmetry
) {
    gpu_getnp4_launch_impl<ghost_width, tetradtype>(
        stream, ex, pex,
        d_X, d_Y, d_Z,
        d_chi, d_trK,
        d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
        d_Axx, d_Axy, d_Axz, d_Ayy, d_Ayz, d_Azz,
        d_Gamxxx, d_Gamxxy, d_Gamxxz,
        d_Gamxyy, d_Gamxyz, d_Gamxzz,
        d_Gamyxx, d_Gamyxy, d_Gamyxz,
        d_Gamyyy, d_Gamyyz, d_Gamyzz,
        d_Gamzxx, d_Gamzxy, d_Gamzxz,
        d_Gamzyy, d_Gamzyz, d_Gamzzz,
        d_Rxx, d_Rxy, d_Rxz, d_Ryy, d_Ryz, d_Rzz,
        d_Rpsi4, d_Ipsi4,
        symmetry
    );
}

template<int GW, int TETRAD>
__global__ void getnp4_ss_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* X, const double* Y, const double* Z,
    const double* phyx, const double* phyy, const double* phyz,
    const double* drhodx,  const double* drhody,  const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx,    const double* dRdy,    const double* dRdz,
    const double* chi, const double* trK,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    const double* Axx, const double* Axy, const double* Axz,
    const double* Ayy, const double* Ayz, const double* Azz,
    const double* Gamxxx, const double* Gamxxy, const double* Gamxxz,
    const double* Gamxyy, const double* Gamxyz, const double* Gamxzz,
    const double* Gamyxx, const double* Gamyxy, const double* Gamyxz,
    const double* Gamyyy, const double* Gamyyz, const double* Gamyzz,
    const double* Gamzxx, const double* Gamzxy, const double* Gamzxz,
    const double* Gamzyy, const double* Gamzyz, const double* Gamzzz,
    const double* Rxx, const double* Rxy, const double* Rxz,
    const double* Ryy, const double* Ryz, const double* Rzz,
    double* Rpsi4, double* Ipsi4,
    int symmetry, int sst
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    size_t idx = i + j * pex0 + k * (size_t)pex0 * pex1;

    const int ex[3]  = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    const double ZEO = 0.0;
    const double ONE = 1.0;
    const double TWO = 2.0;
    const double F1o3 = 1.0 / 3.0;
    const double SYM = 1.0;
    const double ANTI = -1.0;
    const double TINYRR = 1e-14;

    double cx = phyx[idx];
    double cy = phyy[idx];
    double cz = phyz[idx];

    // Load physical metric
    double l_gxx = dxx[idx] + ONE;
    double l_gyy = dyy[idx] + ONE;
    double l_gzz = dzz[idx] + ONE;
    double l_gxy = gxy[idx];
    double l_gxz = gxz[idx];
    double l_gyz = gyz[idx];

    double l_chi = chi[idx];
    double chipn1 = l_chi + ONE;
    double f = ONE / chipn1;

    // Inverse metric
    double detG = l_gxx * l_gyy * l_gzz + l_gxy * l_gyz * l_gxz + l_gxz * l_gxy * l_gyz -
                  l_gxz * l_gyy * l_gxz - l_gxy * l_gxy * l_gzz - l_gxx * l_gyz * l_gyz;

    double l_gupxx =  (l_gyy * l_gzz - l_gyz * l_gyz) / detG;
    double l_gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) / detG;
    double l_gupxz =  (l_gxy * l_gyz - l_gyy * l_gxz) / detG;
    double l_gupyy =  (l_gxx * l_gzz - l_gxz * l_gxz) / detG;
    double l_gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) / detG;
    double l_gupzz =  (l_gxx * l_gyy - l_gxy * l_gxy) / detG;

    // Initialize tetrad vectors U, V, W
    double vx, vy, vz, ux, uy, uz, wx, wy, wz;

    if (fabs(cx) < TINYRR && fabs(cy) < TINYRR && fabs(cz) < TINYRR) {
        vx = TINYRR; vy = TINYRR; vz = TINYRR;
    } else {
        vx = cx; vy = cy; vz = cz;
    }

    if (fabs(cx) < TINYRR && fabs(cy) < TINYRR) {
        ux = -TINYRR; uy = TINYRR; uz = ZEO;
        wx = TINYRR * cz; wy = TINYRR * cz; wz = -TWO * TINYRR * TINYRR;
    } else {
        ux = -cy; uy = cx; uz = ZEO;
        wx = cx * cz; wy = cy * cz; wz = -(cx * cx + cy * cy);
    }

    // Gram-Schmidt orthogonalization (order depends on tetradtype)
    if constexpr (TETRAD == 2) {
        double fxt = vx, fyt = vy, fzt = vz;
        vx = l_gupxx * fxt + l_gupxy * fyt + l_gupxz * fzt;
        vy = l_gupxy * fxt + l_gupyy * fyt + l_gupyz * fzt;
        vz = l_gupxz * fxt + l_gupyz * fyt + l_gupzz * fzt;

        double temp;
        temp = l_gxx*vx*vx + l_gyy*vy*vy + l_gzz*vz*vz + (l_gxy*vx*vy + l_gxz*vx*vz + l_gyz*vy*vz)*TWO;
        temp = sqrt(temp*f);
        vx /= temp; vy /= temp; vz /= temp;

        temp = l_gxx*vx*ux + l_gxy*vx*uy + l_gxz*vx*uz +
               l_gxy*vy*ux + l_gyy*vy*uy + l_gyz*vy*uz +
               l_gxz*vz*ux + l_gyz*vz*uy + l_gzz*vz*uz;
        temp *= f;
        ux -= temp*vx; uy -= temp*vy; uz -= temp*vz;

        temp = l_gxx*ux*ux + l_gyy*uy*uy + l_gzz*uz*uz + (l_gxy*ux*uy + l_gxz*ux*uz + l_gyz*uy*uz)*TWO;
        temp = sqrt(temp*f);
        ux /= temp; uy /= temp; uz /= temp;

        temp = l_gxx*vx*wx + l_gxy*vx*wy + l_gxz*vx*wz +
               l_gxy*vy*wx + l_gyy*vy*wy + l_gyz*vy*wz +
               l_gxz*vz*wx + l_gyz*vz*wy + l_gzz*vz*wz;
        temp *= f;
        wx -= temp*vx; wy -= temp*vy; wz -= temp*vz;

        temp = l_gxx*ux*wx + l_gxy*ux*wy + l_gxz*ux*wz +
               l_gxy*uy*wx + l_gyy*uy*wy + l_gyz*uy*wz +
               l_gxz*uz*wx + l_gyz*uz*wy + l_gzz*uz*wz;
        temp *= f;
        wx -= temp*ux; wy -= temp*uy; wz -= temp*uz;

        temp = l_gxx*wx*wx + l_gyy*wy*wy + l_gzz*wz*wz + (l_gxy*wx*wy + l_gxz*wx*wz + l_gyz*wy*wz)*TWO;
        temp = sqrt(temp*f);
        wx /= temp; wy /= temp; wz /= temp;

    } else if constexpr (TETRAD == 0) {
        double temp;
        temp = l_gxx*vx*vx + l_gyy*vy*vy + l_gzz*vz*vz + (l_gxy*vx*vy + l_gxz*vx*vz + l_gyz*vy*vz)*TWO;
        temp = sqrt(temp*f);
        vx /= temp; vy /= temp; vz /= temp;

        temp = l_gxx*vx*ux + l_gxy*vx*uy + l_gxz*vx*uz +
               l_gxy*vy*ux + l_gyy*vy*uy + l_gyz*vy*uz +
               l_gxz*vz*ux + l_gyz*vz*uy + l_gzz*vz*uz;
        temp *= f;
        ux -= temp*vx; uy -= temp*vy; uz -= temp*vz;

        temp = l_gxx*ux*ux + l_gyy*uy*uy + l_gzz*uz*uz + (l_gxy*ux*uy + l_gxz*ux*uz + l_gyz*uy*uz)*TWO;
        temp = sqrt(temp*f);
        ux /= temp; uy /= temp; uz /= temp;

        temp = l_gxx*vx*wx + l_gxy*vx*wy + l_gxz*vx*wz +
               l_gxy*vy*wx + l_gyy*vy*wy + l_gyz*vy*wz +
               l_gxz*vz*wx + l_gyz*vz*wy + l_gzz*vz*wz;
        temp *= f;
        wx -= temp*vx; wy -= temp*vy; wz -= temp*vz;

        temp = l_gxx*ux*wx + l_gxy*ux*wy + l_gxz*ux*wz +
               l_gxy*uy*wx + l_gyy*uy*wy + l_gyz*uy*wz +
               l_gxz*uz*wx + l_gyz*uz*wy + l_gzz*uz*wz;
        temp *= f;
        wx -= temp*ux; wy -= temp*uy; wz -= temp*uz;

        temp = l_gxx*wx*wx + l_gyy*wy*wy + l_gzz*wz*wz + (l_gxy*wx*wy + l_gxz*wx*wz + l_gyz*wy*wz)*TWO;
        temp = sqrt(temp*f);
        wx /= temp; wy /= temp; wz /= temp;

    } else if constexpr (TETRAD == 1) {
        double temp;

        temp = l_gxx*wx*wx + l_gyy*wy*wy + l_gzz*wz*wz + (l_gxy*wx*wy + l_gxz*wx*wz + l_gyz*wy*wz)*TWO;
        temp = sqrt(temp*f);
        wx /= temp; wy /= temp; wz /= temp;

        temp = l_gxx*wx*ux + l_gxy*wx*uy + l_gxz*wx*uz +
               l_gxy*wy*ux + l_gyy*wy*uy + l_gyz*wy*uz +
               l_gxz*wz*ux + l_gyz*wz*uy + l_gzz*wz*uz;
        temp *= f;
        ux -= temp*wx; uy -= temp*wy; uz -= temp*wz;

        temp = l_gxx*ux*ux + l_gyy*uy*uy + l_gzz*uz*uz + (l_gxy*ux*uy + l_gxz*ux*uz + l_gyz*uy*uz)*TWO;
        temp = sqrt(temp*f);
        ux /= temp; uy /= temp; uz /= temp;

        temp = l_gxx*vx*wx + l_gxy*vx*wy + l_gxz*vx*wz +
               l_gxy*vy*wx + l_gyy*vy*wy + l_gyz*vy*wz +
               l_gxz*vz*wx + l_gyz*vz*wy + l_gzz*vz*wz;
        temp *= f;
        vx -= temp*wx; vy -= temp*wy; vz -= temp*wz;

        temp = l_gxx*ux*vx + l_gxy*ux*vy + l_gxz*ux*vz +
               l_gxy*uy*vx + l_gyy*uy*vy + l_gyz*uy*vz +
               l_gxz*uz*vx + l_gyz*uz*vy + l_gzz*uz*vz;
        temp *= f;
        vx -= temp*ux; vy -= temp*uy; vz -= temp*uz;

        temp = l_gxx*vx*vx + l_gyy*vy*vy + l_gzz*vz*vz + (l_gxy*vx*vy + l_gxz*vx*vz + l_gyz*vy*vz)*TWO;
        temp = sqrt(temp*f);
        vx /= temp; vy /= temp; vz /= temp;
    }

    // Compute derivatives of A_ij, chi, trK using Shell coordinate derivatives
    double Axxx, Axxy, Axxz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Axx, &Axxx, &Axxy, &Axxz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double Axyx, Axyy, Axyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Axy, &Axyx, &Axyy, &Axyz, X, Y, Z, ANTI, ANTI, SYM,  symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double Axzx, Axzy, Axzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Axz, &Axzx, &Axzy, &Axzz, X, Y, Z, ANTI, SYM,  ANTI, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double Ayyx, Ayyy, Ayyz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Ayy, &Ayyx, &Ayyy, &Ayyz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double Ayzx, Ayzy, Ayzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Ayz, &Ayzx, &Ayzy, &Ayzz, X, Y, Z, SYM,  ANTI, ANTI, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double Azzx, Azzy, Azzz;
    d_point_fderivs_shc_gw<GW>(ex, pex, Azz, &Azzx, &Azzy, &Azzz, X, Y, Z, SYM,  SYM,  SYM,  symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double chix, chiy, chiz;
    d_point_fderivs_shc_gw<GW>(ex, pex, chi, &chix, &chiy, &chiz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    double fx_trk, fy_trk, fz_trk;
    d_point_fderivs_shc_gw<GW>(ex, pex, trK, &fx_trk, &fy_trk, &fz_trk, X, Y, Z, SYM, SYM, SYM, symmetry, 0, sst, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz, i, j, k);

    // Load traceless extrinsic curvature at current point
    double l_Axx = Axx[idx];
    double l_Axy = Axy[idx];
    double l_Axz = Axz[idx];
    double l_Ayy = Ayy[idx];
    double l_Ayz = Ayz[idx];
    double l_Azz = Azz[idx];

    // Compute D_k K_ij (covariant derivative corrections using Christoffel symbols)
    Axxx -= (Gamxxx[idx]*l_Axx + Gamyxx[idx]*l_Axy + Gamzxx[idx]*l_Axz)*TWO + chix/chipn1*l_Axx - F1o3*l_gxx*fx_trk;
    Axxy -= (Gamxxy[idx]*l_Axx + Gamyxy[idx]*l_Axy + Gamzxy[idx]*l_Axz)*TWO + chiy/chipn1*l_Axx - F1o3*l_gxx*fy_trk;
    Axxz -= (Gamxxz[idx]*l_Axx + Gamyxz[idx]*l_Axy + Gamzxz[idx]*l_Axz)*TWO + chiz/chipn1*l_Axx - F1o3*l_gxx*fz_trk;

    Ayyx -= (Gamxxy[idx]*l_Axy + Gamyxy[idx]*l_Ayy + Gamzxy[idx]*l_Ayz)*TWO + chix/chipn1*l_Ayy - F1o3*l_gyy*fx_trk;
    Ayyy -= (Gamxyy[idx]*l_Axy + Gamyyy[idx]*l_Ayy + Gamzyy[idx]*l_Ayz)*TWO + chiy/chipn1*l_Ayy - F1o3*l_gyy*fy_trk;
    Ayyz -= (Gamxyz[idx]*l_Axy + Gamyyz[idx]*l_Ayy + Gamzyz[idx]*l_Ayz)*TWO + chiz/chipn1*l_Ayy - F1o3*l_gyy*fz_trk;

    Azzx -= (Gamxxz[idx]*l_Axz + Gamyxz[idx]*l_Ayz + Gamzxz[idx]*l_Azz)*TWO + chix/chipn1*l_Azz - F1o3*l_gzz*fx_trk;
    Azzy -= (Gamxyz[idx]*l_Axz + Gamyyz[idx]*l_Ayz + Gamzyz[idx]*l_Azz)*TWO + chiy/chipn1*l_Azz - F1o3*l_gzz*fy_trk;
    Azzz -= (Gamxzz[idx]*l_Axz + Gamyzz[idx]*l_Ayz + Gamzzz[idx]*l_Azz)*TWO + chiz/chipn1*l_Azz - F1o3*l_gzz*fz_trk;

    Axyx -= (Gamxxy[idx]*l_Axx + Gamyxy[idx]*l_Axy + Gamzxy[idx]*l_Axz + Gamxxx[idx]*l_Axy + Gamyxx[idx]*l_Ayy + Gamzxx[idx]*l_Ayz) + chix/chipn1*l_Axy - F1o3*l_gxy*fx_trk;
    Axyy -= (Gamxyy[idx]*l_Axx + Gamyyy[idx]*l_Axy + Gamzyy[idx]*l_Axz + Gamxxy[idx]*l_Axy + Gamyxy[idx]*l_Ayy + Gamzxy[idx]*l_Ayz) + chiy/chipn1*l_Axy - F1o3*l_gxy*fy_trk;
    Axyz -= (Gamxyz[idx]*l_Axx + Gamyyz[idx]*l_Axy + Gamzyz[idx]*l_Axz + Gamxxz[idx]*l_Axy + Gamyxz[idx]*l_Ayy + Gamzxz[idx]*l_Ayz) + chiz/chipn1*l_Axy - F1o3*l_gxy*fz_trk;

    Axzx -= (Gamxxz[idx]*l_Axx + Gamyxz[idx]*l_Axy + Gamzxz[idx]*l_Axz + Gamxxx[idx]*l_Axz + Gamyxx[idx]*l_Ayz + Gamzxx[idx]*l_Azz) + chix/chipn1*l_Axz - F1o3*l_gxz*fx_trk;
    Axzy -= (Gamxyz[idx]*l_Axx + Gamyyz[idx]*l_Axy + Gamzyz[idx]*l_Axz + Gamxxy[idx]*l_Axz + Gamyxy[idx]*l_Ayz + Gamzxy[idx]*l_Azz) + chiy/chipn1*l_Axz - F1o3*l_gxz*fy_trk;
    Axzz -= (Gamxzz[idx]*l_Axx + Gamyzz[idx]*l_Axy + Gamzzz[idx]*l_Axz + Gamxxz[idx]*l_Axz + Gamyxz[idx]*l_Ayz + Gamzxz[idx]*l_Azz) + chiz/chipn1*l_Axz - F1o3*l_gxz*fz_trk;

    Ayzx -= (Gamxxz[idx]*l_Axy + Gamyxz[idx]*l_Ayy + Gamzxz[idx]*l_Ayz + Gamxxy[idx]*l_Axz + Gamyxy[idx]*l_Ayz + Gamzxy[idx]*l_Azz) + chix/chipn1*l_Ayz - F1o3*l_gyz*fx_trk;
    Ayzy -= (Gamxyz[idx]*l_Axy + Gamyyz[idx]*l_Ayy + Gamzyz[idx]*l_Ayz + Gamxyy[idx]*l_Axz + Gamyyy[idx]*l_Ayz + Gamzyy[idx]*l_Azz) + chiy/chipn1*l_Ayz - F1o3*l_gyz*fy_trk;
    Ayzz -= (Gamxzz[idx]*l_Axy + Gamyzz[idx]*l_Ayy + Gamzzz[idx]*l_Ayz + Gamxyz[idx]*l_Axz + Gamyyz[idx]*l_Ayz + Gamzyz[idx]*l_Azz) + chiz/chipn1*l_Ayz - F1o3*l_gyz*fz_trk;

    // Symmetrize B_ij = v^k (D_k K_ij - D_j K_ik)
    double Bxx = (vy*(Axxy - Axyx) + vz*(Axxz - Axzx)) * f;
    double Byy = (vx*(Ayyx - Axyy) + vz*(Ayyz - Ayzy)) * f;
    double Bzz = (vx*(Azzx - Axzz) + vy*(Azzy - Ayzz)) * f;
    double Bxy = (vx*(Axyx - (Axxy+Axyx)/TWO) + vy*(Axyy-Ayyx)/TWO + vz*(Axyz - (Axzy+Ayzx)/TWO)) * f;
    double Bxz = (vx*(Axzx - (Axxz+Axzx)/TWO) + vy*(Axzy - (Axyz+Ayzx)/TWO) + vz*(Axzz-Azzx)/TWO) * f;
    double Byz = (vx*(Ayzx - (Axyz+Axzy)/TWO) + vy*(Ayzy - (Ayyz+Ayzy)/TWO) + vz*(Ayzz-Azzy)/TWO) * f;

    double l_trK = trK[idx];

    // Physical K_ij = A_ij + (1/3)*trK*g_ij
    double Kxx_c = l_Axx + F1o3 * l_trK * l_gxx;
    double Kxy_c = l_Axy + F1o3 * l_trK * l_gxy;
    double Kxz_c = l_Axz + F1o3 * l_trK * l_gxz;
    double Kyy_c = l_Ayy + F1o3 * l_trK * l_gyy;
    double Kyz_c = l_Ayz + F1o3 * l_trK * l_gyz;
    double Kzz_c = l_Azz + F1o3 * l_trK * l_gzz;

    // Compute E_ij = K_ik K^k_j
    double Exx = l_gupxx*Kxx_c*Kxx_c + l_gupyy*Kxy_c*Kxy_c + l_gupzz*Kxz_c*Kxz_c +
                 TWO*(l_gupxy*Kxx_c*Kxy_c + l_gupxz*Kxx_c*Kxz_c + l_gupyz*Kxy_c*Kxz_c);
    double Eyy = l_gupxx*Kxy_c*Kxy_c + l_gupyy*Kyy_c*Kyy_c + l_gupzz*Kyz_c*Kyz_c +
                 TWO*(l_gupxy*Kxy_c*Kyy_c + l_gupxz*Kxy_c*Kyz_c + l_gupyz*Kyy_c*Kyz_c);
    double Ezz = l_gupxx*Kxz_c*Kxz_c + l_gupyy*Kyz_c*Kyz_c + l_gupzz*Kzz_c*Kzz_c +
                 TWO*(l_gupxy*Kxz_c*Kyz_c + l_gupxz*Kxz_c*Kzz_c + l_gupyz*Kyz_c*Kzz_c);

    double Exy = l_gupxx*Kxx_c*Kxy_c + l_gupyy*Kxy_c*Kyy_c + l_gupzz*Kxz_c*Kyz_c +
                 l_gupxy*(Kxx_c*Kyy_c + Kxy_c*Kxy_c) +
                 l_gupxz*(Kxx_c*Kyz_c + Kxz_c*Kxy_c) +
                 l_gupyz*(Kxy_c*Kyz_c + Kxz_c*Kyy_c);
    double Exz = l_gupxx*Kxx_c*Kxz_c + l_gupyy*Kxy_c*Kyz_c + l_gupzz*Kxz_c*Kzz_c +
                 l_gupxy*(Kxx_c*Kyz_c + Kxy_c*Kxz_c) +
                 l_gupxz*(Kxx_c*Kzz_c + Kxz_c*Kxz_c) +
                 l_gupyz*(Kxy_c*Kzz_c + Kxz_c*Kyz_c);
    double Eyz = l_gupxx*Kxy_c*Kxz_c + l_gupyy*Kyy_c*Kyz_c + l_gupzz*Kyz_c*Kzz_c +
                 l_gupxy*(Kxy_c*Kyz_c + Kyy_c*Kxz_c) +
                 l_gupxz*(Kxy_c*Kzz_c + Kyz_c*Kxz_c) +
                 l_gupyz*(Kyy_c*Kzz_c + Kyz_c*Kyz_c);

    Exx = Rxx[idx] - (Exx - Kxx_c*l_trK)*f - Bxx;
    Exy = Rxy[idx] - (Exy - Kxy_c*l_trK)*f - Bxy;
    Exz = Rxz[idx] - (Exz - Kxz_c*l_trK)*f - Bxz;
    Eyy = Ryy[idx] - (Eyy - Kyy_c*l_trK)*f - Byy;
    Eyz = Ryz[idx] - (Eyz - Kyz_c*l_trK)*f - Byz;
    Ezz = Rzz[idx] - (Ezz - Kzz_c*l_trK)*f - Bzz;

    // Projection tensors: uuww^ij = u^i u^j - w^i w^j, uw^ij = u^i w^j + w^i u^j
    double uuwwxx = ux * ux - wx * wx;
    double uuwwxy = ux * uy - wx * wy;
    double uuwwxz = ux * uz - wx * wz;
    double uuwwyy = uy * uy - wy * wy;
    double uuwwyz = uy * uz - wy * wz;
    double uuwwzz = uz * uz - wz * wz;

    double uwxx = ux * wx + wx * ux;
    double uwxy = ux * wy + wx * uy;
    double uwxz = ux * wz + wx * uz;
    double uwyy = uy * wy + wy * uy;
    double uwyz = uy * wz + wy * uz;
    double uwzz = uz * wz + wz * uz;

    // Combine Psi4
    double l_Rpsi4 = Exx*uuwwxx + Eyy*uuwwyy + Ezz*uuwwzz + (Exy*uuwwxy + Exz*uuwwxz + Eyz*uuwwyz)*TWO;
    double l_Ipsi4 = Exx*uwxx + Eyy*uwyy + Ezz*uwzz + (Exy*uwxy + Exz*uwxz + Eyz*uwyz)*TWO;

    Rpsi4[idx] = -l_Rpsi4 / TWO;
    Ipsi4[idx] = -l_Ipsi4 / TWO;
}

// ---------------------------------------------------------
// Template launch function
// ---------------------------------------------------------
template<int GW, int TETRAD>
void gpu_getnp4_ss_launch_impl(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_phyx, const double* d_phyy, const double* d_phyz,
    const double* d_drhodx, const double* d_drhody, const double* d_drhodz,
    const double* d_dsigmadx, const double* d_dsigmady, const double* d_dsigmadz,
    const double* d_dRdx, const double* d_dRdy, const double* d_dRdz,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    const double* d_Gamxxx, const double* d_Gamxxy, const double* d_Gamxxz,
    const double* d_Gamxyy, const double* d_Gamxyz, const double* d_Gamxzz,
    const double* d_Gamyxx, const double* d_Gamyxy, const double* d_Gamyxz,
    const double* d_Gamyyy, const double* d_Gamyyz, const double* d_Gamyzz,
    const double* d_Gamzxx, const double* d_Gamzxy, const double* d_Gamzxz,
    const double* d_Gamzyy, const double* d_Gamzyz, const double* d_Gamzzz,
    const double* d_Rxx, const double* d_Rxy, const double* d_Rxz,
    const double* d_Ryy, const double* d_Ryz, const double* d_Rzz,
    double* d_Rpsi4, double* d_Ipsi4,
    int symmetry, int sst
) {
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    getnp4_ss_kernel<GW, TETRAD><<<gridDim, blockDim, 0, stream>>>(
        ex[0], ex[1], ex[2],
        pex[0], pex[1], pex[2],
        d_X, d_Y, d_Z,
        d_phyx, d_phyy, d_phyz,
        d_drhodx, d_drhody, d_drhodz,
        d_dsigmadx, d_dsigmady, d_dsigmadz,
        d_dRdx, d_dRdy, d_dRdz,
        d_chi, d_trK,
        d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
        d_Axx, d_Axy, d_Axz, d_Ayy, d_Ayz, d_Azz,
        d_Gamxxx, d_Gamxxy, d_Gamxxz,
        d_Gamxyy, d_Gamxyz, d_Gamxzz,
        d_Gamyxx, d_Gamyxy, d_Gamyxz,
        d_Gamyyy, d_Gamyyz, d_Gamyzz,
        d_Gamzxx, d_Gamzxy, d_Gamzxz,
        d_Gamzyy, d_Gamzyz, d_Gamzzz,
        d_Rxx, d_Rxy, d_Rxz, d_Ryy, d_Ryz, d_Rzz,
        d_Rpsi4, d_Ipsi4,
        symmetry, sst
    );
}

// Backward-compatible wrapper: default to GW=ghost_width, TETRAD=tetradtype via macros
void gpu_getnp4_ss_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_phyx, const double* d_phyy, const double* d_phyz,
    const double* d_drhodx, const double* d_drhody, const double* d_drhodz,
    const double* d_dsigmadx, const double* d_dsigmady, const double* d_dsigmadz,
    const double* d_dRdx, const double* d_dRdy, const double* d_dRdz,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    const double* d_Gamxxx, const double* d_Gamxxy, const double* d_Gamxxz,
    const double* d_Gamxyy, const double* d_Gamxyz, const double* d_Gamxzz,
    const double* d_Gamyxx, const double* d_Gamyxy, const double* d_Gamyxz,
    const double* d_Gamyyy, const double* d_Gamyyz, const double* d_Gamyzz,
    const double* d_Gamzxx, const double* d_Gamzxy, const double* d_Gamzxz,
    const double* d_Gamzyy, const double* d_Gamzyz, const double* d_Gamzzz,
    const double* d_Rxx, const double* d_Rxy, const double* d_Rxz,
    const double* d_Ryy, const double* d_Ryz, const double* d_Rzz,
    double* d_Rpsi4, double* d_Ipsi4,
    int symmetry, int sst
) {
    gpu_getnp4_ss_launch_impl<ghost_width, tetradtype>(
        stream, ex, pex,
        d_X, d_Y, d_Z,
        d_phyx, d_phyy, d_phyz,
        d_drhodx, d_drhody, d_drhodz,
        d_dsigmadx, d_dsigmady, d_dsigmadz,
        d_dRdx, d_dRdy, d_dRdz,
        d_chi, d_trK,
        d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
        d_Axx, d_Axy, d_Axz, d_Ayy, d_Ayz, d_Azz,
        d_Gamxxx, d_Gamxxy, d_Gamxxz,
        d_Gamxyy, d_Gamxyz, d_Gamxzz,
        d_Gamyxx, d_Gamyxy, d_Gamyxz,
        d_Gamyyy, d_Gamyyz, d_Gamyzz,
        d_Gamzxx, d_Gamzxy, d_Gamzxz,
        d_Gamzyy, d_Gamzyz, d_Gamzzz,
        d_Rxx, d_Rxy, d_Rxz, d_Ryy, d_Ryz, d_Rzz,
        d_Rpsi4, d_Ipsi4,
        symmetry, sst
    );
}


// ---------------------------------------------------------
// ABV preprocessing: compute ADM Christoffel + Ricci from BSSN.
// Caller must provide 6 temporary GPU arrays for ADM metric.
// Writes Christoffel symbols into Gam arrays and Ricci into R arrays.
// ---------------------------------------------------------
template<int GW>
void gpu_bssn2adm_ricci_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    double* d_Gamxxx, double* d_Gamxxy, double* d_Gamxxz,
    double* d_Gamxyy, double* d_Gamxyz, double* d_Gamxzz,
    double* d_Gamyxx, double* d_Gamyxy, double* d_Gamyxz,
    double* d_Gamyyy, double* d_Gamyyz, double* d_Gamyzz,
    double* d_Gamzxx, double* d_Gamzxy, double* d_Gamzxz,
    double* d_Gamzyy, double* d_Gamzyz, double* d_Gamzzz,
    double* d_Rxx, double* d_Rxy, double* d_Rxz,
    double* d_Ryy, double* d_Ryz, double* d_Rzz,
    double* d_temp_adm_dxx, double* d_temp_adm_gxy, double* d_temp_adm_gxz,
    double* d_temp_adm_dyy, double* d_temp_adm_gyz, double* d_temp_adm_dzz,
    int symmetry
) {
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    // Phase 1: BSSN -> ADM metric conversion
    bssn2adm_kernel<GW><<<gridDim, blockDim, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_chi, d_trK,
        d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
        d_Axx, d_Axy, d_Axz, d_Ayy, d_Ayz, d_Azz,
        d_temp_adm_dxx, d_temp_adm_gxy, d_temp_adm_gxz,
        d_temp_adm_dyy, d_temp_adm_gyz, d_temp_adm_dzz
    );

    // Phase 2: Compute Christoffel + Ricci from ADM metric temp arrays
    adm_ricci_from_adm_kernel<GW><<<gridDim, blockDim, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_X, d_Y, d_Z,
        d_temp_adm_dxx, d_temp_adm_gxy, d_temp_adm_gxz,
        d_temp_adm_dyy, d_temp_adm_gyz, d_temp_adm_dzz,
        d_Gamxxx, d_Gamxxy, d_Gamxxz,
        d_Gamxyy, d_Gamxyz, d_Gamxzz,
        d_Gamyxx, d_Gamyxy, d_Gamyxz,
        d_Gamyyy, d_Gamyyz, d_Gamyzz,
        d_Gamzxx, d_Gamzxy, d_Gamzxz,
        d_Gamzyy, d_Gamzyz, d_Gamzzz,
        d_Rxx, d_Rxy, d_Rxz, d_Ryy, d_Ryz, d_Rzz,
        symmetry
    );
}
// Non-template wrapper for calling from C++ code
void gpu_bssn2adm_ricci_launch_wrapper(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    double* d_Gamxxx, double* d_Gamxxy, double* d_Gamxxz,
    double* d_Gamxyy, double* d_Gamxyz, double* d_Gamxzz,
    double* d_Gamyxx, double* d_Gamyxy, double* d_Gamyxz,
    double* d_Gamyyy, double* d_Gamyyz, double* d_Gamyzz,
    double* d_Gamzxx, double* d_Gamzxy, double* d_Gamzxz,
    double* d_Gamzyy, double* d_Gamzyz, double* d_Gamzzz,
    double* d_Rxx, double* d_Rxy, double* d_Rxz,
    double* d_Ryy, double* d_Ryz, double* d_Rzz,
    double* d_temp_adm_dxx, double* d_temp_adm_gxy, double* d_temp_adm_gxz,
    double* d_temp_adm_dyy, double* d_temp_adm_gyz, double* d_temp_adm_dzz,
    int symmetry
) {
    gpu_bssn2adm_ricci_launch<ghost_width>(
        stream, ex, pex,
        d_X, d_Y, d_Z,
        d_chi, d_trK,
        d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
        d_Axx, d_Axy, d_Axz, d_Ayy, d_Ayz, d_Azz,
        d_Gamxxx, d_Gamxxy, d_Gamxxz,
        d_Gamxyy, d_Gamxyz, d_Gamxzz,
        d_Gamyxx, d_Gamyxy, d_Gamyxz,
        d_Gamyyy, d_Gamyyz, d_Gamyzz,
        d_Gamzxx, d_Gamzxy, d_Gamzxz,
        d_Gamzyy, d_Gamzyz, d_Gamzzz,
        d_Rxx, d_Rxy, d_Rxz, d_Ryy, d_Ryz, d_Rzz,
        d_temp_adm_dxx, d_temp_adm_gxy, d_temp_adm_gxz,
        d_temp_adm_dyy, d_temp_adm_gyz, d_temp_adm_dzz,
        symmetry
    );
}
