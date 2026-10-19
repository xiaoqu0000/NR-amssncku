#ifndef DIFF_NEW_GPU_CUH
#define DIFF_NEW_GPU_CUH

#include "derivatives.h"

#include "fmisc.h"
#include "fmisc_gpu.cuh"

#ifndef ghost_width
#include "macrodef.fh"
#endif
#include <cmath>
#include <iostream>

#if ghost_width == 2
static constexpr int kDiffCenteredBoundaryRadius = 1;
#elif ghost_width == 3
static constexpr int kDiffCenteredBoundaryRadius = 2;
#elif ghost_width == 4
static constexpr int kDiffCenteredBoundaryRadius = 3;
#elif ghost_width == 5
static constexpr int kDiffCenteredBoundaryRadius = 4;
#else
#error Unsupported ghost_width in diff_new_gpu.cuh; expected ghost_width == 2, 3, 4, or 5.
#endif

// ==========================================
// Device Function: 一阶导数 (4th Order)
// ==========================================
__device__ __forceinline__ void d_fderivs_point(
    const int ex[3], const double* f,
    double* fx, double* fy, double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff,
    int i, int j, int k
) {
    const double ONE = 1.0;
    const double TWO = 2.0;
    const double ZEO = 0.0;
    const int NO_SYMM = 0, EQ_SYMM = 1;

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    const int imax = ex[0] - 1;
    const int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;

    *fx = ZEO;
    *fy = ZEO;
    *fz = ZEO;

    // Fortran 循环范围是 1 到 ex-1，对应 CUDA 0 到 ex-2。
    // 如果 i >= imax (即 i >= ex-1)，直接返回，保持输出为 0。
    if (i >= imax || j >= jmax || k >= kmax) return;

    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -kDiffCenteredBoundaryRadius;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -kDiffCenteredBoundaryRadius;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -kDiffCenteredBoundaryRadius;

    double SoA[3] = {SYM1, SYM2, SYM3};

    const double d2dx = ONE / TWO / dX;
    const double d2dy = ONE / TWO / dY;
    const double d2dz = ONE / TWO / dZ;
#if ghost_width == 4 || ghost_width == 5
    const double d60dx = ONE / 60.0 / dX;
    const double d60dy = ONE / 60.0 / dY;
    const double d60dz = ONE / 60.0 / dZ;
    const double d12dx = ONE / 12.0 / dX;
    const double d12dy = ONE / 12.0 / dY;
    const double d12dz = ONE / 12.0 / dZ;
#endif
#if ghost_width == 5
    const double d840dx = ONE / 840.0 / dX;
    const double d840dy = ONE / 840.0 / dY;
    const double d840dz = ONE / 840.0 / dZ;
#endif

    // Helper lambda for symmetry boundary access
    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b(kDiffCenteredBoundaryRadius, ex[0], ex[1], ex[2], f, ii, jj, kk, SoA[0], SoA[1], SoA[2]);
    };

#if ghost_width == 2
    if (i + 1 <= imax && i - 1 >= imin) {
        *fx = d2dx * (-fh(i-1,j,k) + fh(i+1,j,k));
    }
    if (j + 1 <= jmax && j - 1 >= jmin) {
        *fy = d2dy * (-fh(i,j-1,k) + fh(i,j+1,k));
    }
    if (k + 1 <= kmax && k - 1 >= kmin) {
        *fz = d2dz * (-fh(i,j,k-1) + fh(i,j,k+1));
    }
#elif ghost_width == 3
    const double d12dx = ONE / 12.0 / dX;
    const double d12dy = ONE / 12.0 / dY;
    const double d12dz = ONE / 12.0 / dZ;

    if (i + 2 <= imax && i - 2 >= imin &&
        j + 2 <= jmax && j - 2 >= jmin &&
        k + 2 <= kmax && k - 2 >= kmin) {

        *fx = d12dx * (fh(i-2,j,k) - 8.0*fh(i-1,j,k) + 8.0*fh(i+1,j,k) - fh(i+2,j,k));
        *fy = d12dy * (fh(i,j-2,k) - 8.0*fh(i,j-1,k) + 8.0*fh(i,j+1,k) - fh(i,j+2,k));
        *fz = d12dz * (fh(i,j,k-2) - 8.0*fh(i,j,k-1) + 8.0*fh(i,j,k+1) - fh(i,j,k+2));

    }
    else if (i + 1 <= imax && i - 1 >= imin &&
               j + 1 <= jmax && j - 1 >= jmin &&
               k + 1 <= kmax && k - 1 >= kmin) {

        *fx = d2dx * (-fh(i-1,j,k) + fh(i+1,j,k));
        *fy = d2dy * (-fh(i,j-1,k) + fh(i,j+1,k));
        *fz = d2dz * (-fh(i,j,k-1) + fh(i,j,k+1));
    }
#elif ghost_width == 4
    if (i + 3 <= imax && i - 3 >= imin) {
        *fx = d60dx * (-fh(i-3,j,k) + 9.0*fh(i-2,j,k) - 45.0*fh(i-1,j,k) + 45.0*fh(i+1,j,k) - 9.0*fh(i+2,j,k) + fh(i+3,j,k));
    } else if (i + 2 <= imax && i - 2 >= imin) {
        *fx = d12dx * (fh(i-2,j,k) - 8.0*fh(i-1,j,k) + 8.0*fh(i+1,j,k) - fh(i+2,j,k));
    } else if (i + 1 <= imax && i - 1 >= imin) {
        *fx = d2dx * (-fh(i-1,j,k) + fh(i+1,j,k));
    }

    if (j + 3 <= jmax && j - 3 >= jmin) {
        *fy = d60dy * (-fh(i,j-3,k) + 9.0*fh(i,j-2,k) - 45.0*fh(i,j-1,k) + 45.0*fh(i,j+1,k) - 9.0*fh(i,j+2,k) + fh(i,j+3,k));
    } else if (j + 2 <= jmax && j - 2 >= jmin) {
        *fy = d12dy * (fh(i,j-2,k) - 8.0*fh(i,j-1,k) + 8.0*fh(i,j+1,k) - fh(i,j+2,k));
    } else if (j + 1 <= jmax && j - 1 >= jmin) {
        *fy = d2dy * (-fh(i,j-1,k) + fh(i,j+1,k));
    }

    if (k + 3 <= kmax && k - 3 >= kmin) {
        *fz = d60dz * (-fh(i,j,k-3) + 9.0*fh(i,j,k-2) - 45.0*fh(i,j,k-1) + 45.0*fh(i,j,k+1) - 9.0*fh(i,j,k+2) + fh(i,j,k+3));
    } else if (k + 2 <= kmax && k - 2 >= kmin) {
        *fz = d12dz * (fh(i,j,k-2) - 8.0*fh(i,j,k-1) + 8.0*fh(i,j,k+1) - fh(i,j,k+2));
    } else if (k + 1 <= kmax && k - 1 >= kmin) {
        *fz = d2dz * (-fh(i,j,k-1) + fh(i,j,k+1));
    }
#elif ghost_width == 5
    if (i + 4 <= imax && i - 4 >= imin) {
        *fx = d840dx * (3.0*fh(i-4,j,k) - 32.0*fh(i-3,j,k) + 168.0*fh(i-2,j,k) - 672.0*fh(i-1,j,k) + 672.0*fh(i+1,j,k) - 168.0*fh(i+2,j,k) + 32.0*fh(i+3,j,k) - 3.0*fh(i+4,j,k));
    } else if (i + 3 <= imax && i - 3 >= imin) {
        *fx = d60dx * (-fh(i-3,j,k) + 9.0*fh(i-2,j,k) - 45.0*fh(i-1,j,k) + 45.0*fh(i+1,j,k) - 9.0*fh(i+2,j,k) + fh(i+3,j,k));
    } else if (i + 2 <= imax && i - 2 >= imin) {
        *fx = d12dx * (fh(i-2,j,k) - 8.0*fh(i-1,j,k) + 8.0*fh(i+1,j,k) - fh(i+2,j,k));
    } else if (i + 1 <= imax && i - 1 >= imin) {
        *fx = d2dx * (-fh(i-1,j,k) + fh(i+1,j,k));
    }

    if (j + 4 <= jmax && j - 4 >= jmin) {
        *fy = d840dy * (3.0*fh(i,j-4,k) - 32.0*fh(i,j-3,k) + 168.0*fh(i,j-2,k) - 672.0*fh(i,j-1,k) + 672.0*fh(i,j+1,k) - 168.0*fh(i,j+2,k) + 32.0*fh(i,j+3,k) - 3.0*fh(i,j+4,k));
    } else if (j + 3 <= jmax && j - 3 >= jmin) {
        *fy = d60dy * (-fh(i,j-3,k) + 9.0*fh(i,j-2,k) - 45.0*fh(i,j-1,k) + 45.0*fh(i,j+1,k) - 9.0*fh(i,j+2,k) + fh(i,j+3,k));
    } else if (j + 2 <= jmax && j - 2 >= jmin) {
        *fy = d12dy * (fh(i,j-2,k) - 8.0*fh(i,j-1,k) + 8.0*fh(i,j+1,k) - fh(i,j+2,k));
    } else if (j + 1 <= jmax && j - 1 >= jmin) {
        *fy = d2dy * (-fh(i,j-1,k) + fh(i,j+1,k));
    }

    if (k + 4 <= kmax && k - 4 >= kmin) {
        *fz = d840dz * (3.0*fh(i,j,k-4) - 32.0*fh(i,j,k-3) + 168.0*fh(i,j,k-2) - 672.0*fh(i,j,k-1) + 672.0*fh(i,j,k+1) - 168.0*fh(i,j,k+2) + 32.0*fh(i,j,k+3) - 3.0*fh(i,j,k+4));
    } else if (k + 3 <= kmax && k - 3 >= kmin) {
        *fz = d60dz * (-fh(i,j,k-3) + 9.0*fh(i,j,k-2) - 45.0*fh(i,j,k-1) + 45.0*fh(i,j,k+1) - 9.0*fh(i,j,k+2) + fh(i,j,k+3));
    } else if (k + 2 <= kmax && k - 2 >= kmin) {
        *fz = d12dz * (fh(i,j,k-2) - 8.0*fh(i,j,k-1) + 8.0*fh(i,j,k+1) - fh(i,j,k+2));
    } else if (k + 1 <= kmax && k - 1 >= kmin) {
        *fz = d2dz * (-fh(i,j,k-1) + fh(i,j,k+1));
    }
#endif

    (void)onoff;
}

// ==========================================
// Device Function: 二阶导数 (4th Order)
// ==========================================
__device__ __forceinline__ void d_fdderivs_point(
    const int ex[3], const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff,
    int i, int j, int k
) {
    const double ONE = 1.0;
    const double TWO = 2.0;
    const double F1o4 = 0.25;
    const double ZEO = 0.0;
    const int NO_SYMM = 0, EQ_SYMM = 1;

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    const int imax = ex[0] - 1;
    const int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;

    *fxx = ZEO; *fyy = ZEO; *fzz = ZEO;
    *fxy = ZEO; *fxz = ZEO; *fyz = ZEO;

    if (i >= imax || j >= jmax || k >= kmax) return;

    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -kDiffCenteredBoundaryRadius;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -kDiffCenteredBoundaryRadius;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -kDiffCenteredBoundaryRadius;

    double SoA[3] = {SYM1, SYM2, SYM3};

    const double Sdxdx = ONE / (dX * dX);
    const double Sdydy = ONE / (dY * dY);
    const double Sdzdz = ONE / (dZ * dZ);

    const double Sdxdy = F1o4 / (dX * dY);
    const double Sdxdz = F1o4 / (dX * dZ);
    const double Sdydz = F1o4 / (dY * dZ);

#if ghost_width == 3 || ghost_width == 4 || ghost_width == 5
    const double F1o12 = ONE / 12.0;
    const double F1o144 = ONE / 144.0;
    const double F8 = 8.0;
    const double F16 = 16.0;
    const double F30 = 30.0;
    const double Fdxdx = F1o12 / (dX * dX);
    const double Fdydy = F1o12 / (dY * dY);
    const double Fdzdz = F1o12 / (dZ * dZ);
    const double Fdxdy = F1o144 / (dX * dY);
    const double Fdxdz = F1o144 / (dX * dZ);
    const double Fdydz = F1o144 / (dY * dZ);
#endif
#if ghost_width == 4 || ghost_width == 5
    const double Xdxdx = (ONE / 180.0) / (dX * dX);
    const double Xdydy = (ONE / 180.0) / (dY * dY);
    const double Xdzdz = (ONE / 180.0) / (dZ * dZ);
    const double Xdxdy = (ONE / 3600.0) / (dX * dY);
    const double Xdxdz = (ONE / 3600.0) / (dX * dZ);
    const double Xdydz = (ONE / 3600.0) / (dY * dZ);
#endif
#if ghost_width == 5
    const double Edxdx = (ONE / 5040.0) / (dX * dX);
    const double Edydy = (ONE / 5040.0) / (dY * dY);
    const double Edzdz = (ONE / 5040.0) / (dZ * dZ);
    const double Edxdy = (ONE / 705600.0) / (dX * dY);
    const double Edxdz = (ONE / 705600.0) / (dX * dZ);
    const double Edydz = (ONE / 705600.0) / (dY * dZ);
#endif

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b(kDiffCenteredBoundaryRadius, ex[0], ex[1], ex[2], f, ii, jj, kk, SoA[0], SoA[1], SoA[2]);
    };

#if ghost_width == 2
    if (i + 1 <= imax && i - 1 >= imin) {
        *fxx = Sdxdx * (fh(i-1,j,k) - TWO*fh(i,j,k) + fh(i+1,j,k));
    }
    if (j + 1 <= jmax && j - 1 >= jmin) {
        *fyy = Sdydy * (fh(i,j-1,k) - TWO*fh(i,j,k) + fh(i,j+1,k));
    }
    if (k + 1 <= kmax && k - 1 >= kmin) {
        *fzz = Sdzdz * (fh(i,j,k-1) - TWO*fh(i,j,k) + fh(i,j,k+1));
    }

    if (i + 1 <= imax && i - 1 >= imin &&
        j + 1 <= jmax && j - 1 >= jmin) {
        *fxy = Sdxdy * (fh(i-1,j-1,k) - fh(i+1,j-1,k) - fh(i-1,j+1,k) + fh(i+1,j+1,k));
    }
    if (i + 1 <= imax && i - 1 >= imin &&
        k + 1 <= kmax && k - 1 >= kmin) {
        *fxz = Sdxdz * (fh(i-1,j,k-1) - fh(i+1,j,k-1) - fh(i-1,j,k+1) + fh(i+1,j,k+1));
    }
    if (j + 1 <= jmax && j - 1 >= jmin &&
        k + 1 <= kmax && k - 1 >= kmin) {
        *fyz = Sdydz * (fh(i,j-1,k-1) - fh(i,j+1,k-1) - fh(i,j-1,k+1) + fh(i,j+1,k+1));
    }
#elif ghost_width == 3
    // --- 4th Order Accuracy ---
    if (i + 2 <= imax && i - 2 >= imin &&
        j + 2 <= jmax && j - 2 >= jmin &&
        k + 2 <= kmax && k - 2 >= kmin) {

        *fxx = Fdxdx * (-fh(i-2,j,k) + F16*fh(i-1,j,k) - F30*fh(i,j,k)
                        -fh(i+2,j,k) + F16*fh(i+1,j,k));
        *fyy = Fdydy * (-fh(i,j-2,k) + F16*fh(i,j-1,k) - F30*fh(i,j,k)
                        -fh(i,j+2,k) + F16*fh(i,j+1,k));
        *fzz = Fdzdz * (-fh(i,j,k-2) + F16*fh(i,j,k-1) - F30*fh(i,j,k)
                        -fh(i,j,k+2) + F16*fh(i,j,k+1));

        *fxy = Fdxdy * (    (fh(i-2,j-2,k) - F8*fh(i-1,j-2,k) + F8*fh(i+1,j-2,k) - fh(i+2,j-2,k))
                        -F8*(fh(i-2,j-1,k) - F8*fh(i-1,j-1,k) + F8*fh(i+1,j-1,k) - fh(i+2,j-1,k))
                        +F8*(fh(i-2,j+1,k) - F8*fh(i-1,j+1,k) + F8*fh(i+1,j+1,k) - fh(i+2,j+1,k))
                        -   (fh(i-2,j+2,k) - F8*fh(i-1,j+2,k) + F8*fh(i+1,j+2,k) - fh(i+2,j+2,k)) );
        *fxz = Fdxdz * (    (fh(i-2,j,k-2) - F8*fh(i-1,j,k-2) + F8*fh(i+1,j,k-2) - fh(i+2,j,k-2))
                        -F8*(fh(i-2,j,k-1) - F8*fh(i-1,j,k-1) + F8*fh(i+1,j,k-1) - fh(i+2,j,k-1))
                        +F8*(fh(i-2,j,k+1) - F8*fh(i-1,j,k+1) + F8*fh(i+1,j,k+1) - fh(i+2,j,k+1))
                        -   (fh(i-2,j,k+2) - F8*fh(i-1,j,k+2) + F8*fh(i+1,j,k+2) - fh(i+2,j,k+2)) );
        *fyz = Fdydz * (    (fh(i,j-2,k-2) - F8*fh(i,j-1,k-2) + F8*fh(i,j+1,k-2) - fh(i,j+2,k-2))
                        -F8*(fh(i,j-2,k-1) - F8*fh(i,j-1,k-1) + F8*fh(i,j+1,k-1) - fh(i,j+2,k-1))
                        +F8*(fh(i,j-2,k+1) - F8*fh(i,j-1,k+1) + F8*fh(i,j+1,k+1) - fh(i,j+2,k+1))
                        -   (fh(i,j-2,k+2) - F8*fh(i,j-1,k+2) + F8*fh(i,j+1,k+2) - fh(i,j+2,k+2)) );

    } 
    // --- 2nd Order Accuracy ---
    else if (i + 1 <= imax && i - 1 >= imin &&
               j + 1 <= jmax && j - 1 >= jmin &&
               k + 1 <= kmax && k - 1 >= kmin) {

        *fxx = Sdxdx * (fh(i-1,j,k) - TWO*fh(i,j,k) + fh(i+1,j,k));
        *fyy = Sdydy * (fh(i,j-1,k) - TWO*fh(i,j,k) + fh(i,j+1,k));
        *fzz = Sdzdz * (fh(i,j,k-1) - TWO*fh(i,j,k) + fh(i,j,k+1));

        *fxy = Sdxdy * (fh(i-1,j-1,k) - fh(i+1,j-1,k) - fh(i-1,j+1,k) + fh(i+1,j+1,k));
        *fxz = Sdxdz * (fh(i-1,j,k-1) - fh(i+1,j,k-1) - fh(i-1,j,k+1) + fh(i+1,j,k+1));
        *fyz = Sdydz * (fh(i,j-1,k-1) - fh(i,j+1,k-1) - fh(i,j-1,k+1) + fh(i,j+1,k+1));
    }
#elif ghost_width == 4
    if (i + 3 <= imax && i - 3 >= imin) {
        *fxx = Xdxdx * (2.0*fh(i-3,j,k) - 27.0*fh(i-2,j,k) + 270.0*fh(i-1,j,k) - 490.0*fh(i,j,k) + 270.0*fh(i+1,j,k) - 27.0*fh(i+2,j,k) + 2.0*fh(i+3,j,k));
    } else if (i + 2 <= imax && i - 2 >= imin) {
        *fxx = Fdxdx * (-fh(i-2,j,k) + F16*fh(i-1,j,k) - F30*fh(i,j,k) - fh(i+2,j,k) + F16*fh(i+1,j,k));
    } else if (i + 1 <= imax && i - 1 >= imin) {
        *fxx = Sdxdx * (fh(i-1,j,k) - TWO*fh(i,j,k) + fh(i+1,j,k));
    }

    if (j + 3 <= jmax && j - 3 >= jmin) {
        *fyy = Xdydy * (2.0*fh(i,j-3,k) - 27.0*fh(i,j-2,k) + 270.0*fh(i,j-1,k) - 490.0*fh(i,j,k) + 270.0*fh(i,j+1,k) - 27.0*fh(i,j+2,k) + 2.0*fh(i,j+3,k));
    } else if (j + 2 <= jmax && j - 2 >= jmin) {
        *fyy = Fdydy * (-fh(i,j-2,k) + F16*fh(i,j-1,k) - F30*fh(i,j,k) - fh(i,j+2,k) + F16*fh(i,j+1,k));
    } else if (j + 1 <= jmax && j - 1 >= jmin) {
        *fyy = Sdydy * (fh(i,j-1,k) - TWO*fh(i,j,k) + fh(i,j+1,k));
    }

    if (k + 3 <= kmax && k - 3 >= kmin) {
        *fzz = Xdzdz * (2.0*fh(i,j,k-3) - 27.0*fh(i,j,k-2) + 270.0*fh(i,j,k-1) - 490.0*fh(i,j,k) + 270.0*fh(i,j,k+1) - 27.0*fh(i,j,k+2) + 2.0*fh(i,j,k+3));
    } else if (k + 2 <= kmax && k - 2 >= kmin) {
        *fzz = Fdzdz * (-fh(i,j,k-2) + F16*fh(i,j,k-1) - F30*fh(i,j,k) - fh(i,j,k+2) + F16*fh(i,j,k+1));
    } else if (k + 1 <= kmax && k - 1 >= kmin) {
        *fzz = Sdzdz * (fh(i,j,k-1) - TWO*fh(i,j,k) + fh(i,j,k+1));
    }

    if (i + 3 <= imax && i - 3 >= imin && j + 3 <= jmax && j - 3 >= jmin) {
        *fxy = Xdxdy * (
             (fh(i-3,j-3,k) - 9.0*fh(i-3,j-2,k) + 45.0*fh(i-3,j-1,k) - 45.0*fh(i-3,j+1,k) + 9.0*fh(i-3,j+2,k) - fh(i-3,j+3,k))
          -9.0*(fh(i-2,j-3,k) - 9.0*fh(i-2,j-2,k) + 45.0*fh(i-2,j-1,k) - 45.0*fh(i-2,j+1,k) + 9.0*fh(i-2,j+2,k) - fh(i-2,j+3,k))
         +45.0*(fh(i-1,j-3,k) - 9.0*fh(i-1,j-2,k) + 45.0*fh(i-1,j-1,k) - 45.0*fh(i-1,j+1,k) + 9.0*fh(i-1,j+2,k) - fh(i-1,j+3,k))
         -45.0*(fh(i+1,j-3,k) - 9.0*fh(i+1,j-2,k) + 45.0*fh(i+1,j-1,k) - 45.0*fh(i+1,j+1,k) + 9.0*fh(i+1,j+2,k) - fh(i+1,j+3,k))
          +9.0*(fh(i+2,j-3,k) - 9.0*fh(i+2,j-2,k) + 45.0*fh(i+2,j-1,k) - 45.0*fh(i+2,j+1,k) + 9.0*fh(i+2,j+2,k) - fh(i+2,j+3,k))
             -(fh(i+3,j-3,k) - 9.0*fh(i+3,j-2,k) + 45.0*fh(i+3,j-1,k) - 45.0*fh(i+3,j+1,k) + 9.0*fh(i+3,j+2,k) - fh(i+3,j+3,k)));
    } else if (i + 2 <= imax && i - 2 >= imin && j + 2 <= jmax && j - 2 >= jmin) {
        *fxy = Fdxdy * ((fh(i-2,j-2,k) - F8*fh(i-1,j-2,k) + F8*fh(i+1,j-2,k) - fh(i+2,j-2,k))
                 -F8*(fh(i-2,j-1,k) - F8*fh(i-1,j-1,k) + F8*fh(i+1,j-1,k) - fh(i+2,j-1,k))
                 +F8*(fh(i-2,j+1,k) - F8*fh(i-1,j+1,k) + F8*fh(i+1,j+1,k) - fh(i+2,j+1,k))
                     -(fh(i-2,j+2,k) - F8*fh(i-1,j+2,k) + F8*fh(i+1,j+2,k) - fh(i+2,j+2,k)));
    } else if (i + 1 <= imax && i - 1 >= imin && j + 1 <= jmax && j - 1 >= jmin) {
        *fxy = Sdxdy * (fh(i-1,j-1,k) - fh(i+1,j-1,k) - fh(i-1,j+1,k) + fh(i+1,j+1,k));
    }

    if (i + 3 <= imax && i - 3 >= imin && k + 3 <= kmax && k - 3 >= kmin) {
        *fxz = Xdxdz * (
             (fh(i-3,j,k-3) - 9.0*fh(i-3,j,k-2) + 45.0*fh(i-3,j,k-1) - 45.0*fh(i-3,j,k+1) + 9.0*fh(i-3,j,k+2) - fh(i-3,j,k+3))
          -9.0*(fh(i-2,j,k-3) - 9.0*fh(i-2,j,k-2) + 45.0*fh(i-2,j,k-1) - 45.0*fh(i-2,j,k+1) + 9.0*fh(i-2,j,k+2) - fh(i-2,j,k+3))
         +45.0*(fh(i-1,j,k-3) - 9.0*fh(i-1,j,k-2) + 45.0*fh(i-1,j,k-1) - 45.0*fh(i-1,j,k+1) + 9.0*fh(i-1,j,k+2) - fh(i-1,j,k+3))
         -45.0*(fh(i+1,j,k-3) - 9.0*fh(i+1,j,k-2) + 45.0*fh(i+1,j,k-1) - 45.0*fh(i+1,j,k+1) + 9.0*fh(i+1,j,k+2) - fh(i+1,j,k+3))
          +9.0*(fh(i+2,j,k-3) - 9.0*fh(i+2,j,k-2) + 45.0*fh(i+2,j,k-1) - 45.0*fh(i+2,j,k+1) + 9.0*fh(i+2,j,k+2) - fh(i+2,j,k+3))
             -(fh(i+3,j,k-3) - 9.0*fh(i+3,j,k-2) + 45.0*fh(i+3,j,k-1) - 45.0*fh(i+3,j,k+1) + 9.0*fh(i+3,j,k+2) - fh(i+3,j,k+3)));
    } else if (i + 2 <= imax && i - 2 >= imin && k + 2 <= kmax && k - 2 >= kmin) {
        *fxz = Fdxdz * ((fh(i-2,j,k-2) - F8*fh(i-1,j,k-2) + F8*fh(i+1,j,k-2) - fh(i+2,j,k-2))
                 -F8*(fh(i-2,j,k-1) - F8*fh(i-1,j,k-1) + F8*fh(i+1,j,k-1) - fh(i+2,j,k-1))
                 +F8*(fh(i-2,j,k+1) - F8*fh(i-1,j,k+1) + F8*fh(i+1,j,k+1) - fh(i+2,j,k+1))
                     -(fh(i-2,j,k+2) - F8*fh(i-1,j,k+2) + F8*fh(i+1,j,k+2) - fh(i+2,j,k+2)));
    } else if (i + 1 <= imax && i - 1 >= imin && k + 1 <= kmax && k - 1 >= kmin) {
        *fxz = Sdxdz * (fh(i-1,j,k-1) - fh(i+1,j,k-1) - fh(i-1,j,k+1) + fh(i+1,j,k+1));
    }

    if (j + 3 <= jmax && j - 3 >= jmin && k + 3 <= kmax && k - 3 >= kmin) {
        *fyz = Xdydz * (
             (fh(i,j-3,k-3) - 9.0*fh(i,j-3,k-2) + 45.0*fh(i,j-3,k-1) - 45.0*fh(i,j-3,k+1) + 9.0*fh(i,j-3,k+2) - fh(i,j-3,k+3))
          -9.0*(fh(i,j-2,k-3) - 9.0*fh(i,j-2,k-2) + 45.0*fh(i,j-2,k-1) - 45.0*fh(i,j-2,k+1) + 9.0*fh(i,j-2,k+2) - fh(i,j-2,k+3))
         +45.0*(fh(i,j-1,k-3) - 9.0*fh(i,j-1,k-2) + 45.0*fh(i,j-1,k-1) - 45.0*fh(i,j-1,k+1) + 9.0*fh(i,j-1,k+2) - fh(i,j-1,k+3))
         -45.0*(fh(i,j+1,k-3) - 9.0*fh(i,j+1,k-2) + 45.0*fh(i,j+1,k-1) - 45.0*fh(i,j+1,k+1) + 9.0*fh(i,j+1,k+2) - fh(i,j+1,k+3))
          +9.0*(fh(i,j+2,k-3) - 9.0*fh(i,j+2,k-2) + 45.0*fh(i,j+2,k-1) - 45.0*fh(i,j+2,k+1) + 9.0*fh(i,j+2,k+2) - fh(i,j+2,k+3))
             -(fh(i,j+3,k-3) - 9.0*fh(i,j+3,k-2) + 45.0*fh(i,j+3,k-1) - 45.0*fh(i,j+3,k+1) + 9.0*fh(i,j+3,k+2) - fh(i,j+3,k+3)));
    } else if (j + 2 <= jmax && j - 2 >= jmin && k + 2 <= kmax && k - 2 >= kmin) {
        *fyz = Fdydz * ((fh(i,j-2,k-2) - F8*fh(i,j-1,k-2) + F8*fh(i,j+1,k-2) - fh(i,j+2,k-2))
                 -F8*(fh(i,j-2,k-1) - F8*fh(i,j-1,k-1) + F8*fh(i,j+1,k-1) - fh(i,j+2,k-1))
                 +F8*(fh(i,j-2,k+1) - F8*fh(i,j-1,k+1) + F8*fh(i,j+1,k+1) - fh(i,j+2,k+1))
                     -(fh(i,j-2,k+2) - F8*fh(i,j-1,k+2) + F8*fh(i,j+1,k+2) - fh(i,j+2,k+2)));
    } else if (j + 1 <= jmax && j - 1 >= jmin && k + 1 <= kmax && k - 1 >= kmin) {
        *fyz = Sdydz * (fh(i,j-1,k-1) - fh(i,j+1,k-1) - fh(i,j-1,k+1) + fh(i,j+1,k+1));
    }
#elif ghost_width == 5
    const bool fit_x_4 = (i + 4 <= imax && i - 4 >= imin);
    const bool fit_y_4 = (j + 4 <= jmax && j - 4 >= jmin);
    const bool fit_z_4 = (k + 4 <= kmax && k - 4 >= kmin);
    const bool fit_x_3 = (i + 3 <= imax && i - 3 >= imin);
    const bool fit_y_3 = (j + 3 <= jmax && j - 3 >= jmin);
    const bool fit_z_3 = (k + 3 <= kmax && k - 3 >= kmin);
    const bool fit_x_2 = (i + 2 <= imax && i - 2 >= imin);
    const bool fit_y_2 = (j + 2 <= jmax && j - 2 >= jmin);
    const bool fit_z_2 = (k + 2 <= kmax && k - 2 >= kmin);
    const bool fit_x_1 = (i + 1 <= imax && i - 1 >= imin);
    const bool fit_y_1 = (j + 1 <= jmax && j - 1 >= jmin);
    const bool fit_z_1 = (k + 1 <= kmax && k - 1 >= kmin);

    auto dx8 = [&](int jj, int kk) {
        return 3.0*fh(i-4,jj,kk) - 32.0*fh(i-3,jj,kk) + 168.0*fh(i-2,jj,kk) - 672.0*fh(i-1,jj,kk)
             + 672.0*fh(i+1,jj,kk) - 168.0*fh(i+2,jj,kk) + 32.0*fh(i+3,jj,kk) - 3.0*fh(i+4,jj,kk);
    };
    auto dy8 = [&](int ii, int kk) {
        return 3.0*fh(ii,j-4,kk) - 32.0*fh(ii,j-3,kk) + 168.0*fh(ii,j-2,kk) - 672.0*fh(ii,j-1,kk)
             + 672.0*fh(ii,j+1,kk) - 168.0*fh(ii,j+2,kk) + 32.0*fh(ii,j+3,kk) - 3.0*fh(ii,j+4,kk);
    };
    auto dz8 = [&](int ii, int jj) {
        return 3.0*fh(ii,jj,k-4) - 32.0*fh(ii,jj,k-3) + 168.0*fh(ii,jj,k-2) - 672.0*fh(ii,jj,k-1)
             + 672.0*fh(ii,jj,k+1) - 168.0*fh(ii,jj,k+2) + 32.0*fh(ii,jj,k+3) - 3.0*fh(ii,jj,k+4);
    };
    auto dx6 = [&](int jj, int kk) { return -fh(i-3,jj,kk) + 9.0*fh(i-2,jj,kk) - 45.0*fh(i-1,jj,kk) + 45.0*fh(i+1,jj,kk) - 9.0*fh(i+2,jj,kk) + fh(i+3,jj,kk); };
    auto dy6 = [&](int ii, int kk) { return -fh(ii,j-3,kk) + 9.0*fh(ii,j-2,kk) - 45.0*fh(ii,j-1,kk) + 45.0*fh(ii,j+1,kk) - 9.0*fh(ii,j+2,kk) + fh(ii,j+3,kk); };
    auto dz6 = [&](int ii, int jj) { return -fh(ii,jj,k-3) + 9.0*fh(ii,jj,k-2) - 45.0*fh(ii,jj,k-1) + 45.0*fh(ii,jj,k+1) - 9.0*fh(ii,jj,k+2) + fh(ii,jj,k+3); };
    auto dx4 = [&](int jj, int kk) { return fh(i-2,jj,kk) - 8.0*fh(i-1,jj,kk) + 8.0*fh(i+1,jj,kk) - fh(i+2,jj,kk); };
    auto dy4 = [&](int ii, int kk) { return fh(ii,j-2,kk) - 8.0*fh(ii,j-1,kk) + 8.0*fh(ii,j+1,kk) - fh(ii,j+2,kk); };
    auto dz4 = [&](int ii, int jj) { return fh(ii,jj,k-2) - 8.0*fh(ii,jj,k-1) + 8.0*fh(ii,jj,k+1) - fh(ii,jj,k+2); };

    if (fit_x_4) *fxx = Edxdx * (-9.0*fh(i-4,j,k) + 128.0*fh(i-3,j,k) - 1008.0*fh(i-2,j,k) + 8064.0*fh(i-1,j,k) - 14350.0*fh(i,j,k) + 8064.0*fh(i+1,j,k) - 1008.0*fh(i+2,j,k) + 128.0*fh(i+3,j,k) - 9.0*fh(i+4,j,k));
    else if (fit_x_3) *fxx = Xdxdx * (2.0*fh(i-3,j,k) - 27.0*fh(i-2,j,k) + 270.0*fh(i-1,j,k) - 490.0*fh(i,j,k) + 270.0*fh(i+1,j,k) - 27.0*fh(i+2,j,k) + 2.0*fh(i+3,j,k));
    else if (fit_x_2) *fxx = Fdxdx * (-fh(i-2,j,k) + F16*fh(i-1,j,k) - F30*fh(i,j,k) - fh(i+2,j,k) + F16*fh(i+1,j,k));
    else if (fit_x_1) *fxx = Sdxdx * (fh(i-1,j,k) - TWO*fh(i,j,k) + fh(i+1,j,k));

    if (fit_y_4) *fyy = Edydy * (-9.0*fh(i,j-4,k) + 128.0*fh(i,j-3,k) - 1008.0*fh(i,j-2,k) + 8064.0*fh(i,j-1,k) - 14350.0*fh(i,j,k) + 8064.0*fh(i,j+1,k) - 1008.0*fh(i,j+2,k) + 128.0*fh(i,j+3,k) - 9.0*fh(i,j+4,k));
    else if (fit_y_3) *fyy = Xdydy * (2.0*fh(i,j-3,k) - 27.0*fh(i,j-2,k) + 270.0*fh(i,j-1,k) - 490.0*fh(i,j,k) + 270.0*fh(i,j+1,k) - 27.0*fh(i,j+2,k) + 2.0*fh(i,j+3,k));
    else if (fit_y_2) *fyy = Fdydy * (-fh(i,j-2,k) + F16*fh(i,j-1,k) - F30*fh(i,j,k) - fh(i,j+2,k) + F16*fh(i,j+1,k));
    else if (fit_y_1) *fyy = Sdydy * (fh(i,j-1,k) - TWO*fh(i,j,k) + fh(i,j+1,k));

    if (fit_z_4) *fzz = Edzdz * (-9.0*fh(i,j,k-4) + 128.0*fh(i,j,k-3) - 1008.0*fh(i,j,k-2) + 8064.0*fh(i,j,k-1) - 14350.0*fh(i,j,k) + 8064.0*fh(i,j,k+1) - 1008.0*fh(i,j,k+2) + 128.0*fh(i,j,k+3) - 9.0*fh(i,j,k+4));
    else if (fit_z_3) *fzz = Xdzdz * (2.0*fh(i,j,k-3) - 27.0*fh(i,j,k-2) + 270.0*fh(i,j,k-1) - 490.0*fh(i,j,k) + 270.0*fh(i,j,k+1) - 27.0*fh(i,j,k+2) + 2.0*fh(i,j,k+3));
    else if (fit_z_2) *fzz = Fdzdz * (-fh(i,j,k-2) + F16*fh(i,j,k-1) - F30*fh(i,j,k) - fh(i,j,k+2) + F16*fh(i,j,k+1));
    else if (fit_z_1) *fzz = Sdzdz * (fh(i,j,k-1) - TWO*fh(i,j,k) + fh(i,j,k+1));

    if (fit_x_4 && fit_y_4) *fxy = Edxdy * (3.0*dx8(j-4,k) - 32.0*dx8(j-3,k) + 168.0*dx8(j-2,k) - 672.0*dx8(j-1,k) + 672.0*dx8(j+1,k) - 168.0*dx8(j+2,k) + 32.0*dx8(j+3,k) - 3.0*dx8(j+4,k));
    else if (fit_x_3 && fit_y_3) *fxy = Xdxdy * (-dx6(j-3,k) + 9.0*dx6(j-2,k) - 45.0*dx6(j-1,k) + 45.0*dx6(j+1,k) - 9.0*dx6(j+2,k) + dx6(j+3,k));
    else if (fit_x_2 && fit_y_2) *fxy = Fdxdy * (dy4(i-2,k) - F8*dy4(i-1,k) + F8*dy4(i+1,k) - dy4(i+2,k));
    else if (fit_x_1 && fit_y_1) *fxy = Sdxdy * (fh(i-1,j-1,k) - fh(i+1,j-1,k) - fh(i-1,j+1,k) + fh(i+1,j+1,k));

    if (fit_x_4 && fit_z_4) *fxz = Edxdz * (3.0*dx8(j,k-4) - 32.0*dx8(j,k-3) + 168.0*dx8(j,k-2) - 672.0*dx8(j,k-1) + 672.0*dx8(j,k+1) - 168.0*dx8(j,k+2) + 32.0*dx8(j,k+3) - 3.0*dx8(j,k+4));
    else if (fit_x_3 && fit_z_3) *fxz = Xdxdz * (-dx6(j,k-3) + 9.0*dx6(j,k-2) - 45.0*dx6(j,k-1) + 45.0*dx6(j,k+1) - 9.0*dx6(j,k+2) + dx6(j,k+3));
    else if (fit_x_2 && fit_z_2) *fxz = Fdxdz * (dz4(i-2,j) - F8*dz4(i-1,j) + F8*dz4(i+1,j) - dz4(i+2,j));
    else if (fit_x_1 && fit_z_1) *fxz = Sdxdz * (fh(i-1,j,k-1) - fh(i+1,j,k-1) - fh(i-1,j,k+1) + fh(i+1,j,k+1));

    if (fit_y_4 && fit_z_4) *fyz = Edydz * (3.0*dy8(i,k-4) - 32.0*dy8(i,k-3) + 168.0*dy8(i,k-2) - 672.0*dy8(i,k-1) + 672.0*dy8(i,k+1) - 168.0*dy8(i,k+2) + 32.0*dy8(i,k+3) - 3.0*dy8(i,k+4));
    else if (fit_y_3 && fit_z_3) *fyz = Xdydz * (-dy6(i,k-3) + 9.0*dy6(i,k-2) - 45.0*dy6(i,k-1) + 45.0*dy6(i,k+1) - 9.0*dy6(i,k+2) + dy6(i,k+3));
    else if (fit_y_2 && fit_z_2) *fyz = Fdydz * (dz4(i,j-2) - F8*dz4(i,j-1) + F8*dz4(i,j+1) - dz4(i,j+2));
    else if (fit_y_1 && fit_z_1) *fyz = Sdydz * (fh(i,j-1,k-1) - fh(i,j+1,k-1) - fh(i,j-1,k+1) + fh(i,j+1,k+1));
#endif

    (void)onoff;
}

#endif // DIFF_NEW_GPU_CUH
