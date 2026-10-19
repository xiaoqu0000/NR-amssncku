// diff_newwb_gpu.cuh — wide-band finite differences (one-sided boundary fallback)
// Port of diff_newwb.f90. Key difference from diff_new_gpu.cuh:
//   At boundaries, uses ONE-SIDED forward/backward stencils instead of returning 0.
// Supports ghost_width = 2,3,4,5 (2nd,4th,6th,8th order interiors).
// NOTE: Stencil coefficients are inlined (not via helper functions) to avoid
//       an nvcc optimization bug with if-constexpr + __forceinline__ helpers.

#ifndef DIFF_NEWWB_GPU_CUH
#define DIFF_NEWWB_GPU_CUH

#include "derivatives.h"
#include "fmisc.h"
#include "fmisc_gpu.cuh"
#include "macrodef.fh"
#include <cmath>

// ============================================================================
// First derivatives with one-sided boundary fallback
// Direct function (non-template) to avoid an nvcc template + if-constexpr optimization bug.
// Supports ghost_width = 2,3,4,5 via the compile-time constant `ghost_width`
// and the corresponding `if constexpr` stencil branches below.
// ============================================================================
__device__ void d_fderivs_newwb_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fx, double* fy, double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff,
    int i, int j, int k
) {
    constexpr int GW = ghost_width, ORD = GW - 1;
    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];
    const int NO_SYMM = 0, EQ_SYMM = 1;

    int imax = ex[0] - 1, jmax = ex[1] - 1, kmax = ex[2] - 1;
    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -ORD;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -ORD;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -ORD;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b_pad(ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
                                     f, ii, jj, kk, SYM1, SYM2, SYM3);
    };

    *fx = 0.0; *fy = 0.0; *fz = 0.0;

    // --- x-direction ---
    do {
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin) {
                *fx = (3.0*fh(i-4,j,k) - 32.0*fh(i-3,j,k) + 168.0*fh(i-2,j,k) - 672.0*fh(i-1,j,k)
                    + 672.0*fh(i+1,j,k) - 168.0*fh(i+2,j,k) + 32.0*fh(i+3,j,k) - 3.0*fh(i+4,j,k))
                    / (840.0 * dX);
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (i+3 <= imax && i-3 >= imin) {
                *fx = (-fh(i-3,j,k) + 9.0*fh(i-2,j,k) - 45.0*fh(i-1,j,k) + 45.0*fh(i+1,j,k)
                     - 9.0*fh(i+2,j,k) + fh(i+3,j,k)) / (60.0 * dX);
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (i+2 <= imax && i-2 >= imin) {
                *fx = (fh(i-2,j,k) - 8.0*fh(i-1,j,k) + 8.0*fh(i+1,j,k) - fh(i+2,j,k)) / (12.0 * dX);
                break;
            }
        }
        if (i+1 <= imax && i-1 >= imin) {
            *fx = (-fh(i-1,j,k) + fh(i+1,j,k)) / (2.0 * dX);
            break;
        }
        if (i == imin) {
            *fx = (-fh(i,j,k) + fh(i+1,j,k)) / dX;
        } else if (i == imax) {
            *fx = (-fh(i-1,j,k) + fh(i,j,k)) / dX;
        }
    } while(false);

    // --- y-direction ---
    do {
        if constexpr (GW >= 5) {
            if (j+4 <= jmax && j-4 >= jmin) {
                *fy = (3.0*fh(i,j-4,k) - 32.0*fh(i,j-3,k) + 168.0*fh(i,j-2,k) - 672.0*fh(i,j-1,k)
                    + 672.0*fh(i,j+1,k) - 168.0*fh(i,j+2,k) + 32.0*fh(i,j+3,k) - 3.0*fh(i,j+4,k))
                    / (840.0 * dY);
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (j+3 <= jmax && j-3 >= jmin) {
                *fy = (-fh(i,j-3,k) + 9.0*fh(i,j-2,k) - 45.0*fh(i,j-1,k) + 45.0*fh(i,j+1,k)
                     - 9.0*fh(i,j+2,k) + fh(i,j+3,k)) / (60.0 * dY);
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (j+2 <= jmax && j-2 >= jmin) {
                *fy = (fh(i,j-2,k) - 8.0*fh(i,j-1,k) + 8.0*fh(i,j+1,k) - fh(i,j+2,k)) / (12.0 * dY);
                break;
            }
        }
        if (j+1 <= jmax && j-1 >= jmin) {
            *fy = (-fh(i,j-1,k) + fh(i,j+1,k)) / (2.0 * dY);
            break;
        }
        if (j == jmin) {
            *fy = (-fh(i,j,k) + fh(i,j+1,k)) / dY;
        } else if (j == jmax) {
            *fy = (-fh(i,j-1,k) + fh(i,j,k)) / dY;
        }
    } while(false);

    // --- z-direction ---
    do {
        if constexpr (GW >= 5) {
            if (k+4 <= kmax && k-4 >= kmin) {
                *fz = (3.0*fh(i,j,k-4) - 32.0*fh(i,j,k-3) + 168.0*fh(i,j,k-2) - 672.0*fh(i,j,k-1)
                    + 672.0*fh(i,j,k+1) - 168.0*fh(i,j,k+2) + 32.0*fh(i,j,k+3) - 3.0*fh(i,j,k+4))
                    / (840.0 * dZ);
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (k+3 <= kmax && k-3 >= kmin) {
                *fz = (-fh(i,j,k-3) + 9.0*fh(i,j,k-2) - 45.0*fh(i,j,k-1) + 45.0*fh(i,j,k+1)
                     - 9.0*fh(i,j,k+2) + fh(i,j,k+3)) / (60.0 * dZ);
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (k+2 <= kmax && k-2 >= kmin) {
                *fz = (fh(i,j,k-2) - 8.0*fh(i,j,k-1) + 8.0*fh(i,j,k+1) - fh(i,j,k+2)) / (12.0 * dZ);
                break;
            }
        }
        if (k+1 <= kmax && k-1 >= kmin) {
            *fz = (-fh(i,j,k-1) + fh(i,j,k+1)) / (2.0 * dZ);
            break;
        }
        if (k == kmin) {
            *fz = (-fh(i,j,k) + fh(i,j,k+1)) / dZ;
        } else if (k == kmax) {
            *fz = (-fh(i,j,k-1) + fh(i,j,k)) / dZ;
        }
    } while(false);

    (void)onoff;
}

// ============================================================================
// Template: Second derivatives with one-sided boundary fallback
// ============================================================================
__device__ void d_fdderivs_newwb_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff,
    int i, int j, int k
) {
    constexpr int GW = ghost_width, ORD = GW - 1;
    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];
    const double Sdxdx = 1.0 / (dX * dX);
    const double Sdydy = 1.0 / (dY * dY);
    const double Sdzdz = 1.0 / (dZ * dZ);
    const double F4dxdx = 1.0 / (12.0 * dX * dX);
    const double F4dydy = 1.0 / (12.0 * dY * dY);
    const double F4dzdz = 1.0 / (12.0 * dZ * dZ);
    const double Sdxdy = 0.25 / (dX * dY);
    const double Sdxdz = 0.25 / (dX * dZ);
    const double Sdydz = 0.25 / (dY * dZ);
    const double F4dxdy = 1.0 / (144.0 * dX * dY);
    const double F4dxdz = 1.0 / (144.0 * dX * dZ);
    const double F4dydz = 1.0 / (144.0 * dY * dZ);

    const int NO_SYMM = 0, EQ_SYMM = 1;
    int imax = ex[0] - 1, jmax = ex[1] - 1, kmax = ex[2] - 1;
    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -ORD;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -ORD;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -ORD;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b_pad(ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
                                     f, ii, jj, kk, SYM1, SYM2, SYM3);
    };

    *fxx = 0.0; *fxy = 0.0; *fxz = 0.0;
    *fyy = 0.0; *fyz = 0.0; *fzz = 0.0;

    // --- fxx ---
    do {
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin) {
                *fxx = (-9.0*fh(i-4,j,k) + 128.0*fh(i-3,j,k) - 1008.0*fh(i-2,j,k)
                       + 8064.0*fh(i-1,j,k) - 14350.0*fh(i,j,k) + 8064.0*fh(i+1,j,k)
                       - 1008.0*fh(i+2,j,k) + 128.0*fh(i+3,j,k) - 9.0*fh(i+4,j,k))
                       / (5040.0 * dX * dX);
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (i+3 <= imax && i-3 >= imin) {
                *fxx = (2.0*fh(i-3,j,k) - 27.0*fh(i-2,j,k) + 270.0*fh(i-1,j,k)
                      - 490.0*fh(i,j,k) + 270.0*fh(i+1,j,k) - 27.0*fh(i+2,j,k)
                      + 2.0*fh(i+3,j,k)) / (180.0 * dX * dX);
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (i+2 <= imax && i-2 >= imin) {
                *fxx = F4dxdx * (-fh(i-2,j,k) + 16.0*fh(i-1,j,k) - 30.0*fh(i,j,k)
                                 - fh(i+2,j,k) + 16.0*fh(i+1,j,k));
                break;
            }
        }
        if (i+1 <= imax && i-1 >= imin) {
            *fxx = Sdxdx * (fh(i-1,j,k) - 2.0*fh(i,j,k) + fh(i+1,j,k));
            break;
        }
        if (i == imin)
            *fxx = (fh(i,j,k) - 2.0*fh(i+1,j,k) + fh(i+2,j,k)) / (dX * dX);
        else if (i == imax)
            *fxx = (fh(i-2,j,k) - 2.0*fh(i-1,j,k) + fh(i,j,k)) / (dX * dX);
    } while(false);

    // --- fyy ---
    do {
        if constexpr (GW >= 5) {
            if (j+4 <= jmax && j-4 >= jmin) {
                *fyy = (-9.0*fh(i,j-4,k) + 128.0*fh(i,j-3,k) - 1008.0*fh(i,j-2,k)
                       + 8064.0*fh(i,j-1,k) - 14350.0*fh(i,j,k) + 8064.0*fh(i,j+1,k)
                       - 1008.0*fh(i,j+2,k) + 128.0*fh(i,j+3,k) - 9.0*fh(i,j+4,k))
                       / (5040.0 * dY * dY);
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (j+3 <= jmax && j-3 >= jmin) {
                *fyy = (2.0*fh(i,j-3,k) - 27.0*fh(i,j-2,k) + 270.0*fh(i,j-1,k)
                      - 490.0*fh(i,j,k) + 270.0*fh(i,j+1,k) - 27.0*fh(i,j+2,k)
                      + 2.0*fh(i,j+3,k)) / (180.0 * dY * dY);
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (j+2 <= jmax && j-2 >= jmin) {
                *fyy = F4dydy * (-fh(i,j-2,k) + 16.0*fh(i,j-1,k) - 30.0*fh(i,j,k)
                                 - fh(i,j+2,k) + 16.0*fh(i,j+1,k));
                break;
            }
        }
        if (j+1 <= jmax && j-1 >= jmin) {
            *fyy = Sdydy * (fh(i,j-1,k) - 2.0*fh(i,j,k) + fh(i,j+1,k));
            break;
        }
        if (j == jmin)
            *fyy = (fh(i,j,k) - 2.0*fh(i,j+1,k) + fh(i,j+2,k)) / (dY * dY);
        else if (j == jmax)
            *fyy = (fh(i,j-2,k) - 2.0*fh(i,j-1,k) + fh(i,j,k)) / (dY * dY);
    } while(false);

    // --- fzz ---
    do {
        if constexpr (GW >= 5) {
            if (k+4 <= kmax && k-4 >= kmin) {
                *fzz = (-9.0*fh(i,j,k-4) + 128.0*fh(i,j,k-3) - 1008.0*fh(i,j,k-2)
                       + 8064.0*fh(i,j,k-1) - 14350.0*fh(i,j,k) + 8064.0*fh(i,j,k+1)
                       - 1008.0*fh(i,j,k+2) + 128.0*fh(i,j,k+3) - 9.0*fh(i,j,k+4))
                       / (5040.0 * dZ * dZ);
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (k+3 <= kmax && k-3 >= kmin) {
                *fzz = (2.0*fh(i,j,k-3) - 27.0*fh(i,j,k-2) + 270.0*fh(i,j,k-1)
                      - 490.0*fh(i,j,k) + 270.0*fh(i,j,k+1) - 27.0*fh(i,j,k+2)
                      + 2.0*fh(i,j,k+3)) / (180.0 * dZ * dZ);
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (k+2 <= kmax && k-2 >= kmin) {
                *fzz = F4dzdz * (-fh(i,j,k-2) + 16.0*fh(i,j,k-1) - 30.0*fh(i,j,k)
                                 - fh(i,j,k+2) + 16.0*fh(i,j,k+1));
                break;
            }
        }
        if (k+1 <= kmax && k-1 >= kmin) {
            *fzz = Sdzdz * (fh(i,j,k-1) - 2.0*fh(i,j,k) + fh(i,j,k+1));
            break;
        }
        if (k == kmin)
            *fzz = (fh(i,j,k) - 2.0*fh(i,j,k+1) + fh(i,j,k+2)) / (dZ * dZ);
        else if (k == kmax)
            *fzz = (fh(i,j,k-2) - 2.0*fh(i,j,k-1) + fh(i,j,k)) / (dZ * dZ);
    } while(false);

    // --- fxy ---
    do {
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin && j+4 <= jmax && j-4 >= jmin) {
                const int off8[] = {-4, -3, -2, -1, 1, 2, 3, 4};
                const double c8[] = {3.0, -32.0, 168.0, -672.0, 672.0, -168.0, 32.0, -3.0};
                double s = 0.0;
                #pragma unroll
                for (int jn = 0; jn < 8; jn++) {
                    double row = 0.0;
                    #pragma unroll
                    for (int in = 0; in < 8; in++)
                        row += c8[in] * fh(i+off8[in], j+off8[jn], k);
                    s += c8[jn] * row;
                }
                *fxy = s / (840.0 * 840.0 * dX * dY);
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (i+3 <= imax && i-3 >= imin && j+3 <= jmax && j-3 >= jmin) {
                const int off6[] = {-3, -2, -1, 1, 2, 3};
                const double c6[] = {-1.0, 9.0, -45.0, 45.0, -9.0, 1.0};
                double s = 0.0;
                #pragma unroll
                for (int jn = 0; jn < 6; jn++) {
                    double row = 0.0;
                    #pragma unroll
                    for (int in = 0; in < 6; in++)
                        row += c6[in] * fh(i+off6[in], j+off6[jn], k);
                    s += c6[jn] * row;
                }
                *fxy = s / (60.0 * 60.0 * dX * dY);
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (i+2 <= imax && i-2 >= imin && j+2 <= jmax && j-2 >= jmin) {
                *fxy = F4dxdy * (
                    (fh(i-2,j-2,k) - 8.0*fh(i-1,j-2,k) + 8.0*fh(i+1,j-2,k) - fh(i+2,j-2,k))
                  - 8.0*(fh(i-2,j-1,k) - 8.0*fh(i-1,j-1,k) + 8.0*fh(i+1,j-1,k) - fh(i+2,j-1,k))
                  + 8.0*(fh(i-2,j+1,k) - 8.0*fh(i-1,j+1,k) + 8.0*fh(i+1,j+1,k) - fh(i+2,j+1,k))
                  -     (fh(i-2,j+2,k) - 8.0*fh(i-1,j+2,k) + 8.0*fh(i+1,j+2,k) - fh(i+2,j+2,k)));
                break;
            }
        }
        if (i+1 <= imax && i-1 >= imin && j+1 <= jmax && j-1 >= jmin) {
            *fxy = Sdxdy * (fh(i-1,j-1,k) - fh(i+1,j-1,k) - fh(i-1,j+1,k) + fh(i+1,j+1,k));
            break;
        }
        // One-sided corner/edge fallback (matches Fortran diff_newwb)
        if (i == imin && j+1 <= jmax && j-1 >= jmin)
            *fxy = (fh(i,j-1,k) - fh(i+1,j-1,k) - fh(i,j+1,k) + fh(i+1,j+1,k)) / (2.0 * dX * dY);
        else if (i == imax && j+1 <= jmax && j-1 >= jmin)
            *fxy = (fh(i-1,j-1,k) - fh(i,j-1,k) - fh(i-1,j+1,k) + fh(i,j+1,k)) / (2.0 * dX * dY);
        else if (j == jmin && i+1 <= imax && i-1 >= imin)
            *fxy = (fh(i-1,j,k) - fh(i-1,j+1,k) - fh(i+1,j,k) + fh(i+1,j+1,k)) / (2.0 * dX * dY);
        else if (j == jmax && i+1 <= imax && i-1 >= imin)
            *fxy = (fh(i-1,j-1,k) - fh(i-1,j,k) - fh(i+1,j-1,k) + fh(i+1,j,k)) / (2.0 * dX * dY);
        else if (i == imin && j == jmin)
            *fxy = (fh(i+1,j+1,k) - fh(i+1,j,k) - fh(i,j+1,k) + fh(i,j,k)) / (dX * dY);
        else if (i == imin && j == jmax)
            *fxy = (fh(i+1,j,k) - fh(i+1,j-1,k) - fh(i,j,k) + fh(i,j-1,k)) / (dX * dY);
        else if (i == imax && j == jmin)
            *fxy = (fh(i,j+1,k) - fh(i,j,k) - fh(i-1,j+1,k) + fh(i-1,j,k)) / (dX * dY);
        else if (i == imax && j == jmax)
            *fxy = (fh(i,j,k) - fh(i,j-1,k) - fh(i-1,j,k) + fh(i-1,j-1,k)) / (dX * dY);
    } while(false);

    // --- fxz ---
    do {
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin && k+4 <= kmax && k-4 >= kmin) {
                const int off8[] = {-4, -3, -2, -1, 1, 2, 3, 4};
                const double c8[] = {3.0, -32.0, 168.0, -672.0, 672.0, -168.0, 32.0, -3.0};
                double s = 0.0;
                #pragma unroll
                for (int kn = 0; kn < 8; kn++) {
                    double row = 0.0;
                    #pragma unroll
                    for (int in = 0; in < 8; in++)
                        row += c8[in] * fh(i+off8[in], j, k+off8[kn]);
                    s += c8[kn] * row;
                }
                *fxz = s / (840.0 * 840.0 * dX * dZ);
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (i+3 <= imax && i-3 >= imin && k+3 <= kmax && k-3 >= kmin) {
                const int off6[] = {-3, -2, -1, 1, 2, 3};
                const double c6[] = {-1.0, 9.0, -45.0, 45.0, -9.0, 1.0};
                double s = 0.0;
                #pragma unroll
                for (int kn = 0; kn < 6; kn++) {
                    double row = 0.0;
                    #pragma unroll
                    for (int in = 0; in < 6; in++)
                        row += c6[in] * fh(i+off6[in], j, k+off6[kn]);
                    s += c6[kn] * row;
                }
                *fxz = s / (60.0 * 60.0 * dX * dZ);
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (i+2 <= imax && i-2 >= imin && k+2 <= kmax && k-2 >= kmin) {
                *fxz = F4dxdz * (
                    (fh(i-2,j,k-2) - 8.0*fh(i-1,j,k-2) + 8.0*fh(i+1,j,k-2) - fh(i+2,j,k-2))
                  - 8.0*(fh(i-2,j,k-1) - 8.0*fh(i-1,j,k-1) + 8.0*fh(i+1,j,k-1) - fh(i+2,j,k-1))
                  + 8.0*(fh(i-2,j,k+1) - 8.0*fh(i-1,j,k+1) + 8.0*fh(i+1,j,k+1) - fh(i+2,j,k+1))
                  -     (fh(i-2,j,k+2) - 8.0*fh(i-1,j,k+2) + 8.0*fh(i+1,j,k+2) - fh(i+2,j,k+2)));
                break;
            }
        }
        if (i+1 <= imax && i-1 >= imin && k+1 <= kmax && k-1 >= kmin) {
            *fxz = Sdxdz * (fh(i-1,j,k-1) - fh(i+1,j,k-1) - fh(i-1,j,k+1) + fh(i+1,j,k+1));
            break;
        }
        if (i == imin && k+1 <= kmax && k-1 >= kmin)
            *fxz = (fh(i,j,k-1) - fh(i+1,j,k-1) - fh(i,j,k+1) + fh(i+1,j,k+1)) / (2.0 * dX * dZ);
        else if (i == imax && k+1 <= kmax && k-1 >= kmin)
            *fxz = (fh(i-1,j,k-1) - fh(i,j,k-1) - fh(i-1,j,k+1) + fh(i,j,k+1)) / (2.0 * dX * dZ);
        else if (k == kmin && i+1 <= imax && i-1 >= imin)
            *fxz = (fh(i-1,j,k) - fh(i-1,j,k+1) - fh(i+1,j,k) + fh(i+1,j,k+1)) / (2.0 * dX * dZ);
        else if (k == kmax && i+1 <= imax && i-1 >= imin)
            *fxz = (fh(i-1,j,k-1) - fh(i-1,j,k) - fh(i+1,j,k-1) + fh(i+1,j,k)) / (2.0 * dX * dZ);
        else if (i == imin && k == kmin)
            *fxz = (fh(i+1,j,k+1) - fh(i+1,j,k) - fh(i,j,k+1) + fh(i,j,k)) / (dX * dZ);
        else if (i == imin && k == kmax)
            *fxz = (fh(i+1,j,k) - fh(i+1,j,k-1) - fh(i,j,k) + fh(i,j,k-1)) / (dX * dZ);
        else if (i == imax && k == kmin)
            *fxz = (fh(i,j,k+1) - fh(i,j,k) - fh(i-1,j,k+1) + fh(i-1,j,k)) / (dX * dZ);
        else if (i == imax && k == kmax)
            *fxz = (fh(i,j,k) - fh(i,j,k-1) - fh(i-1,j,k) + fh(i-1,j,k-1)) / (dX * dZ);
    } while(false);

    // --- fyz ---
    do {
        if constexpr (GW >= 5) {
            if (j+4 <= jmax && j-4 >= jmin && k+4 <= kmax && k-4 >= kmin) {
                const int off8[] = {-4, -3, -2, -1, 1, 2, 3, 4};
                const double c8[] = {3.0, -32.0, 168.0, -672.0, 672.0, -168.0, 32.0, -3.0};
                double s = 0.0;
                #pragma unroll
                for (int kn = 0; kn < 8; kn++) {
                    double row = 0.0;
                    #pragma unroll
                    for (int jn = 0; jn < 8; jn++)
                        row += c8[jn] * fh(i, j+off8[jn], k+off8[kn]);
                    s += c8[kn] * row;
                }
                *fyz = s / (840.0 * 840.0 * dY * dZ);
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (j+3 <= jmax && j-3 >= jmin && k+3 <= kmax && k-3 >= kmin) {
                const int off6[] = {-3, -2, -1, 1, 2, 3};
                const double c6[] = {-1.0, 9.0, -45.0, 45.0, -9.0, 1.0};
                double s = 0.0;
                #pragma unroll
                for (int kn = 0; kn < 6; kn++) {
                    double row = 0.0;
                    #pragma unroll
                    for (int jn = 0; jn < 6; jn++)
                        row += c6[jn] * fh(i, j+off6[jn], k+off6[kn]);
                    s += c6[kn] * row;
                }
                *fyz = s / (60.0 * 60.0 * dY * dZ);
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (j+2 <= jmax && j-2 >= jmin && k+2 <= kmax && k-2 >= kmin) {
                *fyz = F4dydz * (
                    (fh(i,j-2,k-2) - 8.0*fh(i,j-1,k-2) + 8.0*fh(i,j+1,k-2) - fh(i,j+2,k-2))
                  - 8.0*(fh(i,j-2,k-1) - 8.0*fh(i,j-1,k-1) + 8.0*fh(i,j+1,k-1) - fh(i,j+2,k-1))
                  + 8.0*(fh(i,j-2,k+1) - 8.0*fh(i,j-1,k+1) + 8.0*fh(i,j+1,k+1) - fh(i,j+2,k+1))
                  -     (fh(i,j-2,k+2) - 8.0*fh(i,j-1,k+2) + 8.0*fh(i,j+1,k+2) - fh(i,j+2,k+2)));
                break;
            }
        }
        if (j+1 <= jmax && j-1 >= jmin && k+1 <= kmax && k-1 >= kmin) {
            *fyz = Sdydz * (fh(i,j-1,k-1) - fh(i,j+1,k-1) - fh(i,j-1,k+1) + fh(i,j+1,k+1));
            break;
        }
        if (j == jmin && k+1 <= kmax && k-1 >= kmin)
            *fyz = (fh(i,j,k-1) - fh(i,j+1,k-1) - fh(i,j,k+1) + fh(i,j+1,k+1)) / (2.0 * dY * dZ);
        else if (j == jmax && k+1 <= kmax && k-1 >= kmin)
            *fyz = (fh(i,j-1,k-1) - fh(i,j,k-1) - fh(i,j-1,k+1) + fh(i,j,k+1)) / (2.0 * dY * dZ);
        else if (k == kmin && j+1 <= jmax && j-1 >= jmin)
            *fyz = (fh(i,j,k) - fh(i,j,k+1) - fh(i,j+1,k) + fh(i,j+1,k+1)) / (2.0 * dY * dZ);
        else if (k == kmax && j+1 <= jmax && j-1 >= jmin)
            *fyz = (fh(i,j-1,k-1) - fh(i,j-1,k) - fh(i,j,k-1) + fh(i,j,k)) / (2.0 * dY * dZ);
        else if (j == jmin && k == kmin)
            *fyz = (fh(i,j+1,k+1) - fh(i,j+1,k) - fh(i,j,k+1) + fh(i,j,k)) / (dY * dZ);
        else if (j == jmin && k == kmax)
            *fyz = (fh(i,j+1,k) - fh(i,j+1,k-1) - fh(i,j,k) + fh(i,j,k-1)) / (dY * dZ);
        else if (j == jmax && k == kmin)
            *fyz = (fh(i,j,k+1) - fh(i,j,k) - fh(i,j-1,k+1) + fh(i,j-1,k)) / (dY * dZ);
        else if (j == jmax && k == kmax)
            *fyz = (fh(i,j,k) - fh(i,j,k-1) - fh(i,j-1,k) + fh(i,j-1,k-1)) / (dY * dZ);
    } while(false);

    (void)onoff;
}

// ============================================================================
// Single-direction first derivatives
// ============================================================================
__device__ void d_fdx_newwb_gw(
    const int ex[3], const int pex[3], const double* f, double* fx,
    const double* X, double SYM1, int symmetry, int onoff,
    int i, int j, int k
) {
    constexpr int GW = ghost_width, ORD = GW - 1;
    const double dX = X[1] - X[0];
    const int EQ_SYMM = 1;
    int imax = ex[0] - 1, imin = 0;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -ORD;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b_pad(ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
                                     f, ii, jj, kk, SYM1, SYM1, SYM1);
    };

    *fx = 0.0;
    do {
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin) {
                *fx = (3.0*fh(i-4,j,k) - 32.0*fh(i-3,j,k) + 168.0*fh(i-2,j,k) - 672.0*fh(i-1,j,k)
                    + 672.0*fh(i+1,j,k) - 168.0*fh(i+2,j,k) + 32.0*fh(i+3,j,k) - 3.0*fh(i+4,j,k))
                    / (840.0 * dX); break;
            }
        }
        if constexpr (GW >= 4) {
            if (i+3 <= imax && i-3 >= imin) {
                *fx = (-fh(i-3,j,k) + 9.0*fh(i-2,j,k) - 45.0*fh(i-1,j,k) + 45.0*fh(i+1,j,k)
                     - 9.0*fh(i+2,j,k) + fh(i+3,j,k)) / (60.0 * dX); break;
            }
        }
        if constexpr (GW >= 3) {
            if (i+2 <= imax && i-2 >= imin) {
                *fx = (fh(i-2,j,k) - 8.0*fh(i-1,j,k) + 8.0*fh(i+1,j,k) - fh(i+2,j,k)) / (12.0 * dX); break;
            }
        }
        if (i+1 <= imax && i-1 >= imin) { *fx = (-fh(i-1,j,k) + fh(i+1,j,k)) / (2.0 * dX); break; }
        if (i == imin) { *fx = (-fh(i,j,k) + fh(i+1,j,k)) / dX; }
        else if (i == imax) { *fx = (-fh(i-1,j,k) + fh(i,j,k)) / dX; }
    } while(false);
    (void)onoff;
}

__device__ void d_fdy_newwb_gw(
    const int ex[3], const int pex[3], const double* f, double* fy,
    const double* Y, double SYM2, int symmetry, int onoff,
    int i, int j, int k
) {
    constexpr int GW = ghost_width, ORD = GW - 1;
    const double dY = Y[1] - Y[0];
    const int EQ_SYMM = 1;
    int jmax = ex[1] - 1, jmin = 0;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -ORD;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b_pad(ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
                                     f, ii, jj, kk, SYM2, SYM2, SYM2);
    };

    *fy = 0.0;
    do {
        if constexpr (GW >= 5) {
            if (j+4 <= jmax && j-4 >= jmin) {
                *fy = (3.0*fh(i,j-4,k) - 32.0*fh(i,j-3,k) + 168.0*fh(i,j-2,k) - 672.0*fh(i,j-1,k)
                    + 672.0*fh(i,j+1,k) - 168.0*fh(i,j+2,k) + 32.0*fh(i,j+3,k) - 3.0*fh(i,j+4,k))
                    / (840.0 * dY); break;
            }
        }
        if constexpr (GW >= 4) {
            if (j+3 <= jmax && j-3 >= jmin) {
                *fy = (-fh(i,j-3,k) + 9.0*fh(i,j-2,k) - 45.0*fh(i,j-1,k) + 45.0*fh(i,j+1,k)
                     - 9.0*fh(i,j+2,k) + fh(i,j+3,k)) / (60.0 * dY); break;
            }
        }
        if constexpr (GW >= 3) {
            if (j+2 <= jmax && j-2 >= jmin) {
                *fy = (fh(i,j-2,k) - 8.0*fh(i,j-1,k) + 8.0*fh(i,j+1,k) - fh(i,j+2,k)) / (12.0 * dY); break;
            }
        }
        if (j+1 <= jmax && j-1 >= jmin) { *fy = (-fh(i,j-1,k) + fh(i,j+1,k)) / (2.0 * dY); break; }
        if (j == jmin) { *fy = (-fh(i,j,k) + fh(i,j+1,k)) / dY; }
        else if (j == jmax) { *fy = (-fh(i,j-1,k) + fh(i,j,k)) / dY; }
    } while(false);
    (void)onoff;
}

__device__ void d_fdz_newwb_gw(
    const int ex[3], const int pex[3], const double* f, double* fz,
    const double* Z, double SYM3, int symmetry, int onoff,
    int i, int j, int k
) {
    constexpr int GW = ghost_width, ORD = GW - 1;
    const double dZ = Z[1] - Z[0];
    const int NO_SYMM = 0;
    int kmax = ex[2] - 1, kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -ORD;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b_pad(ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
                                     f, ii, jj, kk, SYM3, SYM3, SYM3);
    };

    *fz = 0.0;
    do {
        if constexpr (GW >= 5) {
            if (k+4 <= kmax && k-4 >= kmin) {
                *fz = (3.0*fh(i,j,k-4) - 32.0*fh(i,j,k-3) + 168.0*fh(i,j,k-2) - 672.0*fh(i,j,k-1)
                    + 672.0*fh(i,j,k+1) - 168.0*fh(i,j,k+2) + 32.0*fh(i,j,k+3) - 3.0*fh(i,j,k+4))
                    / (840.0 * dZ); break;
            }
        }
        if constexpr (GW >= 4) {
            if (k+3 <= kmax && k-3 >= kmin) {
                *fz = (-fh(i,j,k-3) + 9.0*fh(i,j,k-2) - 45.0*fh(i,j,k-1) + 45.0*fh(i,j,k+1)
                     - 9.0*fh(i,j,k+2) + fh(i,j,k+3)) / (60.0 * dZ); break;
            }
        }
        if constexpr (GW >= 3) {
            if (k+2 <= kmax && k-2 >= kmin) {
                *fz = (fh(i,j,k-2) - 8.0*fh(i,j,k-1) + 8.0*fh(i,j,k+1) - fh(i,j,k+2)) / (12.0 * dZ); break;
            }
        }
        if (k+1 <= kmax && k-1 >= kmin) { *fz = (-fh(i,j,k-1) + fh(i,j,k+1)) / (2.0 * dZ); break; }
        if (k == kmin) { *fz = (-fh(i,j,k) + fh(i,j,k+1)) / dZ; }
        else if (k == kmax) { *fz = (-fh(i,j,k-1) + fh(i,j,k)) / dZ; }
    } while(false);
    (void)onoff;
}

// ============================================================================
// Single second derivative (fddxx)
// ============================================================================
__device__ void d_fddxx_newwb_gw(
    const int ex[3], const int pex[3], const double* f, double* fxx,
    const double* X, double SYM1, int symmetry, int onoff,
    int i, int j, int k
) {
    constexpr int GW = ghost_width, ORD = GW - 1;
    const double dX = X[1] - X[0];
    const double Sdxdx = 1.0/(dX*dX), F4dxdx = 1.0/(12.0*dX*dX);
    const int EQ_SYMM = 1;
    int imax = ex[0] - 1, imin = 0;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -ORD;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b_pad(ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
                                     f, ii, jj, kk, SYM1, SYM1, SYM1);
    };
    *fxx = 0.0;
    do {
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin) {
                *fxx = (-9.0*fh(i-4,j,k)+128.0*fh(i-3,j,k)-1008.0*fh(i-2,j,k)+8064.0*fh(i-1,j,k)
                        -14350.0*fh(i,j,k)+8064.0*fh(i+1,j,k)-1008.0*fh(i+2,j,k)+128.0*fh(i+3,j,k)
                        -9.0*fh(i+4,j,k))/(5040.0*dX*dX); break;
            }
        }
        if constexpr (GW >= 4) {
            if (i+3 <= imax && i-3 >= imin) {
                *fxx = (2.0*fh(i-3,j,k)-27.0*fh(i-2,j,k)+270.0*fh(i-1,j,k)-490.0*fh(i,j,k)
                        +270.0*fh(i+1,j,k)-27.0*fh(i+2,j,k)+2.0*fh(i+3,j,k))/(180.0*dX*dX); break;
            }
        }
        if constexpr (GW >= 3) {
            if (i+2 <= imax && i-2 >= imin) {
                *fxx = F4dxdx*(-fh(i-2,j,k)+16.0*fh(i-1,j,k)-30.0*fh(i,j,k)-fh(i+2,j,k)+16.0*fh(i+1,j,k));
                break;
            }
        }
        if (i+1 <= imax && i-1 >= imin) {
            *fxx = Sdxdx*(fh(i-1,j,k)-2.0*fh(i,j,k)+fh(i+1,j,k)); break;
        }
        if (i == imin) { *fxx = (fh(i,j,k)-2.0*fh(i+1,j,k)+fh(i+2,j,k))/(dX*dX); }
        else if (i == imax) { *fxx = (fh(i-2,j,k)-2.0*fh(i-1,j,k)+fh(i,j,k))/(dX*dX); }
    } while(false);
    (void)onoff;
}

// ============================================================================
// fdderivsdavid (GW=2 only) — strided 2nd derivative
// ============================================================================
#if ghost_width == 2
__device__ void d_fdderivsdavid_newwb(
    const int ex[3], const int pex[3], const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry,
    int i, int j, int k
) {
    const int ORD = 1;
    const double dX = X[1] - X[0], dY = Y[1] - Y[0], dZ = Z[1] - Z[0];
    const int NO_SYMM = 0, EQ_SYMM = 1;
    int imax = ex[0] - 1, jmax = ex[1] - 1, kmax = ex[2] - 1;
    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -ORD;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -ORD;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -ORD;

    *fxx = 0.0; *fyy = 0.0; *fzz = 0.0;
    *fxy = 0.0; *fxz = 0.0; *fyz = 0.0;

    // Fortran fdderivsdavid loops 1..ex-1 (C: 0..ex-2), no one-sided fallback
    if (i >= imax || j >= jmax || k >= kmax) return;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b_pad(ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
                                     f, ii, jj, kk, SYM1, SYM2, SYM3);
    };

    if (i+2 <= imax && i-2 >= imin)
        *fxx = (fh(i-2,j,k) - 2.0*fh(i,j,k) + fh(i+2,j,k)) / (4.0 * dX * dX);
    else if (i+1 <= imax && i-1 >= imin)
        *fxx = (fh(i-1,j,k) - 2.0*fh(i,j,k) + fh(i+1,j,k)) / (dX * dX);

    if (j+2 <= jmax && j-2 >= jmin)
        *fyy = (fh(i,j-2,k) - 2.0*fh(i,j,k) + fh(i,j+2,k)) / (4.0 * dY * dY);
    else if (j+1 <= jmax && j-1 >= jmin)
        *fyy = (fh(i,j-1,k) - 2.0*fh(i,j,k) + fh(i,j+1,k)) / (dY * dY);

    if (k+2 <= kmax && k-2 >= kmin)
        *fzz = (fh(i,j,k-2) - 2.0*fh(i,j,k) + fh(i,j,k+2)) / (4.0 * dZ * dZ);
    else if (k+1 <= kmax && k-1 >= kmin)
        *fzz = (fh(i,j,k-1) - 2.0*fh(i,j,k) + fh(i,j,k+1)) / (dZ * dZ);
}
#endif // ghost_width == 2

// ============================================================================
// Backward-compatible wrappers
// ============================================================================
__device__ __forceinline__ void d_fderivs_newwb_point(
    const int ex[3], const int pex[3], const double* f,
    double* fx, double* fy, double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff,
    int i, int j, int k
) {
    d_fderivs_newwb_gw(ex, pex, f, fx, fy, fz, X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, i, j, k);
}

__device__ __forceinline__ void d_fdderivs_newwb_point(
    const int ex[3], const int pex[3], const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff,
    int i, int j, int k
) {
    d_fdderivs_newwb_gw(ex, pex, f, fxx, fxy, fxz, fyy, fyz, fzz, X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, i, j, k);
}

__device__ __forceinline__ void d_fdx_newwb_point(
    const int ex[3], const int pex[3], const double* f, double* fx,
    const double* X, double SYM1, int symmetry, int onoff,
    int i, int j, int k
) {
    d_fdx_newwb_gw(ex, pex, f, fx, X, SYM1, symmetry, onoff, i, j, k);
}

__device__ __forceinline__ void d_fdy_newwb_point(
    const int ex[3], const int pex[3], const double* f, double* fy,
    const double* Y, double SYM2, int symmetry, int onoff,
    int i, int j, int k
) {
    d_fdy_newwb_gw(ex, pex, f, fy, Y, SYM2, symmetry, onoff, i, j, k);
}

__device__ __forceinline__ void d_fdz_newwb_point(
    const int ex[3], const int pex[3], const double* f, double* fz,
    const double* Z, double SYM3, int symmetry, int onoff,
    int i, int j, int k
) {
    d_fdz_newwb_gw(ex, pex, f, fz, Z, SYM3, symmetry, onoff, i, j, k);
}

__device__ __forceinline__ void d_fddxx_newwb_point(
    const int ex[3], const int pex[3], const double* f, double* fxx,
    const double* X, double SYM1, int symmetry, int onoff,
    int i, int j, int k
) {
    d_fddxx_newwb_gw(ex, pex, f, fxx, X, SYM1, symmetry, onoff, i, j, k);
}

#endif // DIFF_NEWWB_GPU_CUH
