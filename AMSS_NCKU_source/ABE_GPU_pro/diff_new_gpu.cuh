#ifndef DIFF_NEW_GPU_CUH
#define DIFF_NEW_GPU_CUH

#include "derivatives.h"

#include "fmisc.h"
#include "fmisc_gpu.cuh"

#include "macrodef.fh"
#include <cmath>
#include <iostream>

// ============================================================
// Stencil coefficient tables for cross derivatives
// ============================================================
namespace diff_new_detail {

// We use __device__ inline functions instead of static constexpr arrays
// to ensure broad CUDA compatibility

__device__ __forceinline__ void get_fd1_2(double* c, int* o, double& d) {
    d = 2.0;
    c[0] = -1.0; c[1] = 1.0;
    o[0] = -1; o[1] = 1;
}

__device__ __forceinline__ void get_fd1_4(double* c, int* o, double& d) {
    d = 12.0;
    c[0] = 1.0; c[1] = -8.0; c[2] = 8.0; c[3] = -1.0;
    o[0] = -2; o[1] = -1; o[2] = 1; o[3] = 2;
}

__device__ __forceinline__ void get_fd1_6(double* c, int* o, double& d) {
    d = 60.0;
    c[0] = -1.0; c[1] = 9.0; c[2] = -45.0; c[3] = 45.0; c[4] = -9.0; c[5] = 1.0;
    o[0] = -3; o[1] = -2; o[2] = -1; o[3] = 1; o[4] = 2; o[5] = 3;
}

__device__ __forceinline__ void get_fd1_8(double* c, int* o, double& d) {
    d = 840.0;
    c[0] = 3.0; c[1] = -32.0; c[2] = 168.0; c[3] = -672.0;
    c[4] = 672.0; c[5] = -168.0; c[6] = 32.0; c[7] = -3.0;
    o[0] = -4; o[1] = -3; o[2] = -2; o[3] = -1;
    o[4] = 1; o[5] = 2; o[6] = 3; o[7] = 4;
}

// Cross derivative helper: computes outer product of two first-derivative stencils
// Axis1 = inner (x-type), Axis2 = outer (y-type)
// 0=x, 1=y, 2=z
template<int NC, int Axis1, int Axis2, typename FH>
__device__ __forceinline__ double cross_deriv(
    FH fh, int i, int j, int k,
    const double* c, const int* off, double denom2, double inv_h1h2
) {
    double sum = 0.0;
    #pragma unroll
    for (int jn = 0; jn < NC; jn++) {
        double row = 0.0;
        #pragma unroll
        for (int in = 0; in < NC; in++) {
            int ii = i, jj = j, kk = k;
            if constexpr (Axis1 == 0) ii += off[in];
            else if constexpr (Axis1 == 1) jj += off[in];
            else kk += off[in];
            if constexpr (Axis2 == 0) ii += off[jn];
            else if constexpr (Axis2 == 1) jj += off[jn];
            else kk += off[jn];
            row += c[in] * fh(ii, jj, kk);
        }
        sum += c[jn] * row;
    }
    return sum / denom2 * inv_h1h2;
}

} // namespace diff_new_detail

// ============================================================
// Multi-GW first derivatives (device function)
// ============================================================
template<int GW>
__device__ __forceinline__ void d_fderivs_point_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fx, double* fy, double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff,
    int i, int j, int k
) {
    constexpr int ORD = GW - 1;

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    const int imax = ex[0] - 1;
    const int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;

    *fx = 0.0; *fy = 0.0; *fz = 0.0;
    if (i >= imax || j >= jmax || k >= kmax) return;

    // Match Fortran fderivs: symmetry_bd fills (GW-1) ghost layers below first point.
    // Fortran kmin=-1 is 1-indexed (2 layers below k=1), equivalent to C 0-indexed -(GW-1).
    constexpr int sym_ghost = GW - 1;
    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > 0 && fabs(Z[0]) < dZ) kmin = -sym_ghost;
    if (symmetry > 1 && fabs(X[0]) < dX) imin = -sym_ghost;
    if (symmetry > 1 && fabs(Y[0]) < dY) jmin = -sym_ghost;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b_pad(
            ORD,
            ex[0], ex[1], ex[2],
            pex[0], pex[1], pex[2],
            f, ii, jj, kk,
            SYM1, SYM2, SYM3
        );
    };

    if constexpr (GW == 2) {
        const double d2dx = 0.5/dX;
        const double d2dy = 0.5/dY;
        const double d2dz = 0.5/dZ;
        if (i+1 <= imax && i-1 >= imin) {
            *fx = d2dx*(-fh(i-1,j,k) + fh(i+1,j,k));
        }
        if (j+1 <= jmax && j-1 >= jmin) {
            *fy = d2dy*(-fh(i,j-1,k) + fh(i,j+1,k));
        }
        if (k+1 <= kmax && k-1 >= kmin) {
            *fz = d2dz*(-fh(i,j,k-1) + fh(i,j,k+1));
        }
    } else if constexpr (GW == 3) {
        // GW=3: isotropic fallback (matches Fortran diff_new.f90 #else "for bam comparison")
        if (i+2 <= imax && i-2 >= imin &&
            j+2 <= jmax && j-2 >= jmin &&
            k+2 <= kmax && k-2 >= kmin) {
            const double d12dx = 1.0/12.0/dX;
            const double d12dy = 1.0/12.0/dY;
            const double d12dz = 1.0/12.0/dZ;
            *fx = d12dx*(fh(i-2,j,k) - 8*fh(i-1,j,k) + 8*fh(i+1,j,k) - fh(i+2,j,k));
            *fy = d12dy*(fh(i,j-2,k) - 8*fh(i,j-1,k) + 8*fh(i,j+1,k) - fh(i,j+2,k));
            *fz = d12dz*(fh(i,j,k-2) - 8*fh(i,j,k-1) + 8*fh(i,j,k+1) - fh(i,j,k+2));
            return;
        }
        if (i+1 <= imax && i-1 >= imin &&
            j+1 <= jmax && j-1 >= jmin &&
            k+1 <= kmax && k-1 >= kmin) {
            const double d2dx = 0.5/dX;
            const double d2dy = 0.5/dY;
            const double d2dz = 0.5/dZ;
            *fx = d2dx*(-fh(i-1,j,k) + fh(i+1,j,k));
            *fy = d2dy*(-fh(i,j-1,k) + fh(i,j+1,k));
            *fz = d2dz*(-fh(i,j,k-1) + fh(i,j,k+1));
        }
    } else {
        // GW >= 4: per-direction independent stencil order selection
        // ---------- x direction ----------
        do {
            if constexpr (GW >= 5) {
                if (i+4 <= imax && i-4 >= imin) {
                    const double d840dx = 1.0/840.0/dX;
                    *fx = d840dx*(  3*fh(i-4,j,k) - 32*fh(i-3,j,k) + 168*fh(i-2,j,k) - 672*fh(i-1,j,k)
                                  + 672*fh(i+1,j,k) - 168*fh(i+2,j,k) + 32*fh(i+3,j,k) - 3*fh(i+4,j,k));
                    break;
                }
            }
            if (i+3 <= imax && i-3 >= imin) {
                const double d60dx = 1.0/60.0/dX;
                *fx = d60dx*(-fh(i-3,j,k) + 9*fh(i-2,j,k) - 45*fh(i-1,j,k)
                             + 45*fh(i+1,j,k) - 9*fh(i+2,j,k) + fh(i+3,j,k));
                break;
            }
            if (i+2 <= imax && i-2 >= imin) {
                const double d12dx = 1.0/12.0/dX;
                *fx = d12dx*(fh(i-2,j,k) - 8*fh(i-1,j,k) + 8*fh(i+1,j,k) - fh(i+2,j,k));
                break;
            }
            if (i+1 <= imax && i-1 >= imin) {
                const double d2dx = 0.5/dX;
                *fx = d2dx*(-fh(i-1,j,k) + fh(i+1,j,k));
            }
        } while(false);

        // ---------- y direction ----------
        do {
            if constexpr (GW >= 5) {
                if (j+4 <= jmax && j-4 >= jmin) {
                    const double d840dy = 1.0/840.0/dY;
                    *fy = d840dy*(  3*fh(i,j-4,k) - 32*fh(i,j-3,k) + 168*fh(i,j-2,k) - 672*fh(i,j-1,k)
                                  + 672*fh(i,j+1,k) - 168*fh(i,j+2,k) + 32*fh(i,j+3,k) - 3*fh(i,j+4,k));
                    break;
                }
            }
            if (j+3 <= jmax && j-3 >= jmin) {
                const double d60dy = 1.0/60.0/dY;
                *fy = d60dy*(-fh(i,j-3,k) + 9*fh(i,j-2,k) - 45*fh(i,j-1,k)
                             + 45*fh(i,j+1,k) - 9*fh(i,j+2,k) + fh(i,j+3,k));
                break;
            }
            if (j+2 <= jmax && j-2 >= jmin) {
                const double d12dy = 1.0/12.0/dY;
                *fy = d12dy*(fh(i,j-2,k) - 8*fh(i,j-1,k) + 8*fh(i,j+1,k) - fh(i,j+2,k));
                break;
            }
            if (j+1 <= jmax && j-1 >= jmin) {
                const double d2dy = 0.5/dY;
                *fy = d2dy*(-fh(i,j-1,k) + fh(i,j+1,k));
            }
        } while(false);

        // ---------- z direction ----------
        do {
            if constexpr (GW >= 5) {
                if (k+4 <= kmax && k-4 >= kmin) {
                    const double d840dz = 1.0/840.0/dZ;
                    *fz = d840dz*(  3*fh(i,j,k-4) - 32*fh(i,j,k-3) + 168*fh(i,j,k-2) - 672*fh(i,j,k-1)
                                  + 672*fh(i,j,k+1) - 168*fh(i,j,k+2) + 32*fh(i,j,k+3) - 3*fh(i,j,k+4));
                    break;
                }
            }
            if (k+3 <= kmax && k-3 >= kmin) {
                const double d60dz = 1.0/60.0/dZ;
                *fz = d60dz*(-fh(i,j,k-3) + 9*fh(i,j,k-2) - 45*fh(i,j,k-1)
                             + 45*fh(i,j,k+1) - 9*fh(i,j,k+2) + fh(i,j,k+3));
                break;
            }
            if (k+2 <= kmax && k-2 >= kmin) {
                const double d12dz = 1.0/12.0/dZ;
                *fz = d12dz*(fh(i,j,k-2) - 8*fh(i,j,k-1) + 8*fh(i,j,k+1) - fh(i,j,k+2));
                break;
            }
            if (k+1 <= kmax && k-1 >= kmin) {
                const double d2dz = 0.5/dZ;
                *fz = d2dz*(-fh(i,j,k-1) + fh(i,j,k+1));
            }
        } while(false);
    }

    (void)onoff;
}

// ============================================================
// Multi-GW second derivatives (device function)
// ============================================================
template<int GW>
__device__ __forceinline__ void d_fdderivs_point_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff,
    int i, int j, int k
) {
    constexpr int ORD = GW - 1;

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    const int imax = ex[0] - 1;
    const int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;

    *fxx = 0.0; *fyy = 0.0; *fzz = 0.0;
    *fxy = 0.0; *fxz = 0.0; *fyz = 0.0;

    if (i >= imax || j >= jmax || k >= kmax) return;

    // Match Fortran fderivs: symmetry_bd fills (GW-1) ghost layers below first point.
    // Fortran kmin=-1 is 1-indexed (2 layers below k=1), equivalent to C 0-indexed -(GW-1).
    constexpr int sym_ghost_dd = GW - 1;
    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > 0 && fabs(Z[0]) < dZ) kmin = -sym_ghost_dd;
    if (symmetry > 1 && fabs(X[0]) < dX) imin = -sym_ghost_dd;
    if (symmetry > 1 && fabs(Y[0]) < dY) jmin = -sym_ghost_dd;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b_pad(
            ORD,
            ex[0], ex[1], ex[2],
            pex[0], pex[1], pex[2],
            f, ii, jj, kk,
            SYM1, SYM2, SYM3
        );
    };

    if constexpr (GW == 2) {
        const bool fit_x = (i+1 <= imax && i-1 >= imin);
        const bool fit_y = (j+1 <= jmax && j-1 >= jmin);
        const bool fit_z = (k+1 <= kmax && k-1 >= kmin);

        const double Sdxdx = 1.0/(dX*dX);
        const double Sdydy = 1.0/(dY*dY);
        const double Sdzdz = 1.0/(dZ*dZ);
        const double Sdxdy = 0.25/(dX*dY);
        const double Sdxdz = 0.25/(dX*dZ);
        const double Sdydz = 0.25/(dY*dZ);

        if (fit_x) {
            *fxx = Sdxdx*(fh(i-1,j,k) - 2*fh(i,j,k) + fh(i+1,j,k));
        }
        if (fit_y) {
            *fyy = Sdydy*(fh(i,j-1,k) - 2*fh(i,j,k) + fh(i,j+1,k));
        }
        if (fit_z) {
            *fzz = Sdzdz*(fh(i,j,k-1) - 2*fh(i,j,k) + fh(i,j,k+1));
        }
        if (fit_x && fit_y) {
            *fxy = Sdxdy*(fh(i-1,j-1,k) - fh(i+1,j-1,k) - fh(i-1,j+1,k) + fh(i+1,j+1,k));
        }
        if (fit_x && fit_z) {
            *fxz = Sdxdz*(fh(i-1,j,k-1) - fh(i+1,j,k-1) - fh(i-1,j,k+1) + fh(i+1,j,k+1));
        }
        if (fit_y && fit_z) {
            *fyz = Sdydz*(fh(i,j-1,k-1) - fh(i,j+1,k-1) - fh(i,j-1,k+1) + fh(i,j+1,k+1));
        }
    } else if constexpr (GW == 3) {
        // GW <= 3: isotropic fallback (matches Fortran diff_new.f90 behavior)
        // ---------- 4th order (GW >= 3) ----------
        if (i+2 <= imax && i-2 >= imin &&
            j+2 <= jmax && j-2 >= jmin &&
            k+2 <= kmax && k-2 >= kmin) {
            const double Fdxdx = 1.0/12.0/(dX*dX);
            const double Fdydy = 1.0/12.0/(dY*dY);
            const double Fdzdz = 1.0/12.0/(dZ*dZ);
            const double Fdxdy = 1.0/144.0/(dX*dY);
            const double Fdxdz = 1.0/144.0/(dX*dZ);
            const double Fdydz = 1.0/144.0/(dY*dZ);

            *fxx = Fdxdx*(-fh(i-2,j,k)+16*fh(i-1,j,k)-30*fh(i,j,k)
                          -fh(i+2,j,k)+16*fh(i+1,j,k));
            *fyy = Fdydy*(-fh(i,j-2,k)+16*fh(i,j-1,k)-30*fh(i,j,k)
                          -fh(i,j+2,k)+16*fh(i,j+1,k));
            *fzz = Fdzdz*(-fh(i,j,k-2)+16*fh(i,j,k-1)-30*fh(i,j,k)
                          -fh(i,j,k+2)+16*fh(i,j,k+1));

            *fxy = Fdxdy*(   (fh(i-2,j-2,k)- 8*fh(i-1,j-2,k)+ 8*fh(i+1,j-2,k)-fh(i+2,j-2,k))
                        - 8*(fh(i-2,j-1,k)- 8*fh(i-1,j-1,k)+ 8*fh(i+1,j-1,k)-fh(i+2,j-1,k))
                        + 8*(fh(i-2,j+1,k)- 8*fh(i-1,j+1,k)+ 8*fh(i+1,j+1,k)-fh(i+2,j+1,k))
                        -   (fh(i-2,j+2,k)- 8*fh(i-1,j+2,k)+ 8*fh(i+1,j+2,k)-fh(i+2,j+2,k)));
            *fxz = Fdxdz*(   (fh(i-2,j,k-2)- 8*fh(i-1,j,k-2)+ 8*fh(i+1,j,k-2)-fh(i+2,j,k-2))
                        - 8*(fh(i-2,j,k-1)- 8*fh(i-1,j,k-1)+ 8*fh(i+1,j,k-1)-fh(i+2,j,k-1))
                        + 8*(fh(i-2,j,k+1)- 8*fh(i-1,j,k+1)+ 8*fh(i+1,j,k+1)-fh(i+2,j,k+1))
                        -   (fh(i-2,j,k+2)- 8*fh(i-1,j,k+2)+ 8*fh(i+1,j,k+2)-fh(i+2,j,k+2)));
            *fyz = Fdydz*(   (fh(i,j-2,k-2)- 8*fh(i,j-1,k-2)+ 8*fh(i,j+1,k-2)-fh(i,j+2,k-2))
                        - 8*(fh(i,j-2,k-1)- 8*fh(i,j-1,k-1)+ 8*fh(i,j+1,k-1)-fh(i,j+2,k-1))
                        + 8*(fh(i,j-2,k+1)- 8*fh(i,j-1,k+1)+ 8*fh(i,j+1,k+1)-fh(i,j+2,k+1))
                        -   (fh(i,j-2,k+2)- 8*fh(i,j-1,k+2)+ 8*fh(i,j+1,k+2)-fh(i,j+2,k+2)));
            return;
        }
        // ---------- 2nd order (always) ----------
        if (i+1 <= imax && i-1 >= imin &&
            j+1 <= jmax && j-1 >= jmin &&
            k+1 <= kmax && k-1 >= kmin) {
            const double Sdxdx = 1.0/(dX*dX);
            const double Sdydy = 1.0/(dY*dY);
            const double Sdzdz = 1.0/(dZ*dZ);
            const double Sdxdy = 0.25/(dX*dY);
            const double Sdxdz = 0.25/(dX*dZ);
            const double Sdydz = 0.25/(dY*dZ);

            *fxx = Sdxdx*(fh(i-1,j,k) - 2*fh(i,j,k) + fh(i+1,j,k));
            *fyy = Sdydy*(fh(i,j-1,k) - 2*fh(i,j,k) + fh(i,j+1,k));
            *fzz = Sdzdz*(fh(i,j,k-1) - 2*fh(i,j,k) + fh(i,j,k+1));

            *fxy = Sdxdy*(fh(i-1,j-1,k) - fh(i+1,j-1,k) - fh(i-1,j+1,k) + fh(i+1,j+1,k));
            *fxz = Sdxdz*(fh(i-1,j,k-1) - fh(i+1,j,k-1) - fh(i-1,j,k+1) + fh(i+1,j,k+1));
            *fyz = Sdydz*(fh(i,j-1,k-1) - fh(i,j+1,k-1) - fh(i,j-1,k+1) + fh(i,j+1,k+1));
        }
    } else {
        // GW >= 4: per-direction independent stencil order selection
        // ---------- fxx ----------
        do {
            if constexpr (GW >= 5) {
                if (i+4 <= imax && i-4 >= imin) {
                    const double Edxdx = 1.0/5040.0/(dX*dX);
                    *fxx = Edxdx*(- 9*fh(i-4,j,k)+128*fh(i-3,j,k)-1008*fh(i-2,j,k)+8064*fh(i-1,j,k)-14350*fh(i,j,k)
                                  +8064*fh(i+1,j,k)-1008*fh(i+2,j,k)+128*fh(i+3,j,k)- 9*fh(i+4,j,k));
                    break;
                }
            }
            if (i+3 <= imax && i-3 >= imin) {
                const double Xdxdx = 1.0/180.0/(dX*dX);
                *fxx = Xdxdx*(2*fh(i-3,j,k)-27*fh(i-2,j,k)+270*fh(i-1,j,k)-490*fh(i,j,k)
                              +270*fh(i+1,j,k)-27*fh(i+2,j,k)+2*fh(i+3,j,k));
                break;
            }
            if (i+2 <= imax && i-2 >= imin) {
                const double Fdxdx = 1.0/12.0/(dX*dX);
                *fxx = Fdxdx*(-fh(i-2,j,k)+16*fh(i-1,j,k)-30*fh(i,j,k)
                              -fh(i+2,j,k)+16*fh(i+1,j,k));
                break;
            }
            if (i+1 <= imax && i-1 >= imin) {
                const double Sdxdx = 1.0/(dX*dX);
                *fxx = Sdxdx*(fh(i-1,j,k) - 2*fh(i,j,k) + fh(i+1,j,k));
            }
        } while(false);

        // ---------- fyy ----------
        do {
            if constexpr (GW >= 5) {
                if (j+4 <= jmax && j-4 >= jmin) {
                    const double Edydy = 1.0/5040.0/(dY*dY);
                    *fyy = Edydy*(- 9*fh(i,j-4,k)+128*fh(i,j-3,k)-1008*fh(i,j-2,k)+8064*fh(i,j-1,k)-14350*fh(i,j,k)
                                  +8064*fh(i,j+1,k)-1008*fh(i,j+2,k)+128*fh(i,j+3,k)- 9*fh(i,j+4,k));
                    break;
                }
            }
            if (j+3 <= jmax && j-3 >= jmin) {
                const double Xdydy = 1.0/180.0/(dY*dY);
                *fyy = Xdydy*(2*fh(i,j-3,k)-27*fh(i,j-2,k)+270*fh(i,j-1,k)-490*fh(i,j,k)
                              +270*fh(i,j+1,k)-27*fh(i,j+2,k)+2*fh(i,j+3,k));
                break;
            }
            if (j+2 <= jmax && j-2 >= jmin) {
                const double Fdydy = 1.0/12.0/(dY*dY);
                *fyy = Fdydy*(-fh(i,j-2,k)+16*fh(i,j-1,k)-30*fh(i,j,k)
                              -fh(i,j+2,k)+16*fh(i,j+1,k));
                break;
            }
            if (j+1 <= jmax && j-1 >= jmin) {
                const double Sdydy = 1.0/(dY*dY);
                *fyy = Sdydy*(fh(i,j-1,k) - 2*fh(i,j,k) + fh(i,j+1,k));
            }
        } while(false);

        // ---------- fzz ----------
        do {
            if constexpr (GW >= 5) {
                if (k+4 <= kmax && k-4 >= kmin) {
                    const double Edzdz = 1.0/5040.0/(dZ*dZ);
                    *fzz = Edzdz*(- 9*fh(i,j,k-4)+128*fh(i,j,k-3)-1008*fh(i,j,k-2)+8064*fh(i,j,k-1)-14350*fh(i,j,k)
                                  +8064*fh(i,j,k+1)-1008*fh(i,j,k+2)+128*fh(i,j,k+3)- 9*fh(i,j,k+4));
                    break;
                }
            }
            if (k+3 <= kmax && k-3 >= kmin) {
                const double Xdzdz = 1.0/180.0/(dZ*dZ);
                *fzz = Xdzdz*(2*fh(i,j,k-3)-27*fh(i,j,k-2)+270*fh(i,j,k-1)-490*fh(i,j,k)
                              +270*fh(i,j,k+1)-27*fh(i,j,k+2)+2*fh(i,j,k+3));
                break;
            }
            if (k+2 <= kmax && k-2 >= kmin) {
                const double Fdzdz = 1.0/12.0/(dZ*dZ);
                *fzz = Fdzdz*(-fh(i,j,k-2)+16*fh(i,j,k-1)-30*fh(i,j,k)
                              -fh(i,j,k+2)+16*fh(i,j,k+1));
                break;
            }
            if (k+1 <= kmax && k-1 >= kmin) {
                const double Sdzdz = 1.0/(dZ*dZ);
                *fzz = Sdzdz*(fh(i,j,k-1) - 2*fh(i,j,k) + fh(i,j,k+1));
            }
        } while(false);

        // ---------- fxy (cross: x*y) ----------
        do {
            if constexpr (GW >= 5) {
                if (i+4 <= imax && i-4 >= imin && j+4 <= jmax && j-4 >= jmin) {
                    double c[8]; int o[8]; double fd1_denom_;
                    diff_new_detail::get_fd1_8(c, o, fd1_denom_);
                    const double d2 = fd1_denom_ * fd1_denom_;
                    *fxy = diff_new_detail::cross_deriv<8, 0, 1>(fh, i, j, k, c, o, d2, 1.0/(dX*dY));
                    break;
                }
            }
            if (i+3 <= imax && i-3 >= imin && j+3 <= jmax && j-3 >= jmin) {
                double c[6]; int o[6]; double fd1_denom_;
                diff_new_detail::get_fd1_6(c, o, fd1_denom_);
                const double d2 = fd1_denom_ * fd1_denom_;
                *fxy = diff_new_detail::cross_deriv<6, 0, 1>(fh, i, j, k, c, o, d2, 1.0/(dX*dY));
                break;
            }
            if (i+2 <= imax && i-2 >= imin && j+2 <= jmax && j-2 >= jmin) {
                const double Fdxdy = 1.0/144.0/(dX*dY);
                *fxy = Fdxdy*(   (fh(i-2,j-2,k)- 8*fh(i-1,j-2,k)+ 8*fh(i+1,j-2,k)-fh(i+2,j-2,k))
                            - 8*(fh(i-2,j-1,k)- 8*fh(i-1,j-1,k)+ 8*fh(i+1,j-1,k)-fh(i+2,j-1,k))
                            + 8*(fh(i-2,j+1,k)- 8*fh(i-1,j+1,k)+ 8*fh(i+1,j+1,k)-fh(i+2,j+1,k))
                            -   (fh(i-2,j+2,k)- 8*fh(i-1,j+2,k)+ 8*fh(i+1,j+2,k)-fh(i+2,j+2,k)));
                break;
            }
            if (i+1 <= imax && i-1 >= imin && j+1 <= jmax && j-1 >= jmin) {
                const double Sdxdy = 0.25/(dX*dY);
                *fxy = Sdxdy*(fh(i-1,j-1,k) - fh(i+1,j-1,k) - fh(i-1,j+1,k) + fh(i+1,j+1,k));
            }
        } while(false);

        // ---------- fxz (cross: x*z) ----------
        do {
            if constexpr (GW >= 5) {
                if (i+4 <= imax && i-4 >= imin && k+4 <= kmax && k-4 >= kmin) {
                    double c[8]; int o[8]; double fd1_denom_;
                    diff_new_detail::get_fd1_8(c, o, fd1_denom_);
                    const double d2 = fd1_denom_ * fd1_denom_;
                    *fxz = diff_new_detail::cross_deriv<8, 0, 2>(fh, i, j, k, c, o, d2, 1.0/(dX*dZ));
                    break;
                }
            }
            if (i+3 <= imax && i-3 >= imin && k+3 <= kmax && k-3 >= kmin) {
                double c[6]; int o[6]; double fd1_denom_;
                diff_new_detail::get_fd1_6(c, o, fd1_denom_);
                const double d2 = fd1_denom_ * fd1_denom_;
                *fxz = diff_new_detail::cross_deriv<6, 0, 2>(fh, i, j, k, c, o, d2, 1.0/(dX*dZ));
                break;
            }
            if (i+2 <= imax && i-2 >= imin && k+2 <= kmax && k-2 >= kmin) {
                const double Fdxdz = 1.0/144.0/(dX*dZ);
                *fxz = Fdxdz*(   (fh(i-2,j,k-2)- 8*fh(i-1,j,k-2)+ 8*fh(i+1,j,k-2)-fh(i+2,j,k-2))
                            - 8*(fh(i-2,j,k-1)- 8*fh(i-1,j,k-1)+ 8*fh(i+1,j,k-1)-fh(i+2,j,k-1))
                            + 8*(fh(i-2,j,k+1)- 8*fh(i-1,j,k+1)+ 8*fh(i+1,j,k+1)-fh(i+2,j,k+1))
                            -   (fh(i-2,j,k+2)- 8*fh(i-1,j,k+2)+ 8*fh(i+1,j,k+2)-fh(i+2,j,k+2)));
                break;
            }
            if (i+1 <= imax && i-1 >= imin && k+1 <= kmax && k-1 >= kmin) {
                const double Sdxdz = 0.25/(dX*dZ);
                *fxz = Sdxdz*(fh(i-1,j,k-1) - fh(i+1,j,k-1) - fh(i-1,j,k+1) + fh(i+1,j,k+1));
            }
        } while(false);

        // ---------- fyz (cross: y*z) ----------
        do {
            if constexpr (GW >= 5) {
                if (j+4 <= jmax && j-4 >= jmin && k+4 <= kmax && k-4 >= kmin) {
                    double c[8]; int o[8]; double fd1_denom_;
                    diff_new_detail::get_fd1_8(c, o, fd1_denom_);
                    const double d2 = fd1_denom_ * fd1_denom_;
                    *fyz = diff_new_detail::cross_deriv<8, 1, 2>(fh, i, j, k, c, o, d2, 1.0/(dY*dZ));
                    break;
                }
            }
            if (j+3 <= jmax && j-3 >= jmin && k+3 <= kmax && k-3 >= kmin) {
                double c[6]; int o[6]; double fd1_denom_;
                diff_new_detail::get_fd1_6(c, o, fd1_denom_);
                const double d2 = fd1_denom_ * fd1_denom_;
                *fyz = diff_new_detail::cross_deriv<6, 1, 2>(fh, i, j, k, c, o, d2, 1.0/(dY*dZ));
                break;
            }
            if (j+2 <= jmax && j-2 >= jmin && k+2 <= kmax && k-2 >= kmin) {
                const double Fdydz = 1.0/144.0/(dY*dZ);
                *fyz = Fdydz*(   (fh(i,j-2,k-2)- 8*fh(i,j-1,k-2)+ 8*fh(i,j+1,k-2)-fh(i,j+2,k-2))
                            - 8*(fh(i,j-2,k-1)- 8*fh(i,j-1,k-1)+ 8*fh(i,j+1,k-1)-fh(i,j+2,k-1))
                            + 8*(fh(i,j-2,k+1)- 8*fh(i,j-1,k+1)+ 8*fh(i,j+1,k+1)-fh(i,j+2,k+1))
                            -   (fh(i,j-2,k+2)- 8*fh(i,j-1,k+2)+ 8*fh(i,j+1,k+2)-fh(i,j+2,k+2)));
                break;
            }
            if (j+1 <= jmax && j-1 >= jmin && k+1 <= kmax && k-1 >= kmin) {
                const double Sdydz = 0.25/(dY*dZ);
                *fyz = Sdydz*(fh(i,j-1,k-1) - fh(i,j+1,k-1) - fh(i,j-1,k+1) + fh(i,j+1,k+1));
            }
        } while(false);
    }

    (void)onoff;
}

// ============================================================
// Backward-compatible wrappers (default GW=3)
// ============================================================
__device__ __forceinline__ void d_fderivs_point(
    const int ex[3], const int pex[3], const double* f,
    double* fx, double* fy, double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff,
    int i, int j, int k
) {
    d_fderivs_point_gw<ghost_width>(
        ex, pex, f, fx, fy, fz,
        X, Y, Z, SYM1, SYM2, SYM3,
        symmetry, onoff, i, j, k
    );
}

__device__ __forceinline__ void d_fdderivs_point(
    const int ex[3], const int pex[3], const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff,
    int i, int j, int k
) {
    d_fdderivs_point_gw<ghost_width>(
        ex, pex, f, fxx, fxy, fxz, fyy, fyz, fzz,
        X, Y, Z, SYM1, SYM2, SYM3,
        symmetry, onoff, i, j, k
    );
}

#endif // DIFF_NEW_GPU_CUH
