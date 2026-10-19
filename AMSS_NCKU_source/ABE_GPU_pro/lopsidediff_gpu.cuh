#ifndef LOPSIDEDIFF_GPU_CUH
#define LOPSIDEDIFF_GPU_CUH

#ifndef FMISC_GPU_CUH
#include "fmisc_gpu.cuh"
#endif

#ifndef MACRODEF_FH_INCLUDED
#define MACRODEF_FH_INCLUDED
#include "macrodef.fh"
#endif
#include <cmath>

#ifndef IDX3D_PAD
#define IDX3D_PAD(i, j, k, pnx, pny) ((i) + (pnx) * ((j) + (pny) * (k)))
#endif

template<int GW>
__device__ __forceinline__ double d_lopsided_point(
    const int ex[3], const int pex[3], const double* f,
    const double* Sfx, const double* Sfy, const double* Sfz,
    const double* X, const double* Y, const double* Z,
    int symmetry, double SYM1, double SYM2, double SYM3,
    int i, int j, int k
) {
    constexpr int NO_SYMM = 0, EQ_SYMM = 1;
    // Lop-sided upwind advection at order 2*GW uses fh(i-GW..i+GW) at the
    // interior-most stencil; the Fortran reference allocates fh(-GW+1:ex) and
    // call symmetry_bd(GW, ...). Use ORD = GW so d_symmetry_bd_0b_pad's bound
    // check (`i < -ord`) admits the deepest ghost index −GW.
    constexpr int ORD = GW;

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    const int imax = ex[0] - 1;
    const int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;

    if (i >= imax || j >= jmax || k >= kmax) return 0.0;

    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -ORD;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -ORD;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -ORD;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SYM1, SYM2, SYM3
        );
    };

    const int idx = IDX3D_PAD(i, j, k, pex[0], pex[1]);
    const double vx = Sfx[idx];
    const double vy = Sfy[idx];
    const double vz = Sfz[idx];

    double rhs_add = 0.0;

    // ================================================================
    // GW=2: 2nd-order upwind advection
    // ================================================================
    if constexpr (GW == 2) {
        const double d2dx = 0.5 / dX;
        const double d2dy = 0.5 / dY;
        const double d2dz = 0.5 / dZ;

        // --- X direction ---
        if (vx > 0.0) {
            if (i + 2 <= imax && i >= imin) {
                rhs_add += vx * d2dx * (-3.0*fh(i,j,k) + 4.0*fh(i+1,j,k) - fh(i+2,j,k));
            } else if (i + 1 <= imax && i >= imin) {
                rhs_add += vx * d2dx * (-fh(i,j,k) + fh(i+1,j,k));
            }
        } else if (vx < 0.0) {
            if (i - 2 >= imin && i <= imax) {
                rhs_add -= vx * d2dx * (-3.0*fh(i,j,k) + 4.0*fh(i-1,j,k) - fh(i-2,j,k));
            } else if (i - 1 >= imin && i <= imax) {
                rhs_add -= vx * d2dx * (-fh(i,j,k) + fh(i-1,j,k));
            }
        }

        // --- Y direction ---
        if (vy > 0.0) {
            if (j + 2 <= jmax && j >= jmin) {
                rhs_add += vy * d2dy * (-3.0*fh(i,j,k) + 4.0*fh(i,j+1,k) - fh(i,j+2,k));
            } else if (j + 1 <= jmax && j >= jmin) {
                rhs_add += vy * d2dy * (-fh(i,j,k) + fh(i,j+1,k));
            }
        } else if (vy < 0.0) {
            if (j - 2 >= jmin && j <= jmax) {
                rhs_add -= vy * d2dy * (-3.0*fh(i,j,k) + 4.0*fh(i,j-1,k) - fh(i,j-2,k));
            } else if (j - 1 >= jmin && j <= jmax) {
                rhs_add -= vy * d2dy * (-fh(i,j,k) + fh(i,j-1,k));
            }
        }

        // --- Z direction ---
        if (vz > 0.0) {
            if (k + 2 <= kmax && k >= kmin) {
                rhs_add += vz * d2dz * (-3.0*fh(i,j,k) + 4.0*fh(i,j,k+1) - fh(i,j,k+2));
            } else if (k + 1 <= kmax && k >= kmin) {
                rhs_add += vz * d2dz * (-fh(i,j,k) + fh(i,j,k+1));
            }
        } else if (vz < 0.0) {
            if (k - 2 >= kmin && k <= kmax) {
                rhs_add -= vz * d2dz * (-3.0*fh(i,j,k) + 4.0*fh(i,j,k-1) - fh(i,j,k-2));
            } else if (k - 1 >= kmin && k <= kmax) {
                rhs_add -= vz * d2dz * (-fh(i,j,k) + fh(i,j,k-1));
            }
        }

    // ================================================================
    // GW=3: 4th-order upwind advection
    // ================================================================
    } else if constexpr (GW == 3) {
        const double d12dx = 1.0 / 12.0 / dX;
        const double d12dy = 1.0 / 12.0 / dY;
        const double d12dz = 1.0 / 12.0 / dZ;

        // --- X direction ---
        if (vx > 0.0) {
            if (i + 3 <= imax) {
                rhs_add += vx * d12dx * (-3.0*fh(i-1,j,k) - 10.0*fh(i,j,k) + 18.0*fh(i+1,j,k)
                                         - 6.0*fh(i+2,j,k) + fh(i+3,j,k));
            } else if (i + 2 <= imax) {
                rhs_add += vx * d12dx * (fh(i-2,j,k) - 8.0*fh(i-1,j,k) + 8.0*fh(i+1,j,k) - fh(i+2,j,k));
            } else if (i + 1 <= imax) {
                rhs_add -= vx * d12dx * (-3.0*fh(i+1,j,k) - 10.0*fh(i,j,k) + 18.0*fh(i-1,j,k)
                                         - 6.0*fh(i-2,j,k) + fh(i-3,j,k));
            }
        } else if (vx < 0.0) {
            if (i - 3 >= imin) {
                rhs_add -= vx * d12dx * (-3.0*fh(i+1,j,k) - 10.0*fh(i,j,k) + 18.0*fh(i-1,j,k)
                                         - 6.0*fh(i-2,j,k) + fh(i-3,j,k));
            } else if (i - 2 >= imin) {
                rhs_add += vx * d12dx * (fh(i-2,j,k) - 8.0*fh(i-1,j,k) + 8.0*fh(i+1,j,k) - fh(i+2,j,k));
            } else if (i - 1 >= imin) {
                rhs_add += vx * d12dx * (-3.0*fh(i-1,j,k) - 10.0*fh(i,j,k) + 18.0*fh(i+1,j,k)
                                         - 6.0*fh(i+2,j,k) + fh(i+3,j,k));
            }
        }

        // --- Y direction ---
        if (vy > 0.0) {
            if (j + 3 <= jmax) {
                rhs_add += vy * d12dy * (-3.0*fh(i,j-1,k) - 10.0*fh(i,j,k) + 18.0*fh(i,j+1,k)
                                         - 6.0*fh(i,j+2,k) + fh(i,j+3,k));
            } else if (j + 2 <= jmax) {
                rhs_add += vy * d12dy * (fh(i,j-2,k) - 8.0*fh(i,j-1,k) + 8.0*fh(i,j+1,k) - fh(i,j+2,k));
            } else if (j + 1 <= jmax) {
                rhs_add -= vy * d12dy * (-3.0*fh(i,j+1,k) - 10.0*fh(i,j,k) + 18.0*fh(i,j-1,k)
                                         - 6.0*fh(i,j-2,k) + fh(i,j-3,k));
            }
        } else if (vy < 0.0) {
            if (j - 3 >= jmin) {
                rhs_add -= vy * d12dy * (-3.0*fh(i,j+1,k) - 10.0*fh(i,j,k) + 18.0*fh(i,j-1,k)
                                         - 6.0*fh(i,j-2,k) + fh(i,j-3,k));
            } else if (j - 2 >= jmin) {
                rhs_add += vy * d12dy * (fh(i,j-2,k) - 8.0*fh(i,j-1,k) + 8.0*fh(i,j+1,k) - fh(i,j+2,k));
            } else if (j - 1 >= jmin) {
                rhs_add += vy * d12dy * (-3.0*fh(i,j-1,k) - 10.0*fh(i,j,k) + 18.0*fh(i,j+1,k)
                                         - 6.0*fh(i,j+2,k) + fh(i,j+3,k));
            }
        }

        // --- Z direction ---
        if (vz > 0.0) {
            if (k + 3 <= kmax) {
                rhs_add += vz * d12dz * (-3.0*fh(i,j,k-1) - 10.0*fh(i,j,k) + 18.0*fh(i,j,k+1)
                                         - 6.0*fh(i,j,k+2) + fh(i,j,k+3));
            } else if (k + 2 <= kmax) {
                rhs_add += vz * d12dz * (fh(i,j,k-2) - 8.0*fh(i,j,k-1) + 8.0*fh(i,j,k+1) - fh(i,j,k+2));
            } else if (k + 1 <= kmax) {
                rhs_add -= vz * d12dz * (-3.0*fh(i,j,k+1) - 10.0*fh(i,j,k) + 18.0*fh(i,j,k-1)
                                         - 6.0*fh(i,j,k-2) + fh(i,j,k-3));
            }
        } else if (vz < 0.0) {
            if (k - 3 >= kmin) {
                rhs_add -= vz * d12dz * (-3.0*fh(i,j,k+1) - 10.0*fh(i,j,k) + 18.0*fh(i,j,k-1)
                                         - 6.0*fh(i,j,k-2) + fh(i,j,k-3));
            } else if (k - 2 >= kmin) {
                rhs_add += vz * d12dz * (fh(i,j,k-2) - 8.0*fh(i,j,k-1) + 8.0*fh(i,j,k+1) - fh(i,j,k+2));
            } else if (k - 1 >= kmin) {
                rhs_add += vz * d12dz * (-3.0*fh(i,j,k-1) - 10.0*fh(i,j,k) + 18.0*fh(i,j,k+1)
                                         - 6.0*fh(i,j,k+2) + fh(i,j,k+3));
            }
        }

    // ================================================================
    // GW=4: 6th-order upwind advection
    // ================================================================
    } else if constexpr (GW == 4) {
        const double d60dx = 1.0 / 60.0 / dX;
        const double d60dy = 1.0 / 60.0 / dY;
        const double d60dz = 1.0 / 60.0 / dZ;
        const double d12dx = 1.0 / 12.0 / dX;
        const double d12dy = 1.0 / 12.0 / dY;
        const double d12dz = 1.0 / 12.0 / dZ;
        const double d2dx = 0.5 / dX;
        const double d2dy = 0.5 / dY;
        const double d2dz = 0.5 / dZ;

        // --- X direction ---
        if (vx >= 0.0 && i + 4 <= imax && i - 2 >= imin) {
            rhs_add += vx * d60dx * (2.0*fh(i-2,j,k) - 24.0*fh(i-1,j,k) - 35.0*fh(i,j,k)
                                     + 80.0*fh(i+1,j,k) - 30.0*fh(i+2,j,k) + 8.0*fh(i+3,j,k) - fh(i+4,j,k));
        } else if (vx >= 0.0 && i + 5 <= imax && i - 1 >= imin) {
            rhs_add += vx * d60dx * (-10.0*fh(i-1,j,k) - 77.0*fh(i,j,k) + 150.0*fh(i+1,j,k)
                                     - 100.0*fh(i+2,j,k) + 50.0*fh(i+3,j,k) - 15.0*fh(i+4,j,k) + 2.0*fh(i+5,j,k));
        } else if (vx <= 0.0 && i - 4 >= imin && i + 2 <= imax) {
            rhs_add -= vx * d60dx * (2.0*fh(i+2,j,k) - 24.0*fh(i+1,j,k) - 35.0*fh(i,j,k)
                                     + 80.0*fh(i-1,j,k) - 30.0*fh(i-2,j,k) + 8.0*fh(i-3,j,k) - fh(i-4,j,k));
        } else if (vx <= 0.0 && i - 5 >= imin && i + 1 <= imax) {
            rhs_add -= vx * d60dx * (-10.0*fh(i+1,j,k) - 77.0*fh(i,j,k) + 150.0*fh(i-1,j,k)
                                     - 100.0*fh(i-2,j,k) + 50.0*fh(i-3,j,k) - 15.0*fh(i-4,j,k) + 2.0*fh(i-5,j,k));
        } else if (i + 3 <= imax && i - 3 >= imin) {
            rhs_add += vx * d60dx * (-fh(i-3,j,k) + 9.0*fh(i-2,j,k) - 45.0*fh(i-1,j,k)
                                     + 45.0*fh(i+1,j,k) - 9.0*fh(i+2,j,k) + fh(i+3,j,k));
        } else if (i + 2 <= imax && i - 2 >= imin) {
            rhs_add += vx * d12dx * (fh(i-2,j,k) - 8.0*fh(i-1,j,k) + 8.0*fh(i+1,j,k) - fh(i+2,j,k));
        } else if (i + 1 <= imax && i - 1 >= imin) {
            rhs_add += vx * d2dx * (-fh(i-1,j,k) + fh(i+1,j,k));
        }

        // --- Y direction ---
        if (vy >= 0.0 && j + 4 <= jmax && j - 2 >= jmin) {
            rhs_add += vy * d60dy * (2.0*fh(i,j-2,k) - 24.0*fh(i,j-1,k) - 35.0*fh(i,j,k)
                                     + 80.0*fh(i,j+1,k) - 30.0*fh(i,j+2,k) + 8.0*fh(i,j+3,k) - fh(i,j+4,k));
        } else if (vy >= 0.0 && j + 5 <= jmax && j - 1 >= jmin) {
            rhs_add += vy * d60dy * (-10.0*fh(i,j-1,k) - 77.0*fh(i,j,k) + 150.0*fh(i,j+1,k)
                                     - 100.0*fh(i,j+2,k) + 50.0*fh(i,j+3,k) - 15.0*fh(i,j+4,k) + 2.0*fh(i,j+5,k));
        } else if (vy <= 0.0 && j - 4 >= jmin && j + 2 <= jmax) {
            rhs_add -= vy * d60dy * (2.0*fh(i,j+2,k) - 24.0*fh(i,j+1,k) - 35.0*fh(i,j,k)
                                     + 80.0*fh(i,j-1,k) - 30.0*fh(i,j-2,k) + 8.0*fh(i,j-3,k) - fh(i,j-4,k));
        } else if (vy <= 0.0 && j - 5 >= jmin && j + 1 <= jmax) {
            rhs_add -= vy * d60dy * (-10.0*fh(i,j+1,k) - 77.0*fh(i,j,k) + 150.0*fh(i,j-1,k)
                                     - 100.0*fh(i,j-2,k) + 50.0*fh(i,j-3,k) - 15.0*fh(i,j-4,k) + 2.0*fh(i,j-5,k));
        } else if (j + 3 <= jmax && j - 3 >= jmin) {
            rhs_add += vy * d60dy * (-fh(i,j-3,k) + 9.0*fh(i,j-2,k) - 45.0*fh(i,j-1,k)
                                     + 45.0*fh(i,j+1,k) - 9.0*fh(i,j+2,k) + fh(i,j+3,k));
        } else if (j + 2 <= jmax && j - 2 >= jmin) {
            rhs_add += vy * d12dy * (fh(i,j-2,k) - 8.0*fh(i,j-1,k) + 8.0*fh(i,j+1,k) - fh(i,j+2,k));
        } else if (j + 1 <= jmax && j - 1 >= jmin) {
            rhs_add += vy * d2dy * (-fh(i,j-1,k) + fh(i,j+1,k));
        }

        // --- Z direction ---
        if (vz >= 0.0 && k + 4 <= kmax && k - 2 >= kmin) {
            rhs_add += vz * d60dz * (2.0*fh(i,j,k-2) - 24.0*fh(i,j,k-1) - 35.0*fh(i,j,k)
                                     + 80.0*fh(i,j,k+1) - 30.0*fh(i,j,k+2) + 8.0*fh(i,j,k+3) - fh(i,j,k+4));
        } else if (vz >= 0.0 && k + 5 <= kmax && k - 1 >= kmin) {
            rhs_add += vz * d60dz * (-10.0*fh(i,j,k-1) - 77.0*fh(i,j,k) + 150.0*fh(i,j,k+1)
                                     - 100.0*fh(i,j,k+2) + 50.0*fh(i,j,k+3) - 15.0*fh(i,j,k+4) + 2.0*fh(i,j,k+5));
        } else if (vz <= 0.0 && k - 4 >= kmin && k + 2 <= kmax) {
            rhs_add -= vz * d60dz * (2.0*fh(i,j,k+2) - 24.0*fh(i,j,k+1) - 35.0*fh(i,j,k)
                                     + 80.0*fh(i,j,k-1) - 30.0*fh(i,j,k-2) + 8.0*fh(i,j,k-3) - fh(i,j,k-4));
        } else if (vz <= 0.0 && k - 5 >= kmin && k + 1 <= kmax) {
            rhs_add -= vz * d60dz * (-10.0*fh(i,j,k+1) - 77.0*fh(i,j,k) + 150.0*fh(i,j,k-1)
                                     - 100.0*fh(i,j,k-2) + 50.0*fh(i,j,k-3) - 15.0*fh(i,j,k-4) + 2.0*fh(i,j,k-5));
        } else if (k + 3 <= kmax && k - 3 >= kmin) {
            rhs_add += vz * d60dz * (-fh(i,j,k-3) + 9.0*fh(i,j,k-2) - 45.0*fh(i,j,k-1)
                                     + 45.0*fh(i,j,k+1) - 9.0*fh(i,j,k+2) + fh(i,j,k+3));
        } else if (k + 2 <= kmax && k - 2 >= kmin) {
            rhs_add += vz * d12dz * (fh(i,j,k-2) - 8.0*fh(i,j,k-1) + 8.0*fh(i,j,k+1) - fh(i,j,k+2));
        } else if (k + 1 <= kmax && k - 1 >= kmin) {
            rhs_add += vz * d2dz * (-fh(i,j,k-1) + fh(i,j,k+1));
        }

    // ================================================================
    // GW=5: 8th-order upwind advection
    // ================================================================
    } else if constexpr (GW == 5) {
        const double d840dx = 1.0 / 840.0 / dX;
        const double d840dy = 1.0 / 840.0 / dY;
        const double d840dz = 1.0 / 840.0 / dZ;
        const double d60dx = 1.0 / 60.0 / dX;
        const double d60dy = 1.0 / 60.0 / dY;
        const double d60dz = 1.0 / 60.0 / dZ;
        const double d12dx = 1.0 / 12.0 / dX;
        const double d12dy = 1.0 / 12.0 / dY;
        const double d12dz = 1.0 / 12.0 / dZ;
        const double d2dx = 0.5 / dX;
        const double d2dy = 0.5 / dY;
        const double d2dz = 0.5 / dZ;

        // --- X direction ---
        if (vx >= 0.0 && i + 5 <= imax && i - 3 >= imin) {
            rhs_add += vx * d840dx * (-5.0*fh(i-3,j,k) + 60.0*fh(i-2,j,k) - 420.0*fh(i-1,j,k) - 378.0*fh(i,j,k)
                                      + 1050.0*fh(i+1,j,k) - 420.0*fh(i+2,j,k) + 140.0*fh(i+3,j,k)
                                      - 30.0*fh(i+4,j,k) + 3.0*fh(i+5,j,k));
        } else if (vx <= 0.0 && i - 5 >= imin && i + 3 <= imax) {
            rhs_add -= vx * d840dx * (-5.0*fh(i+3,j,k) + 60.0*fh(i+2,j,k) - 420.0*fh(i+1,j,k) - 378.0*fh(i,j,k)
                                      + 1050.0*fh(i-1,j,k) - 420.0*fh(i-2,j,k) + 140.0*fh(i-3,j,k)
                                      - 30.0*fh(i-4,j,k) + 3.0*fh(i-5,j,k));
        } else if (i + 4 <= imax && i - 4 >= imin) {
            rhs_add += vx * d840dx * (3.0*fh(i-4,j,k) - 32.0*fh(i-3,j,k) + 168.0*fh(i-2,j,k) - 672.0*fh(i-1,j,k)
                                      + 672.0*fh(i+1,j,k) - 168.0*fh(i+2,j,k) + 32.0*fh(i+3,j,k) - 3.0*fh(i+4,j,k));
        } else if (i + 3 <= imax && i - 3 >= imin) {
            rhs_add += vx * d60dx * (-fh(i-3,j,k) + 9.0*fh(i-2,j,k) - 45.0*fh(i-1,j,k)
                                     + 45.0*fh(i+1,j,k) - 9.0*fh(i+2,j,k) + fh(i+3,j,k));
        } else if (i + 2 <= imax && i - 2 >= imin) {
            rhs_add += vx * d12dx * (fh(i-2,j,k) - 8.0*fh(i-1,j,k) + 8.0*fh(i+1,j,k) - fh(i+2,j,k));
        } else if (i + 1 <= imax && i - 1 >= imin) {
            rhs_add += vx * d2dx * (-fh(i-1,j,k) + fh(i+1,j,k));
        }

        // --- Y direction ---
        if (vy >= 0.0 && j + 5 <= jmax && j - 3 >= jmin) {
            rhs_add += vy * d840dy * (-5.0*fh(i,j-3,k) + 60.0*fh(i,j-2,k) - 420.0*fh(i,j-1,k) - 378.0*fh(i,j,k)
                                      + 1050.0*fh(i,j+1,k) - 420.0*fh(i,j+2,k) + 140.0*fh(i,j+3,k)
                                      - 30.0*fh(i,j+4,k) + 3.0*fh(i,j+5,k));
        } else if (vy <= 0.0 && j - 5 >= jmin && j + 3 <= jmax) {
            rhs_add -= vy * d840dy * (-5.0*fh(i,j+3,k) + 60.0*fh(i,j+2,k) - 420.0*fh(i,j+1,k) - 378.0*fh(i,j,k)
                                      + 1050.0*fh(i,j-1,k) - 420.0*fh(i,j-2,k) + 140.0*fh(i,j-3,k)
                                      - 30.0*fh(i,j-4,k) + 3.0*fh(i,j-5,k));
        } else if (j + 4 <= jmax && j - 4 >= jmin) {
            rhs_add += vy * d840dy * (3.0*fh(i,j-4,k) - 32.0*fh(i,j-3,k) + 168.0*fh(i,j-2,k) - 672.0*fh(i,j-1,k)
                                      + 672.0*fh(i,j+1,k) - 168.0*fh(i,j+2,k) + 32.0*fh(i,j+3,k) - 3.0*fh(i,j+4,k));
        } else if (j + 3 <= jmax && j - 3 >= jmin) {
            rhs_add += vy * d60dy * (-fh(i,j-3,k) + 9.0*fh(i,j-2,k) - 45.0*fh(i,j-1,k)
                                     + 45.0*fh(i,j+1,k) - 9.0*fh(i,j+2,k) + fh(i,j+3,k));
        } else if (j + 2 <= jmax && j - 2 >= jmin) {
            rhs_add += vy * d12dy * (fh(i,j-2,k) - 8.0*fh(i,j-1,k) + 8.0*fh(i,j+1,k) - fh(i,j+2,k));
        } else if (j + 1 <= jmax && j - 1 >= jmin) {
            rhs_add += vy * d2dy * (-fh(i,j-1,k) + fh(i,j+1,k));
        }

        // --- Z direction ---
        if (vz >= 0.0 && k + 5 <= kmax && k - 3 >= kmin) {
            rhs_add += vz * d840dz * (-5.0*fh(i,j,k-3) + 60.0*fh(i,j,k-2) - 420.0*fh(i,j,k-1) - 378.0*fh(i,j,k)
                                      + 1050.0*fh(i,j,k+1) - 420.0*fh(i,j,k+2) + 140.0*fh(i,j,k+3)
                                      - 30.0*fh(i,j,k+4) + 3.0*fh(i,j,k+5));
        } else if (vz <= 0.0 && k - 5 >= kmin && k + 3 <= kmax) {
            rhs_add -= vz * d840dz * (-5.0*fh(i,j,k+3) + 60.0*fh(i,j,k+2) - 420.0*fh(i,j,k+1) - 378.0*fh(i,j,k)
                                      + 1050.0*fh(i,j,k-1) - 420.0*fh(i,j,k-2) + 140.0*fh(i,j,k-3)
                                      - 30.0*fh(i,j,k-4) + 3.0*fh(i,j,k-5));
        } else if (k + 4 <= kmax && k - 4 >= kmin) {
            rhs_add += vz * d840dz * (3.0*fh(i,j,k-4) - 32.0*fh(i,j,k-3) + 168.0*fh(i,j,k-2) - 672.0*fh(i,j,k-1)
                                      + 672.0*fh(i,j,k+1) - 168.0*fh(i,j,k+2) + 32.0*fh(i,j,k+3) - 3.0*fh(i,j,k+4));
        } else if (k + 3 <= kmax && k - 3 >= kmin) {
            rhs_add += vz * d60dz * (-fh(i,j,k-3) + 9.0*fh(i,j,k-2) - 45.0*fh(i,j,k-1)
                                     + 45.0*fh(i,j,k+1) - 9.0*fh(i,j,k+2) + fh(i,j,k+3));
        } else if (k + 2 <= kmax && k - 2 >= kmin) {
            rhs_add += vz * d12dz * (fh(i,j,k-2) - 8.0*fh(i,j,k-1) + 8.0*fh(i,j,k+1) - fh(i,j,k+2));
        } else if (k + 1 <= kmax && k - 1 >= kmin) {
            rhs_add += vz * d2dz * (-fh(i,j,k-1) + fh(i,j,k+1));
        }
    }

    return rhs_add;
}

#endif // LOPSIDEDIFF_GPU_CUH
