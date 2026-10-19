#ifndef KODISS_SH_GPU_CUH
#define KODISS_SH_GPU_CUH

#ifndef FMISC_GPU_CUH
#include "fmisc_gpu.cuh"
#endif

#ifndef DIFF_NEW_SH_GPU_CUH
#include "diff_new_sh_gpu.cuh"  // for d_symmetry_stbd_0b_pad, d_sh_get_SoA
#endif

#ifndef MACRODEF_FH_INCLUDED
#define MACRODEF_FH_INCLUDED
#include "macrodef.fh"
#endif
#include <cmath>

// d_sh_get_SoA and d_symmetry_stbd_0b_pad come from diff_new_sh_gpu.cuh

// ==========================================================================
// Helper: compute boundary limits for shell KO dissipation (0-based).
// For KO, ORD = GW (stencil half-width equals GW).
// z has NO ghost zone: kmin=0, kmax=ex[2].
// ==========================================================================
__device__ __forceinline__ void d_kodis_sh_bounds(
    int GW,
    const int ex[3],
    const double* X, const double* Y,
    int symmetry, int sst,
    double dX, double dY,
    int &imin, int &jmin, int &imax, int &jmax,
    int &kmin, int &kmax)
{
    constexpr int NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2;

    imax = ex[0] - 1;
    jmax = ex[1] - 1;
    kmax = ex[2] - 1;

    imin = 0;
    jmin = 0;
    kmin = 0;

    if (symmetry == OCTANT) {
        if (fabs(X[0]) < dX) imin = -GW;
        if (fabs(Y[0]) < dY) jmin = -GW;
    } else if (symmetry == EQ_SYMM) {
        if ((sst == 2 || sst == 4) && fabs(Y[0]) < dY) jmin = -GW;
        if ((sst == 3 || sst == 5) && fabs(Y[ex[1] - 1]) < dY) jmax = ex[1] - 1 + GW;
    }
}

// ==========================================================================
// d_kodis_sh_point<GW> -- Main KO dissipation (all 3 directions)
//
// Shell coordinate version using d_symmetry_stbd_0b_pad with 2 SoA factors.
// z direction has NO ghost zone (kmin=0, kmax=ex[2]).
// ORD for KO = GW (unlike derivatives which use GW-1).
// Returns the value to be ADDED to f_rhs.
// ==========================================================================
template<int GW>
__device__ __forceinline__ double d_kodis_sh_point(
    const int ex[3], const int pex[3], const double* f,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW;
    // cof = 2^(2*GW), sign = (-1)^(GW-1)
    // GW=2: cof=16, sign=-1 | GW=3: cof=64, sign=+1
    // GW=4: cof=256, sign=-1 | GW=5: cof=1024, sign=+1
    constexpr double cof = double(1 << (2 * GW));
    constexpr double sign = ((GW - 1) % 2 == 0) ? 1.0 : -1.0;

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    int imin, jmin, imax, jmax, kmin, kmax;
    d_kodis_sh_bounds(GW, ex, X, Y, symmetry, sst, dX, dY,
                      imin, jmin, imax, jmax, kmin, kmax);

    // Boundary accessor: shell version with 2 SoA factors, no z ghost zone
    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SoA0, SoA1
        );
    };

    // Return 0.0 where the stencil doesn't fit in any direction
    if (i - GW < imin || i + GW > imax ||
        j - GW < jmin || j + GW > jmax ||
        k - GW < kmin || k + GW > kmax) {
        return 0.0;
    }

    const double inv_cof_dx = sign * eps / (cof * dX);
    const double inv_cof_dy = sign * eps / (cof * dY);
    const double inv_cof_dz = sign * eps / (cof * dZ);

    double rhs_add = 0.0;

    // ================================================================
    // GW=2: (D+D-)^2, 5-point stencil, coefficients [1,-4,6,-4,1]
    // ================================================================
    if constexpr (GW == 2) {
        double sx = (fh(i-2,j,k) + fh(i+2,j,k))
                  - 4.0 * (fh(i-1,j,k) + fh(i+1,j,k))
                  + 6.0 * fh(i,j,k);
        double sy = (fh(i,j-2,k) + fh(i,j+2,k))
                  - 4.0 * (fh(i,j-1,k) + fh(i,j+1,k))
                  + 6.0 * fh(i,j,k);
        double sz = (fh(i,j,k-2) + fh(i,j,k+2))
                  - 4.0 * (fh(i,j,k-1) + fh(i,j,k+1))
                  + 6.0 * fh(i,j,k);

        rhs_add = inv_cof_dx * sx + inv_cof_dy * sy + inv_cof_dz * sz;

    // ================================================================
    // GW=3: (D+D-)^3, 7-point stencil, coefficients [1,-6,15,-20,15,-6,1]
    // ================================================================
    } else if constexpr (GW == 3) {
        double sx = (fh(i-3,j,k) + fh(i+3,j,k))
                  - 6.0 * (fh(i-2,j,k) + fh(i+2,j,k))
                  + 15.0 * (fh(i-1,j,k) + fh(i+1,j,k))
                  - 20.0 * fh(i,j,k);
        double sy = (fh(i,j-3,k) + fh(i,j+3,k))
                  - 6.0 * (fh(i,j-2,k) + fh(i,j+2,k))
                  + 15.0 * (fh(i,j-1,k) + fh(i,j+1,k))
                  - 20.0 * fh(i,j,k);
        double sz = (fh(i,j,k-3) + fh(i,j,k+3))
                  - 6.0 * (fh(i,j,k-2) + fh(i,j,k+2))
                  + 15.0 * (fh(i,j,k-1) + fh(i,j,k+1))
                  - 20.0 * fh(i,j,k);

        rhs_add = inv_cof_dx * sx + inv_cof_dy * sy + inv_cof_dz * sz;

    // ================================================================
    // GW=4: (D+D-)^4, 9-point stencil, coefficients [1,-8,28,-56,70,-56,28,-8,1]
    // ================================================================
    } else if constexpr (GW == 4) {
        double sx = (fh(i-4,j,k) + fh(i+4,j,k))
                  - 8.0 * (fh(i-3,j,k) + fh(i+3,j,k))
                  + 28.0 * (fh(i-2,j,k) + fh(i+2,j,k))
                  - 56.0 * (fh(i-1,j,k) + fh(i+1,j,k))
                  + 70.0 * fh(i,j,k);
        double sy = (fh(i,j-4,k) + fh(i,j+4,k))
                  - 8.0 * (fh(i,j-3,k) + fh(i,j+3,k))
                  + 28.0 * (fh(i,j-2,k) + fh(i,j+2,k))
                  - 56.0 * (fh(i,j-1,k) + fh(i,j+1,k))
                  + 70.0 * fh(i,j,k);
        double sz = (fh(i,j,k-4) + fh(i,j,k+4))
                  - 8.0 * (fh(i,j,k-3) + fh(i,j,k+3))
                  + 28.0 * (fh(i,j,k-2) + fh(i,j,k+2))
                  - 56.0 * (fh(i,j,k-1) + fh(i,j,k+1))
                  + 70.0 * fh(i,j,k);

        rhs_add = inv_cof_dx * sx + inv_cof_dy * sy + inv_cof_dz * sz;

    // ================================================================
    // GW=5: (D+D-)^5, 11-point stencil
    // coefficients [1,-10,45,-120,210,-252,210,-120,45,-10,1]
    // ================================================================
    } else if constexpr (GW == 5) {
        double sx = (fh(i-5,j,k) + fh(i+5,j,k))
                  - 10.0 * (fh(i-4,j,k) + fh(i+4,j,k))
                  + 45.0 * (fh(i-3,j,k) + fh(i+3,j,k))
                  - 120.0 * (fh(i-2,j,k) + fh(i+2,j,k))
                  + 210.0 * (fh(i-1,j,k) + fh(i+1,j,k))
                  - 252.0 * fh(i,j,k);
        double sy = (fh(i,j-5,k) + fh(i,j+5,k))
                  - 10.0 * (fh(i,j-4,k) + fh(i,j+4,k))
                  + 45.0 * (fh(i,j-3,k) + fh(i,j+3,k))
                  - 120.0 * (fh(i,j-2,k) + fh(i,j+2,k))
                  + 210.0 * (fh(i,j-1,k) + fh(i,j+1,k))
                  - 252.0 * fh(i,j,k);
        double sz = (fh(i,j,k-5) + fh(i,j,k+5))
                  - 10.0 * (fh(i,j,k-4) + fh(i,j,k+4))
                  + 45.0 * (fh(i,j,k-3) + fh(i,j,k+3))
                  - 120.0 * (fh(i,j,k-2) + fh(i,j,k+2))
                  + 210.0 * (fh(i,j,k-1) + fh(i,j,k+1))
                  - 252.0 * fh(i,j,k);

        rhs_add = inv_cof_dx * sx + inv_cof_dy * sy + inv_cof_dz * sz;
    }

    return rhs_add;
}

// ==========================================================================
// d_kodis_shor_point<GW> -- Radial-only (z direction) KO dissipation
//
// Same stencil as d_kodis_sh_point but ONLY in the z (radial) direction.
// Used for all GW values.
// ==========================================================================
template<int GW>
__device__ __forceinline__ double d_kodis_shor_point(
    const int ex[3], const int pex[3], const double* f,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW;
    constexpr double cof = double(1 << (2 * GW));
    constexpr double sign = ((GW - 1) % 2 == 0) ? 1.0 : -1.0;

    const double dZ = Z[1] - Z[0];

    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    // For radial-only, we need imin/jmin/kmin for the fh accessor,
    // but the boundary check is only on k.
    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    int imin, jmin, imax, jmax, kmin, kmax;
    d_kodis_sh_bounds(GW, ex, X, Y, symmetry, sst, dX, dY,
                      imin, jmin, imax, jmax, kmin, kmax);

    // z has no ghost zone: kmin=0, kmax=ex[2]
    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SoA0, SoA1
        );
    };

    // Return 0.0 if stencil doesn't fit in z direction
    if (k - GW < kmin || k + GW > kmax) {
        return 0.0;
    }

    const double inv_cof_dz = sign * eps / (cof * dZ);
    double rhs_add = 0.0;

    if constexpr (GW == 2) {
        double sz = (fh(i,j,k-2) + fh(i,j,k+2))
                  - 4.0 * (fh(i,j,k-1) + fh(i,j,k+1))
                  + 6.0 * fh(i,j,k);
        rhs_add = inv_cof_dz * sz;

    } else if constexpr (GW == 3) {
        double sz = (fh(i,j,k-3) + fh(i,j,k+3))
                  - 6.0 * (fh(i,j,k-2) + fh(i,j,k+2))
                  + 15.0 * (fh(i,j,k-1) + fh(i,j,k+1))
                  - 20.0 * fh(i,j,k);
        rhs_add = inv_cof_dz * sz;

    } else if constexpr (GW == 4) {
        double sz = (fh(i,j,k-4) + fh(i,j,k+4))
                  - 8.0 * (fh(i,j,k-3) + fh(i,j,k+3))
                  + 28.0 * (fh(i,j,k-2) + fh(i,j,k+2))
                  - 56.0 * (fh(i,j,k-1) + fh(i,j,k+1))
                  + 70.0 * fh(i,j,k);
        rhs_add = inv_cof_dz * sz;

    } else if constexpr (GW == 5) {
        double sz = (fh(i,j,k-5) + fh(i,j,k+5))
                  - 10.0 * (fh(i,j,k-4) + fh(i,j,k+4))
                  + 45.0 * (fh(i,j,k-3) + fh(i,j,k+3))
                  - 120.0 * (fh(i,j,k-2) + fh(i,j,k+2))
                  + 210.0 * (fh(i,j,k-1) + fh(i,j,k+1))
                  - 252.0 * fh(i,j,k);
        rhs_add = inv_cof_dz * sz;
    }

    return rhs_add;
}

// ==========================================================================
// d_kodis_shcr_point<GW> -- Constant-r sphere (x,y only) KO dissipation
//
// Only x and y directions. Only exists for GW>=3.
// ==========================================================================
template<int GW>
__device__ __forceinline__ double d_kodis_shcr_point(
    const int ex[3], const int pex[3], const double* f,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst,
    int i, int j, int k)
{
    if constexpr (GW < 3) {
        return 0.0;
    } else {
        constexpr int ORD = GW;
        constexpr double cof = double(1 << (2 * GW));
        constexpr double sign = ((GW - 1) % 2 == 0) ? 1.0 : -1.0;

        const double dX = X[1] - X[0];
        const double dY = Y[1] - Y[0];

        double SoA0, SoA1;
        d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

        int imin, jmin, imax, jmax, kmin, kmax;
        d_kodis_sh_bounds(GW, ex, X, Y, symmetry, sst, dX, dY,
                          imin, jmin, imax, jmax, kmin, kmax);

        const auto fh = [&](int ii, int jj, int kk) -> double {
            return d_symmetry_stbd_0b_pad(
                ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
                f, ii, jj, kk, SoA0, SoA1
            );
        };

        // Return 0.0 if stencil doesn't fit in x or y direction
        if (i - GW < imin || i + GW > imax ||
            j - GW < jmin || j + GW > jmax) {
            return 0.0;
        }

        const double inv_cof_dx = sign * eps / (cof * dX);
        const double inv_cof_dy = sign * eps / (cof * dY);
        double rhs_add = 0.0;

        if constexpr (GW == 3) {
            double sx = (fh(i-3,j,k) + fh(i+3,j,k))
                      - 6.0 * (fh(i-2,j,k) + fh(i+2,j,k))
                      + 15.0 * (fh(i-1,j,k) + fh(i+1,j,k))
                      - 20.0 * fh(i,j,k);
            double sy = (fh(i,j-3,k) + fh(i,j+3,k))
                      - 6.0 * (fh(i,j-2,k) + fh(i,j+2,k))
                      + 15.0 * (fh(i,j-1,k) + fh(i,j+1,k))
                      - 20.0 * fh(i,j,k);
            rhs_add = inv_cof_dx * sx + inv_cof_dy * sy;

        } else if constexpr (GW == 4) {
            double sx = (fh(i-4,j,k) + fh(i+4,j,k))
                      - 8.0 * (fh(i-3,j,k) + fh(i+3,j,k))
                      + 28.0 * (fh(i-2,j,k) + fh(i+2,j,k))
                      - 56.0 * (fh(i-1,j,k) + fh(i+1,j,k))
                      + 70.0 * fh(i,j,k);
            double sy = (fh(i,j-4,k) + fh(i,j+4,k))
                      - 8.0 * (fh(i,j-3,k) + fh(i,j+3,k))
                      + 28.0 * (fh(i,j-2,k) + fh(i,j+2,k))
                      - 56.0 * (fh(i,j-1,k) + fh(i,j+1,k))
                      + 70.0 * fh(i,j,k);
            rhs_add = inv_cof_dx * sx + inv_cof_dy * sy;

        } else if constexpr (GW == 5) {
            double sx = (fh(i-5,j,k) + fh(i+5,j,k))
                      - 10.0 * (fh(i-4,j,k) + fh(i+4,j,k))
                      + 45.0 * (fh(i-3,j,k) + fh(i+3,j,k))
                      - 120.0 * (fh(i-2,j,k) + fh(i+2,j,k))
                      + 210.0 * (fh(i-1,j,k) + fh(i+1,j,k))
                      - 252.0 * fh(i,j,k);
            double sy = (fh(i,j-5,k) + fh(i,j+5,k))
                      - 10.0 * (fh(i,j-4,k) + fh(i,j+4,k))
                      + 45.0 * (fh(i,j-3,k) + fh(i,j+3,k))
                      - 120.0 * (fh(i,j-2,k) + fh(i,j+2,k))
                      + 210.0 * (fh(i,j-1,k) + fh(i,j+1,k))
                      - 252.0 * fh(i,j,k);
            rhs_add = inv_cof_dx * sx + inv_cof_dy * sy;
        }

        return rhs_add;
    }
}

// ==========================================================================
// d_kodis_sh_new_point<GW> -- Near-boundary tangential KO (GW=2 only)
//
// Per-direction independent boundary check. Each direction is checked
// separately, unlike d_kodis_sh_point which requires all directions to fit.
// Only exists for GW=2.
// ==========================================================================
template<int GW>
__device__ __forceinline__ double d_kodis_sh_new_point(
    const int ex[3], const int pex[3], const double* f,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps, int sst,
    int i, int j, int k)
{
    if constexpr (GW != 2) {
        return 0.0;
    } else {
        constexpr int ORD = GW;
        constexpr double cof = 16.0;  // 2^4 for GW=2
        constexpr double sign = -1.0; // GW=2: sign=-1

        const double dX = X[1] - X[0];
        const double dY = Y[1] - Y[0];
        const double dZ = Z[1] - Z[0];

        double SoA0, SoA1;
        d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

        int imin, jmin, imax, jmax, kmin, kmax;
        d_kodis_sh_bounds(GW, ex, X, Y, symmetry, sst, dX, dY,
                          imin, jmin, imax, jmax, kmin, kmax);

        const auto fh = [&](int ii, int jj, int kk) -> double {
            return d_symmetry_stbd_0b_pad(
                ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
                f, ii, jj, kk, SoA0, SoA1
            );
        };

        const double inv_cof_dx = sign * eps / (cof * dX);
        const double inv_cof_dy = sign * eps / (cof * dY);
        const double inv_cof_dz = sign * eps / (cof * dZ);

        double rhs_add = 0.0;

        // x direction: independent boundary check
        if (i - GW >= imin && i + GW <= imax) {
            double sx = (fh(i-2,j,k) + fh(i+2,j,k))
                      - 4.0 * (fh(i-1,j,k) + fh(i+1,j,k))
                      + 6.0 * fh(i,j,k);
            rhs_add += inv_cof_dx * sx;
        }

        // y direction: independent boundary check
        if (j - GW >= jmin && j + GW <= jmax) {
            double sy = (fh(i,j-2,k) + fh(i,j+2,k))
                      - 4.0 * (fh(i,j-1,k) + fh(i,j+1,k))
                      + 6.0 * fh(i,j,k);
            rhs_add += inv_cof_dy * sy;
        }

        // z direction: independent boundary check
        if (k - GW >= kmin && k + GW <= kmax) {
            double sz = (fh(i,j,k-2) + fh(i,j,k+2))
                      - 4.0 * (fh(i,j,k-1) + fh(i,j,k+1))
                      + 6.0 * fh(i,j,k);
            rhs_add += inv_cof_dz * sz;
        }

        return rhs_add;
    }
}

#endif // KODISS_SH_GPU_CUH
