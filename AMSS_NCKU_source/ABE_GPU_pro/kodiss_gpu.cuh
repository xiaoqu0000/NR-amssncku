#ifndef KODISS_GPU_CUH
#define KODISS_GPU_CUH

#ifndef FMISC_GPU_CUH
#include "fmisc_gpu.cuh"
#endif

#ifndef MACRODEF_FH_INCLUDED
#define MACRODEF_FH_INCLUDED
#include "macrodef.fh"
#endif
#include <cmath>

template<int GW>
__device__ __forceinline__ double d_kodis_point(
    const int ex[3], const int pex[3], const double* f,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps,
    int i, int j, int k
) {
    constexpr int NO_SYMM = 0, EQ_SYMM = 1;
    // KO dissipation uses a (D+D-)^GW stencil with half-width = GW (e.g. GW=3 ⇒ 7-pt
    // stencil reaching f(i±3)). At the equatorial symmetry boundary the Fortran
    // reference allocates fh(-GW+1:ex) and call symmetry_bd(GW, ...) (ord=GW). For
    // bit-exact match we must use ORD = GW here so d_symmetry_bd_0b_pad's range
    // check (`i < -ord`) admits the deepest ghost index −GW.
    constexpr int ORD = GW;
    // cof = 2^(2*GW), sign = (-1)^(GW-1)
    // GW=2: cof=16, sign=-1 | GW=3: cof=64, sign=+1
    // GW=4: cof=256, sign=-1 | GW=5: cof=1024, sign=+1
    constexpr double cof = double(1 << (2 * GW));
    constexpr double sign = ((GW - 1) % 2 == 0) ? 1.0 : -1.0;

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    const int imax = ex[0] - 1;
    const int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;

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

    if (i - ORD < imin || i + ORD > imax ||
        j - ORD < jmin || j + ORD > jmax ||
        k - ORD < kmin || k + ORD > kmax) {
        return 0.0;
    }

    double rhs_add = 0.0;

    // ================================================================
    // GW=2: (D+D-)^2, 5-point stencil, coefficients [1,-4,6,-4,1]
    // ================================================================
    if constexpr (GW == 2) {
        const double inv_cof_dx = sign * eps / (cof * dX);
        const double inv_cof_dy = sign * eps / (cof * dY);
        const double inv_cof_dz = sign * eps / (cof * dZ);

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
        const double inv_cof_dx = sign * eps / (cof * dX);
        const double inv_cof_dy = sign * eps / (cof * dY);
        const double inv_cof_dz = sign * eps / (cof * dZ);

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
        const double inv_cof_dx = sign * eps / (cof * dX);
        const double inv_cof_dy = sign * eps / (cof * dY);
        const double inv_cof_dz = sign * eps / (cof * dZ);

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
        const double inv_cof_dx = sign * eps / (cof * dX);
        const double inv_cof_dy = sign * eps / (cof * dY);
        const double inv_cof_dz = sign * eps / (cof * dZ);

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

#endif // KODISS_GPU_CUH
