#ifndef KODISS_GPU_CUH
#define KODISS_GPU_CUH

#include "kodiss.h"

#include "fmisc.h"
#include "fmisc_gpu.cuh"

#include "macrodef.fh"
#include <cmath>

__device__ __forceinline__ double d_kodis_point(
    const int ex[3], const double* f,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, double eps,
    int i, int j, int k // 0-based
) {
    const double SIX = 6.0;
    const int NO_SYMM = 0, EQ_SYMM = 1;

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    const int imax = ex[0] - 1;
    const int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;

    double SoA[3] = {SYM1, SYM2, SYM3};

    double rhs_add = 0.0;

#if ghost_width == 2
    const double FOUR = 4.0;

    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -2;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -2;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -2;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b(2, ex[0], ex[1], ex[2], f, ii, jj, kk, SoA[0], SoA[1], SoA[2]);
    };

    if (i - 2 >= imin && i + 2 <= imax &&
        j - 2 >= jmin && j + 2 <= jmax &&
        k - 2 >= kmin && k + 2 <= kmax) {

        rhs_add = -eps / 16.0 * (
            ((fh(i-2,j,k) + fh(i+2,j,k)) - FOUR * (fh(i-1,j,k) + fh(i+1,j,k)) + SIX * fh(i,j,k)) / dX +
            ((fh(i,j-2,k) + fh(i,j+2,k)) - FOUR * (fh(i,j-1,k) + fh(i,j+1,k)) + SIX * fh(i,j,k)) / dY +
            ((fh(i,j,k-2) + fh(i,j,k+2)) - FOUR * (fh(i,j,k-1) + fh(i,j,k+1)) + SIX * fh(i,j,k)) / dZ
        );
    }
#elif ghost_width == 3
    const double FIT = 15.0;
    const double TWT = 20.0;
    const double cof = 64.0;

    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -3;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -3;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -3;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b(3, ex[0], ex[1], ex[2], f, ii, jj, kk, SoA[0], SoA[1], SoA[2]);
    };

    if (i - 3 >= imin && i + 3 <= imax &&
        j - 3 >= jmin && j + 3 <= jmax &&
        k - 3 >= kmin && k + 3 <= kmax) {

        rhs_add = eps / cof * (
            ((fh(i-3,j,k) + fh(i+3,j,k)) - SIX*(fh(i-2,j,k) + fh(i+2,j,k)) +
             FIT*(fh(i-1,j,k) + fh(i+1,j,k)) - TWT*fh(i,j,k)) / dX +
            ((fh(i,j-3,k) + fh(i,j+3,k)) - SIX*(fh(i,j-2,k) + fh(i,j+2,k)) +
             FIT*(fh(i,j-1,k) + fh(i,j+1,k)) - TWT*fh(i,j,k)) / dY +
            ((fh(i,j,k-3) + fh(i,j,k+3)) - SIX*(fh(i,j,k-2) + fh(i,j,k+2)) +
             FIT*(fh(i,j,k-1) + fh(i,j,k+1)) - TWT*fh(i,j,k)) / dZ
        );
    }
#elif ghost_width == 4
    const double C8 = 8.0;
    const double C28 = 28.0;
    const double C56 = 56.0;
    const double C70 = 70.0;
    const double cof = 256.0;

    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -4;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -4;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -4;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b(4, ex[0], ex[1], ex[2], f, ii, jj, kk, SoA[0], SoA[1], SoA[2]);
    };

    if (i - 4 >= imin && i + 4 <= imax &&
        j - 4 >= jmin && j + 4 <= jmax &&
        k - 4 >= kmin && k + 4 <= kmax) {
        rhs_add = -eps / cof * (
            ((fh(i-4,j,k) + fh(i+4,j,k)) - C8*(fh(i-3,j,k) + fh(i+3,j,k)) +
             C28*(fh(i-2,j,k) + fh(i+2,j,k)) - C56*(fh(i-1,j,k) + fh(i+1,j,k)) +
             C70*fh(i,j,k)) / dX +
            ((fh(i,j-4,k) + fh(i,j+4,k)) - C8*(fh(i,j-3,k) + fh(i,j+3,k)) +
             C28*(fh(i,j-2,k) + fh(i,j+2,k)) - C56*(fh(i,j-1,k) + fh(i,j+1,k)) +
             C70*fh(i,j,k)) / dY +
            ((fh(i,j,k-4) + fh(i,j,k+4)) - C8*(fh(i,j,k-3) + fh(i,j,k+3)) +
             C28*(fh(i,j,k-2) + fh(i,j,k+2)) - C56*(fh(i,j,k-1) + fh(i,j,k+1)) +
            C70*fh(i,j,k)) / dZ
        );
    }
#elif ghost_width == 5
    const double C10 = 10.0;
    const double C45 = 45.0;
    const double C120 = 120.0;
    const double C210 = 210.0;
    const double C252 = 252.0;
    const double cof = 1024.0;

    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -5;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -5;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -5;

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_bd_0b(5, ex[0], ex[1], ex[2], f, ii, jj, kk, SoA[0], SoA[1], SoA[2]);
    };

    if (i - 5 >= imin && i + 5 <= imax &&
        j - 5 >= jmin && j + 5 <= jmax &&
        k - 5 >= kmin && k + 5 <= kmax) {
        rhs_add = eps / cof * (
            ((fh(i-5,j,k) + fh(i+5,j,k)) - C10*(fh(i-4,j,k) + fh(i+4,j,k)) +
             C45*(fh(i-3,j,k) + fh(i+3,j,k)) - C120*(fh(i-2,j,k) + fh(i+2,j,k)) +
             C210*(fh(i-1,j,k) + fh(i+1,j,k)) - C252*fh(i,j,k)) / dX +
            ((fh(i,j-5,k) + fh(i,j+5,k)) - C10*(fh(i,j-4,k) + fh(i,j+4,k)) +
             C45*(fh(i,j-3,k) + fh(i,j+3,k)) - C120*(fh(i,j-2,k) + fh(i,j+2,k)) +
             C210*(fh(i,j-1,k) + fh(i,j+1,k)) - C252*fh(i,j,k)) / dY +
            ((fh(i,j,k-5) + fh(i,j,k+5)) - C10*(fh(i,j,k-4) + fh(i,j,k+4)) +
             C45*(fh(i,j,k-3) + fh(i,j,k+3)) - C120*(fh(i,j,k-2) + fh(i,j,k+2)) +
             C210*(fh(i,j,k-1) + fh(i,j,k+1)) - C252*fh(i,j,k)) / dZ
        );
    }
#else
#error Unsupported ghost_width in d_kodis_point; expected ghost_width == 2, 3, 4, or 5.
#endif

    return rhs_add;
}

#endif // KODISS_GPU_CUH
