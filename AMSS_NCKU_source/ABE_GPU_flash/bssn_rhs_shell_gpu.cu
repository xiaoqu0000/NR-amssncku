#include "bssn_rhs_shell_gpu.cuh"

#include "bssn_rhs.h"
#ifndef ghost_width
#include "macrodef.fh"
#endif

#include <cuda_runtime.h>

#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <iostream>

#if (GAUGE == 6 || GAUGE == 7)
extern int BHN;
extern double Mass[3];
extern double PBH[9];

__constant__ int d_bssn_shell_gauge_bh_count;
__constant__ double d_bssn_shell_gauge_mass[3];
__constant__ double d_bssn_shell_gauge_pbh[9];
#endif

namespace {

constexpr int kGuardBlockSize = 256;
constexpr int kShellDerivativeRadius = ghost_width - 1;
constexpr int kShellKODissRadius = ghost_width;
constexpr double kF1o3 = 1.0 / 3.0;
constexpr double kF2o3 = 2.0 / 3.0;
constexpr double kF3o2 = 1.5;
constexpr double kFF = 0.75;
constexpr double kEta = 2.0;
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kEightPi = 8.0 * kPi;
constexpr double kFourPi = 4.0 * kPi;
constexpr double kF16Pi = 16.0 * kPi;
constexpr double kZ4Ckappa1 = 2.0e-2;
constexpr double kZ4Ckappa2 = 0.0;

struct ShellBounds {
    int imin;
    int jmin;
    int kmin;
    int imax;
    int jmax;
    int kmax;
};

inline int shell_rhs_cell_count(const int ex[3]) {
    return ex[0] * ex[1] * ex[2];
}

inline int shell_rhs_grid_size(int count) {
    return (count + kGuardBlockSize - 1) / kGuardBlockSize;
}

void abort_on_cuda(cudaError_t status, const char* label) {
    if (status == cudaSuccess) {
        return;
    }
    std::cerr << label << ": " << cudaGetErrorString(status) << std::endl;
    std::abort();
}

double* alloc_shell_temp(int count, const char* label) {
    double* ptr = nullptr;
    abort_on_cuda(cudaMalloc(&ptr, static_cast<std::size_t>(count) * sizeof(double)), label);
    return ptr;
}

void free_shell_temp(double* ptr, const char* label) {
    abort_on_cuda(cudaFree(ptr), label);
}

__device__ __forceinline__ int shell_flatten(int i, int j, int k, int ex0, int ex1) {
    return (k * ex1 + j) * ex0 + i;
}

__device__ __forceinline__ void shell_select_soa(
    int sst,
    double sym1,
    double sym2,
    double sym3,
    double& soa0,
    double& soa1
) {
    if (sst == 2 || sst == 3) {
        soa0 = sym2;
        soa1 = sym3;
    } else if (sst == 4 || sst == 5) {
        soa0 = sym1;
        soa1 = sym3;
    } else {
        soa0 = sym1;
        soa1 = sym2;
    }
}

__device__ __forceinline__ double d_shell_symmetry_stbd_0b(
    int ord,
    int ex0,
    int ex1,
    int ex2,
    const double* f,
    int i,
    int j,
    int k,
    double soa0,
    double soa1
) {
    if (i < -ord || i >= ex0 + ord) return 0.0;
    if (j < -ord || j >= ex1 + ord) return 0.0;
    if (k < 0 || k >= ex2) return 0.0;

    double factor = 1.0;
#ifdef Vertex
    if (i < 0) {
        i = -i;
        factor *= soa0;
    } else if (i >= ex0) {
        i = 2 * ex0 - i - 2;
        factor *= soa0;
    }
    if (j < 0) {
        j = -j;
        factor *= soa1;
    } else if (j >= ex1) {
        j = 2 * ex1 - j - 2;
        factor *= soa1;
    }
#else
    if (i < 0) {
        i = -i - 1;
        factor *= soa0;
    } else if (i >= ex0) {
        i = 2 * ex0 - i - 1;
        factor *= soa0;
    }
    if (j < 0) {
        j = -j - 1;
        factor *= soa1;
    } else if (j >= ex1) {
        j = 2 * ex1 - j - 1;
        factor *= soa1;
    }
#endif

    if (i < 0 || i >= ex0) return 0.0;
    if (j < 0 || j >= ex1) return 0.0;
    return f[shell_flatten(i, j, k, ex0, ex1)] * factor;
}

__device__ __forceinline__ ShellBounds make_shell_bounds(
    int ex0,
    int ex1,
    int ex2,
    const double* crho,
    const double* sigma,
    const double* R,
    int symmetry,
    int sst,
    int radius
) {
    (void)R;
    const int no_symm = 0;
    const int eq_symm = 1;
    const int octant = 2;
    const double dX = crho[1] - crho[0];
    const double dY = sigma[1] - sigma[0];

    ShellBounds bounds{0, 0, 0, ex0 - 1, ex1 - 1, ex2 - 1};
    if (symmetry == octant) {
        if (fabs(crho[0]) < dX) bounds.imin = -radius;
        if (fabs(sigma[0]) < dY) bounds.jmin = -radius;
    } else if (symmetry == eq_symm) {
        if ((sst == 2 || sst == 4) && fabs(sigma[0]) < dY) {
            bounds.jmin = -radius;
        }
        if ((sst == 3 || sst == 5) && fabs(sigma[ex1 - 1] - 0.0) < dY) {
            bounds.jmax = ex1 - 1 + radius;
        }
    } else if (symmetry != no_symm) {
        bounds.imin = 0;
        bounds.jmin = 0;
    }
    return bounds;
}

__device__ __forceinline__ bool axis_fits(int center, int radius, int lower, int upper) {
    return center - radius >= lower && center + radius <= upper;
}

__device__ __forceinline__ double shell_axis_spacing(
    int axis,
    const double* crho,
    const double* sigma,
    const double* R
) {
    if (axis == 0) return crho[1] - crho[0];
    if (axis == 1) return sigma[1] - sigma[0];
    return R[1] - R[0];
}

__device__ __forceinline__ bool shell_axis_fits(
    int axis,
    int i,
    int j,
    int k,
    int radius,
    const ShellBounds& bounds
) {
    if (axis == 0) return axis_fits(i, radius, bounds.imin, bounds.imax);
    if (axis == 1) return axis_fits(j, radius, bounds.jmin, bounds.jmax);
    return axis_fits(k, radius, bounds.kmin, bounds.kmax);
}

__device__ __forceinline__ double shell_fetch_axis_offset(
    int ex0,
    int ex1,
    int ex2,
    const double* f,
    int i,
    int j,
    int k,
    int axis,
    int offset,
    int ord,
    double soa0,
    double soa1
) {
    if (axis == 0) i += offset;
    else if (axis == 1) j += offset;
    else k += offset;
    return d_shell_symmetry_stbd_0b(ord, ex0, ex1, ex2, f, i, j, k, soa0, soa1);
}

__device__ __forceinline__ double shell_fetch_two_axis_offset(
    int ex0,
    int ex1,
    int ex2,
    const double* f,
    int i,
    int j,
    int k,
    int axis_a,
    int offset_a,
    int axis_b,
    int offset_b,
    int ord,
    double soa0,
    double soa1
) {
    if (axis_a == 0) i += offset_a;
    else if (axis_a == 1) j += offset_a;
    else k += offset_a;

    if (axis_b == 0) i += offset_b;
    else if (axis_b == 1) j += offset_b;
    else k += offset_b;

    return d_shell_symmetry_stbd_0b(ord, ex0, ex1, ex2, f, i, j, k, soa0, soa1);
}

__device__ __forceinline__ int first_deriv_coeffs(
    int radius,
    int* offsets,
    double* coeffs,
    double& denom
) {
    if (radius == 4) {
        const int local_offsets[8] = {-4, -3, -2, -1, 1, 2, 3, 4};
        const double local_coeffs[8] = {3.0, -32.0, 168.0, -672.0, 672.0, -168.0, 32.0, -3.0};
        for (int n = 0; n < 8; ++n) { offsets[n] = local_offsets[n]; coeffs[n] = local_coeffs[n]; }
        denom = 840.0;
        return 8;
    }
    if (radius == 3) {
        const int local_offsets[6] = {-3, -2, -1, 1, 2, 3};
        const double local_coeffs[6] = {-1.0, 9.0, -45.0, 45.0, -9.0, 1.0};
        for (int n = 0; n < 6; ++n) { offsets[n] = local_offsets[n]; coeffs[n] = local_coeffs[n]; }
        denom = 60.0;
        return 6;
    }
    if (radius == 2) {
        const int local_offsets[4] = {-2, -1, 1, 2};
        const double local_coeffs[4] = {1.0, -8.0, 8.0, -1.0};
        for (int n = 0; n < 4; ++n) { offsets[n] = local_offsets[n]; coeffs[n] = local_coeffs[n]; }
        denom = 12.0;
        return 4;
    }
    offsets[0] = -1;
    offsets[1] = 1;
    coeffs[0] = -1.0;
    coeffs[1] = 1.0;
    denom = 2.0;
    return 2;
}

__device__ __forceinline__ int second_deriv_coeffs(
    int radius,
    int* offsets,
    double* coeffs,
    double& denom
) {
    if (radius == 4) {
        const int local_offsets[9] = {-4, -3, -2, -1, 0, 1, 2, 3, 4};
        const double local_coeffs[9] = {-9.0, 128.0, -1008.0, 8064.0, -14350.0, 8064.0, -1008.0, 128.0, -9.0};
        for (int n = 0; n < 9; ++n) { offsets[n] = local_offsets[n]; coeffs[n] = local_coeffs[n]; }
        denom = 5040.0;
        return 9;
    }
    if (radius == 3) {
        const int local_offsets[7] = {-3, -2, -1, 0, 1, 2, 3};
        const double local_coeffs[7] = {2.0, -27.0, 270.0, -490.0, 270.0, -27.0, 2.0};
        for (int n = 0; n < 7; ++n) { offsets[n] = local_offsets[n]; coeffs[n] = local_coeffs[n]; }
        denom = 180.0;
        return 7;
    }
    if (radius == 2) {
        const int local_offsets[5] = {-2, -1, 0, 1, 2};
        const double local_coeffs[5] = {-1.0, 16.0, -30.0, 16.0, -1.0};
        for (int n = 0; n < 5; ++n) { offsets[n] = local_offsets[n]; coeffs[n] = local_coeffs[n]; }
        denom = 12.0;
        return 5;
    }
    offsets[0] = -1;
    offsets[1] = 0;
    offsets[2] = 1;
    coeffs[0] = 1.0;
    coeffs[1] = -2.0;
    coeffs[2] = 1.0;
    denom = 1.0;
    return 3;
}

__device__ __forceinline__ double d_shell_first_axis(
    int ex0,
    int ex1,
    int ex2,
    const double* f,
    const double* crho,
    const double* sigma,
    const double* R,
    int i,
    int j,
    int k,
    int axis,
    int radius,
    double soa0,
    double soa1
) {
    int offsets[8];
    double coeffs[8];
    double denom = 1.0;
    const int count = first_deriv_coeffs(radius, offsets, coeffs, denom);
    double sum = 0.0;
    for (int n = 0; n < count; ++n) {
        sum += coeffs[n] * shell_fetch_axis_offset(
            ex0, ex1, ex2, f, i, j, k, axis, offsets[n],
            radius, soa0, soa1);
    }
    return sum / (denom * shell_axis_spacing(axis, crho, sigma, R));
}

__device__ __forceinline__ double d_shell_second_axis(
    int ex0,
    int ex1,
    int ex2,
    const double* f,
    const double* crho,
    const double* sigma,
    const double* R,
    int i,
    int j,
    int k,
    int axis,
    int radius,
    double soa0,
    double soa1
) {
    int offsets[9];
    double coeffs[9];
    double denom = 1.0;
    const int count = second_deriv_coeffs(radius, offsets, coeffs, denom);
    double sum = 0.0;
    for (int n = 0; n < count; ++n) {
        sum += coeffs[n] * shell_fetch_axis_offset(
            ex0, ex1, ex2, f, i, j, k, axis, offsets[n],
            radius, soa0, soa1);
    }
    const double spacing = shell_axis_spacing(axis, crho, sigma, R);
    return sum / (denom * spacing * spacing);
}

__device__ __forceinline__ double d_shell_mixed_axis(
    int ex0,
    int ex1,
    int ex2,
    const double* f,
    const double* crho,
    const double* sigma,
    const double* R,
    int i,
    int j,
    int k,
    int axis_a,
    int axis_b,
    int radius,
    double soa0,
    double soa1
) {
    int offsets_a[8];
    int offsets_b[8];
    double coeffs_a[8];
    double coeffs_b[8];
    double denom_a = 1.0;
    double denom_b = 1.0;
    const int count_a = first_deriv_coeffs(radius, offsets_a, coeffs_a, denom_a);
    const int count_b = first_deriv_coeffs(radius, offsets_b, coeffs_b, denom_b);
    double sum = 0.0;
    for (int a = 0; a < count_a; ++a) {
        for (int b = 0; b < count_b; ++b) {
            sum += coeffs_a[a] * coeffs_b[b] * shell_fetch_two_axis_offset(
                ex0, ex1, ex2, f, i, j, k,
                axis_a, offsets_a[a],
                axis_b, offsets_b[b],
                radius, soa0, soa1);
        }
    }
    const double spacing_a = shell_axis_spacing(axis_a, crho, sigma, R);
    const double spacing_b = shell_axis_spacing(axis_b, crho, sigma, R);
    return sum / (denom_a * denom_b * spacing_a * spacing_b);
}

__device__ __forceinline__ int best_axis_radius(
    int axis,
    int i,
    int j,
    int k,
    const ShellBounds& bounds,
    int max_radius
) {
    for (int radius = max_radius; radius >= 1; --radius) {
        if (shell_axis_fits(axis, i, j, k, radius, bounds)) {
            return radius;
        }
    }
    return 0;
}

__device__ __forceinline__ int best_mixed_radius(
    int axis_a,
    int axis_b,
    int i,
    int j,
    int k,
    const ShellBounds& bounds,
    int max_radius
) {
    for (int radius = max_radius; radius >= 1; --radius) {
        if (shell_axis_fits(axis_a, i, j, k, radius, bounds) &&
            shell_axis_fits(axis_b, i, j, k, radius, bounds)) {
            return radius;
        }
    }
    return 0;
}

__device__ __forceinline__ void d_fderivs_sh_point(
    const int ex[3],
    const double* f,
    double* fx,
    double* fy,
    double* fz,
    const double* crho,
    const double* sigma,
    const double* R,
    double sym1,
    double sym2,
    double sym3,
    int symmetry,
    int sst,
    int i,
    int j,
    int k
) {
    *fx = 0.0;
    *fy = 0.0;
    *fz = 0.0;

    double soa0 = 1.0;
    double soa1 = 1.0;
    shell_select_soa(sst, sym1, sym2, sym3, soa0, soa1);
    const ShellBounds bounds = make_shell_bounds(
        ex[0], ex[1], ex[2], crho, sigma, R, symmetry, sst, kShellDerivativeRadius);

    const int rx = best_axis_radius(0, i, j, k, bounds, kShellDerivativeRadius);
    const int ry = best_axis_radius(1, i, j, k, bounds, kShellDerivativeRadius);
    const int rz = best_axis_radius(2, i, j, k, bounds, kShellDerivativeRadius);
    if (rx > 0) *fx = d_shell_first_axis(ex[0], ex[1], ex[2], f, crho, sigma, R, i, j, k, 0, rx, soa0, soa1);
    if (ry > 0) *fy = d_shell_first_axis(ex[0], ex[1], ex[2], f, crho, sigma, R, i, j, k, 1, ry, soa0, soa1);
    if (rz > 0) *fz = d_shell_first_axis(ex[0], ex[1], ex[2], f, crho, sigma, R, i, j, k, 2, rz, soa0, soa1);
}

__device__ __forceinline__ void d_fdderivs_sh_point(
    const int ex[3],
    const double* f,
    double* fxx,
    double* fxy,
    double* fxz,
    double* fyy,
    double* fyz,
    double* fzz,
    const double* crho,
    const double* sigma,
    const double* R,
    double sym1,
    double sym2,
    double sym3,
    int symmetry,
    int sst,
    int i,
    int j,
    int k
) {
    *fxx = 0.0;
    *fxy = 0.0;
    *fxz = 0.0;
    *fyy = 0.0;
    *fyz = 0.0;
    *fzz = 0.0;

    double soa0 = 1.0;
    double soa1 = 1.0;
    shell_select_soa(sst, sym1, sym2, sym3, soa0, soa1);
    const ShellBounds bounds = make_shell_bounds(
        ex[0], ex[1], ex[2], crho, sigma, R, symmetry, sst, kShellDerivativeRadius);

    const int rx = best_axis_radius(0, i, j, k, bounds, kShellDerivativeRadius);
    const int ry = best_axis_radius(1, i, j, k, bounds, kShellDerivativeRadius);
    const int rz = best_axis_radius(2, i, j, k, bounds, kShellDerivativeRadius);
    if (rx > 0) *fxx = d_shell_second_axis(ex[0], ex[1], ex[2], f, crho, sigma, R, i, j, k, 0, rx, soa0, soa1);
    if (ry > 0) *fyy = d_shell_second_axis(ex[0], ex[1], ex[2], f, crho, sigma, R, i, j, k, 1, ry, soa0, soa1);
    if (rz > 0) *fzz = d_shell_second_axis(ex[0], ex[1], ex[2], f, crho, sigma, R, i, j, k, 2, rz, soa0, soa1);

    const int rxy = best_mixed_radius(0, 1, i, j, k, bounds, kShellDerivativeRadius);
    const int rxz = best_mixed_radius(0, 2, i, j, k, bounds, kShellDerivativeRadius);
    const int ryz = best_mixed_radius(1, 2, i, j, k, bounds, kShellDerivativeRadius);
    if (rxy > 0) *fxy = d_shell_mixed_axis(ex[0], ex[1], ex[2], f, crho, sigma, R, i, j, k, 0, 1, rxy, soa0, soa1);
    if (rxz > 0) *fxz = d_shell_mixed_axis(ex[0], ex[1], ex[2], f, crho, sigma, R, i, j, k, 0, 2, rxz, soa0, soa1);
    if (ryz > 0) *fyz = d_shell_mixed_axis(ex[0], ex[1], ex[2], f, crho, sigma, R, i, j, k, 1, 2, ryz, soa0, soa1);
}

__device__ __forceinline__ void d_fderivs_shc_point(
    const int ex[3],
    const double* f,
    double* fx,
    double* fy,
    double* fz,
    const double* crho,
    const double* sigma,
    const double* R,
    const double* drhodx,
    const double* drhody,
    const double* drhodz,
    const double* dsigmadx,
    const double* dsigmady,
    const double* dsigmadz,
    const double* dRdx,
    const double* dRdy,
    const double* dRdz,
    double sym1,
    double sym2,
    double sym3,
    int symmetry,
    int sst,
    int i,
    int j,
    int k
) {
    double gx = 0.0;
    double gy = 0.0;
    double gz = 0.0;
    d_fderivs_sh_point(ex, f, &gx, &gy, &gz, crho, sigma, R, sym1, sym2, sym3, symmetry, sst, i, j, k);

    const int idx = shell_flatten(i, j, k, ex[0], ex[1]);
    *fx = dRdx[idx] * gz + drhodx[idx] * gx + dsigmadx[idx] * gy;
    *fy = dRdy[idx] * gz + drhody[idx] * gx + dsigmady[idx] * gy;
    *fz = dRdz[idx] * gz + drhodz[idx] * gx + dsigmadz[idx] * gy;
}

__device__ __forceinline__ void d_fdderivs_shc_point(
    const int ex[3],
    const double* f,
    double* fxx,
    double* fxy,
    double* fxz,
    double* fyy,
    double* fyz,
    double* fzz,
    const double* crho,
    const double* sigma,
    const double* R,
    const double* drhodx,
    const double* drhody,
    const double* drhodz,
    const double* dsigmadx,
    const double* dsigmady,
    const double* dsigmadz,
    const double* dRdx,
    const double* dRdy,
    const double* dRdz,
    const double* drhodxx,
    const double* drhodxy,
    const double* drhodxz,
    const double* drhodyy,
    const double* drhodyz,
    const double* drhodzz,
    const double* dsigmadxx,
    const double* dsigmadxy,
    const double* dsigmadxz,
    const double* dsigmadyy,
    const double* dsigmadyz,
    const double* dsigmadzz,
    const double* dRdxx,
    const double* dRdxy,
    const double* dRdxz,
    const double* dRdyy,
    const double* dRdyz,
    const double* dRdzz,
    double sym1,
    double sym2,
    double sym3,
    int symmetry,
    int sst,
    int i,
    int j,
    int k
) {
    double gx = 0.0;
    double gy = 0.0;
    double gz = 0.0;
    double gxx = 0.0;
    double gxy = 0.0;
    double gxz = 0.0;
    double gyy = 0.0;
    double gyz = 0.0;
    double gzz = 0.0;
    d_fderivs_sh_point(ex, f, &gx, &gy, &gz, crho, sigma, R, sym1, sym2, sym3, symmetry, sst, i, j, k);
    d_fdderivs_sh_point(ex, f, &gxx, &gxy, &gxz, &gyy, &gyz, &gzz,
                        crho, sigma, R, sym1, sym2, sym3, symmetry, sst, i, j, k);

    const int idx = shell_flatten(i, j, k, ex[0], ex[1]);
    *fxx = dRdxx[idx] * gz + drhodxx[idx] * gx + dsigmadxx[idx] * gy
         + dRdx[idx] * dRdx[idx] * gzz + drhodx[idx] * drhodx[idx] * gxx + dsigmadx[idx] * dsigmadx[idx] * gyy
         + 2.0 * (dRdx[idx] * drhodx[idx] * gxz + dRdx[idx] * dsigmadx[idx] * gyz + drhodx[idx] * dsigmadx[idx] * gxy);
    *fyy = dRdyy[idx] * gz + drhodyy[idx] * gx + dsigmadyy[idx] * gy
         + dRdy[idx] * dRdy[idx] * gzz + drhody[idx] * drhody[idx] * gxx + dsigmady[idx] * dsigmady[idx] * gyy
         + 2.0 * (dRdy[idx] * drhody[idx] * gxz + dRdy[idx] * dsigmady[idx] * gyz + drhody[idx] * dsigmady[idx] * gxy);
    *fzz = dRdzz[idx] * gz + drhodzz[idx] * gx + dsigmadzz[idx] * gy
         + dRdz[idx] * dRdz[idx] * gzz + drhodz[idx] * drhodz[idx] * gxx + dsigmadz[idx] * dsigmadz[idx] * gyy
         + 2.0 * (dRdz[idx] * drhodz[idx] * gxz + dRdz[idx] * dsigmadz[idx] * gyz + drhodz[idx] * dsigmadz[idx] * gxy);
    *fxy = dRdxy[idx] * gz + drhodxy[idx] * gx + dsigmadxy[idx] * gy
         + dRdx[idx] * drhody[idx] * gxz + dRdx[idx] * dsigmady[idx] * gyz + drhodx[idx] * dsigmady[idx] * gxy
         + dRdy[idx] * drhodx[idx] * gxz + dRdy[idx] * dsigmadx[idx] * gyz + drhody[idx] * dsigmadx[idx] * gxy
         + dRdx[idx] * dRdy[idx] * gzz + drhodx[idx] * drhody[idx] * gxx + dsigmadx[idx] * dsigmady[idx] * gyy;
    *fxz = dRdxz[idx] * gz + drhodxz[idx] * gx + dsigmadxz[idx] * gy
         + dRdx[idx] * drhodz[idx] * gxz + dRdx[idx] * dsigmadz[idx] * gyz + drhodx[idx] * dsigmadz[idx] * gxy
         + dRdz[idx] * drhodx[idx] * gxz + dRdz[idx] * dsigmadx[idx] * gyz + drhodz[idx] * dsigmadx[idx] * gxy
         + dRdx[idx] * dRdz[idx] * gzz + drhodx[idx] * drhodz[idx] * gxx + dsigmadx[idx] * dsigmadz[idx] * gyy;
    *fyz = dRdyz[idx] * gz + drhodyz[idx] * gx + dsigmadyz[idx] * gy
         + dRdz[idx] * drhody[idx] * gxz + dRdz[idx] * dsigmady[idx] * gyz + drhodz[idx] * dsigmady[idx] * gxy
         + dRdy[idx] * drhodz[idx] * gxz + dRdy[idx] * dsigmadz[idx] * gyz + drhody[idx] * dsigmadz[idx] * gxy
         + dRdz[idx] * dRdy[idx] * gzz + drhodz[idx] * drhody[idx] * gxx + dsigmadz[idx] * dsigmady[idx] * gyy;
}

__device__ __forceinline__ double d_shell_contract_second(
    double gupxx,
    double gupxy,
    double gupxz,
    double gupyy,
    double gupyz,
    double gupzz,
    double fxx,
    double fxy,
    double fxz,
    double fyy,
    double fyz,
    double fzz
) {
    return gupxx * fxx + gupyy * fyy + gupzz * fzz
        + 2.0 * (gupxy * fxy + gupxz * fxz + gupyz * fyz);
}

__device__ __forceinline__ int kodis_coeffs(int radius, int* offsets, double* coeffs, double& scale) {
    if (radius == 5) {
        const int local_offsets[11] = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5};
        const double local_coeffs[11] = {1.0, -10.0, 45.0, -120.0, 210.0, -252.0, 210.0, -120.0, 45.0, -10.0, 1.0};
        for (int n = 0; n < 11; ++n) { offsets[n] = local_offsets[n]; coeffs[n] = local_coeffs[n]; }
        scale = 1.0 / 1024.0;
        return 11;
    }
    if (radius == 4) {
        const int local_offsets[9] = {-4, -3, -2, -1, 0, 1, 2, 3, 4};
        const double local_coeffs[9] = {1.0, -8.0, 28.0, -56.0, 70.0, -56.0, 28.0, -8.0, 1.0};
        for (int n = 0; n < 9; ++n) { offsets[n] = local_offsets[n]; coeffs[n] = local_coeffs[n]; }
        scale = -1.0 / 256.0;
        return 9;
    }
    if (radius == 3) {
        const int local_offsets[7] = {-3, -2, -1, 0, 1, 2, 3};
        const double local_coeffs[7] = {1.0, -6.0, 15.0, -20.0, 15.0, -6.0, 1.0};
        for (int n = 0; n < 7; ++n) { offsets[n] = local_offsets[n]; coeffs[n] = local_coeffs[n]; }
        scale = 1.0 / 64.0;
        return 7;
    }
    const int local_offsets[5] = {-2, -1, 0, 1, 2};
    const double local_coeffs[5] = {1.0, -4.0, 6.0, -4.0, 1.0};
    for (int n = 0; n < 5; ++n) { offsets[n] = local_offsets[n]; coeffs[n] = local_coeffs[n]; }
    scale = -1.0 / 16.0;
    return 5;
}

__device__ __forceinline__ double d_kodis_sh_point(
    const int ex[3],
    const double* f,
    const double* crho,
    const double* sigma,
    const double* R,
    double sym1,
    double sym2,
    double sym3,
    int symmetry,
    double eps,
    int sst,
    int i,
    int j,
    int k
) {
    double soa0 = 1.0;
    double soa1 = 1.0;
    shell_select_soa(sst, sym1, sym2, sym3, soa0, soa1);
    const ShellBounds bounds = make_shell_bounds(
        ex[0], ex[1], ex[2], crho, sigma, R, symmetry, sst, kShellKODissRadius);
    if (!shell_axis_fits(0, i, j, k, kShellKODissRadius, bounds) ||
        !shell_axis_fits(1, i, j, k, kShellKODissRadius, bounds) ||
        !shell_axis_fits(2, i, j, k, kShellKODissRadius, bounds)) {
        return 0.0;
    }

    int offsets[11];
    double coeffs[11];
    double scale = 0.0;
    const int count = kodis_coeffs(kShellKODissRadius, offsets, coeffs, scale);
    double sx = 0.0;
    double sy = 0.0;
    double sz = 0.0;
    for (int n = 0; n < count; ++n) {
        sx += coeffs[n] * shell_fetch_axis_offset(
            ex[0], ex[1], ex[2], f, i, j, k, 0, offsets[n], kShellKODissRadius, soa0, soa1);
        sy += coeffs[n] * shell_fetch_axis_offset(
            ex[0], ex[1], ex[2], f, i, j, k, 1, offsets[n], kShellKODissRadius, soa0, soa1);
        sz += coeffs[n] * shell_fetch_axis_offset(
            ex[0], ex[1], ex[2], f, i, j, k, 2, offsets[n], kShellKODissRadius, soa0, soa1);
    }
    return eps * scale * (
        sx / shell_axis_spacing(0, crho, sigma, R) +
        sy / shell_axis_spacing(1, crho, sigma, R) +
        sz / shell_axis_spacing(2, crho, sigma, R));
}

__device__ __forceinline__ void d_apply_transport_ko_rhs(
    const int ex[3],
    const double* crho,
    const double* sigma,
    const double* R,
    const double* drhodx,
    const double* drhody,
    const double* drhodz,
    const double* dsigmadx,
    const double* dsigmady,
    const double* dsigmadz,
    const double* dRdx,
    const double* dRdy,
    const double* dRdz,
    const double* field,
    const double* betax,
    const double* betay,
    const double* betaz,
    double* rhs,
    double sym1,
    double sym2,
    double sym3,
    int symmetry,
    double eps,
    int sst,
    int i,
    int j,
    int k
) {
    const int idx = shell_flatten(i, j, k, ex[0], ex[1]);
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    d_fderivs_shc_point(
        ex, field, &fx, &fy, &fz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        sym1, sym2, sym3, symmetry, sst, i, j, k);

    double adv = betax[idx] * fx;
    adv += betay[idx] * fy;
    adv += betaz[idx] * fz;
    rhs[idx] += adv;

    if (eps > 0.0) {
        rhs[idx] += d_kodis_sh_point(
            ex, field, crho, sigma, R,
            sym1, sym2, sym3, symmetry, eps, sst, i, j, k);
    }
}

__global__ void bssn_shell_derivative_probe_kernel(
    int ex0,
    int ex1,
    int ex2,
    const double* crho,
    const double* sigma,
    const double* R,
    const double* f,
    const double* drhodx,
    const double* drhody,
    const double* drhodz,
    const double* dsigmadx,
    const double* dsigmady,
    const double* dsigmadz,
    const double* dRdx,
    const double* dRdy,
    const double* dRdz,
    const double* drhodxx,
    const double* drhodxy,
    const double* drhodxz,
    const double* drhodyy,
    const double* drhodyz,
    const double* drhodzz,
    const double* dsigmadxx,
    const double* dsigmadxy,
    const double* dsigmadxz,
    const double* dsigmadyy,
    const double* dsigmadyz,
    const double* dsigmadzz,
    const double* dRdxx,
    const double* dRdxy,
    const double* dRdxz,
    const double* dRdyy,
    const double* dRdyz,
    const double* dRdzz,
    double* fx,
    double* fy,
    double* fz,
    double* fxx,
    double* fxy,
    double* fxz,
    double* fyy,
    double* fyz,
    double* fzz,
    int symmetry,
    int sst,
    double sym1,
    double sym2,
    double sym3
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int cell_count = ex0 * ex1 * ex2;
    if (idx >= cell_count) {
        return;
    }
    const int i = idx % ex0;
    const int j = (idx / ex0) % ex1;
    const int k = idx / (ex0 * ex1);
    const int ex[3] = {ex0, ex1, ex2};

    d_fderivs_shc_point(
        ex, f, &fx[idx], &fy[idx], &fz[idx],
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        sym1, sym2, sym3, symmetry, sst, i, j, k);
    d_fdderivs_shc_point(
        ex, f, &fxx[idx], &fxy[idx], &fxz[idx], &fyy[idx], &fyz[idx], &fzz[idx],
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        sym1, sym2, sym3, symmetry, sst, i, j, k);
}

__global__ void bssn_shell_metric_rhs_probe_kernel(
    int ex0,
    int ex1,
    int ex2,
    const double* crho,
    const double* sigma,
    const double* R,
    const double* drhodx,
    const double* drhody,
    const double* drhodz,
    const double* dsigmadx,
    const double* dsigmady,
    const double* dsigmadz,
    const double* dRdx,
    const double* dRdy,
    const double* dRdz,
    const double* chi,
    const double* trK,
    const double* dxx,
    const double* gxy,
    const double* gxz,
    const double* dyy,
    const double* gyz,
    const double* dzz,
    const double* Axx,
    const double* Axy,
    const double* Axz,
    const double* Ayy,
    const double* Ayz,
    const double* Azz,
    const double* Lap,
    const double* betax,
    const double* betay,
    const double* betaz,
    double* chi_rhs,
    double* gxx_rhs,
    double* gxy_rhs,
    double* gxz_rhs,
    double* gyy_rhs,
    double* gyz_rhs,
    double* gzz_rhs,
    int symmetry,
    int sst
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int cell_count = ex0 * ex1 * ex2;
    if (idx >= cell_count) {
        return;
    }
    const int i = idx % ex0;
    const int j = (idx / ex0) % ex1;
    const int k = idx / (ex0 * ex1);
    const int ex[3] = {ex0, ex1, ex2};

    double betaxx = 0.0;
    double betaxy = 0.0;
    double betaxz = 0.0;
    double betayx = 0.0;
    double betayy = 0.0;
    double betayz = 0.0;
    double betazx = 0.0;
    double betazy = 0.0;
    double betazz = 0.0;
    d_fderivs_shc_point(
        ex, betax, &betaxx, &betaxy, &betaxz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        -1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, betay, &betayx, &betayy, &betayz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, -1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, betaz, &betazx, &betazy, &betazz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, -1.0, symmetry, sst, i, j, k);

    const double div_beta = betaxx + betayy + betazz;
    const double alpn1 = Lap[idx] + 1.0;
    const double chin1 = chi[idx] + 1.0;
    const double gxx = dxx[idx] + 1.0;
    const double gyy = dyy[idx] + 1.0;
    const double gzz = dzz[idx] + 1.0;

    chi_rhs[idx] = kF2o3 * chin1 * (alpn1 * trK[idx] - div_beta);
    gxx_rhs[idx] = -2.0 * alpn1 * Axx[idx] - kF2o3 * gxx * div_beta
        + 2.0 * (gxx * betaxx + gxy[idx] * betayx + gxz[idx] * betazx);
    gyy_rhs[idx] = -2.0 * alpn1 * Ayy[idx] - kF2o3 * gyy * div_beta
        + 2.0 * (gxy[idx] * betaxy + gyy * betayy + gyz[idx] * betazy);
    gzz_rhs[idx] = -2.0 * alpn1 * Azz[idx] - kF2o3 * gzz * div_beta
        + 2.0 * (gxz[idx] * betaxz + gyz[idx] * betayz + gzz * betazz);
    gxy_rhs[idx] = -2.0 * alpn1 * Axy[idx] + kF1o3 * gxy[idx] * div_beta
        + gxx * betaxy + gxz[idx] * betazy
        + gyy * betayx + gyz[idx] * betazx
        - gxy[idx] * betazz;
    gyz_rhs[idx] = -2.0 * alpn1 * Ayz[idx] + kF1o3 * gyz[idx] * div_beta
        + gxy[idx] * betaxz + gyy * betayz
        + gxz[idx] * betaxy + gzz * betazy
        - gyz[idx] * betaxx;
    gxz_rhs[idx] = -2.0 * alpn1 * Axz[idx] + kF1o3 * gxz[idx] * div_beta
        + gxx * betaxz + gxy[idx] * betayz
        + gyz[idx] * betayx + gzz * betazx
        - gxz[idx] * betayy;
}

__global__ void bssn_shell_connection_probe_kernel(
    int ex0,
    int ex1,
    int ex2,
    const double* crho,
    const double* sigma,
    const double* R,
    const double* drhodx,
    const double* drhody,
    const double* drhodz,
    const double* dsigmadx,
    const double* dsigmady,
    const double* dsigmadz,
    const double* dRdx,
    const double* dRdy,
    const double* dRdz,
    const double* dxx,
    const double* gxy,
    const double* gxz,
    const double* dyy,
    const double* gyz,
    const double* dzz,
    const double* Axx,
    const double* Axy,
    const double* Axz,
    const double* Ayy,
    const double* Ayz,
    const double* Azz,
    double* Gamxxx,
    double* Gamxxy,
    double* Gamxxz,
    double* Gamxyy,
    double* Gamxyz,
    double* Gamxzz,
    double* Gamyxx,
    double* Gamyxy,
    double* Gamyxz,
    double* Gamyyy,
    double* Gamyyz,
    double* Gamyzz,
    double* Gamzxx,
    double* Gamzxy,
    double* Gamzxz,
    double* Gamzyy,
    double* Gamzyz,
    double* Gamzzz,
    double* Rxx,
    double* Rxy,
    double* Rxz,
    double* Ryy,
    double* Ryz,
    double* Rzz,
    int symmetry,
    int sst
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int cell_count = ex0 * ex1 * ex2;
    if (idx >= cell_count) {
        return;
    }
    const int i = idx % ex0;
    const int j = (idx / ex0) % ex1;
    const int k = idx / (ex0 * ex1);
    const int ex[3] = {ex0, ex1, ex2};

    double gxxx = 0.0;
    double gxxy = 0.0;
    double gxxz = 0.0;
    double gxyx = 0.0;
    double gxyy = 0.0;
    double gxyz = 0.0;
    double gxzx = 0.0;
    double gxzy = 0.0;
    double gxzz = 0.0;
    double gyyx = 0.0;
    double gyyy = 0.0;
    double gyyz = 0.0;
    double gyzx = 0.0;
    double gyzy = 0.0;
    double gyzz = 0.0;
    double gzzx = 0.0;
    double gzzy = 0.0;
    double gzzz = 0.0;

    d_fderivs_shc_point(
        ex, dxx, &gxxx, &gxxy, &gxxz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, gxy, &gxyx, &gxyy, &gxyz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        -1.0, -1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, gxz, &gxzx, &gxzy, &gxzz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        -1.0, 1.0, -1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, dyy, &gyyx, &gyyy, &gyyz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, gyz, &gyzx, &gyzy, &gyzz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, -1.0, -1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, dzz, &gzzx, &gzzy, &gzzz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, 1.0, symmetry, sst, i, j, k);

    const double l_gxx = dxx[idx] + 1.0;
    const double l_gyy = dyy[idx] + 1.0;
    const double l_gzz = dzz[idx] + 1.0;
    const double l_gxy = gxy[idx];
    const double l_gxz = gxz[idx];
    const double l_gyz = gyz[idx];
    const double det = l_gxx * l_gyy * l_gzz
        + l_gxy * l_gyz * l_gxz + l_gxz * l_gxy * l_gyz
        - l_gxz * l_gyy * l_gxz - l_gxy * l_gxy * l_gzz
        - l_gxx * l_gyz * l_gyz;
    const double gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) / det;
    const double gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) / det;
    const double gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) / det;
    const double gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) / det;
    const double gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) / det;
    const double gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) / det;

    const double cxxx = 0.5 * (gupxx * gxxx + gupxy * (2.0 * gxyx - gxxy) + gupxz * (2.0 * gxzx - gxxz));
    const double cyxx = 0.5 * (gupxy * gxxx + gupyy * (2.0 * gxyx - gxxy) + gupyz * (2.0 * gxzx - gxxz));
    const double czxx = 0.5 * (gupxz * gxxx + gupyz * (2.0 * gxyx - gxxy) + gupzz * (2.0 * gxzx - gxxz));
    const double cxyy = 0.5 * (gupxx * (2.0 * gxyy - gyyx) + gupxy * gyyy + gupxz * (2.0 * gyzy - gyyz));
    const double cyyy = 0.5 * (gupxy * (2.0 * gxyy - gyyx) + gupyy * gyyy + gupyz * (2.0 * gyzy - gyyz));
    const double czyy = 0.5 * (gupxz * (2.0 * gxyy - gyyx) + gupyz * gyyy + gupzz * (2.0 * gyzy - gyyz));
    const double cxzz = 0.5 * (gupxx * (2.0 * gxzz - gzzx) + gupxy * (2.0 * gyzz - gzzy) + gupxz * gzzz);
    const double cyzz = 0.5 * (gupxy * (2.0 * gxzz - gzzx) + gupyy * (2.0 * gyzz - gzzy) + gupyz * gzzz);
    const double czzz = 0.5 * (gupxz * (2.0 * gxzz - gzzx) + gupyz * (2.0 * gyzz - gzzy) + gupzz * gzzz);
    const double cxxy = 0.5 * (gupxx * gxxy + gupxy * gyyx + gupxz * (gxzy + gyzx - gxyz));
    const double cyxy = 0.5 * (gupxy * gxxy + gupyy * gyyx + gupyz * (gxzy + gyzx - gxyz));
    const double czxy = 0.5 * (gupxz * gxxy + gupyz * gyyx + gupzz * (gxzy + gyzx - gxyz));
    const double cxxz = 0.5 * (gupxx * gxxz + gupxy * (gxyz + gyzx - gxzy) + gupxz * gzzx);
    const double cyxz = 0.5 * (gupxy * gxxz + gupyy * (gxyz + gyzx - gxzy) + gupyz * gzzx);
    const double czxz = 0.5 * (gupxz * gxxz + gupyz * (gxyz + gyzx - gxzy) + gupzz * gzzx);
    const double cxyz = 0.5 * (gupxx * (gxyz + gxzy - gyzx) + gupxy * gyyz + gupxz * gzzy);
    const double cyyz = 0.5 * (gupxy * (gxyz + gxzy - gyzx) + gupyy * gyyz + gupyz * gzzy);
    const double czyz = 0.5 * (gupxz * (gxyz + gxzy - gyzx) + gupyz * gyyz + gupzz * gzzy);

    Gamxxx[idx] = cxxx;
    Gamyxx[idx] = cyxx;
    Gamzxx[idx] = czxx;
    Gamxyy[idx] = cxyy;
    Gamyyy[idx] = cyyy;
    Gamzyy[idx] = czyy;
    Gamxzz[idx] = cxzz;
    Gamyzz[idx] = cyzz;
    Gamzzz[idx] = czzz;
    Gamxxy[idx] = cxxy;
    Gamyxy[idx] = cyxy;
    Gamzxy[idx] = czxy;
    Gamxxz[idx] = cxxz;
    Gamyxz[idx] = cyxz;
    Gamzxz[idx] = czxz;
    Gamxyz[idx] = cxyz;
    Gamyyz[idx] = cyyz;
    Gamzyz[idx] = czyz;

    Rxx[idx] = gupxx * gupxx * Axx[idx] + gupxy * gupxy * Ayy[idx] + gupxz * gupxz * Azz[idx]
        + 2.0 * (gupxx * gupxy * Axy[idx] + gupxx * gupxz * Axz[idx] + gupxy * gupxz * Ayz[idx]);
    Ryy[idx] = gupxy * gupxy * Axx[idx] + gupyy * gupyy * Ayy[idx] + gupyz * gupyz * Azz[idx]
        + 2.0 * (gupxy * gupyy * Axy[idx] + gupxy * gupyz * Axz[idx] + gupyy * gupyz * Ayz[idx]);
    Rzz[idx] = gupxz * gupxz * Axx[idx] + gupyz * gupyz * Ayy[idx] + gupzz * gupzz * Azz[idx]
        + 2.0 * (gupxz * gupyz * Axy[idx] + gupxz * gupzz * Axz[idx] + gupyz * gupzz * Ayz[idx]);
    Rxy[idx] = gupxx * gupxy * Axx[idx] + gupxy * gupyy * Ayy[idx] + gupxz * gupyz * Azz[idx]
        + (gupxx * gupyy + gupxy * gupxy) * Axy[idx]
        + (gupxx * gupyz + gupxz * gupxy) * Axz[idx]
        + (gupxy * gupyz + gupxz * gupyy) * Ayz[idx];
    Rxz[idx] = gupxx * gupxz * Axx[idx] + gupxy * gupyz * Ayy[idx] + gupxz * gupzz * Azz[idx]
        + (gupxx * gupyz + gupxy * gupxz) * Axy[idx]
        + (gupxx * gupzz + gupxz * gupxz) * Axz[idx]
        + (gupxy * gupzz + gupxz * gupyz) * Ayz[idx];
    Ryz[idx] = gupxy * gupxz * Axx[idx] + gupyy * gupyz * Ayy[idx] + gupyz * gupzz * Azz[idx]
        + (gupxy * gupyz + gupyy * gupxz) * Axy[idx]
        + (gupxy * gupzz + gupyz * gupxz) * Axz[idx]
        + (gupyy * gupzz + gupyz * gupyz) * Ayz[idx];
}

__global__ void bssn_shell_ricci_probe_kernel(
    int ex0,
    int ex1,
    int ex2,
    const double* crho,
    const double* sigma,
    const double* R,
    const double* drhodx,
    const double* drhody,
    const double* drhodz,
    const double* dsigmadx,
    const double* dsigmady,
    const double* dsigmadz,
    const double* dRdx,
    const double* dRdy,
    const double* dRdz,
    const double* drhodxx,
    const double* drhodxy,
    const double* drhodxz,
    const double* drhodyy,
    const double* drhodyz,
    const double* drhodzz,
    const double* dsigmadxx,
    const double* dsigmadxy,
    const double* dsigmadxz,
    const double* dsigmadyy,
    const double* dsigmadyz,
    const double* dsigmadzz,
    const double* dRdxx,
    const double* dRdxy,
    const double* dRdxz,
    const double* dRdyy,
    const double* dRdyz,
    const double* dRdzz,
    const double* dxx,
    const double* gxy,
    const double* gxz,
    const double* dyy,
    const double* gyz,
    const double* dzz,
    const double* Gamx,
    const double* Gamy,
    const double* Gamz,
    const double* Gamxxx,
    const double* Gamxxy,
    const double* Gamxxz,
    const double* Gamxyy,
    const double* Gamxyz,
    const double* Gamxzz,
    const double* Gamyxx,
    const double* Gamyxy,
    const double* Gamyxz,
    const double* Gamyyy,
    const double* Gamyyz,
    const double* Gamyzz,
    const double* Gamzxx,
    const double* Gamzxy,
    const double* Gamzxz,
    const double* Gamzyy,
    const double* Gamzyz,
    const double* Gamzzz,
    double* Ricci_xx,
    double* Ricci_xy,
    double* Ricci_xz,
    double* Ricci_yy,
    double* Ricci_yz,
    double* Ricci_zz,
    int symmetry,
    int sst
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int cell_count = ex0 * ex1 * ex2;
    if (idx >= cell_count) {
        return;
    }
    const int i = idx % ex0;
    const int j = (idx / ex0) % ex1;
    const int k = idx / (ex0 * ex1);
    const int ex[3] = {ex0, ex1, ex2};

    double dxx_xx = 0.0, dxx_xy = 0.0, dxx_xz = 0.0, dxx_yy = 0.0, dxx_yz = 0.0, dxx_zz = 0.0;
    double gxy_xx = 0.0, gxy_xy = 0.0, gxy_xz = 0.0, gxy_yy = 0.0, gxy_yz = 0.0, gxy_zz = 0.0;
    double gxz_xx = 0.0, gxz_xy = 0.0, gxz_xz = 0.0, gxz_yy = 0.0, gxz_yz = 0.0, gxz_zz = 0.0;
    double dyy_xx = 0.0, dyy_xy = 0.0, dyy_xz = 0.0, dyy_yy = 0.0, dyy_yz = 0.0, dyy_zz = 0.0;
    double gyz_xx = 0.0, gyz_xy = 0.0, gyz_xz = 0.0, gyz_yy = 0.0, gyz_yz = 0.0, gyz_zz = 0.0;
    double dzz_xx = 0.0, dzz_xy = 0.0, dzz_xz = 0.0, dzz_yy = 0.0, dzz_yz = 0.0, dzz_zz = 0.0;

    d_fdderivs_shc_point(ex, dxx, &dxx_xx, &dxx_xy, &dxx_xz, &dxx_yy, &dxx_yz, &dxx_zz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, 1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fdderivs_shc_point(ex, dyy, &dyy_xx, &dyy_xy, &dyy_xz, &dyy_yy, &dyy_yz, &dyy_zz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, 1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fdderivs_shc_point(ex, dzz, &dzz_xx, &dzz_xy, &dzz_xz, &dzz_yy, &dzz_yz, &dzz_zz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, 1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fdderivs_shc_point(ex, gxy, &gxy_xx, &gxy_xy, &gxy_xz, &gxy_yy, &gxy_yz, &gxy_zz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, -1.0, -1.0, 1.0, symmetry, sst, i, j, k);
    d_fdderivs_shc_point(ex, gxz, &gxz_xx, &gxz_xy, &gxz_xz, &gxz_yy, &gxz_yz, &gxz_zz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, -1.0, 1.0, -1.0, symmetry, sst, i, j, k);
    d_fdderivs_shc_point(ex, gyz, &gyz_xx, &gyz_xy, &gyz_xz, &gyz_yy, &gyz_yz, &gyz_zz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz, 1.0, -1.0, -1.0, symmetry, sst, i, j, k);

    double Gamxx = 0.0, Gamxy = 0.0, Gamxz = 0.0;
    double Gamyx = 0.0, Gamyy = 0.0, Gamyz = 0.0;
    double Gamzx = 0.0, Gamzy = 0.0, Gamzz = 0.0;
    d_fderivs_shc_point(ex, Gamx, &Gamxx, &Gamxy, &Gamxz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, -1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, Gamy, &Gamyx, &Gamyy, &Gamyz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, -1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, Gamz, &Gamzx, &Gamzy, &Gamzz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, 1.0, -1.0, symmetry, sst, i, j, k);

    const double l_gxx = dxx[idx] + 1.0;
    const double l_gyy = dyy[idx] + 1.0;
    const double l_gzz = dzz[idx] + 1.0;
    const double l_gxy = gxy[idx];
    const double l_gxz = gxz[idx];
    const double l_gyz = gyz[idx];
    const double det = l_gxx * l_gyy * l_gzz
        + l_gxy * l_gyz * l_gxz + l_gxz * l_gxy * l_gyz
        - l_gxz * l_gyy * l_gxz - l_gxy * l_gxy * l_gzz
        - l_gxx * l_gyz * l_gyz;
    const double gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) / det;
    const double gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) / det;
    const double gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) / det;
    const double gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) / det;
    const double gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) / det;
    const double gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) / det;

    const double l_Gamxxx = Gamxxx[idx], l_Gamxxy = Gamxxy[idx], l_Gamxxz = Gamxxz[idx];
    const double l_Gamxyy = Gamxyy[idx], l_Gamxyz = Gamxyz[idx], l_Gamxzz = Gamxzz[idx];
    const double l_Gamyxx = Gamyxx[idx], l_Gamyxy = Gamyxy[idx], l_Gamyxz = Gamyxz[idx];
    const double l_Gamyyy = Gamyyy[idx], l_Gamyyz = Gamyyz[idx], l_Gamyzz = Gamyzz[idx];
    const double l_Gamzxx = Gamzxx[idx], l_Gamzxy = Gamzxy[idx], l_Gamzxz = Gamzxz[idx];
    const double l_Gamzyy = Gamzyy[idx], l_Gamzyz = Gamzyz[idx], l_Gamzzz = Gamzzz[idx];

    const double Gamxa = gupxx * l_Gamxxx + gupyy * l_Gamxyy + gupzz * l_Gamxzz
        + 2.0 * (gupxy * l_Gamxxy + gupxz * l_Gamxxz + gupyz * l_Gamxyz);
    const double Gamya = gupxx * l_Gamyxx + gupyy * l_Gamyyy + gupzz * l_Gamyzz
        + 2.0 * (gupxy * l_Gamyxy + gupxz * l_Gamyxz + gupyz * l_Gamyyz);
    const double Gamza = gupxx * l_Gamzxx + gupyy * l_Gamzyy + gupzz * l_Gamzzz
        + 2.0 * (gupxy * l_Gamzxy + gupxz * l_Gamzxz + gupyz * l_Gamzyz);

    const double gxxx = l_gxx * l_Gamxxx + l_gxy * l_Gamyxx + l_gxz * l_Gamzxx;
    const double gxyx = l_gxx * l_Gamxxy + l_gxy * l_Gamyxy + l_gxz * l_Gamzxy;
    const double gxzx = l_gxx * l_Gamxxz + l_gxy * l_Gamyxz + l_gxz * l_Gamzxz;
    const double gyyx = l_gxx * l_Gamxyy + l_gxy * l_Gamyyy + l_gxz * l_Gamzyy;
    const double gyzx = l_gxx * l_Gamxyz + l_gxy * l_Gamyyz + l_gxz * l_Gamzyz;
    const double gzzx = l_gxx * l_Gamxzz + l_gxy * l_Gamyzz + l_gxz * l_Gamzzz;
    const double gxxy = l_gxy * l_Gamxxx + l_gyy * l_Gamyxx + l_gyz * l_Gamzxx;
    const double gxyy = l_gxy * l_Gamxxy + l_gyy * l_Gamyxy + l_gyz * l_Gamzxy;
    const double gxzy = l_gxy * l_Gamxxz + l_gyy * l_Gamyxz + l_gyz * l_Gamzxz;
    const double gyyy = l_gxy * l_Gamxyy + l_gyy * l_Gamyyy + l_gyz * l_Gamzyy;
    const double gyzy = l_gxy * l_Gamxyz + l_gyy * l_Gamyyz + l_gyz * l_Gamzyz;
    const double gzzy = l_gxy * l_Gamxzz + l_gyy * l_Gamyzz + l_gyz * l_Gamzzz;
    const double gxxz = l_gxz * l_Gamxxx + l_gyz * l_Gamyxx + l_gzz * l_Gamzxx;
    const double gxyz = l_gxz * l_Gamxxy + l_gyz * l_Gamyxy + l_gzz * l_Gamzxy;
    const double gxzz = l_gxz * l_Gamxxz + l_gyz * l_Gamyxz + l_gzz * l_Gamzxz;
    const double gyyz = l_gxz * l_Gamxyy + l_gyz * l_Gamyyy + l_gzz * l_Gamzyy;
    const double gyzz = l_gxz * l_Gamxyz + l_gyz * l_Gamyyz + l_gzz * l_Gamzyz;
    const double gzzz = l_gxz * l_Gamxzz + l_gyz * l_Gamyzz + l_gzz * l_Gamzzz;

    double Rxxv = d_shell_contract_second(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        dxx_xx, dxx_xy, dxx_xz, dxx_yy, dxx_yz, dxx_zz);
    double Ryyv = d_shell_contract_second(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        dyy_xx, dyy_xy, dyy_xz, dyy_yy, dyy_yz, dyy_zz);
    double Rzzv = d_shell_contract_second(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        dzz_xx, dzz_xy, dzz_xz, dzz_yy, dzz_yz, dzz_zz);
    double Rxyv = d_shell_contract_second(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        gxy_xx, gxy_xy, gxy_xz, gxy_yy, gxy_yz, gxy_zz);
    double Rxzv = d_shell_contract_second(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        gxz_xx, gxz_xy, gxz_xz, gxz_yy, gxz_yz, gxz_zz);
    double Ryzv = d_shell_contract_second(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        gyz_xx, gyz_xy, gyz_xz, gyz_yy, gyz_yz, gyz_zz);

    Rxxv = -0.5 * Rxxv + l_gxx * Gamxx + l_gxy * Gamyx + l_gxz * Gamzx
        + Gamxa * gxxx + Gamya * gxyx + Gamza * gxzx
        + gupxx * (2.0 * (l_Gamxxx * gxxx + l_Gamyxx * gxyx + l_Gamzxx * gxzx)
                   + l_Gamxxx * gxxx + l_Gamyxx * gxxy + l_Gamzxx * gxxz)
        + gupxy * (2.0 * (l_Gamxxx * gxyx + l_Gamyxx * gyyx + l_Gamzxx * gyzx
                          + l_Gamxxy * gxxx + l_Gamyxy * gxyx + l_Gamzxy * gxzx)
                   + l_Gamxxy * gxxx + l_Gamyxy * gxxy + l_Gamzxy * gxxz
                   + l_Gamxxx * gxyx + l_Gamyxx * gxyy + l_Gamzxx * gxyz)
        + gupxz * (2.0 * (l_Gamxxx * gxzx + l_Gamyxx * gyzx + l_Gamzxx * gzzx
                          + l_Gamxxz * gxxx + l_Gamyxz * gxyx + l_Gamzxz * gxzx)
                   + l_Gamxxz * gxxx + l_Gamyxz * gxxy + l_Gamzxz * gxxz
                   + l_Gamxxx * gxzx + l_Gamyxx * gxzy + l_Gamzxx * gxzz)
        + gupyy * (2.0 * (l_Gamxxy * gxyx + l_Gamyxy * gyyx + l_Gamzxy * gyzx)
                   + l_Gamxxy * gxyx + l_Gamyxy * gxyy + l_Gamzxy * gxyz)
        + gupyz * (2.0 * (l_Gamxxy * gxzx + l_Gamyxy * gyzx + l_Gamzxy * gzzx
                          + l_Gamxxz * gxyx + l_Gamyxz * gyyx + l_Gamzxz * gyzx)
                   + l_Gamxxz * gxyx + l_Gamyxz * gxyy + l_Gamzxz * gxyz
                   + l_Gamxxy * gxzx + l_Gamyxy * gxzy + l_Gamzxy * gxzz)
        + gupzz * (2.0 * (l_Gamxxz * gxzx + l_Gamyxz * gyzx + l_Gamzxz * gzzx)
                   + l_Gamxxz * gxzx + l_Gamyxz * gxzy + l_Gamzxz * gxzz);

    Ryyv = -0.5 * Ryyv + l_gxy * Gamxy + l_gyy * Gamyy + l_gyz * Gamzy
        + Gamxa * gxyy + Gamya * gyyy + Gamza * gyzy
        + gupxx * (2.0 * (l_Gamxxy * gxxy + l_Gamyxy * gxyy + l_Gamzxy * gxzy)
                   + l_Gamxxy * gxyx + l_Gamyxy * gxyy + l_Gamzxy * gxyz)
        + gupxy * (2.0 * (l_Gamxxy * gxyy + l_Gamyxy * gyyy + l_Gamzxy * gyzy
                          + l_Gamxyy * gxxy + l_Gamyyy * gxyy + l_Gamzyy * gxzy)
                   + l_Gamxyy * gxyx + l_Gamyyy * gxyy + l_Gamzyy * gxyz
                   + l_Gamxxy * gyyx + l_Gamyxy * gyyy + l_Gamzxy * gyyz)
        + gupxz * (2.0 * (l_Gamxxy * gxzy + l_Gamyxy * gyzy + l_Gamzxy * gzzy
                          + l_Gamxyz * gxxy + l_Gamyyz * gxyy + l_Gamzyz * gxzy)
                   + l_Gamxyz * gxyx + l_Gamyyz * gxyy + l_Gamzyz * gxyz
                   + l_Gamxxy * gyzx + l_Gamyxy * gyzy + l_Gamzxy * gyzz)
        + gupyy * (2.0 * (l_Gamxyy * gxyy + l_Gamyyy * gyyy + l_Gamzyy * gyzy)
                   + l_Gamxyy * gyyx + l_Gamyyy * gyyy + l_Gamzyy * gyyz)
        + gupyz * (2.0 * (l_Gamxyy * gxzy + l_Gamyyy * gyzy + l_Gamzyy * gzzy
                          + l_Gamxyz * gxyy + l_Gamyyz * gyyy + l_Gamzyz * gyzy)
                   + l_Gamxyz * gyyx + l_Gamyyz * gyyy + l_Gamzyz * gyyz
                   + l_Gamxyy * gyzx + l_Gamyyy * gyzy + l_Gamzyy * gyzz)
        + gupzz * (2.0 * (l_Gamxyz * gxzy + l_Gamyyz * gyzy + l_Gamzyz * gzzy)
                   + l_Gamxyz * gyzx + l_Gamyyz * gyzy + l_Gamzyz * gyzz);

    Rzzv = -0.5 * Rzzv + l_gxz * Gamxz + l_gyz * Gamyz + l_gzz * Gamzz
        + Gamxa * gxzz + Gamya * gyzz + Gamza * gzzz
        + gupxx * (2.0 * (l_Gamxxz * gxxz + l_Gamyxz * gxyz + l_Gamzxz * gxzz)
                   + l_Gamxxz * gxzx + l_Gamyxz * gxzy + l_Gamzxz * gxzz)
        + gupxy * (2.0 * (l_Gamxxz * gxyz + l_Gamyxz * gyyz + l_Gamzxz * gyzz
                          + l_Gamxyz * gxxz + l_Gamyyz * gxyz + l_Gamzyz * gxzz)
                   + l_Gamxyz * gxzx + l_Gamyyz * gxzy + l_Gamzyz * gxzz
                   + l_Gamxxz * gyzx + l_Gamyxz * gyzy + l_Gamzxz * gyzz)
        + gupxz * (2.0 * (l_Gamxxz * gxzz + l_Gamyxz * gyzz + l_Gamzxz * gzzz
                          + l_Gamxzz * gxxz + l_Gamyzz * gxyz + l_Gamzzz * gxzz)
                   + l_Gamxzz * gxzx + l_Gamyzz * gxzy + l_Gamzzz * gxzz
                   + l_Gamxxz * gzzx + l_Gamyxz * gzzy + l_Gamzxz * gzzz)
        + gupyy * (2.0 * (l_Gamxyz * gxyz + l_Gamyyz * gyyz + l_Gamzyz * gyzz)
                   + l_Gamxyz * gyzx + l_Gamyyz * gyzy + l_Gamzyz * gyzz)
        + gupyz * (2.0 * (l_Gamxyz * gxzz + l_Gamyyz * gyzz + l_Gamzyz * gzzz
                          + l_Gamxzz * gxyz + l_Gamyzz * gyyz + l_Gamzzz * gyzz)
                   + l_Gamxzz * gyzx + l_Gamyzz * gyzy + l_Gamzzz * gyzz
                   + l_Gamxyz * gzzx + l_Gamyyz * gzzy + l_Gamzyz * gzzz)
        + gupzz * (2.0 * (l_Gamxzz * gxzz + l_Gamyzz * gyzz + l_Gamzzz * gzzz)
                   + l_Gamxzz * gzzx + l_Gamyzz * gzzy + l_Gamzzz * gzzz);

    Rxyv = 0.5 * (-Rxyv + l_gxx * Gamxy + l_gxy * Gamyy + l_gxz * Gamzy
        + l_gxy * Gamxx + l_gyy * Gamyx + l_gyz * Gamzx
        + Gamxa * gxyx + Gamya * gyyx + Gamza * gyzx
        + Gamxa * gxxy + Gamya * gxyy + Gamza * gxzy)
        + gupxx * (l_Gamxxx * gxxy + l_Gamyxx * gxyy + l_Gamzxx * gxzy
                   + l_Gamxxy * gxxx + l_Gamyxy * gxyx + l_Gamzxy * gxzx
                   + l_Gamxxx * gxyx + l_Gamyxx * gxyy + l_Gamzxx * gxyz)
        + gupxy * (l_Gamxxx * gxyy + l_Gamyxx * gyyy + l_Gamzxx * gyzy
                   + l_Gamxxy * gxyx + l_Gamyxy * gyyx + l_Gamzxy * gyzx
                   + l_Gamxxy * gxyx + l_Gamyxy * gxyy + l_Gamzxy * gxyz
                   + l_Gamxxy * gxxy + l_Gamyxy * gxyy + l_Gamzxy * gxzy
                   + l_Gamxyy * gxxx + l_Gamyyy * gxyx + l_Gamzyy * gxzx
                   + l_Gamxxx * gyyx + l_Gamyxx * gyyy + l_Gamzxx * gyyz)
        + gupxz * (l_Gamxxx * gxzy + l_Gamyxx * gyzy + l_Gamzxx * gzzy
                   + l_Gamxxy * gxzx + l_Gamyxy * gyzx + l_Gamzxy * gzzx
                   + l_Gamxxz * gxyx + l_Gamyxz * gxyy + l_Gamzxz * gxyz
                   + l_Gamxxz * gxxy + l_Gamyxz * gxyy + l_Gamzxz * gxzy
                   + l_Gamxyz * gxxx + l_Gamyyz * gxyx + l_Gamzyz * gxzx
                   + l_Gamxxx * gyzx + l_Gamyxx * gyzy + l_Gamzxx * gyzz)
        + gupyy * (l_Gamxxy * gxyy + l_Gamyxy * gyyy + l_Gamzxy * gyzy
                   + l_Gamxyy * gxyx + l_Gamyyy * gyyx + l_Gamzyy * gyzx
                   + l_Gamxxy * gyyx + l_Gamyxy * gyyy + l_Gamzxy * gyyz)
        + gupyz * (l_Gamxxy * gxzy + l_Gamyxy * gyzy + l_Gamzxy * gzzy
                   + l_Gamxyy * gxzx + l_Gamyyy * gyzx + l_Gamzyy * gzzx
                   + l_Gamxxz * gyyx + l_Gamyxz * gyyy + l_Gamzxz * gyyz
                   + l_Gamxxz * gxyy + l_Gamyxz * gyyy + l_Gamzxz * gyzy
                   + l_Gamxyz * gxyx + l_Gamyyz * gyyx + l_Gamzyz * gyzx
                   + l_Gamxxy * gyzx + l_Gamyxy * gyzy + l_Gamzxy * gyzz)
        + gupzz * (l_Gamxxz * gxzy + l_Gamyxz * gyzy + l_Gamzxz * gzzy
                   + l_Gamxyz * gxzx + l_Gamyyz * gyzx + l_Gamzyz * gzzx
                   + l_Gamxxz * gyzx + l_Gamyxz * gyzy + l_Gamzxz * gyzz);

    Rxzv = 0.5 * (-Rxzv + l_gxx * Gamxz + l_gxy * Gamyz + l_gxz * Gamzz
        + l_gxz * Gamxx + l_gyz * Gamyx + l_gzz * Gamzx
        + Gamxa * gxzx + Gamya * gyzx + Gamza * gzzx
        + Gamxa * gxxz + Gamya * gxyz + Gamza * gxzz)
        + gupxx * (l_Gamxxx * gxxz + l_Gamyxx * gxyz + l_Gamzxx * gxzz
                   + l_Gamxxz * gxxx + l_Gamyxz * gxyx + l_Gamzxz * gxzx
                   + l_Gamxxx * gxzx + l_Gamyxx * gxzy + l_Gamzxx * gxzz)
        + gupxy * (l_Gamxxx * gxyz + l_Gamyxx * gyyz + l_Gamzxx * gyzz
                   + l_Gamxxz * gxyx + l_Gamyxz * gyyx + l_Gamzxz * gyzx
                   + l_Gamxxy * gxzx + l_Gamyxy * gxzy + l_Gamzxy * gxzz
                   + l_Gamxxy * gxxz + l_Gamyxy * gxyz + l_Gamzxy * gxzz
                   + l_Gamxyz * gxxx + l_Gamyyz * gxyx + l_Gamzyz * gxzx
                   + l_Gamxxx * gyzx + l_Gamyxx * gyzy + l_Gamzxx * gyzz)
        + gupxz * (l_Gamxxx * gxzz + l_Gamyxx * gyzz + l_Gamzxx * gzzz
                   + l_Gamxxz * gxzx + l_Gamyxz * gyzx + l_Gamzxz * gzzx
                   + l_Gamxxz * gxzx + l_Gamyxz * gxzy + l_Gamzxz * gxzz
                   + l_Gamxxz * gxxz + l_Gamyxz * gxyz + l_Gamzxz * gxzz
                   + l_Gamxzz * gxxx + l_Gamyzz * gxyx + l_Gamzzz * gxzx
                   + l_Gamxxx * gzzx + l_Gamyxx * gzzy + l_Gamzxx * gzzz)
        + gupyy * (l_Gamxxy * gxyz + l_Gamyxy * gyyz + l_Gamzxy * gyzz
                   + l_Gamxyz * gxyx + l_Gamyyz * gyyx + l_Gamzyz * gyzx
                   + l_Gamxxy * gyzx + l_Gamyxy * gyzy + l_Gamzxy * gyzz)
        + gupyz * (l_Gamxxy * gxzz + l_Gamyxy * gyzz + l_Gamzxy * gzzz
                   + l_Gamxyz * gxzx + l_Gamyyz * gyzx + l_Gamzyz * gzzx
                   + l_Gamxxz * gyzx + l_Gamyxz * gyzy + l_Gamzxz * gyzz
                   + l_Gamxxz * gxyz + l_Gamyxz * gyyz + l_Gamzxz * gyzz
                   + l_Gamxzz * gxyx + l_Gamyzz * gyyx + l_Gamzzz * gyzx
                   + l_Gamxxy * gzzx + l_Gamyxy * gzzy + l_Gamzxy * gzzz)
        + gupzz * (l_Gamxxz * gxzz + l_Gamyxz * gyzz + l_Gamzxz * gzzz
                   + l_Gamxzz * gxzx + l_Gamyzz * gyzx + l_Gamzzz * gzzx
                   + l_Gamxxz * gzzx + l_Gamyxz * gzzy + l_Gamzxz * gzzz);

    Ryzv = 0.5 * (-Ryzv + l_gxy * Gamxz + l_gyy * Gamyz + l_gyz * Gamzz
        + l_gxz * Gamxy + l_gyz * Gamyy + l_gzz * Gamzy
        + Gamxa * gxzy + Gamya * gyzy + Gamza * gzzy
        + Gamxa * gxyz + Gamya * gyyz + Gamza * gyzz)
        + gupxx * (l_Gamxxy * gxxz + l_Gamyxy * gxyz + l_Gamzxy * gxzz
                   + l_Gamxxz * gxxy + l_Gamyxz * gxyy + l_Gamzxz * gxzy
                   + l_Gamxxy * gxzx + l_Gamyxy * gxzy + l_Gamzxy * gxzz)
        + gupxy * (l_Gamxxy * gxyz + l_Gamyxy * gyyz + l_Gamzxy * gyzz
                   + l_Gamxxz * gxyy + l_Gamyxz * gyyy + l_Gamzxz * gyzy
                   + l_Gamxyy * gxzx + l_Gamyyy * gxzy + l_Gamzyy * gxzz
                   + l_Gamxyy * gxxz + l_Gamyyy * gxyz + l_Gamzyy * gxzz
                   + l_Gamxyz * gxxy + l_Gamyyz * gxyy + l_Gamzyz * gxzy
                   + l_Gamxxy * gyzx + l_Gamyxy * gyzy + l_Gamzxy * gyzz)
        + gupxz * (l_Gamxxy * gxzz + l_Gamyxy * gyzz + l_Gamzxy * gzzz
                   + l_Gamxxz * gxzy + l_Gamyxz * gyzy + l_Gamzxz * gzzy
                   + l_Gamxyz * gxzx + l_Gamyyz * gxzy + l_Gamzyz * gxzz
                   + l_Gamxyz * gxxz + l_Gamyyz * gxyz + l_Gamzyz * gxzz
                   + l_Gamxzz * gxxy + l_Gamyzz * gxyy + l_Gamzzz * gxzy
                   + l_Gamxxy * gzzx + l_Gamyxy * gzzy + l_Gamzxy * gzzz)
        + gupyy * (l_Gamxyy * gxyz + l_Gamyyy * gyyz + l_Gamzyy * gyzz
                   + l_Gamxyz * gxyy + l_Gamyyz * gyyy + l_Gamzyz * gyzy
                   + l_Gamxyy * gyzx + l_Gamyyy * gyzy + l_Gamzyy * gyzz)
        + gupyz * (l_Gamxyy * gxzz + l_Gamyyy * gyzz + l_Gamzyy * gzzz
                   + l_Gamxyz * gxzy + l_Gamyyz * gyzy + l_Gamzyz * gzzy
                   + l_Gamxyz * gyzx + l_Gamyyz * gyzy + l_Gamzyz * gyzz
                   + l_Gamxyz * gxyz + l_Gamyyz * gyyz + l_Gamzyz * gyzz
                   + l_Gamxzz * gxyy + l_Gamyzz * gyyy + l_Gamzzz * gyzy
                   + l_Gamxyy * gzzx + l_Gamyyy * gzzy + l_Gamzyy * gzzz)
        + gupzz * (l_Gamxyz * gxzz + l_Gamyyz * gyzz + l_Gamzyz * gzzz
                   + l_Gamxzz * gxzy + l_Gamyzz * gyzy + l_Gamzzz * gzzy
                   + l_Gamxyz * gzzx + l_Gamyyz * gzzy + l_Gamzyz * gzzz);

    Ricci_xx[idx] = Rxxv;
    Ricci_xy[idx] = Rxyv;
    Ricci_xz[idx] = Rxzv;
    Ricci_yy[idx] = Ryyv;
    Ricci_yz[idx] = Ryzv;
    Ricci_zz[idx] = Rzzv;
}

__global__ void bssn_shell_gauge_rhs_probe_kernel(
    int ex0,
    int ex1,
    int ex2,
    const double* crho,
    const double* sigma,
    const double* R,
    const double* x,
    const double* y,
    const double* z,
    const double* drhodx,
    const double* drhody,
    const double* drhodz,
    const double* dsigmadx,
    const double* dsigmady,
    const double* dsigmadz,
    const double* dRdx,
    const double* dRdy,
    const double* dRdz,
    const double* drhodxx,
    const double* drhodxy,
    const double* drhodxz,
    const double* drhodyy,
    const double* drhodyz,
    const double* drhodzz,
    const double* dsigmadxx,
    const double* dsigmadxy,
    const double* dsigmadxz,
    const double* dsigmadyy,
    const double* dsigmadyz,
    const double* dsigmadzz,
    const double* dRdxx,
    const double* dRdxy,
    const double* dRdxz,
    const double* dRdyy,
    const double* dRdyz,
    const double* dRdzz,
    const double* chi,
    const double* trK,
    const double* dxx,
    const double* gxy,
    const double* gxz,
    const double* dyy,
    const double* gyz,
    const double* dzz,
    const double* Gamx,
    const double* Gamy,
    const double* Gamz,
    const double* Lap,
    const double* betax,
    const double* betay,
    const double* betaz,
    const double* dtSfx,
    const double* dtSfy,
    const double* dtSfz,
    const double* Sx,
    const double* Sy,
    const double* Sz,
    const double* Gamxxx,
    const double* Gamxxy,
    const double* Gamxxz,
    const double* Gamxyy,
    const double* Gamxyz,
    const double* Gamxzz,
    const double* Gamyxx,
    const double* Gamyxy,
    const double* Gamyxz,
    const double* Gamyyy,
    const double* Gamyyz,
    const double* Gamyzz,
    const double* Gamzxx,
    const double* Gamzxy,
    const double* Gamzxz,
    const double* Gamzyy,
    const double* Gamzyz,
    const double* Gamzzz,
    const double* Rxx,
    const double* Rxy,
    const double* Rxz,
    const double* Ryy,
    const double* Ryz,
    const double* Rzz,
    double* Gamx_rhs,
    double* Gamy_rhs,
    double* Gamz_rhs,
    double* Lap_rhs,
    double* betax_rhs,
    double* betay_rhs,
    double* betaz_rhs,
    double* dtSfx_rhs,
    double* dtSfy_rhs,
    double* dtSfz_rhs,
    int symmetry,
    int sst
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int cell_count = ex0 * ex1 * ex2;
    if (idx >= cell_count) {
        return;
    }
    const int i = idx % ex0;
    const int j = (idx / ex0) % ex1;
    const int k = idx / (ex0 * ex1);
    const int ex[3] = {ex0, ex1, ex2};

    double betaxx = 0.0, betaxy = 0.0, betaxz = 0.0;
    double betayx = 0.0, betayy = 0.0, betayz = 0.0;
    double betazx = 0.0, betazy = 0.0, betazz = 0.0;
    d_fderivs_shc_point(
        ex, betax, &betaxx, &betaxy, &betaxz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        -1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, betay, &betayx, &betayy, &betayz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, -1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, betaz, &betazx, &betazy, &betazz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, -1.0, symmetry, sst, i, j, k);

    double bx_gxxx = 0.0, bx_gxyx = 0.0, bx_gxzx = 0.0, bx_gyyx = 0.0, bx_gyzx = 0.0, bx_gzzx = 0.0;
    double by_gxxy = 0.0, by_gxyy = 0.0, by_gxzy = 0.0, by_gyyy = 0.0, by_gyzy = 0.0, by_gzzy = 0.0;
    double bz_gxxz = 0.0, bz_gxyz = 0.0, bz_gxzz = 0.0, bz_gyyz = 0.0, bz_gyzz = 0.0, bz_gzzz = 0.0;
    d_fdderivs_shc_point(
        ex, betax, &bx_gxxx, &bx_gxyx, &bx_gxzx, &bx_gyyx, &bx_gyzx, &bx_gzzx,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        -1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fdderivs_shc_point(
        ex, betay, &by_gxxy, &by_gxyy, &by_gxzy, &by_gyyy, &by_gyzy, &by_gzzy,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        1.0, -1.0, 1.0, symmetry, sst, i, j, k);
    d_fdderivs_shc_point(
        ex, betaz, &bz_gxxz, &bz_gxyz, &bz_gxzz, &bz_gyyz, &bz_gyzz, &bz_gzzz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        1.0, 1.0, -1.0, symmetry, sst, i, j, k);

    double Lapx = 0.0, Lapy = 0.0, Lapz = 0.0;
    double chix = 0.0, chiy = 0.0, chiz = 0.0;
    double Kx = 0.0, Ky = 0.0, Kz = 0.0;
    d_fderivs_shc_point(
        ex, Lap, &Lapx, &Lapy, &Lapz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, chi, &chix, &chiy, &chiz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, trK, &Kx, &Ky, &Kz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, 1.0, symmetry, sst, i, j, k);

    const double l_gxx = dxx[idx] + 1.0;
    const double l_gyy = dyy[idx] + 1.0;
    const double l_gzz = dzz[idx] + 1.0;
    const double l_gxy = gxy[idx];
    const double l_gxz = gxz[idx];
    const double l_gyz = gyz[idx];
    const double det = l_gxx * l_gyy * l_gzz
        + l_gxy * l_gyz * l_gxz + l_gxz * l_gxy * l_gyz
        - l_gxz * l_gyy * l_gxz - l_gxy * l_gxy * l_gzz
        - l_gxx * l_gyz * l_gyz;
    const double gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) / det;
    const double gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) / det;
    const double gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) / det;
    const double gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) / det;
    const double gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) / det;
    const double gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) / det;

    const double alpn1 = Lap[idx] + 1.0;
    const double chin1 = chi[idx] + 1.0;
    const double inv_chin1 = 1.0 / chin1;
    const double div_beta = betaxx + betayy + betazz;

    const double Aup_xx = Rxx[idx];
    const double Aup_xy = Rxy[idx];
    const double Aup_xz = Rxz[idx];
    const double Aup_yy = Ryy[idx];
    const double Aup_yz = Ryz[idx];
    const double Aup_zz = Rzz[idx];
    const double l_Gamxxx = Gamxxx[idx], l_Gamxxy = Gamxxy[idx], l_Gamxxz = Gamxxz[idx];
    const double l_Gamxyy = Gamxyy[idx], l_Gamxyz = Gamxyz[idx], l_Gamxzz = Gamxzz[idx];
    const double l_Gamyxx = Gamyxx[idx], l_Gamyxy = Gamyxy[idx], l_Gamyxz = Gamyxz[idx];
    const double l_Gamyyy = Gamyyy[idx], l_Gamyyz = Gamyyz[idx], l_Gamyzz = Gamyzz[idx];
    const double l_Gamzxx = Gamzxx[idx], l_Gamzxy = Gamzxy[idx], l_Gamzxz = Gamzxz[idx];
    const double l_Gamzyy = Gamzyy[idx], l_Gamzyz = Gamzyz[idx], l_Gamzzz = Gamzzz[idx];

    double val_Gamx_rhs = -2.0 * (Lapx * Aup_xx + Lapy * Aup_xy + Lapz * Aup_xz)
        + 2.0 * alpn1 * (
            -kF3o2 * inv_chin1 * (chix * Aup_xx + chiy * Aup_xy + chiz * Aup_xz)
            - gupxx * (kF2o3 * Kx + kEightPi * Sx[idx])
            - gupxy * (kF2o3 * Ky + kEightPi * Sy[idx])
            - gupxz * (kF2o3 * Kz + kEightPi * Sz[idx])
            + l_Gamxxx * Aup_xx + l_Gamxyy * Aup_yy + l_Gamxzz * Aup_zz
            + 2.0 * (l_Gamxxy * Aup_xy + l_Gamxxz * Aup_xz + l_Gamxyz * Aup_yz));
    double val_Gamy_rhs = -2.0 * (Lapx * Aup_xy + Lapy * Aup_yy + Lapz * Aup_yz)
        + 2.0 * alpn1 * (
            -kF3o2 * inv_chin1 * (chix * Aup_xy + chiy * Aup_yy + chiz * Aup_yz)
            - gupxy * (kF2o3 * Kx + kEightPi * Sx[idx])
            - gupyy * (kF2o3 * Ky + kEightPi * Sy[idx])
            - gupyz * (kF2o3 * Kz + kEightPi * Sz[idx])
            + l_Gamyxx * Aup_xx + l_Gamyyy * Aup_yy + l_Gamyzz * Aup_zz
            + 2.0 * (l_Gamyxy * Aup_xy + l_Gamyxz * Aup_xz + l_Gamyyz * Aup_yz));
    double val_Gamz_rhs = -2.0 * (Lapx * Aup_xz + Lapy * Aup_yz + Lapz * Aup_zz)
        + 2.0 * alpn1 * (
            -kF3o2 * inv_chin1 * (chix * Aup_xz + chiy * Aup_yz + chiz * Aup_zz)
            - gupxz * (kF2o3 * Kx + kEightPi * Sx[idx])
            - gupyz * (kF2o3 * Ky + kEightPi * Sy[idx])
            - gupzz * (kF2o3 * Kz + kEightPi * Sz[idx])
            + l_Gamzxx * Aup_xx + l_Gamzyy * Aup_yy + l_Gamzzz * Aup_zz
            + 2.0 * (l_Gamzxy * Aup_xy + l_Gamzxz * Aup_xz + l_Gamzyz * Aup_yz));

    const double fxx_beta = bx_gxxx + by_gxyy + bz_gxzz;
    const double fxy_beta = bx_gxyx + by_gyyy + bz_gyzz;
    const double fxz_beta = bx_gxzx + by_gyzy + bz_gzzz;
    const double Gamxa = gupxx * l_Gamxxx + gupyy * l_Gamxyy + gupzz * l_Gamxzz
        + 2.0 * (gupxy * l_Gamxxy + gupxz * l_Gamxxz + gupyz * l_Gamxyz);
    const double Gamya = gupxx * l_Gamyxx + gupyy * l_Gamyyy + gupzz * l_Gamyzz
        + 2.0 * (gupxy * l_Gamyxy + gupxz * l_Gamyxz + gupyz * l_Gamyyz);
    const double Gamza = gupxx * l_Gamzxx + gupyy * l_Gamzyy + gupzz * l_Gamzzz
        + 2.0 * (gupxy * l_Gamzxy + gupxz * l_Gamzxz + gupyz * l_Gamzyz);

    val_Gamx_rhs += kF2o3 * Gamxa * div_beta
        - Gamxa * betaxx - Gamya * betaxy - Gamza * betaxz
        + kF1o3 * (gupxx * fxx_beta + gupxy * fxy_beta + gupxz * fxz_beta)
        + gupxx * bx_gxxx + gupyy * bx_gyyx + gupzz * bx_gzzx
        + 2.0 * (gupxy * bx_gxyx + gupxz * bx_gxzx + gupyz * bx_gyzx);
    val_Gamy_rhs += kF2o3 * Gamya * div_beta
        - Gamxa * betayx - Gamya * betayy - Gamza * betayz
        + kF1o3 * (gupxy * fxx_beta + gupyy * fxy_beta + gupyz * fxz_beta)
        + gupxx * by_gxxy + gupyy * by_gyyy + gupzz * by_gzzy
        + 2.0 * (gupxy * by_gxyy + gupxz * by_gxzy + gupyz * by_gyzy);
    val_Gamz_rhs += kF2o3 * Gamza * div_beta
        - Gamxa * betazx - Gamya * betazy - Gamza * betazz
        + kF1o3 * (gupxz * fxx_beta + gupyz * fxy_beta + gupzz * fxz_beta)
        + gupxx * bz_gxxz + gupyy * bz_gyyz + gupzz * bz_gzzz
        + 2.0 * (gupxy * bz_gxyz + gupxz * bz_gxzz + gupyz * bz_gyzz);

    Lap_rhs[idx] = -2.0 * alpn1 * trK[idx];
#if (GAUGE == 0)
    betax_rhs[idx] = kFF * dtSfx[idx];
    betay_rhs[idx] = kFF * dtSfy[idx];
    betaz_rhs[idx] = kFF * dtSfz[idx];
    dtSfx_rhs[idx] = val_Gamx_rhs - kEta * dtSfx[idx];
    dtSfy_rhs[idx] = val_Gamy_rhs - kEta * dtSfy[idx];
    dtSfz_rhs[idx] = val_Gamz_rhs - kEta * dtSfz[idx];
#elif (GAUGE == 1)
    betax_rhs[idx] = Gamx[idx] - kEta * betax[idx];
    betay_rhs[idx] = Gamy[idx] - kEta * betay[idx];
    betaz_rhs[idx] = Gamz[idx] - kEta * betaz[idx];
    dtSfx_rhs[idx] = 0.0;
    dtSfy_rhs[idx] = 0.0;
    dtSfz_rhs[idx] = 0.0;
#elif (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5)
    double reta = gupxx * chix * chix + gupyy * chiy * chiy + gupzz * chiz * chiz
        + 2.0 * (gupxy * chix * chiy + gupxz * chix * chiz + gupyz * chiy * chiz);
#if (GAUGE == 2 || GAUGE == 4)
    const double reta_denom = 1.0 - sqrt(chin1);
#else
    const double reta_denom = 1.0 - chin1;
#endif
    reta = 1.31 / 2.0 * sqrt(reta / chin1) / (reta_denom * reta_denom);
#if (GAUGE == 2 || GAUGE == 3)
    betax_rhs[idx] = kFF * dtSfx[idx];
    betay_rhs[idx] = kFF * dtSfy[idx];
    betaz_rhs[idx] = kFF * dtSfz[idx];
    dtSfx_rhs[idx] = val_Gamx_rhs - reta * dtSfx[idx];
    dtSfy_rhs[idx] = val_Gamy_rhs - reta * dtSfy[idx];
    dtSfz_rhs[idx] = val_Gamz_rhs - reta * dtSfz[idx];
#else
    betax_rhs[idx] = kFF * Gamx[idx] - reta * betax[idx];
    betay_rhs[idx] = kFF * Gamy[idx] - reta * betay[idx];
    betaz_rhs[idx] = kFF * Gamz[idx] - reta * betaz[idx];
    dtSfx_rhs[idx] = 0.0;
    dtSfy_rhs[idx] = 0.0;
    dtSfz_rhs[idx] = 0.0;
#endif
#elif (GAUGE == 6 || GAUGE == 7)
    double reta = kEta;
    if (d_bssn_shell_gauge_bh_count == 2) {
        const double total_mass = d_bssn_shell_gauge_mass[0] + d_bssn_shell_gauge_mass[1];
        const double A = 2.0 / total_mass;
        const double w1 = 12.0;
        const double w2 = w1;
        const double C1 = 1.0 / d_bssn_shell_gauge_mass[0] - A;
        const double C2 = 1.0 / d_bssn_shell_gauge_mass[1] - A;
        const double xp = x[idx];
        const double yp = y[idx];
        const double zp = z[idx];
        const double dx12 = d_bssn_shell_gauge_pbh[0] - d_bssn_shell_gauge_pbh[3];
        const double dy12 = d_bssn_shell_gauge_pbh[1] - d_bssn_shell_gauge_pbh[4];
        const double dz12 = d_bssn_shell_gauge_pbh[2] - d_bssn_shell_gauge_pbh[5];
        const double sep2 = dx12 * dx12 + dy12 * dy12 + dz12 * dz12;
        const double dx1 = d_bssn_shell_gauge_pbh[0] - xp;
        const double dy1 = d_bssn_shell_gauge_pbh[1] - yp;
        const double dz1 = d_bssn_shell_gauge_pbh[2] - zp;
        const double dx2 = d_bssn_shell_gauge_pbh[3] - xp;
        const double dy2 = d_bssn_shell_gauge_pbh[4] - yp;
        const double dz2 = d_bssn_shell_gauge_pbh[5] - zp;
        const double r1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1) / sep2;
        const double r2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2) / sep2;
#if (GAUGE == 6)
        reta = A + C1 / (1.0 + w1 * r1) + C2 / (1.0 + w2 * r2);
#else
        reta = A + C1 * exp(-w1 * r1) + C2 * exp(-w2 * r2);
#endif
    }
    betax_rhs[idx] = kFF * dtSfx[idx];
    betay_rhs[idx] = kFF * dtSfy[idx];
    betaz_rhs[idx] = kFF * dtSfz[idx];
    dtSfx_rhs[idx] = val_Gamx_rhs - reta * dtSfx[idx];
    dtSfy_rhs[idx] = val_Gamy_rhs - reta * dtSfy[idx];
    dtSfz_rhs[idx] = val_Gamz_rhs - reta * dtSfz[idx];
#else
#error Unsupported Shell-Patch gauge RHS probe; expected GAUGE == 0..7.
#endif
    Gamx_rhs[idx] = val_Gamx_rhs;
    Gamy_rhs[idx] = val_Gamy_rhs;
    Gamz_rhs[idx] = val_Gamz_rhs;
}

__global__ void bssn_shell_curvature_rhs_probe_kernel(
    int ex0,
    int ex1,
    int ex2,
    const double* crho,
    const double* sigma,
    const double* R,
    const double* drhodx,
    const double* drhody,
    const double* drhodz,
    const double* dsigmadx,
    const double* dsigmady,
    const double* dsigmadz,
    const double* dRdx,
    const double* dRdy,
    const double* dRdz,
    const double* drhodxx,
    const double* drhodxy,
    const double* drhodxz,
    const double* drhodyy,
    const double* drhodyz,
    const double* drhodzz,
    const double* dsigmadxx,
    const double* dsigmadxy,
    const double* dsigmadxz,
    const double* dsigmadyy,
    const double* dsigmadyz,
    const double* dsigmadzz,
    const double* dRdxx,
    const double* dRdxy,
    const double* dRdxz,
    const double* dRdyy,
    const double* dRdyz,
    const double* dRdzz,
    const double* chi,
    const double* trK,
    const double* dxx,
    const double* gxy,
    const double* gxz,
    const double* dyy,
    const double* gyz,
    const double* dzz,
    const double* Axx,
    const double* Axy,
    const double* Axz,
    const double* Ayy,
    const double* Ayz,
    const double* Azz,
    const double* Lap,
    const double* betax,
    const double* betay,
    const double* betaz,
    const double* rho,
    const double* Sxx,
    const double* Sxy,
    const double* Sxz,
    const double* Syy,
    const double* Syz,
    const double* Szz,
    const double* Gamxxx,
    const double* Gamxxy,
    const double* Gamxxz,
    const double* Gamxyy,
    const double* Gamxyz,
    const double* Gamxzz,
    const double* Gamyxx,
    const double* Gamyxy,
    const double* Gamyxz,
    const double* Gamyyy,
    const double* Gamyyz,
    const double* Gamyzz,
    const double* Gamzxx,
    const double* Gamzxy,
    const double* Gamzxz,
    const double* Gamzyy,
    const double* Gamzyz,
    const double* Gamzzz,
    const double* Ricci_xx,
    const double* Ricci_xy,
    const double* Ricci_xz,
    const double* Ricci_yy,
    const double* Ricci_yz,
    const double* Ricci_zz,
    double* phys_Rxx,
    double* phys_Rxy,
    double* phys_Rxz,
    double* phys_Ryy,
    double* phys_Ryz,
    double* phys_Rzz,
    double* trK_rhs,
    double* Axx_rhs,
    double* Axy_rhs,
    double* Axz_rhs,
    double* Ayy_rhs,
    double* Ayz_rhs,
    double* Azz_rhs,
    int symmetry,
    int sst
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int cell_count = ex0 * ex1 * ex2;
    if (idx >= cell_count) {
        return;
    }
    const int i = idx % ex0;
    const int j = (idx / ex0) % ex1;
    const int k = idx / (ex0 * ex1);
    const int ex[3] = {ex0, ex1, ex2};

    double betaxx = 0.0, betaxy = 0.0, betaxz = 0.0;
    double betayx = 0.0, betayy = 0.0, betayz = 0.0;
    double betazx = 0.0, betazy = 0.0, betazz = 0.0;
    d_fderivs_shc_point(
        ex, betax, &betaxx, &betaxy, &betaxz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        -1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, betay, &betayx, &betayy, &betayz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, -1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(
        ex, betaz, &betazx, &betazy, &betazz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, -1.0, symmetry, sst, i, j, k);

    double Lapx = 0.0, Lapy = 0.0, Lapz = 0.0;
    double Lapxx = 0.0, Lapxy = 0.0, Lapxz = 0.0;
    double Lapyy = 0.0, Lapyz = 0.0, Lapzz = 0.0;
    d_fderivs_shc_point(
        ex, Lap, &Lapx, &Lapy, &Lapz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fdderivs_shc_point(
        ex, Lap, &Lapxx, &Lapxy, &Lapxz, &Lapyy, &Lapyz, &Lapzz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        1.0, 1.0, 1.0, symmetry, sst, i, j, k);

    double chix = 0.0, chiy = 0.0, chiz = 0.0;
    double chixx = 0.0, chixy = 0.0, chixz = 0.0;
    double chiyy = 0.0, chiyz = 0.0, chizz = 0.0;
    d_fderivs_shc_point(
        ex, chi, &chix, &chiy, &chiz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fdderivs_shc_point(
        ex, chi, &chixx, &chixy, &chixz, &chiyy, &chiyz, &chizz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        1.0, 1.0, 1.0, symmetry, sst, i, j, k);

    const double l_gxx = dxx[idx] + 1.0;
    const double l_gyy = dyy[idx] + 1.0;
    const double l_gzz = dzz[idx] + 1.0;
    const double l_gxy = gxy[idx];
    const double l_gxz = gxz[idx];
    const double l_gyz = gyz[idx];
    const double det = l_gxx * l_gyy * l_gzz
        + l_gxy * l_gyz * l_gxz + l_gxz * l_gxy * l_gyz
        - l_gxz * l_gyy * l_gxz - l_gxy * l_gxy * l_gzz
        - l_gxx * l_gyz * l_gyz;
    const double gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) / det;
    const double gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) / det;
    const double gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) / det;
    const double gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) / det;
    const double gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) / det;
    const double gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) / det;

    const double alpn1 = Lap[idx] + 1.0;
    const double chin1 = chi[idx] + 1.0;
    const double inv_chin1 = 1.0 / chin1;
    const double div_beta = betaxx + betayy + betazz;

    const double l_Gamxxx = Gamxxx[idx], l_Gamxxy = Gamxxy[idx], l_Gamxxz = Gamxxz[idx];
    const double l_Gamxyy = Gamxyy[idx], l_Gamxyz = Gamxyz[idx], l_Gamxzz = Gamxzz[idx];
    const double l_Gamyxx = Gamyxx[idx], l_Gamyxy = Gamyxy[idx], l_Gamyxz = Gamyxz[idx];
    const double l_Gamyyy = Gamyyy[idx], l_Gamyyz = Gamyyz[idx], l_Gamyzz = Gamyzz[idx];
    const double l_Gamzxx = Gamzxx[idx], l_Gamzxy = Gamzxy[idx], l_Gamzxz = Gamzxz[idx];
    const double l_Gamzyy = Gamzyy[idx], l_Gamzyz = Gamzyz[idx], l_Gamzzz = Gamzzz[idx];

    double cxx = chixx - l_Gamxxx * chix - l_Gamyxx * chiy - l_Gamzxx * chiz;
    double cxy = chixy - l_Gamxxy * chix - l_Gamyxy * chiy - l_Gamzxy * chiz;
    double cxz = chixz - l_Gamxxz * chix - l_Gamyxz * chiy - l_Gamzxz * chiz;
    double cyy = chiyy - l_Gamxyy * chix - l_Gamyyy * chiy - l_Gamzyy * chiz;
    double cyz = chiyz - l_Gamxyz * chix - l_Gamyyz * chiy - l_Gamzyz * chiz;
    double czz = chizz - l_Gamxzz * chix - l_Gamyzz * chiy - l_Gamzzz * chiz;
    const double chi_scalar = gupxx * (cxx - kF3o2 * inv_chin1 * chix * chix)
        + gupyy * (cyy - kF3o2 * inv_chin1 * chiy * chiy)
        + gupzz * (czz - kF3o2 * inv_chin1 * chiz * chiz)
        + 2.0 * (gupxy * (cxy - kF3o2 * inv_chin1 * chix * chiy)
               + gupxz * (cxz - kF3o2 * inv_chin1 * chix * chiz)
               + gupyz * (cyz - kF3o2 * inv_chin1 * chiy * chiz));

    const double Rphys_xx = Ricci_xx[idx] + (cxx - 0.5 * chix * chix * inv_chin1 + l_gxx * chi_scalar) * inv_chin1 * 0.5;
    const double Rphys_yy = Ricci_yy[idx] + (cyy - 0.5 * chiy * chiy * inv_chin1 + l_gyy * chi_scalar) * inv_chin1 * 0.5;
    const double Rphys_zz = Ricci_zz[idx] + (czz - 0.5 * chiz * chiz * inv_chin1 + l_gzz * chi_scalar) * inv_chin1 * 0.5;
    const double Rphys_xy = Ricci_xy[idx] + (cxy - 0.5 * chix * chiy * inv_chin1 + l_gxy * chi_scalar) * inv_chin1 * 0.5;
    const double Rphys_xz = Ricci_xz[idx] + (cxz - 0.5 * chix * chiz * inv_chin1 + l_gxz * chi_scalar) * inv_chin1 * 0.5;
    const double Rphys_yz = Ricci_yz[idx] + (cyz - 0.5 * chiy * chiz * inv_chin1 + l_gyz * chi_scalar) * inv_chin1 * 0.5;

    const double gx_phy = (gupxx * chix + gupxy * chiy + gupxz * chiz) * inv_chin1;
    const double gy_phy = (gupxy * chix + gupyy * chiy + gupyz * chiz) * inv_chin1;
    const double gz_phy = (gupxz * chix + gupyz * chiy + gupzz * chiz) * inv_chin1;
    const double chi_dot_lap = gx_phy * Lapx + gy_phy * Lapy + gz_phy * Lapz;

    const double D2_lap_xx = Lapxx - l_Gamxxx * Lapx - l_Gamyxx * Lapy - l_Gamzxx * Lapz;
    const double D2_lap_xy = Lapxy - l_Gamxxy * Lapx - l_Gamyxy * Lapy - l_Gamzxy * Lapz;
    const double D2_lap_xz = Lapxz - l_Gamxxz * Lapx - l_Gamyxz * Lapy - l_Gamzxz * Lapz;
    const double D2_lap_yy = Lapyy - l_Gamxyy * Lapx - l_Gamyyy * Lapy - l_Gamzyy * Lapz;
    const double D2_lap_yz = Lapyz - l_Gamxyz * Lapx - l_Gamyyz * Lapy - l_Gamzyz * Lapz;
    const double D2_lap_zz = Lapzz - l_Gamxzz * Lapx - l_Gamyzz * Lapy - l_Gamzzz * Lapz;

    const double lap_phy_xx = D2_lap_xx + chix * Lapx * inv_chin1 - 0.5 * l_gxx * chi_dot_lap;
    const double lap_phy_yy = D2_lap_yy + chiy * Lapy * inv_chin1 - 0.5 * l_gyy * chi_dot_lap;
    const double lap_phy_zz = D2_lap_zz + chiz * Lapz * inv_chin1 - 0.5 * l_gzz * chi_dot_lap;
    const double lap_phy_xy = D2_lap_xy + 0.5 * (chix * Lapy + chiy * Lapx) * inv_chin1 - 0.5 * l_gxy * chi_dot_lap;
    const double lap_phy_xz = D2_lap_xz + 0.5 * (chix * Lapz + chiz * Lapx) * inv_chin1 - 0.5 * l_gxz * chi_dot_lap;
    const double lap_phy_yz = D2_lap_yz + 0.5 * (chiy * Lapz + chiz * Lapy) * inv_chin1 - 0.5 * l_gyz * chi_dot_lap;
    const double lap_trace = gupxx * lap_phy_xx + gupyy * lap_phy_yy + gupzz * lap_phy_zz
        + 2.0 * (gupxy * lap_phy_xy + gupxz * lap_phy_xz + gupyz * lap_phy_yz);

    const double src_xx = alpn1 * (Rphys_xx - kEightPi * Sxx[idx]) - lap_phy_xx;
    const double src_xy = alpn1 * (Rphys_xy - kEightPi * Sxy[idx]) - lap_phy_xy;
    const double src_xz = alpn1 * (Rphys_xz - kEightPi * Sxz[idx]) - lap_phy_xz;
    const double src_yy = alpn1 * (Rphys_yy - kEightPi * Syy[idx]) - lap_phy_yy;
    const double src_yz = alpn1 * (Rphys_yz - kEightPi * Syz[idx]) - lap_phy_yz;
    const double src_zz = alpn1 * (Rphys_zz - kEightPi * Szz[idx]) - lap_phy_zz;
    const double trace_free = kF1o3 * (gupxx * src_xx + gupyy * src_yy + gupzz * src_zz
        + 2.0 * (gupxy * src_xy + gupxz * src_xz + gupyz * src_yz));

    const double l_Axx = Axx[idx], l_Axy = Axy[idx], l_Axz = Axz[idx];
    const double l_Ayy = Ayy[idx], l_Ayz = Ayz[idx], l_Azz = Azz[idx];
    const double term_xx = gupxx * l_Axx * l_Axx + gupyy * l_Axy * l_Axy + gupzz * l_Axz * l_Axz
        + 2.0 * (gupxy * l_Axx * l_Axy + gupxz * l_Axx * l_Axz + gupyz * l_Axy * l_Axz);
    const double term_yy = gupxx * l_Axy * l_Axy + gupyy * l_Ayy * l_Ayy + gupzz * l_Ayz * l_Ayz
        + 2.0 * (gupxy * l_Axy * l_Ayy + gupxz * l_Axy * l_Ayz + gupyz * l_Ayy * l_Ayz);
    const double term_zz = gupxx * l_Axz * l_Axz + gupyy * l_Ayz * l_Ayz + gupzz * l_Azz * l_Azz
        + 2.0 * (gupxy * l_Axz * l_Ayz + gupxz * l_Axz * l_Azz + gupyz * l_Ayz * l_Azz);
    const double term_xy = gupxx * l_Axx * l_Axy + gupyy * l_Axy * l_Ayy + gupzz * l_Axz * l_Ayz
        + gupxy * (l_Axx * l_Ayy + l_Axy * l_Axy)
        + gupxz * (l_Axx * l_Ayz + l_Axz * l_Axy)
        + gupyz * (l_Axy * l_Ayz + l_Axz * l_Ayy);
    const double term_xz = gupxx * l_Axx * l_Axz + gupyy * l_Axy * l_Ayz + gupzz * l_Axz * l_Azz
        + gupxy * (l_Axx * l_Ayz + l_Axy * l_Axz)
        + gupxz * (l_Axx * l_Azz + l_Axz * l_Axz)
        + gupyz * (l_Axy * l_Azz + l_Axz * l_Ayz);
    const double term_yz = gupxx * l_Axy * l_Axz + gupyy * l_Ayy * l_Ayz + gupzz * l_Ayz * l_Azz
        + gupxy * (l_Axy * l_Ayz + l_Ayy * l_Axz)
        + gupxz * (l_Axy * l_Azz + l_Ayz * l_Axz)
        + gupyz * (l_Ayy * l_Azz + l_Ayz * l_Ayz);

    const double l_trK = trK[idx];
    const double trA2 = gupxx * term_xx + gupyy * term_yy + gupzz * term_zz
        + 2.0 * (gupxy * term_xy + gupxz * term_xz + gupyz * term_yz);
    const double S_trace = chin1 * (gupxx * Sxx[idx] + gupyy * Syy[idx] + gupzz * Szz[idx]
        + 2.0 * (gupxy * Sxy[idx] + gupxz * Sxz[idx] + gupyz * Syz[idx]));

    phys_Rxx[idx] = Rphys_xx;
    phys_Rxy[idx] = Rphys_xy;
    phys_Rxz[idx] = Rphys_xz;
    phys_Ryy[idx] = Rphys_yy;
    phys_Ryz[idx] = Rphys_yz;
    phys_Rzz[idx] = Rphys_zz;
    Axx_rhs[idx] = chin1 * (src_xx - l_gxx * trace_free) + alpn1 * (l_trK * l_Axx - 2.0 * term_xx)
        + 2.0 * (l_Axx * betaxx + l_Axy * betayx + l_Axz * betazx) - kF2o3 * l_Axx * div_beta;
    Ayy_rhs[idx] = chin1 * (src_yy - l_gyy * trace_free) + alpn1 * (l_trK * l_Ayy - 2.0 * term_yy)
        + 2.0 * (l_Axy * betaxy + l_Ayy * betayy + l_Ayz * betazy) - kF2o3 * l_Ayy * div_beta;
    Azz_rhs[idx] = chin1 * (src_zz - l_gzz * trace_free) + alpn1 * (l_trK * l_Azz - 2.0 * term_zz)
        + 2.0 * (l_Axz * betaxz + l_Ayz * betayz + l_Azz * betazz) - kF2o3 * l_Azz * div_beta;
    Axy_rhs[idx] = chin1 * (src_xy - l_gxy * trace_free) + alpn1 * (l_trK * l_Axy - 2.0 * term_xy)
        + l_Axx * betaxy + l_Axz * betazy
        + l_Ayy * betayx + l_Ayz * betazx
        + kF1o3 * l_Axy * div_beta - l_Axy * betazz;
    Ayz_rhs[idx] = chin1 * (src_yz - l_gyz * trace_free) + alpn1 * (l_trK * l_Ayz - 2.0 * term_yz)
        + l_Axy * betaxz + l_Ayy * betayz
        + l_Axz * betaxy + l_Azz * betazy
        + kF1o3 * l_Ayz * div_beta - l_Ayz * betaxx;
    Axz_rhs[idx] = chin1 * (src_xz - l_gxz * trace_free) + alpn1 * (l_trK * l_Axz - 2.0 * term_xz)
        + l_Axx * betaxz + l_Axy * betayz
        + l_Ayz * betayx + l_Azz * betazx
        + kF1o3 * l_Axz * div_beta - l_Axz * betayy;
    trK_rhs[idx] = -chin1 * lap_trace
        + alpn1 * (kF1o3 * l_trK * l_trK + trA2 + 0.5 * kEightPi * (rho[idx] + S_trace));
}

__global__ void bssn_shell_kodis_probe_kernel(
    int ex0,
    int ex1,
    int ex2,
    const double* crho,
    const double* sigma,
    const double* R,
    const double* f,
    double* f_rhs,
    int symmetry,
    double eps,
    int sst,
    double sym1,
    double sym2,
    double sym3
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int cell_count = ex0 * ex1 * ex2;
    if (idx >= cell_count) {
        return;
    }
    const int i = idx % ex0;
    const int j = (idx / ex0) % ex1;
    const int k = idx / (ex0 * ex1);
    const int ex[3] = {ex0, ex1, ex2};
    f_rhs[idx] += d_kodis_sh_point(
        ex, f, crho, sigma, R, sym1, sym2, sym3, symmetry, eps, sst, i, j, k);
}

__global__ void bssn_shell_transport_ko_probe_kernel(
    int ex0,
    int ex1,
    int ex2,
    const double* crho,
    const double* sigma,
    const double* R,
    const double* drhodx,
    const double* drhody,
    const double* drhodz,
    const double* dsigmadx,
    const double* dsigmady,
    const double* dsigmadz,
    const double* dRdx,
    const double* dRdy,
    const double* dRdz,
    const double* chi,
    const double* trK,
    const double* dxx,
    const double* gxy,
    const double* gxz,
    const double* dyy,
    const double* gyz,
    const double* dzz,
    const double* Axx,
    const double* Axy,
    const double* Axz,
    const double* Ayy,
    const double* Ayz,
    const double* Azz,
    const double* Gamx,
    const double* Gamy,
    const double* Gamz,
    const double* Lap,
    const double* betax,
    const double* betay,
    const double* betaz,
    const double* dtSfx,
    const double* dtSfy,
    const double* dtSfz,
    double* chi_rhs,
    double* trK_rhs,
    double* gxx_rhs,
    double* gxy_rhs,
    double* gxz_rhs,
    double* gyy_rhs,
    double* gyz_rhs,
    double* gzz_rhs,
    double* Axx_rhs,
    double* Axy_rhs,
    double* Axz_rhs,
    double* Ayy_rhs,
    double* Ayz_rhs,
    double* Azz_rhs,
    double* Gamx_rhs,
    double* Gamy_rhs,
    double* Gamz_rhs,
    double* Lap_rhs,
    double* betax_rhs,
    double* betay_rhs,
    double* betaz_rhs,
    double* dtSfx_rhs,
    double* dtSfy_rhs,
    double* dtSfz_rhs,
    int symmetry,
    double eps,
    int sst
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int cell_count = ex0 * ex1 * ex2;
    if (idx >= cell_count) {
        return;
    }
    const int i = idx % ex0;
    const int j = (idx / ex0) % ex1;
    const int k = idx / (ex0 * ex1);
    const int ex[3] = {ex0, ex1, ex2};

#define APPLY_TRANSPORT_KO(field, rhs, s1, s2, s3) \
    d_apply_transport_ko_rhs( \
        ex, crho, sigma, R, \
        drhodx, drhody, drhodz, \
        dsigmadx, dsigmady, dsigmadz, \
        dRdx, dRdy, dRdz, \
        field, betax, betay, betaz, rhs, \
        s1, s2, s3, symmetry, eps, sst, i, j, k)

    APPLY_TRANSPORT_KO(chi, chi_rhs, 1.0, 1.0, 1.0);
    APPLY_TRANSPORT_KO(trK, trK_rhs, 1.0, 1.0, 1.0);
    APPLY_TRANSPORT_KO(dxx, gxx_rhs, 1.0, 1.0, 1.0);
    APPLY_TRANSPORT_KO(gxy, gxy_rhs, -1.0, -1.0, 1.0);
    APPLY_TRANSPORT_KO(gxz, gxz_rhs, -1.0, 1.0, -1.0);
    APPLY_TRANSPORT_KO(dyy, gyy_rhs, 1.0, 1.0, 1.0);
    APPLY_TRANSPORT_KO(gyz, gyz_rhs, 1.0, -1.0, -1.0);
    APPLY_TRANSPORT_KO(dzz, gzz_rhs, 1.0, 1.0, 1.0);
    APPLY_TRANSPORT_KO(Axx, Axx_rhs, 1.0, 1.0, 1.0);
    APPLY_TRANSPORT_KO(Axy, Axy_rhs, -1.0, -1.0, 1.0);
    APPLY_TRANSPORT_KO(Axz, Axz_rhs, -1.0, 1.0, -1.0);
    APPLY_TRANSPORT_KO(Ayy, Ayy_rhs, 1.0, 1.0, 1.0);
    APPLY_TRANSPORT_KO(Ayz, Ayz_rhs, 1.0, -1.0, -1.0);
    APPLY_TRANSPORT_KO(Azz, Azz_rhs, 1.0, 1.0, 1.0);
    APPLY_TRANSPORT_KO(Gamx, Gamx_rhs, -1.0, 1.0, 1.0);
    APPLY_TRANSPORT_KO(Gamy, Gamy_rhs, 1.0, -1.0, 1.0);
    APPLY_TRANSPORT_KO(Gamz, Gamz_rhs, 1.0, 1.0, -1.0);
    APPLY_TRANSPORT_KO(Lap, Lap_rhs, 1.0, 1.0, 1.0);

#if (GAUGE == 0 || GAUGE == 1 || GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5 || GAUGE == 6 || GAUGE == 7)
    APPLY_TRANSPORT_KO(betax, betax_rhs, -1.0, 1.0, 1.0);
    APPLY_TRANSPORT_KO(betay, betay_rhs, 1.0, -1.0, 1.0);
    APPLY_TRANSPORT_KO(betaz, betaz_rhs, 1.0, 1.0, -1.0);
#endif

#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)
    APPLY_TRANSPORT_KO(dtSfx, dtSfx_rhs, -1.0, 1.0, 1.0);
    APPLY_TRANSPORT_KO(dtSfy, dtSfy_rhs, 1.0, -1.0, 1.0);
    APPLY_TRANSPORT_KO(dtSfz, dtSfz_rhs, 1.0, 1.0, -1.0);
#else
    (void)dtSfx;
    (void)dtSfy;
    (void)dtSfz;
    (void)dtSfx_rhs;
    (void)dtSfy_rhs;
    (void)dtSfz_rhs;
#endif

#undef APPLY_TRANSPORT_KO
}

__global__ void bssn_shell_constraint_probe_kernel(
    int ex0,
    int ex1,
    int ex2,
    const double* crho,
    const double* sigma,
    const double* R,
    const double* drhodx,
    const double* drhody,
    const double* drhodz,
    const double* dsigmadx,
    const double* dsigmady,
    const double* dsigmadz,
    const double* dRdx,
    const double* dRdy,
    const double* dRdz,
    const double* chi,
    const double* trK,
    const double* dxx,
    const double* gxy,
    const double* gxz,
    const double* dyy,
    const double* gyz,
    const double* dzz,
    const double* Axx,
    const double* Axy,
    const double* Axz,
    const double* Ayy,
    const double* Ayz,
    const double* Azz,
    const double* Gamx,
    const double* Gamy,
    const double* Gamz,
    const double* rho,
    const double* Sx,
    const double* Sy,
    const double* Sz,
    const double* Gamxxx,
    const double* Gamxxy,
    const double* Gamxxz,
    const double* Gamxyy,
    const double* Gamxyz,
    const double* Gamxzz,
    const double* Gamyxx,
    const double* Gamyxy,
    const double* Gamyxz,
    const double* Gamyyy,
    const double* Gamyyz,
    const double* Gamyzz,
    const double* Gamzxx,
    const double* Gamzxy,
    const double* Gamzxz,
    const double* Gamzyy,
    const double* Gamzyz,
    const double* Gamzzz,
    const double* phys_Rxx,
    const double* phys_Rxy,
    const double* phys_Rxz,
    const double* phys_Ryy,
    const double* phys_Ryz,
    const double* phys_Rzz,
    double* ham_Res,
    double* movx_Res,
    double* movy_Res,
    double* movz_Res,
    double* Gmx_Res,
    double* Gmy_Res,
    double* Gmz_Res,
    int symmetry,
    int sst
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int cell_count = ex0 * ex1 * ex2;
    if (idx >= cell_count) {
        return;
    }
    const int i = idx % ex0;
    const int j = (idx / ex0) % ex1;
    const int k = idx / (ex0 * ex1);
    const int ex[3] = {ex0, ex1, ex2};

    double gxxx = 0.0, gxxy = 0.0, gxxz = 0.0;
    double gxyx = 0.0, gxyy = 0.0, gxyz = 0.0;
    double gxzx = 0.0, gxzy = 0.0, gxzz = 0.0;
    double gyyx = 0.0, gyyy = 0.0, gyyz = 0.0;
    double gyzx = 0.0, gyzy = 0.0, gyzz = 0.0;
    double gzzx = 0.0, gzzy = 0.0, gzzz = 0.0;
    d_fderivs_shc_point(ex, dxx, &gxxx, &gxxy, &gxxz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, gxy, &gxyx, &gxyy, &gxyz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, -1.0, -1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, gxz, &gxzx, &gxzy, &gxzz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, -1.0, 1.0, -1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, dyy, &gyyx, &gyyy, &gyyz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, gyz, &gyzx, &gyzy, &gyzz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, -1.0, -1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, dzz, &gzzx, &gzzy, &gzzz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, 1.0, 1.0, symmetry, sst, i, j, k);

    const double l_gxx = dxx[idx] + 1.0;
    const double l_gyy = dyy[idx] + 1.0;
    const double l_gzz = dzz[idx] + 1.0;
    const double l_gxy = gxy[idx];
    const double l_gxz = gxz[idx];
    const double l_gyz = gyz[idx];
    const double det = l_gxx * l_gyy * l_gzz
        + l_gxy * l_gyz * l_gxz + l_gxz * l_gxy * l_gyz
        - l_gxz * l_gyy * l_gxz - l_gxy * l_gxy * l_gzz
        - l_gxx * l_gyz * l_gyz;
    const double gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) / det;
    const double gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) / det;
    const double gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) / det;
    const double gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) / det;
    const double gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) / det;
    const double gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) / det;

    Gmx_Res[idx] = Gamx[idx] - (
        gupxx * (gupxx * gxxx + gupxy * gxyx + gupxz * gxzx)
        + gupxy * (gupxx * gxyx + gupxy * gyyx + gupxz * gyzx)
        + gupxz * (gupxx * gxzx + gupxy * gyzx + gupxz * gzzx)
        + gupxx * (gupxy * gxxy + gupyy * gxyy + gupyz * gxzy)
        + gupxy * (gupxy * gxyy + gupyy * gyyy + gupyz * gyzy)
        + gupxz * (gupxy * gxzy + gupyy * gyzy + gupyz * gzzy)
        + gupxx * (gupxz * gxxz + gupyz * gxyz + gupzz * gxzz)
        + gupxy * (gupxz * gxyz + gupyz * gyyz + gupzz * gyzz)
        + gupxz * (gupxz * gxzz + gupyz * gyzz + gupzz * gzzz));
    Gmy_Res[idx] = Gamy[idx] - (
        gupxx * (gupxy * gxxx + gupyy * gxyx + gupyz * gxzx)
        + gupxy * (gupxy * gxyx + gupyy * gyyx + gupyz * gyzx)
        + gupxz * (gupxy * gxzx + gupyy * gyzx + gupyz * gzzx)
        + gupxy * (gupxy * gxxy + gupyy * gxyy + gupyz * gxzy)
        + gupyy * (gupxy * gxyy + gupyy * gyyy + gupyz * gyzy)
        + gupyz * (gupxy * gxzy + gupyy * gyzy + gupyz * gzzy)
        + gupxy * (gupxz * gxxz + gupyz * gxyz + gupzz * gxzz)
        + gupyy * (gupxz * gxyz + gupyz * gyyz + gupzz * gyzz)
        + gupyz * (gupxz * gxzz + gupyz * gyzz + gupzz * gzzz));
    Gmz_Res[idx] = Gamz[idx] - (
        gupxx * (gupxz * gxxx + gupyz * gxyx + gupzz * gxzx)
        + gupxy * (gupxz * gxyx + gupyz * gyyx + gupzz * gyzx)
        + gupxz * (gupxz * gxzx + gupyz * gyzx + gupzz * gzzx)
        + gupxy * (gupxz * gxxy + gupyz * gxyy + gupzz * gxzy)
        + gupyy * (gupxz * gxyy + gupyz * gyyy + gupzz * gyzy)
        + gupyz * (gupxz * gxzy + gupyz * gyzy + gupzz * gzzy)
        + gupxz * (gupxz * gxxz + gupyz * gxyz + gupzz * gxzz)
        + gupyz * (gupxz * gxyz + gupyz * gyyz + gupzz * gyzz)
        + gupzz * (gupxz * gxzz + gupyz * gyzz + gupzz * gzzz));

    const double chin1 = chi[idx] + 1.0;
    const double l_Axx = Axx[idx], l_Axy = Axy[idx], l_Axz = Axz[idx];
    const double l_Ayy = Ayy[idx], l_Ayz = Ayz[idx], l_Azz = Azz[idx];
    const double term_xx = gupxx * l_Axx * l_Axx + gupyy * l_Axy * l_Axy + gupzz * l_Axz * l_Axz
        + 2.0 * (gupxy * l_Axx * l_Axy + gupxz * l_Axx * l_Axz + gupyz * l_Axy * l_Axz);
    const double term_xy = gupxx * l_Axx * l_Axy + gupyy * l_Axy * l_Ayy + gupzz * l_Axz * l_Ayz
        + gupxy * (l_Axx * l_Ayy + l_Axy * l_Axy)
        + gupxz * (l_Axx * l_Ayz + l_Axz * l_Axy)
        + gupyz * (l_Axy * l_Ayz + l_Axz * l_Ayy);
    const double term_xz = gupxx * l_Axx * l_Axz + gupyy * l_Axy * l_Ayz + gupzz * l_Axz * l_Azz
        + gupxy * (l_Axx * l_Ayz + l_Axy * l_Axz)
        + gupxz * (l_Axx * l_Azz + l_Axz * l_Axz)
        + gupyz * (l_Axy * l_Azz + l_Axz * l_Ayz);
    const double term_yy = gupxx * l_Axy * l_Axy + gupyy * l_Ayy * l_Ayy + gupzz * l_Ayz * l_Ayz
        + 2.0 * (gupxy * l_Axy * l_Ayy + gupxz * l_Axy * l_Ayz + gupyz * l_Ayy * l_Ayz);
    const double term_yz = gupxx * l_Axy * l_Axz + gupyy * l_Ayy * l_Ayz + gupzz * l_Ayz * l_Azz
        + gupxy * (l_Axy * l_Ayz + l_Ayy * l_Axz)
        + gupxz * (l_Axy * l_Azz + l_Ayz * l_Axz)
        + gupyz * (l_Ayy * l_Azz + l_Ayz * l_Ayz);
    const double term_zz = gupxx * l_Axz * l_Axz + gupyy * l_Ayz * l_Ayz + gupzz * l_Azz * l_Azz
        + 2.0 * (gupxy * l_Axz * l_Ayz + gupxz * l_Axz * l_Azz + gupyz * l_Ayz * l_Azz);
    const double trA2 = gupxx * term_xx + gupyy * term_yy + gupzz * term_zz
        + 2.0 * (gupxy * term_xy + gupxz * term_xz + gupyz * term_yz);
    double trR = gupxx * phys_Rxx[idx] + gupyy * phys_Ryy[idx] + gupzz * phys_Rzz[idx]
        + 2.0 * (gupxy * phys_Rxy[idx] + gupxz * phys_Rxz[idx] + gupyz * phys_Ryz[idx]);
    ham_Res[idx] = chin1 * trR + kF2o3 * trK[idx] * trK[idx] - trA2 - 2.0 * kEightPi * rho[idx];

    double chix = 0.0, chiy = 0.0, chiz = 0.0;
    double Kx = 0.0, Ky = 0.0, Kz = 0.0;
    d_fderivs_shc_point(ex, chi, &chix, &chiy, &chiz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, trK, &Kx, &Ky, &Kz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, 1.0, 1.0, symmetry, sst, i, j, k);

    double Axxx = 0.0, Axxy = 0.0, Axxz = 0.0;
    double Axyx = 0.0, Axyy = 0.0, Axyz = 0.0;
    double Axzx = 0.0, Axzy = 0.0, Axzz = 0.0;
    double Ayyx = 0.0, Ayyy = 0.0, Ayyz = 0.0;
    double Ayzx = 0.0, Ayzy = 0.0, Ayzz = 0.0;
    double Azzx = 0.0, Azzy = 0.0, Azzz = 0.0;
    d_fderivs_shc_point(ex, Axx, &Axxx, &Axxy, &Axxz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, Axy, &Axyx, &Axyy, &Axyz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, -1.0, -1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, Axz, &Axzx, &Axzy, &Axzz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, -1.0, 1.0, -1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, Ayy, &Ayyx, &Ayyy, &Ayyz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, 1.0, 1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, Ayz, &Ayzx, &Ayzy, &Ayzz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, -1.0, -1.0, symmetry, sst, i, j, k);
    d_fderivs_shc_point(ex, Azz, &Azzx, &Azzy, &Azzz,
        crho, sigma, R, drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz, 1.0, 1.0, 1.0, symmetry, sst, i, j, k);

    const double l_Gamxxx = Gamxxx[idx], l_Gamxxy = Gamxxy[idx], l_Gamxxz = Gamxxz[idx];
    const double l_Gamxyy = Gamxyy[idx], l_Gamxyz = Gamxyz[idx], l_Gamxzz = Gamxzz[idx];
    const double l_Gamyxx = Gamyxx[idx], l_Gamyxy = Gamyxy[idx], l_Gamyxz = Gamyxz[idx];
    const double l_Gamyyy = Gamyyy[idx], l_Gamyyz = Gamyyz[idx], l_Gamyzz = Gamyzz[idx];
    const double l_Gamzxx = Gamzxx[idx], l_Gamzxy = Gamzxy[idx], l_Gamzxz = Gamzxz[idx];
    const double l_Gamzyy = Gamzyy[idx], l_Gamzyz = Gamzyz[idx], l_Gamzzz = Gamzzz[idx];
    const double inv_chin1 = 1.0 / chin1;

    gxxx = Axxx - (l_Gamxxx * l_Axx + l_Gamyxx * l_Axy + l_Gamzxx * l_Axz
                 + l_Gamxxx * l_Axx + l_Gamyxx * l_Axy + l_Gamzxx * l_Axz) - chix * l_Axx * inv_chin1;
    gxyx = Axyx - (l_Gamxxy * l_Axx + l_Gamyxy * l_Axy + l_Gamzxy * l_Axz
                 + l_Gamxxx * l_Axy + l_Gamyxx * l_Ayy + l_Gamzxx * l_Ayz) - chix * l_Axy * inv_chin1;
    gxzx = Axzx - (l_Gamxxz * l_Axx + l_Gamyxz * l_Axy + l_Gamzxz * l_Axz
                 + l_Gamxxx * l_Axz + l_Gamyxx * l_Ayz + l_Gamzxx * l_Azz) - chix * l_Axz * inv_chin1;
    gyyx = Ayyx - (l_Gamxxy * l_Axy + l_Gamyxy * l_Ayy + l_Gamzxy * l_Ayz
                 + l_Gamxxy * l_Axy + l_Gamyxy * l_Ayy + l_Gamzxy * l_Ayz) - chix * l_Ayy * inv_chin1;
    gyzx = Ayzx - (l_Gamxxz * l_Axy + l_Gamyxz * l_Ayy + l_Gamzxz * l_Ayz
                 + l_Gamxxy * l_Axz + l_Gamyxy * l_Ayz + l_Gamzxy * l_Azz) - chix * l_Ayz * inv_chin1;
    gzzx = Azzx - (l_Gamxxz * l_Axz + l_Gamyxz * l_Ayz + l_Gamzxz * l_Azz
                 + l_Gamxxz * l_Axz + l_Gamyxz * l_Ayz + l_Gamzxz * l_Azz) - chix * l_Azz * inv_chin1;
    gxxy = Axxy - (l_Gamxxy * l_Axx + l_Gamyxy * l_Axy + l_Gamzxy * l_Axz
                 + l_Gamxxy * l_Axx + l_Gamyxy * l_Axy + l_Gamzxy * l_Axz) - chiy * l_Axx * inv_chin1;
    gxyy = Axyy - (l_Gamxyy * l_Axx + l_Gamyyy * l_Axy + l_Gamzyy * l_Axz
                 + l_Gamxxy * l_Axy + l_Gamyxy * l_Ayy + l_Gamzxy * l_Ayz) - chiy * l_Axy * inv_chin1;
    gxzy = Axzy - (l_Gamxyz * l_Axx + l_Gamyyz * l_Axy + l_Gamzyz * l_Axz
                 + l_Gamxxy * l_Axz + l_Gamyxy * l_Ayz + l_Gamzxy * l_Azz) - chiy * l_Axz * inv_chin1;
    gyyy = Ayyy - (l_Gamxyy * l_Axy + l_Gamyyy * l_Ayy + l_Gamzyy * l_Ayz
                 + l_Gamxyy * l_Axy + l_Gamyyy * l_Ayy + l_Gamzyy * l_Ayz) - chiy * l_Ayy * inv_chin1;
    gyzy = Ayzy - (l_Gamxyz * l_Axy + l_Gamyyz * l_Ayy + l_Gamzyz * l_Ayz
                 + l_Gamxyy * l_Axz + l_Gamyyy * l_Ayz + l_Gamzyy * l_Azz) - chiy * l_Ayz * inv_chin1;
    gzzy = Azzy - (l_Gamxyz * l_Axz + l_Gamyyz * l_Ayz + l_Gamzyz * l_Azz
                 + l_Gamxyz * l_Axz + l_Gamyyz * l_Ayz + l_Gamzyz * l_Azz) - chiy * l_Azz * inv_chin1;
    gxxz = Axxz - (l_Gamxxz * l_Axx + l_Gamyxz * l_Axy + l_Gamzxz * l_Axz
                 + l_Gamxxz * l_Axx + l_Gamyxz * l_Axy + l_Gamzxz * l_Axz) - chiz * l_Axx * inv_chin1;
    gxyz = Axyz - (l_Gamxyz * l_Axx + l_Gamyyz * l_Axy + l_Gamzyz * l_Axz
                 + l_Gamxxz * l_Axy + l_Gamyxz * l_Ayy + l_Gamzxz * l_Ayz) - chiz * l_Axy * inv_chin1;
    gxzz = Axzz - (l_Gamxzz * l_Axx + l_Gamyzz * l_Axy + l_Gamzzz * l_Axz
                 + l_Gamxxz * l_Axz + l_Gamyxz * l_Ayz + l_Gamzxz * l_Azz) - chiz * l_Axz * inv_chin1;
    gyyz = Ayyz - (l_Gamxyz * l_Axy + l_Gamyyz * l_Ayy + l_Gamzyz * l_Ayz
                 + l_Gamxyz * l_Axy + l_Gamyyz * l_Ayy + l_Gamzyz * l_Ayz) - chiz * l_Ayy * inv_chin1;
    gyzz = Ayzz - (l_Gamxzz * l_Axy + l_Gamyzz * l_Ayy + l_Gamzzz * l_Ayz
                 + l_Gamxyz * l_Axz + l_Gamyyz * l_Ayz + l_Gamzyz * l_Azz) - chiz * l_Ayz * inv_chin1;
    gzzz = Azzz - (l_Gamxzz * l_Axz + l_Gamyzz * l_Ayz + l_Gamzzz * l_Azz
                 + l_Gamxzz * l_Axz + l_Gamyzz * l_Ayz + l_Gamzzz * l_Azz) - chiz * l_Azz * inv_chin1;

    movx_Res[idx] = gupxx * gxxx + gupyy * gxyy + gupzz * gxzz
        + gupxy * gxyx + gupxz * gxzx + gupyz * gxzy
        + gupxy * gxxy + gupxz * gxxz + gupyz * gxyz
        - kF2o3 * Kx - kEightPi * Sx[idx];
    movy_Res[idx] = gupxx * gxyx + gupyy * gyyy + gupzz * gyzz
        + gupxy * gyyx + gupxz * gyzx + gupyz * gyzy
        + gupxy * gxyy + gupxz * gxyz + gupyz * gyyz
        - kF2o3 * Ky - kEightPi * Sy[idx];
    movz_Res[idx] = gupxx * gxzx + gupyy * gyzy + gupzz * gzzz
        + gupxy * gyzx + gupxz * gzzx + gupyz * gzzy
        + gupxy * gxzy + gupxz * gxzz + gupyz * gyzz
        - kF2o3 * Kz - kEightPi * Sz[idx];
}

// Z4C corrections on a shell block: adds the theta (TZ) terms to the BSSN
// shell RHS, mirroring apply_z4c_algebraic_corrections_kernel /
// apply_z4c_gauge_corrections_kernel / apply_z4c_transport_damping_kernel
// from bssn_rhs_gpu.cu. Must run after the curvature and constraint probes
// and before the transport/KO probe, so that the algebraic part reads the
// transport-free trK_rhs exactly like the Cartesian Z4C path does.
__global__ void apply_z4c_shell_corrections_kernel(
    int ex0, int ex1, int ex2,
    const double* crho, const double* sigma, const double* R,
    const double* drhodx, const double* drhody, const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx, const double* dRdy, const double* dRdz,
    const double* chi, const double* trK, const double* TZ,
    const double* gxx, const double* gxy, const double* gxz,
    const double* gyy, const double* gyz, const double* gzz,
    const double* Axx, const double* Axy, const double* Axz,
    const double* Ayy, const double* Ayz, const double* Azz,
    const double* Lap,
    const double* betax, const double* betay, const double* betaz,
    const double* rho, const double* Sxx, const double* Sxy, const double* Sxz,
    const double* Syy, const double* Syz, const double* Szz,
    const double* ham_Res, const double* Gmx_Res, const double* Gmy_Res, const double* Gmz_Res,
    double* chi_rhs, double* trK_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs,
    double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* TZ_rhs,
    int symmetry, double eps, int sst
) {
    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    const int count = ex0 * ex1 * ex2;
    if (point_id >= count) return;
    const int i = point_id % ex0;
    const int j = (point_id / ex0) % ex1;
    const int k = point_id / (ex0 * ex1);
    const int ex[3] = {ex0, ex1, ex2};
    const int idx = shell_flatten(i, j, k, ex0, ex1);

    // ---- algebraic corrections (pointwise, as on Cartesian blocks) ----
    const double val_chi = chi[idx];
    const double chin1 = val_chi + 1.0;
    const double inv_chin1 = 1.0 / chin1;
    const double val_trK = trK[idx];
    const double val_TZ = TZ[idx];
    const double trKd = val_trK + 2.0 * val_TZ;
    const double alpn1 = Lap[idx] + 1.0;

    const double l_gxx = gxx[idx] + 1.0;
    const double l_gxy = gxy[idx];
    const double l_gxz = gxz[idx];
    const double l_gyy = gyy[idx] + 1.0;
    const double l_gyz = gyz[idx];
    const double l_gzz = gzz[idx] + 1.0;

    const double l_Axx = Axx[idx];
    const double l_Axy = Axy[idx];
    const double l_Axz = Axz[idx];
    const double l_Ayy = Ayy[idx];
    const double l_Ayz = Ayz[idx];
    const double l_Azz = Azz[idx];

    const double det = l_gxx * (l_gyy * l_gzz - l_gyz * l_gyz)
                     - l_gxy * (l_gxy * l_gzz - l_gyz * l_gxz)
                     + l_gxz * (l_gxy * l_gyz - l_gyy * l_gxz);
    const double inv_det = 1.0 / det;
    const double gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) * inv_det;
    const double gupxy = (l_gxz * l_gyz - l_gxy * l_gzz) * inv_det;
    const double gupxz = (l_gxy * l_gyz - l_gxz * l_gyy) * inv_det;
    const double gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) * inv_det;
    const double gupyz = (l_gxz * l_gxy - l_gxx * l_gyz) * inv_det;
    const double gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) * inv_det;

    // Aup = gup * A * gup
    const double MAxx = gupxx * l_Axx + gupxy * l_Axy + gupxz * l_Axz;
    const double MAxy = gupxx * l_Axy + gupxy * l_Ayy + gupxz * l_Ayz;
    const double MAxz = gupxx * l_Axz + gupxy * l_Ayz + gupxz * l_Azz;
    const double MAyy = gupxy * l_Axy + gupyy * l_Ayy + gupyz * l_Ayz;
    const double MAyz = gupxy * l_Axz + gupyy * l_Ayz + gupyz * l_Azz;
    const double MAzz = gupxz * l_Axz + gupyz * l_Ayz + gupzz * l_Azz;
    const double Aupxx = MAxx * gupxx + MAxy * gupxy + MAxz * gupxz;
    const double Aupxy = MAxx * gupxy + MAxy * gupyy + MAxz * gupyz;
    const double Aupxz = MAxx * gupxz + MAxy * gupyz + MAxz * gupzz;
    const double Aupyy = MAyy * gupyy + MAyz * gupyz;
    const double Aupyz = MAyy * gupyz + MAyz * gupzz;
    const double Aupzz = MAzz * gupzz;

    const double trA2 = l_Axx * Aupxx + l_Ayy * Aupyy + l_Azz * Aupzz
                      + 2.0 * (l_Axy * Aupxy + l_Axz * Aupxz + l_Ayz * Aupyz);
    const double trace_S = gupxx * Sxx[idx] + gupyy * Syy[idx] + gupzz * Szz[idx]
                         + 2.0 * (gupxy * Sxy[idx] + gupxz * Sxz[idx] + gupyz * Syz[idx]);
    const double S = chin1 * trace_S;

    const double d2_lap = (-trK_rhs[idx]
        + alpn1 * (kF1o3 * val_trK * val_trK + trA2 + kFourPi * (rho[idx] + S))) * inv_chin1;
    const double bssn_trace_aux = kF2o3 * val_trK * val_trK - trA2 - kF16Pi * rho[idx] + kEightPi * S;
    const double f_trace_bssn = -kF1o3 * (d2_lap + alpn1 * inv_chin1 * bssn_trace_aux);

    const double ham_bssn = ham_Res[idx];
    const double ricci_trace = (ham_bssn - kF2o3 * val_trK * val_trK + trA2 + kF16Pi * rho[idx]) * inv_chin1;
    const double f_trace_z4c = kF1o3 * (alpn1 * ricci_trace - (kEightPi * alpn1 * trace_S + d2_lap));
    const double trace_delta = chin1 * (f_trace_bssn - f_trace_z4c);
    const double trKd_delta = alpn1 * (trKd - val_trK);

    chi_rhs[idx] += kF2o3 * chin1 * alpn1 * (trKd - val_trK);

    Axx_rhs[idx] += trace_delta * l_gxx + trKd_delta * l_Axx;
    Ayy_rhs[idx] += trace_delta * l_gyy + trKd_delta * l_Ayy;
    Azz_rhs[idx] += trace_delta * l_gzz + trKd_delta * l_Azz;
    Axy_rhs[idx] += trace_delta * l_gxy + trKd_delta * l_Axy;
    Ayz_rhs[idx] += trace_delta * l_gyz + trKd_delta * l_Ayz;
    Axz_rhs[idx] += trace_delta * l_gxz + trKd_delta * l_Axz;

    trK_rhs[idx] += alpn1 * kF1o3 * (trKd * trKd - val_trK * val_trK);

    const double ham_z4c = ham_bssn + kF2o3 * (trKd * trKd - val_trK * val_trK);
    TZ_rhs[idx] = 0.5 * alpn1 * ham_z4c;

    // ---- gauge corrections: TZ gradient via curvilinear chain rule ----
    double TZx = 0.0;
    double TZy = 0.0;
    double TZz = 0.0;
    d_fderivs_shc_point(
        ex, TZ, &TZx, &TZy, &TZz,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        1.0, 1.0, 1.0, symmetry, sst, i, j, k);

    const double delta_x = -kF2o3 * alpn1 * (gupxx * TZx + gupxy * TZy + gupxz * TZz);
    const double delta_y = -kF2o3 * alpn1 * (gupxy * TZx + gupyy * TZy + gupyz * TZz);
    const double delta_z = -kF2o3 * alpn1 * (gupxz * TZx + gupyz * TZy + gupzz * TZz);

    Gamx_rhs[idx] += delta_x;
    Gamy_rhs[idx] += delta_y;
    Gamz_rhs[idx] += delta_z;

#if (GAUGE == 0)
    dtSfx_rhs[idx] += delta_x;
    dtSfy_rhs[idx] += delta_y;
    dtSfz_rhs[idx] += delta_z;
#else
#error Unsupported Z4C GPU shell gauge; the shell Z4C corrections only implement GAUGE == 0.
#endif

    // ---- transport and damping for TZ ----
    TZ_rhs[idx] += betax[idx] * TZx + betay[idx] * TZy + betaz[idx] * TZz;
    if (eps > 0.0) {
        TZ_rhs[idx] += d_kodis_sh_point(
            ex, TZ, crho, sigma, R,
            1.0, 1.0, 1.0, symmetry, eps, sst, i, j, k);
    }

    TZ_rhs[idx] -= alpn1 * (2.0 + kZ4Ckappa2) * kZ4Ckappa1 * val_TZ;
    trK_rhs[idx] += alpn1 * kZ4Ckappa1 * (1.0 - kZ4Ckappa2) * val_TZ;
    Gamx_rhs[idx] -= 2.0 * alpn1 * kZ4Ckappa1 * Gmx_Res[idx];
    Gamy_rhs[idx] -= 2.0 * alpn1 * kZ4Ckappa1 * Gmy_Res[idx];
    Gamz_rhs[idx] -= 2.0 * alpn1 * kZ4Ckappa1 * Gmz_Res[idx];
}

}  // namespace

void gpu_apply_z4c_shell_corrections_launch(
    cudaStream_t &stream,
    int* ex,
    double* crho, double* sigma, double* R,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* chi, double* trK, double* TZ,
    double* gxx, double* gxy, double* gxz,
    double* gyy, double* gyz, double* gzz,
    double* Axx, double* Axy, double* Axz,
    double* Ayy, double* Ayz, double* Azz,
    double* Lap,
    double* betax, double* betay, double* betaz,
    double* rho, double* Sxx, double* Sxy, double* Sxz,
    double* Syy, double* Syz, double* Szz,
    double* ham_Res, double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    double* chi_rhs, double* trK_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs,
    double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* TZ_rhs,
    int symmetry, double eps, int sst
) {
    const int count = shell_rhs_cell_count(ex);
    apply_z4c_shell_corrections_kernel<<<shell_rhs_grid_size(count), kGuardBlockSize, 0, stream>>>(
        ex[0], ex[1], ex[2],
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        chi, trK, TZ,
        gxx, gxy, gxz,
        gyy, gyz, gzz,
        Axx, Axy, Axz,
        Ayy, Ayz, Azz,
        Lap,
        betax, betay, betaz,
        rho, Sxx, Sxy, Sxz,
        Syy, Syz, Szz,
        ham_Res, Gmx_Res, Gmy_Res, Gmz_Res,
        chi_rhs, trK_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs,
        Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        TZ_rhs,
        symmetry, eps, sst);
    abort_on_cuda(cudaGetLastError(), "apply_z4c_shell_corrections_kernel");
}

void gpu_bssn_shell_derivative_probe_launch(
    cudaStream_t &stream,
    int* ex,
    double* crho, double* sigma, double* R,
    double* f,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* drhodxx, double* drhodxy, double* drhodxz,
    double* drhodyy, double* drhodyz, double* drhodzz,
    double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
    double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
    double* dRdxx, double* dRdxy, double* dRdxz,
    double* dRdyy, double* dRdyz, double* dRdzz,
    double* fx, double* fy, double* fz,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    int symmetry, int lev, int sst,
    double sym1, double sym2, double sym3
) {
    (void)lev;
    const int count = shell_rhs_cell_count(ex);
    bssn_shell_derivative_probe_kernel<<<shell_rhs_grid_size(count), kGuardBlockSize, 0, stream>>>(
        ex[0], ex[1], ex[2],
        crho, sigma, R,
        f,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        fx, fy, fz,
        fxx, fxy, fxz,
        fyy, fyz, fzz,
        symmetry, sst, sym1, sym2, sym3);
    abort_on_cuda(cudaGetLastError(), "bssn_shell_derivative_probe_kernel");
}

void gpu_bssn_shell_metric_rhs_probe_launch(
    cudaStream_t &stream,
    int* ex,
    double* crho, double* sigma, double* R,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz,
    double* Ayy, double* Ayz, double* Azz,
    double* Lap,
    double* betax, double* betay, double* betaz,
    double* chi_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs,
    double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    int symmetry, int lev, int sst
) {
    (void)lev;
    const int count = shell_rhs_cell_count(ex);
    bssn_shell_metric_rhs_probe_kernel<<<shell_rhs_grid_size(count), kGuardBlockSize, 0, stream>>>(
        ex[0], ex[1], ex[2],
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        chi, trK,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Axx, Axy, Axz,
        Ayy, Ayz, Azz,
        Lap,
        betax, betay, betaz,
        chi_rhs,
        gxx_rhs, gxy_rhs, gxz_rhs,
        gyy_rhs, gyz_rhs, gzz_rhs,
        symmetry, sst);
    abort_on_cuda(cudaGetLastError(), "bssn_shell_metric_rhs_probe_kernel");
}

void gpu_bssn_shell_connection_probe_launch(
    cudaStream_t &stream,
    int* ex,
    double* crho, double* sigma, double* R,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz,
    double* Ayy, double* Ayz, double* Azz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    int symmetry, int lev, int sst
) {
    (void)lev;
    const int count = shell_rhs_cell_count(ex);
    bssn_shell_connection_probe_kernel<<<shell_rhs_grid_size(count), kGuardBlockSize, 0, stream>>>(
        ex[0], ex[1], ex[2],
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Axx, Axy, Axz,
        Ayy, Ayz, Azz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz,
        Ryy, Ryz, Rzz,
        symmetry, sst);
    abort_on_cuda(cudaGetLastError(), "bssn_shell_connection_probe_kernel");
}

void gpu_bssn_shell_ricci_probe_launch(
    cudaStream_t &stream,
    int* ex,
    double* crho, double* sigma, double* R,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* drhodxx, double* drhodxy, double* drhodxz,
    double* drhodyy, double* drhodyz, double* drhodzz,
    double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
    double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
    double* dRdxx, double* dRdxy, double* dRdxz,
    double* dRdyy, double* dRdyz, double* dRdzz,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Gamx, double* Gamy, double* Gamz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Ricci_xx, double* Ricci_xy, double* Ricci_xz,
    double* Ricci_yy, double* Ricci_yz, double* Ricci_zz,
    int symmetry, int lev, int sst
) {
    (void)lev;
    const int count = shell_rhs_cell_count(ex);
    bssn_shell_ricci_probe_kernel<<<shell_rhs_grid_size(count), kGuardBlockSize, 0, stream>>>(
        ex[0], ex[1], ex[2],
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Gamx, Gamy, Gamz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz,
        Ricci_xx, Ricci_xy, Ricci_xz,
        Ricci_yy, Ricci_yz, Ricci_zz,
        symmetry, sst);
    abort_on_cuda(cudaGetLastError(), "bssn_shell_ricci_probe_kernel");
}

void gpu_bssn_shell_gauge_rhs_probe_launch(
    cudaStream_t &stream,
    int* ex,
    double* crho, double* sigma, double* R,
    double* x, double* y, double* z,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* drhodxx, double* drhodxy, double* drhodxz,
    double* drhodyy, double* drhodyz, double* drhodzz,
    double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
    double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
    double* dRdxx, double* dRdxy, double* dRdxz,
    double* dRdyy, double* dRdyz, double* dRdzz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Gamx, double* Gamy, double* Gamz,
    double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* Sx, double* Sy, double* Sz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs,
    double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    int symmetry, int lev, int sst
) {
    (void)lev;
    const int count = shell_rhs_cell_count(ex);
#if (GAUGE == 6 || GAUGE == 7)
    cudaMemcpyToSymbolAsync(d_bssn_shell_gauge_bh_count, &BHN, sizeof(int), 0, cudaMemcpyHostToDevice, stream);
    cudaMemcpyToSymbolAsync(d_bssn_shell_gauge_mass, Mass, 3 * sizeof(double), 0, cudaMemcpyHostToDevice, stream);
    cudaMemcpyToSymbolAsync(d_bssn_shell_gauge_pbh, PBH, 9 * sizeof(double), 0, cudaMemcpyHostToDevice, stream);
#endif
    bssn_shell_gauge_rhs_probe_kernel<<<shell_rhs_grid_size(count), kGuardBlockSize, 0, stream>>>(
        ex[0], ex[1], ex[2],
        crho, sigma, R,
        x, y, z,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        chi, trK,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Gamx, Gamy, Gamz,
        Lap,
        betax, betay, betaz,
        dtSfx, dtSfy, dtSfz,
        Sx, Sy, Sz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz,
        Ryy, Ryz, Rzz,
        Gamx_rhs, Gamy_rhs, Gamz_rhs,
        Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        symmetry, sst);
    abort_on_cuda(cudaGetLastError(), "bssn_shell_gauge_rhs_probe_kernel");
}

void gpu_bssn_shell_curvature_rhs_probe_launch(
    cudaStream_t &stream,
    int* ex,
    double* crho, double* sigma, double* R,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* drhodxx, double* drhodxy, double* drhodxz,
    double* drhodyy, double* drhodyz, double* drhodzz,
    double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
    double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
    double* dRdxx, double* dRdxy, double* dRdxz,
    double* dRdyy, double* dRdyz, double* dRdzz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz,
    double* Ayy, double* Ayz, double* Azz,
    double* Lap,
    double* betax, double* betay, double* betaz,
    double* rho,
    double* Sxx, double* Sxy, double* Sxz,
    double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Ricci_xx, double* Ricci_xy, double* Ricci_xz,
    double* Ricci_yy, double* Ricci_yz, double* Ricci_zz,
    double* phys_Rxx, double* phys_Rxy, double* phys_Rxz,
    double* phys_Ryy, double* phys_Ryz, double* phys_Rzz,
    double* trK_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs,
    double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    int symmetry, int lev, int sst
) {
    (void)lev;
    const int count = shell_rhs_cell_count(ex);
    bssn_shell_curvature_rhs_probe_kernel<<<shell_rhs_grid_size(count), kGuardBlockSize, 0, stream>>>(
        ex[0], ex[1], ex[2],
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        chi, trK,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Axx, Axy, Axz,
        Ayy, Ayz, Azz,
        Lap,
        betax, betay, betaz,
        rho,
        Sxx, Sxy, Sxz,
        Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz,
        Ricci_xx, Ricci_xy, Ricci_xz,
        Ricci_yy, Ricci_yz, Ricci_zz,
        phys_Rxx, phys_Rxy, phys_Rxz,
        phys_Ryy, phys_Ryz, phys_Rzz,
        trK_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs,
        Ayy_rhs, Ayz_rhs, Azz_rhs,
        symmetry, sst);
    abort_on_cuda(cudaGetLastError(), "bssn_shell_curvature_rhs_probe_kernel");
}

void gpu_bssn_shell_transport_ko_probe_launch(
    cudaStream_t &stream,
    int* ex,
    double* crho, double* sigma, double* R,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz,
    double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz,
    double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs,
    double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs,
    double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs,
    double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    int symmetry, int lev, double eps, int sst
) {
    (void)lev;
    const int count = shell_rhs_cell_count(ex);
    bssn_shell_transport_ko_probe_kernel<<<shell_rhs_grid_size(count), kGuardBlockSize, 0, stream>>>(
        ex[0], ex[1], ex[2],
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        chi, trK,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Axx, Axy, Axz,
        Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz,
        Lap,
        betax, betay, betaz,
        dtSfx, dtSfy, dtSfz,
        chi_rhs, trK_rhs,
        gxx_rhs, gxy_rhs, gxz_rhs,
        gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs,
        Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs,
        Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        symmetry, eps, sst);
    abort_on_cuda(cudaGetLastError(), "bssn_shell_transport_ko_probe_kernel");
}

void gpu_bssn_shell_constraint_probe_launch(
    cudaStream_t &stream,
    int* ex,
    double* crho, double* sigma, double* R,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz,
    double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* phys_Rxx, double* phys_Rxy, double* phys_Rxz,
    double* phys_Ryy, double* phys_Ryz, double* phys_Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, int sst
) {
    (void)lev;
    const int count = shell_rhs_cell_count(ex);
    bssn_shell_constraint_probe_kernel<<<shell_rhs_grid_size(count), kGuardBlockSize, 0, stream>>>(
        ex[0], ex[1], ex[2],
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        chi, trK,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Axx, Axy, Axz,
        Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz,
        rho, Sx, Sy, Sz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz,
        phys_Rxx, phys_Rxy, phys_Rxz,
        phys_Ryy, phys_Ryz, phys_Rzz,
        ham_Res, movx_Res, movy_Res, movz_Res,
        Gmx_Res, Gmy_Res, Gmz_Res,
        symmetry, sst);
    abort_on_cuda(cudaGetLastError(), "bssn_shell_constraint_probe_kernel");
}

void gpu_bssn_shell_kodis_probe_launch(
    cudaStream_t &stream,
    int* ex,
    double* crho, double* sigma, double* R,
    double* f, double* f_rhs,
    int symmetry, double eps, int sst,
    double sym1, double sym2, double sym3
) {
    const int count = shell_rhs_cell_count(ex);
    bssn_shell_kodis_probe_kernel<<<shell_rhs_grid_size(count), kGuardBlockSize, 0, stream>>>(
        ex[0], ex[1], ex[2],
        crho, sigma, R,
        f, f_rhs,
        symmetry, eps, sst, sym1, sym2, sym3);
    abort_on_cuda(cudaGetLastError(), "bssn_shell_kodis_probe_kernel");
}

void gpu_bssn_shell_full_rhs_probe_launch(
    cudaStream_t &stream,
    int* ex,
    double* crho, double* sigma, double* R,
    double* x, double* y, double* z,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* drhodxx, double* drhodxy, double* drhodxz,
    double* drhodyy, double* drhodyz, double* drhodzz,
    double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
    double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
    double* dRdxx, double* dRdxy, double* dRdxz,
    double* dRdyy, double* dRdyz, double* dRdzz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz,
    double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz,
    double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs,
    double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs,
    double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs,
    double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz,
    double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int sst, int co,
    double* TZ, double* TZ_rhs
) {
    const int count = shell_rhs_cell_count(ex);
    double* Aupxx = alloc_shell_temp(count, "cudaMalloc shell full RHS Aupxx");
    double* Aupxy = alloc_shell_temp(count, "cudaMalloc shell full RHS Aupxy");
    double* Aupxz = alloc_shell_temp(count, "cudaMalloc shell full RHS Aupxz");
    double* Aupyy = alloc_shell_temp(count, "cudaMalloc shell full RHS Aupyy");
    double* Aupyz = alloc_shell_temp(count, "cudaMalloc shell full RHS Aupyz");
    double* Aupzz = alloc_shell_temp(count, "cudaMalloc shell full RHS Aupzz");
    double* Ricci_xx = alloc_shell_temp(count, "cudaMalloc shell full RHS Ricci_xx");
    double* Ricci_xy = alloc_shell_temp(count, "cudaMalloc shell full RHS Ricci_xy");
    double* Ricci_xz = alloc_shell_temp(count, "cudaMalloc shell full RHS Ricci_xz");
    double* Ricci_yy = alloc_shell_temp(count, "cudaMalloc shell full RHS Ricci_yy");
    double* Ricci_yz = alloc_shell_temp(count, "cudaMalloc shell full RHS Ricci_yz");
    double* Ricci_zz = alloc_shell_temp(count, "cudaMalloc shell full RHS Ricci_zz");

    gpu_bssn_shell_metric_rhs_probe_launch(
        stream, ex,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        chi, trK,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Axx, Axy, Axz,
        Ayy, Ayz, Azz,
        Lap,
        betax, betay, betaz,
        chi_rhs,
        gxx_rhs, gxy_rhs, gxz_rhs,
        gyy_rhs, gyz_rhs, gzz_rhs,
        symmetry, lev, sst);

    gpu_bssn_shell_connection_probe_launch(
        stream, ex,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Axx, Axy, Axz,
        Ayy, Ayz, Azz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz,
        Aupxx, Aupxy, Aupxz,
        Aupyy, Aupyz, Aupzz,
        symmetry, lev, sst);

    gpu_bssn_shell_ricci_probe_launch(
        stream, ex,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Gamx, Gamy, Gamz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz,
        Ricci_xx, Ricci_xy, Ricci_xz,
        Ricci_yy, Ricci_yz, Ricci_zz,
        symmetry, lev, sst);

    gpu_bssn_shell_gauge_rhs_probe_launch(
        stream, ex,
        crho, sigma, R,
        x, y, z,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        chi, trK,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Gamx, Gamy, Gamz,
        Lap,
        betax, betay, betaz,
        dtSfx, dtSfy, dtSfz,
        Sx, Sy, Sz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz,
        Aupxx, Aupxy, Aupxz,
        Aupyy, Aupyz, Aupzz,
        Gamx_rhs, Gamy_rhs, Gamz_rhs,
        Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        symmetry, lev, sst);

    gpu_bssn_shell_curvature_rhs_probe_launch(
        stream, ex,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        chi, trK,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Axx, Axy, Axz,
        Ayy, Ayz, Azz,
        Lap,
        betax, betay, betaz,
        rho,
        Sxx, Sxy, Sxz,
        Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz,
        Ricci_xx, Ricci_xy, Ricci_xz,
        Ricci_yy, Ricci_yz, Ricci_zz,
        Rxx, Rxy, Rxz,
        Ryy, Ryz, Rzz,
        trK_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs,
        Ayy_rhs, Ayz_rhs, Azz_rhs,
        symmetry, lev, sst);

    // Z4C corrections need fresh constraint residuals in every RK stage
    // (the Gamma-damping and TZ sources use them), not only in the predictor.
    if (co == 0 || (TZ != 0 && TZ_rhs != 0)) {
        gpu_bssn_shell_constraint_probe_launch(
            stream, ex,
            crho, sigma, R,
            drhodx, drhody, drhodz,
            dsigmadx, dsigmady, dsigmadz,
            dRdx, dRdy, dRdz,
            chi, trK,
            dxx, gxy, gxz,
            dyy, gyz, dzz,
            Axx, Axy, Axz,
            Ayy, Ayz, Azz,
            Gamx, Gamy, Gamz,
            rho, Sx, Sy, Sz,
            Gamxxx, Gamxxy, Gamxxz,
            Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz,
            Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz,
            Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz,
            Ryy, Ryz, Rzz,
            ham_Res, movx_Res, movy_Res, movz_Res,
            Gmx_Res, Gmy_Res, Gmz_Res,
            symmetry, lev, sst);
    }

#if (ABEtype == 2)
    // Insert the Z4C theta corrections between the constraint probe and the
    // transport/KO probe so the algebraic part reads the transport-free
    // trK_rhs, matching the Cartesian gpu_compute_rhs_z4c_launch ordering.
    if (TZ != 0 && TZ_rhs != 0) {
        gpu_apply_z4c_shell_corrections_launch(
            stream, ex,
            crho, sigma, R,
            drhodx, drhody, drhodz,
            dsigmadx, dsigmady, dsigmadz,
            dRdx, dRdy, dRdz,
            chi, trK, TZ,
            dxx, gxy, gxz,
            dyy, gyz, dzz,
            Axx, Axy, Axz,
            Ayy, Ayz, Azz,
            Lap,
            betax, betay, betaz,
            rho, Sxx, Sxy, Sxz,
            Syy, Syz, Szz,
            ham_Res, Gmx_Res, Gmy_Res, Gmz_Res,
            chi_rhs, trK_rhs,
            Axx_rhs, Axy_rhs, Axz_rhs,
            Ayy_rhs, Ayz_rhs, Azz_rhs,
            Gamx_rhs, Gamy_rhs, Gamz_rhs,
            dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
            TZ_rhs,
            symmetry, eps, sst);
    }
#endif

    gpu_bssn_shell_transport_ko_probe_launch(
        stream, ex,
        crho, sigma, R,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        chi, trK,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Axx, Axy, Axz,
        Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz,
        Lap,
        betax, betay, betaz,
        dtSfx, dtSfy, dtSfz,
        chi_rhs, trK_rhs,
        gxx_rhs, gxy_rhs, gxz_rhs,
        gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs,
        Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs,
        Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        symmetry, lev, eps, sst);

    free_shell_temp(Aupxx, "cudaFree shell full RHS Aupxx");
    free_shell_temp(Aupxy, "cudaFree shell full RHS Aupxy");
    free_shell_temp(Aupxz, "cudaFree shell full RHS Aupxz");
    free_shell_temp(Aupyy, "cudaFree shell full RHS Aupyy");
    free_shell_temp(Aupyz, "cudaFree shell full RHS Aupyz");
    free_shell_temp(Aupzz, "cudaFree shell full RHS Aupzz");
    free_shell_temp(Ricci_xx, "cudaFree shell full RHS Ricci_xx");
    free_shell_temp(Ricci_xy, "cudaFree shell full RHS Ricci_xy");
    free_shell_temp(Ricci_xz, "cudaFree shell full RHS Ricci_xz");
    free_shell_temp(Ricci_yy, "cudaFree shell full RHS Ricci_yy");
    free_shell_temp(Ricci_yz, "cudaFree shell full RHS Ricci_yz");
    free_shell_temp(Ricci_zz, "cudaFree shell full RHS Ricci_zz");
}

void gpu_compute_rhs_bssn_ss_launch(
    cudaStream_t &stream,
    int* ex, double T,
    double* crho, double* sigma, double* R,
    double* x, double* y, double* z,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* drhodxx, double* drhodxy, double* drhodxz,
    double* drhodyy, double* drhodyz, double* drhodzz,
    double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
    double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
    double* dRdxx, double* dRdxy, double* dRdxz,
    double* dRdyy, double* dRdyz, double* dRdzz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz,
    double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz,
    double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz,
    double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs,
    double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs,
    double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs,
    double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz,
    double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz,
    double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int sst, int co,
    double* TZ, double* TZ_rhs
) {
    (void)T;
    gpu_bssn_shell_full_rhs_probe_launch(
        stream,
        ex,
        crho, sigma, R,
        x, y, z,
        drhodx, drhody, drhodz,
        dsigmadx, dsigmady, dsigmadz,
        dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz,
        drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz,
        dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz,
        dRdyy, dRdyz, dRdzz,
        chi, trK,
        dxx, gxy, gxz,
        dyy, gyz, dzz,
        Axx, Axy, Axz,
        Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz,
        Lap,
        betax, betay, betaz,
        dtSfx, dtSfy, dtSfz,
        chi_rhs, trK_rhs,
        gxx_rhs, gxy_rhs, gxz_rhs,
        gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs,
        Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs,
        Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        rho, Sx, Sy, Sz,
        Sxx, Sxy, Sxz,
        Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz,
        Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz,
        Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz,
        Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz,
        Ryy, Ryz, Rzz,
        ham_Res, movx_Res, movy_Res, movz_Res,
        Gmx_Res, Gmy_Res, Gmz_Res,
        symmetry, lev, eps, sst, co,
        TZ, TZ_rhs
    );
    abort_on_cuda(cudaGetLastError(), "gpu_bssn_shell_full_rhs_probe_launch via shell wrapper");
}
