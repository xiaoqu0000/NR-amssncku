#include "macrodef.h"
#include "sommerfeld_rout.h"
#include <cuda_runtime.h>
#include <math.h>

#include "gpu_manager.h"
#include "fmisc_gpu.cuh"

constexpr double ZEO = 0.0;
constexpr double ONE = 1.0;
constexpr double TWO = 2.0;
constexpr int NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2, CORRECTSTEP = 1;

__device__ __forceinline__ bool is_sommerfeld_boundary(
    int i, int j, int k,
    int ex0, int ex1, int ex2,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    int Symmetry
) {
    double dX = X[1] - X[0];
    double dY = Y[1] - Y[0];
    double dZ = Z[1] - Z[0];

    if (i == ex0 && fabs(X[ex0-1] - xmax) < dX) return true;
    if (j == ex1 && fabs(Y[ex1-1] - ymax) < dY) return true;
    if (k == ex2 && fabs(Z[ex2-1] - zmax) < dZ) return true;

    if (i == 1 && fabs(X[0] - xmin) < dX && !(Symmetry == OCTANT && fabs(xmin) < dX/2.0)) return true;
    if (j == 1 && fabs(Y[0] - ymin) < dY && !(Symmetry == OCTANT && fabs(ymin) < dY/2.0)) return true;
    if (k == 1 && fabs(Z[0] - zmin) < dZ && !(Symmetry > NO_SYMM && fabs(zmin) < dZ/2.0)) return true;

    return false;
}

__global__ void sommerfeld_rout_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double inv_dX, double inv_dY, double inv_dZ,
    bool bound_x_min, bool bound_x_max,
    bool bound_y_min, bool bound_y_max,
    bool bound_z_min, bool bound_z_max,
    double dT,
    const double* __restrict__ chi0, const double* __restrict__ Lap0,
    const double* __restrict__ f0, double* __restrict__ f,
    const double SYM1, const double SYM2, const double SYM3,
    int Symmetry,
    int precor
) {
    int i0 = blockIdx.x * blockDim.x + threadIdx.x;
    int j0 = blockIdx.y * blockDim.y + threadIdx.y;
    int k0 = blockIdx.z * blockDim.z + threadIdx.z;

    if (i0 >= ex0 || j0 >= ex1 || k0 >= ex2) return;
    (void)pex2;

    bool is_bound = false;
    if ((i0 == 0 && bound_x_min) || (i0 == ex0 - 1 && bound_x_max)) is_bound = true;
    if ((j0 == 0 && bound_y_min) || (j0 == ex1 - 1 && bound_y_max)) is_bound = true;
    if ((k0 == 0 && bound_z_min) || (k0 == ex2 - 1 && bound_z_max)) is_bound = true;

    if (!is_bound) return;

    size_t idx = (size_t)i0 + (size_t)pex0 * ((size_t)j0 + (size_t)pex1 * (size_t)k0);

    if (precor == 1) {
        f[idx] = f0[idx];
        return;
    }

    double R_sq = X[i0] * X[i0] + Y[j0] * Y[j0] + Z[k0] * Z[k0];
    double r = (Lap0[idx] + 1.0) * sqrt(1.0 + chi0[idx]) * dT * rsqrt(R_sq);
    double fac = 1.0 - r;

    double cx_orig[3] = {r * X[i0] * inv_dX, r * Y[j0] * inv_dY, r * Z[k0] * inv_dZ};

    int i = i0 + 1; int j = j0 + 1; int k = k0 + 1;
    int cxB_1b[3];
    cxB_1b[0] = (cx_orig[0] > 0.0) ? i - (int)trunc(cx_orig[0]) - 3 : i - (int)trunc(cx_orig[0]) - 2;
    cxB_1b[1] = (cx_orig[1] > 0.0) ? j - (int)trunc(cx_orig[1]) - 3 : j - (int)trunc(cx_orig[1]) - 2;
    cxB_1b[2] = (cx_orig[2] > 0.0) ? k - (int)trunc(cx_orig[2]) - 3 : k - (int)trunc(cx_orig[2]) - 2;

    double cx_eval[3];
    for (int m = 0; m < 3; ++m) {
        cx_eval[m] = (cx_orig[m] > 0.0) ? trunc(cx_orig[m]) - cx_orig[m] + 3.0 : trunc(cx_orig[m]) - cx_orig[m] + 2.0;
    }

    int cxT_1b[3];
    for (int m = 0; m < 3; ++m) cxT_1b[m] = cxB_1b[m] + 5;

    if (Symmetry == 0 && cxB_1b[2] < 1) { cx_eval[2] += (cxB_1b[2] - 1); cxT_1b[2] -= (cxB_1b[2] - 1); cxB_1b[2] = 1; }
    if (Symmetry < 2 && cxB_1b[1] < 1)  { cx_eval[1] += (cxB_1b[1] - 1); cxT_1b[1] -= (cxB_1b[1] - 1); cxB_1b[1] = 1; }
    if (Symmetry < 2 && cxB_1b[0] < 1)  { cx_eval[0] += (cxB_1b[0] - 1); cxT_1b[0] -= (cxB_1b[0] - 1); cxB_1b[0] = 1; }

    int ext[3] = {ex0, ex1, ex2};
    for (int m = 0; m < 3; ++m) {
        if (cxT_1b[m] > ext[m]) {
            cx_eval[m] += (cxT_1b[m] - ext[m]); cxB_1b[m] -= (cxT_1b[m] - ext[m]); cxT_1b[m] = ext[m];
        }
    }

    int cxB_0b[3] = {cxB_1b[0] - 1, cxB_1b[1] - 1, cxB_1b[2] - 1};

    double raw_wx[MAX_ORDN], raw_wy[MAX_ORDN], raw_wz[MAX_ORDN];
    get_lagrange_weights(cx_eval[0], raw_wx);
    get_lagrange_weights(cx_eval[1], raw_wy);
    get_lagrange_weights(cx_eval[2], raw_wz);

    const double SoA[3] = {SYM1, SYM2, SYM3};
    int idx_x[MAX_ORDN], idx_y[MAX_ORDN], idx_z[MAX_ORDN];
    double wx[MAX_ORDN], wy[MAX_ORDN], wz[MAX_ORDN];

    for (int m = 0; m < MAX_ORDN; ++m) {
        int g = cxB_0b[0] + m; double sign = 1.0;
#ifdef Vertex
        if (g < 0) { g = -g; sign *= SoA[0]; }
#else
        if (g < 0) { g = -g - 1; sign *= SoA[0]; }
#endif
        if (g >= ex0) { g = 0; sign = 0.0; }
        idx_x[m] = g; wx[m] = raw_wx[m] * sign;
    }
    for (int m = 0; m < MAX_ORDN; ++m) {
        int g = cxB_0b[1] + m; double sign = 1.0;
#ifdef Vertex
        if (g < 0) { g = -g; sign *= SoA[1]; }
#else
        if (g < 0) { g = -g - 1; sign *= SoA[1]; }
#endif
        if (g >= ex1) { g = 0; sign = 0.0; }
        idx_y[m] = g; wy[m] = raw_wy[m] * sign;
    }
    for (int m = 0; m < MAX_ORDN; ++m) {
        int g = cxB_0b[2] + m; double sign = 1.0;
#ifdef Vertex
        if (g < 0) { g = -g; sign *= SoA[2]; }
#else
        if (g < 0) { g = -g - 1; sign *= SoA[2]; }
#endif
        if (g >= ex2) { g = 0; sign = 0.0; }
        idx_z[m] = g; wz[m] = raw_wz[m] * sign;
    }

    double interp_val = 0.0;

    #pragma unroll
    for (int kk = 0; kk < MAX_ORDN; ++kk) {
        double w_z_val = wz[kk];
        size_t z_offset = (size_t)idx_z[kk] * (size_t)pex1;
        #pragma unroll
        for (int jj = 0; jj < MAX_ORDN; ++jj) {
            double w_zy_val = w_z_val * wy[jj];
            size_t zy_offset = (z_offset + (size_t)idx_y[jj]) * (size_t)pex0;
            #pragma unroll
            for (int ii = 0; ii < MAX_ORDN; ++ii) {
                interp_val += w_zy_val * wx[ii] * __ldg(&f0[zy_offset + (size_t)idx_x[ii]]);
            }
        }
    }

    f[idx] = interp_val * fac;
}

void gpu_sommerfeld_rout_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double dT, const double* d_chi0, const double* d_Lap0,
    const double* d_f0, double* d_f, const double SoA[3],
    int Symmetry, int precor
) {
    double dX = (xmax - xmin) / (ex[0] - 1);
    double dY = (ymax - ymin) / (ex[1] - 1);
    double dZ = (zmax - zmin) / (ex[2] - 1);

    double inv_dX = 1.0 / dX;
    double inv_dY = 1.0 / dY;
    double inv_dZ = 1.0 / dZ;

    bool bound_x_min = !(Symmetry == 2 && fabs(xmin) < dX / 2.0);
    bool bound_y_min = !(Symmetry == 2 && fabs(ymin) < dY / 2.0);
    bool bound_z_min = !(Symmetry > 0 && fabs(zmin) < dZ / 2.0);
    bool bound_x_max = true;
    bool bound_y_max = true;
    bool bound_z_max = true;

    dim3 block(8, 8, 4);
    dim3 grid((ex[0] + block.x - 1) / block.x,
              (ex[1] + block.y - 1) / block.y,
              (ex[2] + block.z - 1) / block.z);

    sommerfeld_rout_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_X, d_Y, d_Z,
        inv_dX, inv_dY, inv_dZ,
        bound_x_min, bound_x_max, bound_y_min, bound_y_max, bound_z_min, bound_z_max,
        dT, d_chi0, d_Lap0, d_f0, d_f, SoA[0], SoA[1], SoA[2], Symmetry, precor
    );
}

__global__ void sommerfeld_routbam_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double* __restrict__ f_rhs,
    const double* __restrict__ f0,
    double velocity,
    const double SYM1, const double SYM2, const double SYM3,
    int Symmetry
) {
    int gx = blockIdx.x * blockDim.x + threadIdx.x;
    int gy = blockIdx.y * blockDim.y + threadIdx.y;
    int gz = blockIdx.z * blockDim.z + threadIdx.z;

    if (gx >= ex0 || gy >= ex1 || gz >= ex2) return;
    (void)pex2;

    int i = gx + 1;
    int j = gy + 1;
    int k = gz + 1;

    if (!is_sommerfeld_boundary(i, j, k, ex0, ex1, ex2, X, Y, Z, xmin, ymin, zmin, xmax, ymax, zmax, Symmetry)) {
        return;
    }

    size_t idx = (size_t)gx + (size_t)pex0 * ((size_t)gy + (size_t)pex1 * (size_t)gz);
    int ext[3] = {ex0, ex1, ex2};
    int pext[3] = {pex0, pex1, pex2};
    const double SoA[3] = {SYM1, SYM2, SYM3};

    double dX = X[1] - X[0];
    double dY = Y[1] - Y[0];
    double dZ = Z[1] - Z[0];

    double d2dx = 0.5 / dX;
    double d2dy = 0.5 / dY;
    double d2dz = 0.5 / dZ;

    int imin_0b = 0, jmin_0b = 0, kmin_0b = 0;
    if (Symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin_0b = -1;
    if (Symmetry > EQ_SYMM && fabs(X[0]) < dX) imin_0b = -1;
    if (Symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin_0b = -1;

    int imax_0b = ex0 - 1;
    int jmax_0b = ex1 - 1;
    int kmax_0b = ex2 - 1;

    double R_sq = X[gx] * X[gx] + Y[gy] * Y[gy] + Z[gz] * Z[gz];
    double R_inv = rsqrt(R_sq);

    double wx = velocity * X[gx] * R_inv;
    double wy = velocity * Y[gy] * R_inv;
    double wz = velocity * Z[gz] * R_inv;

    double f_curr = f0[idx];
    double fx = 0.0, fy = 0.0, fz = 0.0;

    int sx = (wx > 0.0) ? -1 : ((wx < 0.0) ? 1 : 0);
    if (sx != 0) {
        int g1, g2;
        double c0, c1, c2;
        bool can_pref = (sx == -1) ? (gx - 2 >= imin_0b) : (gx + 2 <= imax_0b);
        bool can_cent = (gx - 1 >= imin_0b) && (gx + 1 <= imax_0b);

        if (can_pref) {
            g1 = gx + sx; g2 = gx + 2 * sx;
            c0 = -3.0 * sx; c1 = 4.0 * sx; c2 = -1.0 * sx;
        } else if (can_cent) {
            g1 = gx + sx; g2 = gx - sx;
            c0 = 0.0; c1 = 1.0 * sx; c2 = -1.0 * sx;
        } else {
            g1 = gx - sx; g2 = gx - 2 * sx;
            c0 = 3.0 * sx; c1 = -4.0 * sx; c2 = 1.0 * sx;
        }

        double f1 = fetch_f_sym_0b_pad(f0, ext, pext, g1, gy, gz, SoA);
        double f2 = fetch_f_sym_0b_pad(f0, ext, pext, g2, gy, gz, SoA);
        fx = d2dx * (c0 * f_curr + c1 * f1 + c2 * f2);
    }

    int sy = (wy > 0.0) ? -1 : ((wy < 0.0) ? 1 : 0);
    if (sy != 0) {
        int g1, g2;
        double c0, c1, c2;
        bool can_pref = (sy == -1) ? (gy - 2 >= jmin_0b) : (gy + 2 <= jmax_0b);
        bool can_cent = (gy - 1 >= jmin_0b) && (gy + 1 <= jmax_0b);

        if (can_pref) {
            g1 = gy + sy; g2 = gy + 2 * sy;
            c0 = -3.0 * sy; c1 = 4.0 * sy; c2 = -1.0 * sy;
        } else if (can_cent) {
            g1 = gy + sy; g2 = gy - sy;
            c0 = 0.0; c1 = 1.0 * sy; c2 = -1.0 * sy;
        } else {
            g1 = gy - sy; g2 = gy - 2 * sy;
            c0 = 3.0 * sy; c1 = -4.0 * sy; c2 = 1.0 * sy;
        }
        double f1 = fetch_f_sym_0b_pad(f0, ext, pext, gx, g1, gz, SoA);
        double f2 = fetch_f_sym_0b_pad(f0, ext, pext, gx, g2, gz, SoA);
        fy = d2dy * (c0 * f_curr + c1 * f1 + c2 * f2);
    }

    int sz = (wz > 0.0) ? -1 : ((wz < 0.0) ? 1 : 0);
    if (sz != 0) {
        int g1, g2;
        double c0, c1, c2;
        bool can_pref = (sz == -1) ? (gz - 2 >= kmin_0b) : (gz + 2 <= kmax_0b);
        bool can_cent = (gz - 1 >= kmin_0b) && (gz + 1 <= kmax_0b);

        if (can_pref) {
            g1 = gz + sz; g2 = gz + 2 * sz;
            c0 = -3.0 * sz; c1 = 4.0 * sz; c2 = -1.0 * sz;
        } else if (can_cent) {
            g1 = gz + sz; g2 = gz - sz;
            c0 = 0.0; c1 = 1.0 * sz; c2 = -1.0 * sz;
        } else {
            g1 = gz - sz; g2 = gz - 2 * sz;
            c0 = 3.0 * sz; c1 = -4.0 * sz; c2 = 1.0 * sz;
        }
        double f1 = fetch_f_sym_0b_pad(f0, ext, pext, gx, gy, g1, SoA);
        double f2 = fetch_f_sym_0b_pad(f0, ext, pext, gx, gy, g2, SoA);
        fz = d2dz * (c0 * f_curr + c1 * f1 + c2 * f2);
    }

    f_rhs[idx] = -velocity * (fx * X[gx] + fy * Y[gy] + fz * Z[gz] + f_curr) * R_inv;
}

void gpu_sommerfeld_routbam_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double* d_f_rhs, const double* d_f0,
    double velocity, const double SoA[3], int Symmetry
) {
    dim3 block(8, 8, 4);
    dim3 grid((ex[0] + block.x - 1) / block.x,
              (ex[1] + block.y - 1) / block.y,
              (ex[2] + block.z - 1) / block.z);

    sommerfeld_routbam_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_X, d_Y, d_Z, xmin, ymin, zmin, xmax, ymax, zmax,
        d_f_rhs, d_f0, velocity, SoA[0], SoA[1], SoA[2], Symmetry
    );
}

__global__ void fused_sommerfeld_rout_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double inv_dX, double inv_dY, double inv_dZ,
    bool bound_x_min, bool bound_x_max, 
    bool bound_y_min, bool bound_y_max, 
    bool bound_z_min, bool bound_z_max,
    double dT,
    const double* __restrict__ chi0, const double* __restrict__ Lap0,
    FusedSommerfeldVars vars,
    int Symmetry, int precor
) {
    int i0 = blockIdx.x * blockDim.x + threadIdx.x;
    int j0 = blockIdx.y * blockDim.y + threadIdx.y;
    int k0 = blockIdx.z * blockDim.z + threadIdx.z;

    if (i0 >= ex0 || j0 >= ex1 || k0 >= ex2) return;

    bool is_bound = false;
    if ((i0 == 0 && bound_x_min) || (i0 == ex0 - 1 && bound_x_max)) is_bound = true;
    if ((j0 == 0 && bound_y_min) || (j0 == ex1 - 1 && bound_y_max)) is_bound = true;
    if ((k0 == 0 && bound_z_min) || (k0 == ex2 - 1 && bound_z_max)) is_bound = true;

    if (!is_bound) return;

    int ex_idx = i0 + pex0 * (j0 + pex1 * k0);

    if (precor == 1) {
        for (int v = 0; v < vars.num_vars; ++v) {
            vars.f_out[v][ex_idx] = vars.f0[v][ex_idx];
        }
        return;
    }

    // 几何计算提取到外层，所有变量复用
    double R_sq = X[i0] * X[i0] + Y[j0] * Y[j0] + Z[k0] * Z[k0];
    double r = (Lap0[ex_idx] + 1.0) * sqrt(1.0 + chi0[ex_idx]) * dT * rsqrt(R_sq);
    double fac = 1.0 - r;

    double cx_orig[3] = {r * X[i0] * inv_dX, r * Y[j0] * inv_dY, r * Z[k0] * inv_dZ};
    
    int i = i0 + 1; int j = j0 + 1; int k = k0 + 1;
    int cxB_1b[3];
    cxB_1b[0] = (cx_orig[0] > 0.0) ? i - (int)trunc(cx_orig[0]) - 3 : i - (int)trunc(cx_orig[0]) - 2;
    cxB_1b[1] = (cx_orig[1] > 0.0) ? j - (int)trunc(cx_orig[1]) - 3 : j - (int)trunc(cx_orig[1]) - 2;
    cxB_1b[2] = (cx_orig[2] > 0.0) ? k - (int)trunc(cx_orig[2]) - 3 : k - (int)trunc(cx_orig[2]) - 2;

    double cx_eval[3];
    for (int m = 0; m < 3; ++m) {
        cx_eval[m] = (cx_orig[m] > 0.0) ? trunc(cx_orig[m]) - cx_orig[m] + 3.0 : trunc(cx_orig[m]) - cx_orig[m] + 2.0;
    }

    int cxT_1b[3];
    for (int m = 0; m < 3; ++m) cxT_1b[m] = cxB_1b[m] + 5;

    if (Symmetry == 0 && cxB_1b[2] < 1) { cx_eval[2] += (cxB_1b[2] - 1); cxT_1b[2] -= (cxB_1b[2] - 1); cxB_1b[2] = 1; }
    if (Symmetry < 2 && cxB_1b[1] < 1)  { cx_eval[1] += (cxB_1b[1] - 1); cxT_1b[1] -= (cxB_1b[1] - 1); cxB_1b[1] = 1; }
    if (Symmetry < 2 && cxB_1b[0] < 1)  { cx_eval[0] += (cxB_1b[0] - 1); cxT_1b[0] -= (cxB_1b[0] - 1); cxB_1b[0] = 1; }

    int ext[3] = {ex0, ex1, ex2};
    for (int m = 0; m < 3; ++m) {
        if (cxT_1b[m] > ext[m]) {
            cx_eval[m] += (cxT_1b[m] - ext[m]); cxB_1b[m] -= (cxT_1b[m] - ext[m]); cxT_1b[m] = ext[m];
        }
    }

    int cxB_0b[3] = {cxB_1b[0] - 1, cxB_1b[1] - 1, cxB_1b[2] - 1};

    double raw_wx[MAX_ORDN], raw_wy[MAX_ORDN], raw_wz[MAX_ORDN];
    get_lagrange_weights(cx_eval[0], raw_wx);
    get_lagrange_weights(cx_eval[1], raw_wy);
    get_lagrange_weights(cx_eval[2], raw_wz);

    // 开始变量循环
    for (int v = 0; v < vars.num_vars; ++v) {
        double SoA[3] = {vars.SoA0[v], vars.SoA1[v], vars.SoA2[v]};
        int idx_x[MAX_ORDN], idx_y[MAX_ORDN], idx_z[MAX_ORDN];
        double wx[MAX_ORDN], wy[MAX_ORDN], wz[MAX_ORDN];

        for (int m = 0; m < MAX_ORDN; ++m) {
            int g = cxB_0b[0] + m; double sign = 1.0;
#ifdef Vertex
            if (g < 0) { g = -g; sign *= SoA[0]; }
#else
            if (g < 0) { g = -g - 1; sign *= SoA[0]; }
#endif
            if (g >= ex0) { g = 0; sign = 0.0; }
            idx_x[m] = g; wx[m] = raw_wx[m] * sign;
        }
        for (int m = 0; m < MAX_ORDN; ++m) {
            int g = cxB_0b[1] + m; double sign = 1.0;
#ifdef Vertex
            if (g < 0) { g = -g; sign *= SoA[1]; }
#else
            if (g < 0) { g = -g - 1; sign *= SoA[1]; }
#endif
            if (g >= ex1) { g = 0; sign = 0.0; }
            idx_y[m] = g; wy[m] = raw_wy[m] * sign;
        }
        for (int m = 0; m < MAX_ORDN; ++m) {
            int g = cxB_0b[2] + m; double sign = 1.0;
#ifdef Vertex
            if (g < 0) { g = -g; sign *= SoA[2]; }
#else
            if (g < 0) { g = -g - 1; sign *= SoA[2]; }
#endif
            if (g >= ex2) { g = 0; sign = 0.0; }
            idx_z[m] = g; wz[m] = raw_wz[m] * sign;
        }

        double interp_val = 0.0;
        const double* f0_v = vars.f0[v];

        #pragma unroll
        for (int kk = 0; kk < MAX_ORDN; ++kk) {
            double w_z_val = wz[kk];
            size_t z_offset = (size_t)idx_z[kk] * (size_t)pex1;
            #pragma unroll
            for (int jj = 0; jj < MAX_ORDN; ++jj) {
                double w_zy_val = w_z_val * wy[jj];
                size_t zy_offset = (z_offset + (size_t)idx_y[jj]) * (size_t)pex0;
                #pragma unroll
                for (int ii = 0; ii < MAX_ORDN; ++ii) {
                    interp_val += w_zy_val * wx[ii] * __ldg(&f0_v[zy_offset + (size_t)idx_x[ii]]);
                }
            }
        }
        vars.f_out[v][ex_idx] = interp_val * fac;
    }
}

void gpu_fused_sommerfeld_rout_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double dT, const double* d_chi0, const double* d_Lap0,
    FusedSommerfeldVars vars,
    int Symmetry, int precor
) {
    if (vars.num_vars == 0) return;
    double dX = (xmax - xmin) / (ex[0] - 1);
    double dY = (ymax - ymin) / (ex[1] - 1);
    double dZ = (zmax - zmin) / (ex[2] - 1);

    bool bound_x_min = !(Symmetry == 2 && fabs(xmin) < dX / 2.0);
    bool bound_y_min = !(Symmetry == 2 && fabs(ymin) < dY / 2.0);
    bool bound_z_min = !(Symmetry > 0 && fabs(zmin) < dZ / 2.0);

    dim3 block(8, 8, 4);
    dim3 grid((ex[0] + block.x - 1) / block.x, 
              (ex[1] + block.y - 1) / block.y, 
              (ex[2] + block.z - 1) / block.z);

    fused_sommerfeld_rout_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], d_X, d_Y, d_Z,
        1.0/dX, 1.0/dY, 1.0/dZ,
        bound_x_min, true, bound_y_min, true, bound_z_min, true,
        dT, d_chi0, d_Lap0, vars, Symmetry, precor
    );
}

// 同理改造 sommerfeld_routbam_kernel
__global__ void fused_sommerfeld_routbam_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    FusedSommerfeldVars vars,
    int Symmetry
) {
    int gx = blockIdx.x * blockDim.x + threadIdx.x;
    int gy = blockIdx.y * blockDim.y + threadIdx.y;
    int gz = blockIdx.z * blockDim.z + threadIdx.z;

    if (gx >= ex0 || gy >= ex1 || gz >= ex2) return;

    if (!is_sommerfeld_boundary(gx+1, gy+1, gz+1, ex0, ex1, ex2, X, Y, Z, xmin, ymin, zmin, xmax, ymax, zmax, Symmetry)) {
        return;
    }

    int ex_idx = gx + pex0 * (gy + pex1 * gz);
    int ext[3] = {ex0, ex1, ex2};
    int pext[3] = {pex0, pex1, pex2};
    (void)pex2;

    double d2dx = 0.5 / (X[1] - X[0]);
    double d2dy = 0.5 / (Y[1] - Y[0]);
    double d2dz = 0.5 / (Z[1] - Z[0]);

    int imin_0b = 0, jmin_0b = 0, kmin_0b = 0;
    if (Symmetry > NO_SYMM && fabs(Z[0]) < (Z[1]-Z[0])) kmin_0b = -1;
    if (Symmetry > EQ_SYMM && fabs(X[0]) < (X[1]-X[0])) imin_0b = -1;
    if (Symmetry > EQ_SYMM && fabs(Y[0]) < (Y[1]-Y[0])) jmin_0b = -1;

    int imax_0b = ex0 - 1; int jmax_0b = ex1 - 1; int kmax_0b = ex2 - 1;

    double R_sq = X[gx] * X[gx] + Y[gy] * Y[gy] + Z[gz] * Z[gz];
    double R_inv = rsqrt(R_sq);

    // 开始变量循环
    for (int v = 0; v < vars.num_vars; ++v) {
        double velocity = vars.propspeed[v];
        double SoA[3] = {vars.SoA0[v], vars.SoA1[v], vars.SoA2[v]};
        const double* f0 = vars.f0[v];

        double wx = velocity * X[gx] * R_inv;
        double wy = velocity * Y[gy] * R_inv;
        double wz = velocity * Z[gz] * R_inv;

        double f_curr = f0[ex_idx]; 
        double fx = 0.0, fy = 0.0, fz = 0.0;

        // ... 保留内部的 sx, sy, sz 判断与差分逻辑，仅需替换 f0 抓取 ...
        int sx = (wx > 0.0) ? -1 : ((wx < 0.0) ? 1 : 0);
        if (sx != 0) {
            int g1, g2; double c0, c1, c2;
            bool can_pref = (sx == -1) ? (gx - 2 >= imin_0b) : (gx + 2 <= imax_0b);
            bool can_cent = (gx - 1 >= imin_0b) && (gx + 1 <= imax_0b);
            if (can_pref) { g1 = gx + sx; g2 = gx + 2 * sx; c0 = -3.0 * sx; c1 = 4.0 * sx; c2 = -1.0 * sx; } 
            else if (can_cent) { g1 = gx + sx; g2 = gx - sx; c0 = 0.0; c1 = 1.0 * sx; c2 = -1.0 * sx; } 
            else { g1 = gx - sx; g2 = gx - 2 * sx; c0 = 3.0 * sx; c1 = -4.0 * sx; c2 = 1.0 * sx; }
            fx = d2dx * (c0 * f_curr + c1 * fetch_f_sym_0b_pad(f0, ext, pext, g1, gy, gz, SoA) + c2 * fetch_f_sym_0b_pad(f0, ext, pext, g2, gy, gz, SoA));
        }
        
        int sy = (wy > 0.0) ? -1 : ((wy < 0.0) ? 1 : 0);
        if (sy != 0) {
            int g1, g2; double c0, c1, c2;
            bool can_pref = (sy == -1) ? (gy - 2 >= jmin_0b) : (gy + 2 <= jmax_0b);
            bool can_cent = (gy - 1 >= jmin_0b) && (gy + 1 <= jmax_0b);
            if (can_pref) { g1 = gy + sy; g2 = gy + 2 * sy; c0 = -3.0 * sy; c1 = 4.0 * sy; c2 = -1.0 * sy; } 
            else if (can_cent) { g1 = gy + sy; g2 = gy - sy; c0 = 0.0; c1 = 1.0 * sy; c2 = -1.0 * sy; } 
            else { g1 = gy - sy; g2 = gy - 2 * sy; c0 = 3.0 * sy; c1 = -4.0 * sy; c2 = 1.0 * sy; }
            fy = d2dy * (c0 * f_curr + c1 * fetch_f_sym_0b_pad(f0, ext, pext, gx, g1, gz, SoA) + c2 * fetch_f_sym_0b_pad(f0, ext, pext, gx, g2, gz, SoA));
        }

        int sz = (wz > 0.0) ? -1 : ((wz < 0.0) ? 1 : 0);
        if (sz != 0) {
            int g1, g2; double c0, c1, c2;
            bool can_pref = (sz == -1) ? (gz - 2 >= kmin_0b) : (gz + 2 <= kmax_0b);
            bool can_cent = (gz - 1 >= kmin_0b) && (gz + 1 <= kmax_0b);
            if (can_pref) { g1 = gz + sz; g2 = gz + 2 * sz; c0 = -3.0 * sz; c1 = 4.0 * sz; c2 = -1.0 * sz; } 
            else if (can_cent) { g1 = gz + sz; g2 = gz - sz; c0 = 0.0; c1 = 1.0 * sz; c2 = -1.0 * sz; } 
            else { g1 = gz - sz; g2 = gz - 2 * sz; c0 = 3.0 * sz; c1 = -4.0 * sz; c2 = 1.0 * sz; }
            fz = d2dz * (c0 * f_curr + c1 * fetch_f_sym_0b_pad(f0, ext, pext, gx, gy, g1, SoA) + c2 * fetch_f_sym_0b_pad(f0, ext, pext, gx, gy, g2, SoA));
        }

        vars.f_out[v][ex_idx] = -velocity * (fx * X[gx] + fy * Y[gy] + fz * Z[gz] + f_curr) * R_inv;
    }
}

void gpu_fused_sommerfeld_routbam_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    FusedSommerfeldVars vars, int Symmetry
) {
    if (vars.num_vars == 0) return;
    dim3 block(8, 8, 4);
    dim3 grid((ex[0] + block.x - 1) / block.x, 
              (ex[1] + block.y - 1) / block.y, 
              (ex[2] + block.z - 1) / block.z);

    fused_sommerfeld_routbam_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], d_X, d_Y, d_Z, xmin, ymin, zmin, xmax, ymax, zmax,
        vars, Symmetry
    );
}

__global__ void sommerfeld_routbam_ss_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double* __restrict__ f_rhs,
    const double* __restrict__ f0,
    double velocity,
    const double SYM1, const double SYM2, const double SYM3,
    int Symmetry
) {
    int gx = blockIdx.x * blockDim.x + threadIdx.x;
    int gy = blockIdx.y * blockDim.y + threadIdx.y;
    int gz = blockIdx.z * blockDim.z + threadIdx.z;

    if (gx >= ex0 || gy >= ex1 || gz >= ex2) return;
    (void)pex2;
    (void)X; (void)Y; (void)xmin; (void)ymin; (void)xmax; (void)ymax;
    (void)SYM1; (void)SYM2; (void)SYM3; (void)Symmetry;

    double dZ = Z[1] - Z[0];

    bool is_zmax = (gz + 1 == ex2 && fabs(Z[gz] - zmax) < dZ);
    bool is_zmin = (gz == 0 && fabs(Z[gz] - zmin) < dZ);
    if (!is_zmax && !is_zmin) return;

    double d2dz = 0.5 / dZ;

    int kmin_0b = 0;
    int kmax_0b = ex2 - 1;

    size_t idx = (size_t)gx + (size_t)pex0 * ((size_t)gy + (size_t)pex1 * (size_t)gz);

    if (is_zmax) {
        double f_curr = f0[idx];
        double fz = 0.0;

        if (velocity > 0.0) {
            if (gz - 2 >= kmin_0b) {
                fz = d2dz * (3.0 * f_curr - 4.0 * f0[idx - pex0*pex1] + f0[idx - 2*pex0*pex1]);
            } else if (gz - 1 >= kmin_0b) {
                fz = d2dz * (-f0[idx - pex0*pex1] + f0[idx + pex0*pex1]);
            } else {
                fz = d2dz * (-f0[idx + 2*pex0*pex1] + 4.0 * f0[idx + pex0*pex1] - 3.0 * f_curr);
            }
        } else if (velocity < 0.0) {
            if (gz + 2 <= kmax_0b) {
                fz = d2dz * (-f0[idx + 2*pex0*pex1] + 4.0 * f0[idx + pex0*pex1] - 3.0 * f_curr);
            } else if (gz + 1 <= kmax_0b) {
                fz = d2dz * (-f0[idx - pex0*pex1] + f0[idx + pex0*pex1]);
            } else {
                fz = d2dz * (3.0 * f_curr - 4.0 * f0[idx - pex0*pex1] + f0[idx - 2*pex0*pex1]);
            }
        }

        f_rhs[idx] = -velocity * (fz + f_curr / Z[gz]);
    }

    if (is_zmin) {
        f_rhs[idx] = 0.0;
    }
}

void gpu_sommerfeld_routbam_ss_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double* d_f_rhs, const double* d_f0,
    double velocity, const double SoA[3], int Symmetry
) {
    dim3 block(8, 8, 4);
    dim3 grid((ex[0] + block.x - 1) / block.x,
              (ex[1] + block.y - 1) / block.y,
              (ex[2] + block.z - 1) / block.z);

    sommerfeld_routbam_ss_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_X, d_Y, d_Z, xmin, ymin, zmin, xmax, ymax, zmax,
        d_f_rhs, d_f0, velocity, SoA[0], SoA[1], SoA[2], Symmetry
    );
}

__global__ void fused_sommerfeld_routbam_ss_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    FusedSommerfeldVars vars,
    int Symmetry
) {
    int gx = blockIdx.x * blockDim.x + threadIdx.x;
    int gy = blockIdx.y * blockDim.y + threadIdx.y;
    int gz = blockIdx.z * blockDim.z + threadIdx.z;

    if (gx >= ex0 || gy >= ex1 || gz >= ex2) return;
    (void)pex2;
    (void)X; (void)Y; (void)xmin; (void)ymin; (void)xmax; (void)ymax;
    (void)Symmetry;

    double dZ = Z[1] - Z[0];

    bool is_zmax = (gz + 1 == ex2 && fabs(Z[gz] - zmax) < dZ);
    bool is_zmin = (gz == 0 && fabs(Z[gz] - zmin) < dZ);
    if (!is_zmax && !is_zmin) return;

    double d2dz = 0.5 / dZ;

    int kmin_0b = 0;
    int kmax_0b = ex2 - 1;

    size_t ex_idx = (size_t)gx + (size_t)pex0 * ((size_t)gy + (size_t)pex1 * (size_t)gz);

    if (is_zmax) {
        double z_val = Z[gz];
        size_t sz_step = (size_t)pex0 * (size_t)pex1;

        for (int v = 0; v < vars.num_vars; ++v) {
            double velocity = vars.propspeed[v];
            const double* f0_v = vars.f0[v];
            double f_curr = f0_v[ex_idx];
            double fz = 0.0;

            if (velocity > 0.0) {
                if (gz - 2 >= kmin_0b) {
                    fz = d2dz * (3.0 * f_curr - 4.0 * f0_v[ex_idx - sz_step] + f0_v[ex_idx - 2*sz_step]);
                } else if (gz - 1 >= kmin_0b) {
                    fz = d2dz * (-f0_v[ex_idx - sz_step] + f0_v[ex_idx + sz_step]);
                } else {
                    fz = d2dz * (-f0_v[ex_idx + 2*sz_step] + 4.0 * f0_v[ex_idx + sz_step] - 3.0 * f_curr);
                }
            } else if (velocity < 0.0) {
                if (gz + 2 <= kmax_0b) {
                    fz = d2dz * (-f0_v[ex_idx + 2*sz_step] + 4.0 * f0_v[ex_idx + sz_step] - 3.0 * f_curr);
                } else if (gz + 1 <= kmax_0b) {
                    fz = d2dz * (-f0_v[ex_idx - sz_step] + f0_v[ex_idx + sz_step]);
                } else {
                    fz = d2dz * (3.0 * f_curr - 4.0 * f0_v[ex_idx - sz_step] + f0_v[ex_idx - 2*sz_step]);
                }
            }

            vars.f_out[v][ex_idx] = -velocity * (fz + f_curr / z_val);
        }
    }

    if (is_zmin) {
        for (int v = 0; v < vars.num_vars; ++v) {
            vars.f_out[v][ex_idx] = 0.0;
        }
    }
}

void gpu_fused_sommerfeld_routbam_ss_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    FusedSommerfeldVars vars, int Symmetry
) {
    if (vars.num_vars == 0) return;
    dim3 block(8, 8, 4);
    dim3 grid((ex[0] + block.x - 1) / block.x,
              (ex[1] + block.y - 1) / block.y,
              (ex[2] + block.z - 1) / block.z);

    fused_sommerfeld_routbam_ss_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], d_X, d_Y, d_Z,
        xmin, ymin, zmin, xmax, ymax, zmax,
        vars, Symmetry
    );
}

// ===== falloff_ss GPU kernel =====
// Applies power-law falloff f ~ r^(-n) at the outer radial boundary
// of a Shell patch. Only the outermost k-plane (k = extc2-1) is modified.
// Formula: f(k) = f(k-1) * ((R_avg/(n*dR) - 1) / (R_avg/(n*dR) + 1))
//   where R_avg = (R[k] + R[k-1]) / 2, dR = R[1] - R[0]
__global__ void falloff_ss_kernel(
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1,
    double* f, const double* R,
    int n, double zmax
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= extc0 || j >= extc1) return;

    int k = extc2 - 1;
    double dR = R[1] - R[0];
    if (fabs(R[k] - zmax) >= dR) return;

    int idx    = i + j * pextc0 +  k      * pextc0 * pextc1;
    int idx_km1 = i + j * pextc0 + (k - 1) * pextc0 * pextc1;

    double R_avg = 0.5 * (R[k] + R[k - 1]);
    double ratio = R_avg / (n * dR);
    double factor = (ratio - 1.0) / (ratio + 1.0);
    f[idx] = f[idx_km1] * factor;
}

void gpu_falloff_ss_launch(
    cudaStream_t &stream,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1,
    double* d_f, const double* d_R,
    int n, double zmax
) {
    dim3 block(16, 16);
    dim3 grid((extc0 + block.x - 1) / block.x,
              (extc1 + block.y - 1) / block.y);

    falloff_ss_kernel<<<grid, block, 0, stream>>>(
        extc0, extc1, extc2, pextc0, pextc1, d_f, d_R, n, zmax
    );
}