#ifndef FMISC_GPU_CUH
#define FMISC_GPU_CUH

#include <cuda_runtime.h>
#include <math.h>

#include <iostream>
#include "gpu_manager.h"

#ifdef ghost_width
#define MAX_ORDN (((6) > (2 * ghost_width)) ? (6) : (2 * ghost_width))
#else
#define MAX_ORDN 10
#endif
#ifndef GPU_DEBUG_PRINT
#define GPU_DEBUG_PRINT 0
#endif

__device__ __forceinline__ double d_symmetry_bd_0b(
    int ord, 
    int extc0, int extc1, int extc2,
    const double* __restrict__ func,
    int i, int j, int k, // 0-based coordinates
    double SoA0, double SoA1, double SoA2
) {
    if (i < -ord || i >= extc0) return 0.0;
    if (j < -ord || j >= extc1) return 0.0;
    if (k < -ord || k >= extc2) return 0.0;

    double factor = 1.0;

#ifdef Vertex
    if (i < 0) { i = -i; factor *= SoA0; }
    if (j < 0) { j = -j; factor *= SoA1; }
    if (k < 0) { k = -k; factor *= SoA2; }
#else
    if (i < 0) { i = -i - 1; factor *= SoA0; }
    if (j < 0) { j = -j - 1; factor *= SoA1; }
    if (k < 0) { k = -k - 1; factor *= SoA2; }
#endif

    if (i >= extc0 || j >= extc1 || k >= extc2) return 0.0;

    return func[(k * extc1 + j) * extc0 + i] * factor;
}

__device__ __forceinline__ double d_symmetry_bd_0b_pad(
    int ord, 
    int extc0, int extc1, int extc2,    // original shape
    int pextc0, int pextc1, int pextc2, // padding shape
    const double* __restrict__ func,
    int i, int j, int k, // 0-based logical coordinates
    double SoA0, double SoA1, double SoA2
) {
    if (i < -ord || i >= extc0 || j < -ord || j >= extc1 || k < -ord || k >= extc2) {
        return 0.0;
    }

    double factor = 1.0;

#ifdef Vertex
    factor *= (i < 0) ? SoA0 : 1.0;
    i = (i < 0) ? (-i) : i;
    factor *= (j < 0) ? SoA1 : 1.0;
    j = (j < 0) ? (-j) : j;
    factor *= (k < 0) ? SoA2 : 1.0;
    k = (k < 0) ? (-k) : k;
#else
    factor *= (i < 0) ? SoA0 : 1.0;
    i = (i < 0) ? (-i - 1) : i;
    factor *= (j < 0) ? SoA1 : 1.0;
    j = (j < 0) ? (-j - 1) : j;
    factor *= (k < 0) ? SoA2 : 1.0;
    k = (k < 0) ? (-k - 1) : k;
#endif

    return func[(k * pextc1 + j) * pextc0 + i] * factor;
}

__device__ __forceinline__ void gpu_stop() {
#if GPU_STRICT_STOP
    asm("trap;");
#endif
}

__device__ __forceinline__ double X_at_1b(const double* X, int i1b) {
	return X[i1b - 1];
}
__device__ __forceinline__ double f_at_1b(const double* f, const int ex[3], int i1b, int j1b, int k1b) {
	return f[((k1b - 1) * ex[1] + (j1b - 1)) * ex[0] + (i1b - 1)];
}

__device__ __forceinline__ double f_at_1b_pad(
    const double* f, const int ex[3], const int pex[3], int i1b, int j1b, int k1b
) {
    (void)ex; // 物理边界检查由上层函数负责
    return f[((k1b - 1) * pex[1] + (j1b - 1)) * pex[0] + (i1b - 1)];
}

__device__ __forceinline__ void get_lagrange_weights(double x, double* w) {
    constexpr int LAGRANGE_PTS = 6;
    const double inv_den[LAGRANGE_PTS] = {
        -1.0 / 120.0,  1.0 / 24.0, -1.0 / 12.0,
         1.0 / 12.0,  -1.0 / 24.0,  1.0 / 120.0
    };

    #pragma unroll
    for (int i = 0; i < LAGRANGE_PTS; ++i) {
        double num = 1.0;
        #pragma unroll
        for (int j = 0; j < LAGRANGE_PTS; ++j) {
            if (i != j) num *= (x - (double)j);
        }
        w[i] = num * inv_den[i];
    }
}

__device__ __forceinline__ double fetch_f_sym_0b(
    const double* __restrict__ f, const int ex[3],
    int gx, int gy, int gz,
    const double SoA[3]
) {
    double sign = 1.0;

#ifdef Vertex
    // Vertex: symmetry plane at index 0. Ghost -1 reflects to 1, -2 to 2, etc.
    if (gx < 0) { gx = -gx; sign *= SoA[0]; }
    if (gy < 0) { gy = -gy; sign *= SoA[1]; }
    if (gz < 0) { gz = -gz; sign *= SoA[2]; }
#else
    // Cell: symmetry plane between indices. Ghost -1 reflects to 0, -2 to 1, etc.
    if (gx < 0) { gx = -gx - 1; sign *= SoA[0]; }
    if (gy < 0) { gy = -gy - 1; sign *= SoA[1]; }
    if (gz < 0) { gz = -gz - 1; sign *= SoA[2]; }
#endif

    if (gx >= ex[0] || gy >= ex[1] || gz >= ex[2]) {
        return 0.0;
    }

    return sign * f[(gz * ex[1] + gy) * ex[0] + gx];
}

__device__ __forceinline__ double fetch_f_sym_0b_pad(
    const double* __restrict__ f, const int ex[3], const int pex[3],
    int gx, int gy, int gz, const double SoA[3]
) {
    double sign = 1.0;
#ifdef Vertex
    // Vertex: symmetry plane at index 0. Ghost -1 reflects to 1, -2 to 2, etc.
    if (gx < 0) { gx = -gx; sign *= SoA[0]; }
    if (gy < 0) { gy = -gy; sign *= SoA[1]; }
    if (gz < 0) { gz = -gz; sign *= SoA[2]; }
#else
    // Cell: symmetry plane between indices. Ghost -1 reflects to 0, -2 to 1, etc.
    if (gx < 0) { gx = -gx - 1; sign *= SoA[0]; }
    if (gy < 0) { gy = -gy - 1; sign *= SoA[1]; }
    if (gz < 0) { gz = -gz - 1; sign *= SoA[2]; }
#endif

    if (gx >= ex[0] || gy >= ex[1] || gz >= ex[2]) return 0.0;
    return sign * f[(gz * pex[1] + gy) * pex[0] + gx];
}

__device__ __forceinline__ double fetch_f_sym(
    const double* f, const int ex[3],
    int gx, int gy, int gz,
    const double SoA[3])
{
    double sign = 1.0;

#ifdef Vertex
    // Vertex (1-based): symmetry plane at index 1. Ghost 0 reflects to 2, -1 to 3, etc.
    if (gx < 1) { gx = 2 - gx; sign *= SoA[0]; }
    if (gy < 1) { gy = 2 - gy; sign *= SoA[1]; }
    if (gz < 1) { gz = 2 - gz; sign *= SoA[2]; }
#else
    // Cell (1-based): symmetry plane between index 0 and 1. Ghost 0 reflects to 1, -1 to 2, etc.
    if (gx < 1) { gx = 1 - gx; sign *= SoA[0]; }
    if (gy < 1) { gy = 1 - gy; sign *= SoA[1]; }
    if (gz < 1) { gz = 1 - gz; sign *= SoA[2]; }
#endif

    if (gx > ex[0] || gy > ex[1] || gz > ex[2]) {
        return NAN;
    }

    return sign * f_at_1b(f, ex, gx, gy, gz);
}

__device__ __forceinline__ double fetch_f_sym_pad(
    const double* f, const int ex[3], const int pex[3],
    int gx, int gy, int gz,
    const double SoA[3]
) {
    double sign = 1.0;

#ifdef Vertex
    // Vertex (1-based): symmetry plane at index 1. Ghost 0 reflects to 2, -1 to 3, etc.
    if (gx < 1) { gx = 2 - gx; sign *= SoA[0]; }
    if (gy < 1) { gy = 2 - gy; sign *= SoA[1]; }
    if (gz < 1) { gz = 2 - gz; sign *= SoA[2]; }
#else
    // Cell (1-based): symmetry plane between index 0 and 1. Ghost 0 reflects to 1, -1 to 2, etc.
    if (gx < 1) { gx = 1 - gx; sign *= SoA[0]; }
    if (gy < 1) { gy = 1 - gy; sign *= SoA[1]; }
    if (gz < 1) { gz = 1 - gz; sign *= SoA[2]; }
#endif

    if (gx > ex[0] || gy > ex[1] || gz > ex[2]) {
        return NAN;
    }

    return sign * f_at_1b_pad(f, ex, pex, gx, gy, gz);
}

// ---------------------------------------------------------
// 计算标准均匀网格 [0, 1, ..., ordn-1] 上的拉格朗日插值权重
// ---------------------------------------------------------
__device__ __forceinline__ void compute_lagrange_weights(double x, int ordn, double* w) {
    if (ordn <= 0 || ordn > MAX_ORDN) return;
    for (int i = 0; i < ordn; ++i) {
        double num = 1.0;
        double den = 1.0;
        for (int j = 0; j < ordn; ++j) {
            if (i != j) {
                num *= (x - (double)j);
                den *= (double)(i - j);
            }
        }
        w[i] = num / den;
    }
}

__device__ __forceinline__ double polint_device(
    const double* __restrict__ xa,
    const double* __restrict__ ya,
    double x,
    double* dy,
    int ordn
) {
    double c[MAX_ORDN], d[MAX_ORDN], ho[MAX_ORDN];
    int ns = 0;
    if (ordn <= 0 || ordn > MAX_ORDN) {
        *dy = NAN;
        gpu_stop();
        return NAN;
    }
    double dif = fabs(x - xa[0]);

    for (int m = 0; m < ordn; ++m) {
        double dift = fabs(x - xa[m]);
        if (dift < dif) {
            ns = m;
            dif = dift;
        }
        c[m] = ya[m];
        d[m] = ya[m];
        ho[m] = xa[m] - x;
    }

    double y = ya[ns--];
    *dy = 0.0;
    for (int m = 1; m < ordn; ++m) {
        for (int i = 0; i < ordn - m; ++i) {
            double den = ho[i] - ho[i + m];
            if (den == 0.0) {
                *dy = NAN;
                return NAN;
            }
            den = (c[i + 1] - d[i]) / den;
            d[i] = ho[i + m] * den;
            c[i] = ho[i] * den;
        }
        if (2 * (ns + 1) < ordn - m) {
            *dy = c[ns + 1];
        } else {
            *dy = d[ns];
            --ns;
        }
        y += *dy;
    }
    return y;
}

__device__ __forceinline__ double polin2_device(
    const double* __restrict__ x1a,
    const double* __restrict__ x2a,
    const double* __restrict__ ya,
    double x1, double x2,
    double* dy,
    int ordn
) {
    double ymtmp[MAX_ORDN], yntmp[MAX_ORDN];
    if (ordn <= 0 || ordn > MAX_ORDN) {
        *dy = NAN;
        gpu_stop();
        return NAN;
    }
    for (int i = 0; i < ordn; ++i) {
        for (int j = 0; j < ordn; ++j) yntmp[j] = ya[i * ordn + j];
        ymtmp[i] = polint_device(x2a, yntmp, x2, dy, ordn);
    }
    return polint_device(x1a, ymtmp, x1, dy, ordn);
}

__device__ __forceinline__ double polin3_device(
    const double* __restrict__ x1a,
    const double* __restrict__ x2a,
    const double* __restrict__ x3a,
    const double* __restrict__ ya,
    double x1, double x2, double x3,
    double* dy,
    int ordn
) {
    double yatmp[MAX_ORDN * MAX_ORDN], ymtmp[MAX_ORDN], yntmp[MAX_ORDN], yqtmp[MAX_ORDN];
    if (ordn <= 0 || ordn > MAX_ORDN) {
        *dy = NAN;
        gpu_stop();
        return NAN;
    }
    for (int i = 0; i < ordn; ++i) {
        for (int j = 0; j < ordn; ++j) {
            for (int k = 0; k < ordn; ++k) yqtmp[k] = ya[(i * ordn + j) * ordn + k];
            yatmp[i * ordn + j] = polint_device(x3a, yqtmp, x3, dy, ordn);
        }
        for (int j = 0; j < ordn; ++j) yntmp[j] = yatmp[i * ordn + j];
        ymtmp[i] = polint_device(x2a, yntmp, x2, dy, ordn);
    }
    return polint_device(x1a, ymtmp, x1, dy, ordn);
}

__device__ __forceinline__ double ffact_device(int n) {
    if (n < 0) return NAN;
    double out = 1.0;
    for (int i = 1; i <= n; ++i) out *= (double)i;
    return out;
}

__device__ __forceinline__ double wigner_d_device(int l, int m, int s, double costheta) {
    int c1 = max(0, m - s);
    int c2 = min(l + m, l - s);
    double vv = 0.0;
    double sinht = sqrt((1.0 - costheta) * 0.5);
    double cosht = sqrt((1.0 + costheta) * 0.5);
    int sign0 = ((c1 / 2) * 2 == c1) ? 1 : -1;

    for (int t = c1; t <= c2; ++t) {
        int sign = ((t - c1) & 1) ? -sign0 : sign0;
        double denom = ffact_device(l + m - t) * ffact_device(l - s - t) *
                       ffact_device(t) * ffact_device(t + s - m);
        vv += sign * pow(cosht, 2 * l + m - s - 2 * t) *
              pow(sinht, 2 * t + s - m) / denom;
    }

    return vv * sqrt(ffact_device(l + m) * ffact_device(l - m) *
                     ffact_device(l + s) * ffact_device(l - s));
}

__device__ __forceinline__ void yslm_device(
    int s, int l, int m, double theta, double phi,
    double* real_out, double* imag_out
) {
    double pi = acos(-1.0);
    double rp = wigner_d_device(l, m, s, cos(theta)) * sqrt((2.0 * l + 1.0) / (4.0 * pi));
    if ((s / 2) * 2 != s) rp = -rp;
    *real_out = rp * cos(m * phi);
    *imag_out = rp * sin(m * phi);
}

__device__ __forceinline__ double f_at_0b_pad(
    const double* f, const int pex[3], int i, int j, int k
) {
    return f[(k * pex[1] + j) * pex[0] + i];
}

__device__ __forceinline__ int reflect_cell_index(int idx, int n, double* factor, double soa) {
    if (idx < 0) {
        *factor *= soa;
#ifdef Vertex
        return -idx;
#else
        return -idx - 1;
#endif
    }
    if (idx >= n) {
        *factor *= soa;
#ifdef Vertex
        return 2 * n - idx - 2;
#else
        return 2 * n - idx - 1;
#endif
    }
    return idx;
}

__device__ __forceinline__ double fetch_f_tbd_0b_pad(
    const double* __restrict__ f, const int ex[3], const int pex[3],
    int gx, int gy, int gz, const double SoA[3]
) {
    if (gx < -MAX_ORDN + 1 || gx > ex[0] + MAX_ORDN - 1) return 0.0;
    if (gy < -MAX_ORDN + 1 || gy > ex[1] + MAX_ORDN - 1) return 0.0;
    if (gz < -MAX_ORDN + 1 || gz > ex[2] + MAX_ORDN - 1) return 0.0;

    double factor = 1.0;
    int ix = reflect_cell_index(gx, ex[0], &factor, SoA[0]);
    int iy = reflect_cell_index(gy, ex[1], &factor, SoA[1]);
    int iz = reflect_cell_index(gz, ex[2], &factor, SoA[2]);

    if (ix < 0 || ix >= ex[0] || iy < 0 || iy >= ex[1] || iz < 0 || iz >= ex[2]) {
        return 0.0;
    }
    return factor * f_at_0b_pad(f, pex, ix, iy, iz);
}

__device__ __forceinline__ void shell_soa3(int sst, const double SoA[3], double SoAh[3]) {
    if (sst == 2 || sst == 3) {
        SoAh[0] = SoA[1];
        SoAh[1] = SoA[2];
        SoAh[2] = 0.0;
    } else if (sst == 4 || sst == 5) {
        SoAh[0] = SoA[0];
        SoAh[1] = SoA[2];
        SoAh[2] = 0.0;
    } else {
        SoAh[0] = SoA[0];
        SoAh[1] = SoA[1];
        SoAh[2] = (sst == -1) ? SoA[2] : 0.0;
    }
}

__device__ __forceinline__ double shell_soa1(int sst, int dumyd, const double SoA[3]) {
    if (dumyd == 1) {
        if (sst == 2 || sst == 3) return SoA[1];
        return SoA[0];
    }
    if (sst == 0 || sst == 1) return SoA[1];
    return SoA[2];
}

__device__ __forceinline__ double interp3_indexed_device(
    const int ex[3], const int pex[3],
    const double* __restrict__ f,
    const int inds[3], const double* __restrict__ coef,
    int ordn, const double SoA[3], int sst
) {
    if (ordn <= 0 || ordn > MAX_ORDN) {
        gpu_stop();
        return NAN;
    }
    double SoAh[3];
    shell_soa3(sst, SoA, SoAh);

    double val = 0.0;
    for (int k = 0; k < ordn; ++k) {
        int gz = inds[2] + k;
        double wz = coef[2 * ordn + k];
        for (int j = 0; j < ordn; ++j) {
            int gy = inds[1] + j;
            double wy = coef[ordn + j];
            for (int i = 0; i < ordn; ++i) {
                int gx = inds[0] + i;
                val += coef[i] * wy * wz * fetch_f_tbd_0b_pad(f, ex, pex, gx, gy, gz, SoAh);
            }
        }
    }
    return val;
}

__device__ __forceinline__ double interp2_indexed_device(
    const int ex[3], const int pex[3],
    const double* __restrict__ f,
    const int inds[3], const double* __restrict__ coef,
    int ordn, const double SoA[3], int sst
) {
    if (ordn <= 0 || ordn > MAX_ORDN) {
        gpu_stop();
        return NAN;
    }
    double SoAh[3];
    shell_soa3(sst, SoA, SoAh);

    double val = 0.0;
    int gz = inds[2];
    for (int j = 0; j < ordn; ++j) {
        int gy = inds[1] + j;
        double wy = coef[ordn + j];
        for (int i = 0; i < ordn; ++i) {
            int gx = inds[0] + i;
            val += coef[i] * wy * fetch_f_tbd_0b_pad(f, ex, pex, gx, gy, gz, SoAh);
        }
    }
    return val;
}

__device__ __forceinline__ double interp1_indexed_device(
    const int ex[3], const int pex[3],
    const double* __restrict__ f,
    const int inds[3], const double* __restrict__ coef,
    int ordn, const double SoA[3], int sst, int dumyd
) {
    if (ordn <= 0 || ordn > MAX_ORDN) {
        gpu_stop();
        return NAN;
    }
    double axis_soa = shell_soa1(sst, dumyd, SoA);
    double val = 0.0;

    for (int m = 0; m < ordn; ++m) {
        int gx = inds[0];
        int gy = inds[1];
        int gz = inds[2];
        double factor = 1.0;

        if (dumyd == 1) {
            gx = reflect_cell_index(inds[0] + m, ex[0], &factor, axis_soa);
        } else {
            gy = reflect_cell_index(inds[0] + m, ex[1], &factor, axis_soa);
            gx = inds[1];
        }

        if (gx >= 0 && gx < ex[0] && gy >= 0 && gy < ex[1] && gz >= 0 && gz < ex[2]) {
            val += coef[m] * factor * f_at_0b_pad(f, pex, gx, gy, gz);
        }
    }
    return val;
}

// ---------------------------------------------------------
// 优化后的核心设备端插值函数
// ---------------------------------------------------------
__device__ __forceinline__ void global_interp_device(
    const int* ex, const int* pex,
    const double* X, const double* Y, const double* Z,
    const double* f, double* f_int,
    double x1, double y1, double z1,
    int ORDN, const double* SoA, int symmetry
) {
    if (ORDN <= 0 || ORDN > MAX_ORDN) { f_int[0] = NAN; gpu_stop(); return; }

    const int NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2;
    (void)EQUATORIAL;

    int imin = 1, jmin = 1, kmin = 1;

    double dX = X_at_1b(X, imin + 1) - X_at_1b(X, imin);
    double dY = X_at_1b(Y, jmin + 1) - X_at_1b(Y, jmin);
    double dZ = X_at_1b(Z, kmin + 1) - X_at_1b(Z, kmin);

    int cxI[3];
    cxI[0] = (int)((x1 - X_at_1b(X, 1)) / dX + 0.4) + 1;
    cxI[1] = (int)((y1 - X_at_1b(Y, 1)) / dY + 0.4) + 1;
    cxI[2] = (int)((z1 - X_at_1b(Z, 1)) / dZ + 0.4) + 1;

    int cxB[3], cxT[3], cmin[3], cmax[3];
    for (int m = 0; m < 3; ++m) {
        cxB[m] = cxI[m] - ORDN / 2 + 1;
        cxT[m] = cxB[m] + ORDN - 1;
        cmin[m] = 1;
        cmax[m] = ex[m];
    }

#ifdef Vertex
    int sym_shift = -ORDN / 2 + 2;
#else
    int sym_shift = -ORDN / 2 + 1;
#endif

    if (symmetry == OCTANT && fabs(X_at_1b(X, 1)) < dX) cmin[0] = sym_shift;
    if (symmetry == OCTANT && fabs(X_at_1b(Y, 1)) < dY) cmin[1] = sym_shift;
    if (symmetry != NO_SYMM && fabs(X_at_1b(Z, 1)) < dZ) cmin[2] = sym_shift;

    for (int m = 0; m < 3; ++m) {
        if (cxB[m] < cmin[m]) { cxB[m] = cmin[m]; cxT[m] = cxB[m] + ORDN - 1; }
        if (cxT[m] > cmax[m]) { cxT[m] = cmax[m]; cxB[m] = cxT[m] + 1 - ORDN; }
    }

#ifdef Vertex
    int ref_shift = 2;
#else
    int ref_shift = 1;
#endif

    double cx[3];
    cx[0] = (cxB[0] > 0) ? (x1 - X_at_1b(X, cxB[0])) / dX : (x1 + X_at_1b(X, ref_shift - cxB[0])) / dX;
    cx[1] = (cxB[1] > 0) ? (y1 - X_at_1b(Y, cxB[1])) / dY : (y1 + X_at_1b(Y, ref_shift - cxB[1])) / dY;
    cx[2] = (cxB[2] > 0) ? (z1 - X_at_1b(Z, cxB[2])) / dZ : (z1 + X_at_1b(Z, ref_shift - cxB[2])) / dZ;

    double wx[MAX_ORDN], wy[MAX_ORDN], wz[MAX_ORDN];
    compute_lagrange_weights(cx[0], ORDN, wx);
    compute_lagrange_weights(cx[1], ORDN, wy);
    compute_lagrange_weights(cx[2], ORDN, wz);

    double val = 0.0;
    for (int k = 0; k < ORDN; ++k) {
        int gz = cxB[2] + k;
        for (int j = 0; j < ORDN; ++j) {
            int gy = cxB[1] + j;
            for (int i = 0; i < ORDN; ++i) {
                int gx = cxB[0] + i;

                double f_val = fetch_f_sym_pad(f, ex, pex, gx, gy, gz, SoA);

                if (isnan(f_val)) {
#if GPU_DEBUG_PRINT
                    printf("global_interp error at global idx: %d %d %d\n", gx, gy, gz);
#endif
                    f_int[0] = NAN;
                    gpu_stop();
                    return;
                }

                val += wx[i] * wy[j] * wz[k] * f_val;
            }
        }
    }

    f_int[0] = val;
}

// backward-compatible wrapper: legacy callers (without pex)
__device__ __forceinline__ void global_interp_device(
    const int* ex,
    const double* X, const double* Y, const double* Z,
    const double* f, double* f_int,
    double x1, double y1, double z1,
    int ORDN, const double* SoA, int symmetry
) {
    int pex_local[3] = {ex[0], ex[1], ex[2]};
    global_interp_device(ex, pex_local, X, Y, Z, f, f_int, x1, y1, z1, ORDN, SoA, symmetry);
}

__device__ __forceinline__ double d_symmetry_bd_1b(
    int ord, const int extc[3], const double* func,
    int i1b, int j1b, int k1b, const double SoA[3]
) {
    if (i1b < -ord + 1 || i1b > extc[0]) return 0.0;
    if (j1b < -ord + 1 || j1b > extc[1]) return 0.0;
    if (k1b < -ord + 1 || k1b > extc[2]) return 0.0;

    int ii = i1b, jj = j1b, kk = k1b;
    double factor = 1.0;

#ifdef Vertex
    // Vertex: symmetry plane is ON a grid point (index 1). Ghost 0 reflects to 2, -1 to 3.
    if (ii <= 0) { ii = 2 - ii; factor *= SoA[0]; }
    if (jj <= 0) { jj = 2 - jj; factor *= SoA[1]; }
    if (kk <= 0) { kk = 2 - kk; factor *= SoA[2]; }
#else
    // Cell: symmetry plane is BETWEEN grid points (0 and 1). Ghost 0 reflects to 1, -1 to 2.
    if (ii <= 0) { ii = 1 - ii; factor *= SoA[0]; }
    if (jj <= 0) { jj = 1 - jj; factor *= SoA[1]; }
    if (kk <= 0) { kk = 1 - kk; factor *= SoA[2]; }
#endif

    if (ii < 1 || ii > extc[0]) return 0.0;
    if (jj < 1 || jj > extc[1]) return 0.0;
    if (kk < 1 || kk > extc[2]) return 0.0;

    return f_at_1b(func, extc, ii, jj, kk) * factor;
}

__device__ __forceinline__ double d_symmetry_bd_1b_pad(
    int ord, const int extc[3], const int pextc[3], const double* func,
    int i1b, int j1b, int k1b, const double SoA[3]
) {
    if (i1b < -ord + 1 || i1b > extc[0]) return 0.0;
    if (j1b < -ord + 1 || j1b > extc[1]) return 0.0;
    if (k1b < -ord + 1 || k1b > extc[2]) return 0.0;

    int ii = i1b, jj = j1b, kk = k1b;
    double factor = 1.0;

#ifdef Vertex
    // Vertex: symmetry plane is ON a grid point (index 1). Ghost 0 reflects to 2, -1 to 3.
    if (ii <= 0) { ii = 2 - ii; factor *= SoA[0]; }
    if (jj <= 0) { jj = 2 - jj; factor *= SoA[1]; }
    if (kk <= 0) { kk = 2 - kk; factor *= SoA[2]; }
#else
    // Cell: symmetry plane is BETWEEN grid points (0 and 1). Ghost 0 reflects to 1, -1 to 2.
    if (ii <= 0) { ii = 1 - ii; factor *= SoA[0]; }
    if (jj <= 0) { jj = 1 - jj; factor *= SoA[1]; }
    if (kk <= 0) { kk = 1 - kk; factor *= SoA[2]; }
#endif

    if (ii < 1 || ii > extc[0]) return 0.0;
    if (jj < 1 || jj > extc[1]) return 0.0;
    if (kk < 1 || kk > extc[2]) return 0.0;

    return f_at_1b_pad(func, extc, pextc, ii, jj, kk) * factor;
}

// global_interp_ss_device — shell patch version with sst parameter
// sst modifies SoA permutation and stencil bounds (cmin/cmax)
__device__ __forceinline__ void global_interp_ss_device(
    const int* ex, const int* pex,
    const double* X, const double* Y, const double* Z,
    const double* f, double* f_int,
    double x1, double y1, double z1,
    int ORDN, const double* SoA, int symmetry, int sst
) {
    if (ORDN <= 0 || ORDN > MAX_ORDN) { f_int[0] = NAN; gpu_stop(); return; }

    double dX = X_at_1b(X, 2) - X_at_1b(X, 1);
    double dY = X_at_1b(Y, 2) - X_at_1b(Y, 1);
    double dZ = X_at_1b(Z, 2) - X_at_1b(Z, 1);

    // Fortran's default-REAL 0.4 is rounded to float before double promotion.
    constexpr double fortran_shift = static_cast<double>(0.4f);
    int cxI[3];
    cxI[0] = (int)((x1 - X_at_1b(X, 1)) / dX + fortran_shift) + 1;
    cxI[1] = (int)((y1 - X_at_1b(Y, 1)) / dY + fortran_shift) + 1;
    cxI[2] = (int)((z1 - X_at_1b(Z, 1)) / dZ + fortran_shift) + 1;

    int cxB[3], cxT[3], cmin[3], cmax[3];
    for (int m = 0; m < 3; ++m) {
        cxB[m] = cxI[m] - ORDN / 2 + 1;
        cxT[m] = cxB[m] + ORDN - 1;
        cmin[m] = 1;
        cmax[m] = ex[m];
    }

    // [修复点 1]：必须根据 sst 生成 SoAh
    double SoAh[3];
    shell_soa3(sst, SoA, SoAh);

    if (sst == -1) {
        cmin[0] = -ORDN + 1; cmin[1] = -ORDN + 1; cmin[2] = -ORDN + 1;
        cmax[0] = ex[0]; cmax[1] = ex[1]; cmax[2] = ex[2];
    } else if (sst == 0 || sst == 1 || sst == 2 || sst == 3 || sst == 4 || sst == 5) {
        cmin[0] = -ORDN + 1; cmin[1] = -ORDN + 1;
        cmax[0] = ex[0] + ORDN; cmax[1] = ex[1] + ORDN;
    }

    for (int m = 0; m < 3; ++m) {
        if (cxB[m] < cmin[m]) { cxB[m] = cmin[m]; cxT[m] = cxB[m] + ORDN - 1; }
        if (cxT[m] > cmax[m]) { cxT[m] = cmax[m]; cxB[m] = cxT[m] + 1 - ORDN; }
    }

    // [修复点 2]：忠于 Fortran，直接代数平移计算 cx，不能带 if 条件
    double cx[3];
    cx[0] = (x1 - X_at_1b(X, 1)) / dX - cxB[0] + 1.0;
    cx[1] = (y1 - X_at_1b(Y, 1)) / dY - cxB[1] + 1.0;
    cx[2] = (z1 - X_at_1b(Z, 1)) / dZ - cxB[2] + 1.0;

    double wx[MAX_ORDN], wy[MAX_ORDN], wz[MAX_ORDN];
    compute_lagrange_weights(cx[0], ORDN, wx);
    compute_lagrange_weights(cx[1], ORDN, wy);
    compute_lagrange_weights(cx[2], ORDN, wz);

    double val = 0.0;
    for (int k = 0; k < ORDN; ++k) {
        int gz = cxB[2] + k; // 1-based index
        for (int j = 0; j < ORDN; ++j) {
            int gy = cxB[1] + j; // 1-based index
            for (int i = 0; i < ORDN; ++i) {
                int gx = cxB[0] + i; // 1-based index
                
                // [修复点 3]：必须使用支持双侧边界的 fetch_f_tbd_0b_pad 函数，并传入转化为 0-based 的坐标
                double f_val = fetch_f_tbd_0b_pad(f, ex, pex, gx - 1, gy - 1, gz - 1, SoAh);
                
                if (isnan(f_val)) {
                    f_int[0] = NAN;
                    gpu_stop();
                    return;
                }
                val += wx[i] * wy[j] * wz[k] * f_val;
            }
        }
    }

    f_int[0] = val;
}

// 2D shell interpolation device function (for NullShellPatch Interp_Points_2D)
// Interpolates at fixed z-index indZ (1-based), only varies in x,y
__device__ __forceinline__ void global_interp_ss_2d_device(
    const int* ex, const int* pex,
    const double* X, const double* Y,
    const double* f, double* f_int,
    double x1, double y1,
    int ORDN, const double* SoA, int symmetry, int sst, int indZ
) {
    if (ORDN <= 0 || ORDN > MAX_ORDN) { f_int[0] = NAN; gpu_stop(); return; }

    double dX = X_at_1b(X, 2) - X_at_1b(X, 1);
    double dY = X_at_1b(Y, 2) - X_at_1b(Y, 1);

    int cxI[2];
    cxI[0] = (int)((x1 - X_at_1b(X, 1)) / dX + 0.4) + 1;
    cxI[1] = (int)((y1 - X_at_1b(Y, 1)) / dY + 0.4) + 1;

    int cxB[2], cxT[2], cmin[2], cmax[2];
    for (int m = 0; m < 2; ++m) {
        cxB[m] = cxI[m] - ORDN / 2 + 1;
        cxT[m] = cxB[m] + ORDN - 1;
        cmin[m] = 1;
        cmax[m] = ex[m];
    }

    // [修复点] 同样的 SoAh 变换
    double SoAh[3];
    shell_soa3(sst, SoA, SoAh);

    if (sst == -1) {
        cmin[0] = -ORDN + 1; cmin[1] = -ORDN + 1;
    } else if (sst == 0 || sst == 1 || sst == 2 || sst == 3 || sst == 4 || sst == 5) {
        cmin[0] = -ORDN + 1; cmin[1] = -ORDN + 1;
        cmax[0] = ex[0] + ORDN; cmax[1] = ex[1] + ORDN;
    }

    for (int m = 0; m < 2; ++m) {
        if (cxB[m] < cmin[m]) { cxB[m] = cmin[m]; cxT[m] = cxB[m] + ORDN - 1; }
        if (cxT[m] > cmax[m]) { cxT[m] = cmax[m]; cxB[m] = cxT[m] + 1 - ORDN; }
    }

    // [修复点] 拉平坐标计算
    double cx[2];
    cx[0] = (x1 - X_at_1b(X, 1)) / dX - cxB[0] + 1.0;
    cx[1] = (y1 - X_at_1b(Y, 1)) / dY - cxB[1] + 1.0;

    double wx[MAX_ORDN], wy[MAX_ORDN];
    compute_lagrange_weights(cx[0], ORDN, wx);
    compute_lagrange_weights(cx[1], ORDN, wy);

    double val = 0.0;
    for (int j = 0; j < ORDN; ++j) {
        int gy = cxB[1] + j; // 1-based
        for (int i = 0; i < ORDN; ++i) {
            int gx = cxB[0] + i; // 1-based
            
            // [修复点] 必须用 tbd 的 0-based 获取函数
            double f_val = fetch_f_tbd_0b_pad(f, ex, pex, gx - 1, gy - 1, indZ - 1, SoAh);
            
            if (isnan(f_val)) {
                f_int[0] = NAN;
                gpu_stop();
                return;
            }
            val += wx[i] * wy[j] * f_val;
        }
    }

    f_int[0] = val;
}

__device__ __forceinline__ double warpReduceSum(double val) {
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {
        val += __shfl_down_sync(0xffffffff, val, offset);
    }
    return val;
}

__device__ __forceinline__ int warpReduceSumInt(int val) {
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {
        val += __shfl_down_sync(0xffffffff, val, offset);
    }
    return val;
}

#endif // FMISC_GPU_CUH
