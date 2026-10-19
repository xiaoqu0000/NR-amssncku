#ifndef FMISC_GPU_CUH
#define FMISC_GPU_CUH

#include <cuda_runtime.h>
#include <math.h>

#include <iostream>
#include "gpu_runtime_coordinator.h"
#ifndef ghost_width
#include "macrodef.fh"
#endif

#define MAX_ORDN 10
#ifndef GPU_DEBUG_PRINT
#define GPU_DEBUG_PRINT 0
#endif

template <typename Real>
__device__ __forceinline__ Real d_symmetry_bd_0b_t(
    int ord, 
    int extc0, int extc1, int extc2,
    const Real* __restrict__ func,
    int i, int j, int k, // 0-based coordinates
    double SoA0, double SoA1, double SoA2
) {
    if (i < -ord || i >= extc0) return Real(0.0);
    if (j < -ord || j >= extc1) return Real(0.0);
    if (k < -ord || k >= extc2) return Real(0.0);

    Real factor = Real(1.0);

#ifdef Vertex
    if (i < 0) { i = -i;     factor *= static_cast<Real>(SoA0); }
    if (j < 0) { j = -j;     factor *= static_cast<Real>(SoA1); }
    if (k < 0) { k = -k;     factor *= static_cast<Real>(SoA2); }
#else
    if (i < 0) { i = -i - 1; factor *= static_cast<Real>(SoA0); }
    if (j < 0) { j = -j - 1; factor *= static_cast<Real>(SoA1); }
    if (k < 0) { k = -k - 1; factor *= static_cast<Real>(SoA2); }
#endif

    if (i >= extc0 || j >= extc1 || k >= extc2) return Real(0.0);

    return func[(k * extc1 + j) * extc0 + i] * factor;
}

__device__ __forceinline__ double d_symmetry_bd_0b(
    int ord,
    int extc0, int extc1, int extc2,
    const double* __restrict__ func,
    int i, int j, int k, // 0-based coordinates
    double SoA0, double SoA1, double SoA2
) {
    return d_symmetry_bd_0b_t(
        ord, extc0, extc1, extc2, func, i, j, k, SoA0, SoA1, SoA2
    );
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

__device__ __forceinline__ void polint(const double* xa, const double* ya, double x, double& y, double& dy, int ordn) {
	double c[MAX_ORDN], d[MAX_ORDN], den[MAX_ORDN], ho[MAX_ORDN];
	int ns = 1;
	double dif = fabs(x - xa[0]);
	for (int m = 0; m < ordn; ++m) {
		c[m] = ya[m];
		d[m] = ya[m];
		ho[m] = xa[m] - x;
		double dift = fabs(x - xa[m]);
		if (dift < dif) { ns = m + 1; dif = dift; }
	}
	y = ya[ns - 1];
	ns = ns - 1;
	for (int m = 1; m < ordn; ++m) {
		for (int i = 0; i < ordn - m; ++i) {
			den[i] = ho[i] - ho[i + m];
			if (den[i] == 0.0) {
#if GPU_DEBUG_PRINT
                printf("failure in polint for point %f\n", x);
                printf("with input points: ");
                for (int t = 0; t < ordn; ++t) printf("%f ", xa[t]);
                printf("\n");
#endif
				y = NAN; dy = NAN; gpu_stop(); return;
			}
			den[i] = (c[i + 1] - d[i]) / den[i];
			d[i] = ho[i + m] * den[i];
			c[i] = ho[i] * den[i];
		}
		if (2 * ns < (ordn - m)) {
			dy = c[ns];
		} else {
			dy = d[ns - 1];
			ns = ns - 1;
		}
		y = y + dy;
	}
}

__device__ __forceinline__ void d_polin3_1b(
	const double* x1a, const double* x2a, const double* x3a,
	const double* ya, double x1, double x2, double x3,
	double& y, double& dy, int ordn
) {
	double yatmp[MAX_ORDN * MAX_ORDN];
	double ymtmp[MAX_ORDN];
	double yntmp[MAX_ORDN];
	double yqtmp[MAX_ORDN];

	for (int i = 0; i < ordn; ++i) {
		for (int j = 0; j < ordn; ++j) {
			for (int k = 0; k < ordn; ++k) {
				yqtmp[k] = ya[(k * ordn + j) * ordn + i];
			}
			polint(x3a, yqtmp, x3, yatmp[j * ordn + i], dy, ordn);
		}
		for (int j = 0; j < ordn; ++j) yntmp[j] = yatmp[j * ordn + i];
		polint(x2a, yntmp, x2, ymtmp[i], dy, ordn);
	}
	polint(x1a, ymtmp, x1, y, dy, ordn);
}

__device__ __forceinline__ bool d_decide3d(
	const int ex[3], const double* f, const double* fpi,
	const int cxB[3], const int cxT[3], const double SoA[3],
	double* ya, int ordn, int Symmetry
) {
	(void)fpi;
	(void)Symmetry;
	bool gont = false;
	int fmin1[3], fmin2[3], fmax1[3], fmax2[3];

	for (int m = 0; m < 3; ++m) {
		if (!(abs(cxB[m]) >= 0)) gont = true;
		if (!(abs(cxT[m]) >= 0)) gont = true;
		fmin1[m] = max(1, cxB[m]);
		fmax1[m] = cxT[m];
		fmin2[m] = cxB[m];
		fmax2[m] = min(0, cxT[m]);
		if ((fmin1[m] <= fmax1[m]) && (fmin1[m] < 1 || fmax1[m] > ex[m])) gont = true;
#ifdef Vertex
		if ((fmin2[m] <= fmax2[m]) && (2 - fmax2[m] < 1 || 2 - fmin2[m] > ex[m])) gont = true;
#else
		if ((fmin2[m] <= fmax2[m]) && (1 - fmax2[m] < 1 || 1 - fmin2[m] > ex[m])) gont = true;
#endif
	}
	if (gont) {
#if GPU_DEBUG_PRINT
        printf("error in decide3d\n");
        printf("cxB: %d %d %d, cxT: %d %d %d, ex: %d %d %d\n",
               cxB[0], cxB[1], cxB[2], cxT[0], cxT[1], cxT[2], ex[0], ex[1], ex[2]);
        printf("fmin1: %d %d %d, fmax1: %d %d %d\n",
               fmin1[0], fmin1[1], fmin1[2], fmax1[0], fmax1[1], fmax1[2]);
        printf("fmin2: %d %d %d, fmax2: %d %d %d\n",
               fmin2[0], fmin2[1], fmin2[2], fmax2[0], fmax2[1], fmax2[2]);
#endif
		return true;
	}

	auto idx = [&](int i, int j, int k) {
		return ((k - cxB[2]) * ordn + (j - cxB[1])) * ordn + (i - cxB[0]);
	};

	for (int k = fmin1[2]; k <= fmax1[2]; ++k) {
		for (int j = fmin1[1]; j <= fmax1[1]; ++j) {
			for (int i = fmin1[0]; i <= fmax1[0]; ++i)
				ya[idx(i, j, k)] = f_at_1b(f, ex, i, j, k);
			for (int i = fmin2[0]; i <= fmax2[0]; ++i)
#ifdef Vertex
				ya[idx(i, j, k)] = f_at_1b(f, ex, 2 - i, j, k) * SoA[0];
#else
				ya[idx(i, j, k)] = f_at_1b(f, ex, 1 - i, j, k) * SoA[0];
#endif
		}
		for (int j = fmin2[1]; j <= fmax2[1]; ++j) {
			for (int i = fmin1[0]; i <= fmax1[0]; ++i)
#ifdef Vertex
				ya[idx(i, j, k)] = f_at_1b(f, ex, i, 2 - j, k) * SoA[1];
#else
				ya[idx(i, j, k)] = f_at_1b(f, ex, i, 1 - j, k) * SoA[1];
#endif
			for (int i = fmin2[0]; i <= fmax2[0]; ++i)
#ifdef Vertex
				ya[idx(i, j, k)] = f_at_1b(f, ex, 2 - i, 2 - j, k) * SoA[0] * SoA[1];
#else
				ya[idx(i, j, k)] = f_at_1b(f, ex, 1 - i, 1 - j, k) * SoA[0] * SoA[1];
#endif
		}
	}

	for (int k = fmin2[2]; k <= fmax2[2]; ++k) {
		for (int j = fmin1[1]; j <= fmax1[1]; ++j) {
			for (int i = fmin1[0]; i <= fmax1[0]; ++i)
#ifdef Vertex
				ya[idx(i, j, k)] = f_at_1b(f, ex, i, j, 2 - k) * SoA[2];
#else
				ya[idx(i, j, k)] = f_at_1b(f, ex, i, j, 1 - k) * SoA[2];
#endif
			for (int i = fmin2[0]; i <= fmax2[0]; ++i)
#ifdef Vertex
				ya[idx(i, j, k)] = f_at_1b(f, ex, 2 - i, j, 2 - k) * SoA[0] * SoA[2];
#else
				ya[idx(i, j, k)] = f_at_1b(f, ex, 1 - i, j, 1 - k) * SoA[0] * SoA[2];
#endif
		}
		for (int j = fmin2[1]; j <= fmax2[1]; ++j) {
			for (int i = fmin1[0]; i <= fmax1[0]; ++i)
#ifdef Vertex
				ya[idx(i, j, k)] = f_at_1b(f, ex, i, 2 - j, 2 - k) * SoA[1] * SoA[2];
#else
				ya[idx(i, j, k)] = f_at_1b(f, ex, i, 1 - j, 1 - k) * SoA[1] * SoA[2];
#endif
			for (int i = fmin2[0]; i <= fmax2[0]; ++i)
#ifdef Vertex
				ya[idx(i, j, k)] = f_at_1b(f, ex, 2 - i, 2 - j, 2 - k) * SoA[0] * SoA[1] * SoA[2];
#else
				ya[idx(i, j, k)] = f_at_1b(f, ex, 1 - i, 1 - j, 1 - k) * SoA[0] * SoA[1] * SoA[2];
#endif
		}
	}
	return false;
}

// __device__ __forceinline__ void global_interp_device(
// 	const int* ex, const double* X, const double* Y, const double* Z,
// 	const double* f, double* f_int,
// 	double x1, double y1, double z1,
// 	int ORDN, const double* SoA, int symmetry
// ) {
// 	if (ORDN > MAX_ORDN) { f_int[0] = NAN; gpu_stop(); return; }

// 	const int NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2;
// 	int imin = 1, jmin = 1, kmin = 1;

// 	double dX = X_at_1b(X, imin + 1) - X_at_1b(X, imin);
// 	double dY = X_at_1b(Y, jmin + 1) - X_at_1b(Y, jmin);
// 	double dZ = X_at_1b(Z, kmin + 1) - X_at_1b(Z, kmin);

// 	double x1a[MAX_ORDN];
// 	for (int j = 0; j < ORDN; ++j) x1a[j] = (double)j;

// 	int cxI[3];
// 	cxI[0] = (int)((x1 - X_at_1b(X, 1)) / dX + 0.4) + 1;
// 	cxI[1] = (int)((y1 - X_at_1b(Y, 1)) / dY + 0.4) + 1;
// 	cxI[2] = (int)((z1 - X_at_1b(Z, 1)) / dZ + 0.4) + 1;

// 	int cxB[3], cxT[3], cmin[3], cmax[3];
// 	for (int m = 0; m < 3; ++m) {
// 		cxB[m] = cxI[m] - ORDN / 2 + 1;
// 		cxT[m] = cxB[m] + ORDN - 1;
// 		cmin[m] = 1;
// 		cmax[m] = ex[m];
// 	}
// 	if (symmetry == OCTANT && fabs(X_at_1b(X, 1)) < dX) cmin[0] = -ORDN / 2 + 1;
// 	if (symmetry == OCTANT && fabs(X_at_1b(Y, 1)) < dY) cmin[1] = -ORDN / 2 + 1;
// 	if (symmetry != NO_SYMM && fabs(X_at_1b(Z, 1)) < dZ) cmin[2] = -ORDN / 2 + 1;

// 	for (int m = 0; m < 3; ++m) {
// 		if (cxB[m] < cmin[m]) { cxB[m] = cmin[m]; cxT[m] = cxB[m] + ORDN - 1; }
// 		if (cxT[m] > cmax[m]) { cxT[m] = cmax[m]; cxB[m] = cxT[m] + 1 - ORDN; }
// 	}

// 	double cx[3];
// 	cx[0] = (cxB[0] > 0) ? (x1 - X_at_1b(X, cxB[0])) / dX : (x1 + X_at_1b(X, 1 - cxB[0])) / dX;
// 	cx[1] = (cxB[1] > 0) ? (y1 - X_at_1b(Y, cxB[1])) / dY : (y1 + X_at_1b(Y, 1 - cxB[1])) / dY;
// 	cx[2] = (cxB[2] > 0) ? (z1 - X_at_1b(Z, cxB[2])) / dZ : (z1 + X_at_1b(Z, 1 - cxB[2])) / dZ;

// 	double ya[MAX_ORDN * MAX_ORDN * MAX_ORDN];
// 	if (d_decide3d(ex, f, f, cxB, cxT, SoA, ya, ORDN, symmetry)) {
// #if GPU_DEBUG_PRINT
//         printf("global_interp position: %f %f %f\n", x1, y1, z1);
//         printf("data range: %f %f %f %f %f %f\n",
//                X_at_1b(X, 1), X_at_1b(X, ex[0]),
//                X_at_1b(Y, 1), X_at_1b(Y, ex[1]),
//                X_at_1b(Z, 1), X_at_1b(Z, ex[2]));
// #endif
// 		f_int[0] = NAN;
// 		gpu_stop();
// 		return;
// 	}

// 	double ddy = 0.0;
// 	d_polin3_1b(x1a, x1a, x1a, ya, cx[0], cx[1], cx[2], f_int[0], ddy, ORDN);
// }

__device__ __forceinline__ double fetch_f_sym(
    const double* f, const int ex[3], 
    int gx, int gy, int gz, 
    const double SoA[3]) 
{
    double sign = 1.0;

#ifdef Vertex
    if (gx < 1) { gx = 2 - gx; sign *= SoA[0]; }
    if (gy < 1) { gy = 2 - gy; sign *= SoA[1]; }
    if (gz < 1) { gz = 2 - gz; sign *= SoA[2]; }
#else
    if (gx < 1) { gx = 1 - gx; sign *= SoA[0]; }
    if (gy < 1) { gy = 1 - gy; sign *= SoA[1]; }
    if (gz < 1) { gz = 1 - gz; sign *= SoA[2]; }
#endif

    if (gx > ex[0] || gy > ex[1] || gz > ex[2]) {
        return NAN; 
    }
    
    return sign * f_at_1b(f, ex, gx, gy, gz);
}

__device__ __forceinline__ void compute_lagrange_weights(double x, int ordn, double* w) {
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

__device__ __forceinline__ double uniform_spacing_1b(const double* axis) {
	return X_at_1b(axis, 2) - X_at_1b(axis, 1);
}

__device__ __forceinline__ int nearest_cell_1b(double coord, const double* axis, double spacing) {
	return (int)((coord - X_at_1b(axis, 1)) / spacing + 0.4) + 1;
}

__device__ __forceinline__ int interp_lower_bound(int symmetry, double axis_origin, double spacing, int axis, int ordn) {
	const int NO_SYMM = 0;
	const int OCTANT = 2;
#ifdef Vertex
	const int symmetry_lower = -ordn / 2 + 2;
#else
	const int symmetry_lower = -ordn / 2 + 1;
#endif
	if (axis < 2 && symmetry == OCTANT && fabs(axis_origin) < spacing) return symmetry_lower;
	if (axis == 2 && symmetry != NO_SYMM && fabs(axis_origin) < spacing) return symmetry_lower;
	return 1;
}

__device__ __forceinline__ void build_interp_window(
	const int* ex,
	const double* X,
	const double* Y,
	const double* Z,
	double x1,
	double y1,
	double z1,
	int ORDN,
	int symmetry,
	int cxB[3],
	int cxT[3]
) {
	const double dX = uniform_spacing_1b(X);
	const double dY = uniform_spacing_1b(Y);
	const double dZ = uniform_spacing_1b(Z);
	const int centers[3] = {
		nearest_cell_1b(x1, X, dX),
		nearest_cell_1b(y1, Y, dY),
		nearest_cell_1b(z1, Z, dZ)
	};
	const int lower_bounds[3] = {
		interp_lower_bound(symmetry, X_at_1b(X, 1), dX, 0, ORDN),
		interp_lower_bound(symmetry, X_at_1b(Y, 1), dY, 1, ORDN),
		interp_lower_bound(symmetry, X_at_1b(Z, 1), dZ, 2, ORDN)
	};

	for (int axis = 0; axis < 3; ++axis) {
		cxB[axis] = centers[axis] - ORDN / 2 + 1;
		cxT[axis] = cxB[axis] + ORDN - 1;
		if (cxB[axis] < lower_bounds[axis]) {
			cxB[axis] = lower_bounds[axis];
			cxT[axis] = cxB[axis] + ORDN - 1;
		}
		if (cxT[axis] > ex[axis]) {
			cxT[axis] = ex[axis];
			cxB[axis] = cxT[axis] + 1 - ORDN;
		}
	}
}

__device__ __forceinline__ void compute_local_coords(
	const double* X,
	const double* Y,
	const double* Z,
	const int cxB[3],
	double x1,
	double y1,
	double z1,
	double cx[3]
) {
	const double dX = uniform_spacing_1b(X);
	const double dY = uniform_spacing_1b(Y);
	const double dZ = uniform_spacing_1b(Z);
#ifdef Vertex
	cx[0] = (cxB[0] > 0) ? (x1 - X_at_1b(X, cxB[0])) / dX : (x1 + X_at_1b(X, 2 - cxB[0])) / dX;
	cx[1] = (cxB[1] > 0) ? (y1 - X_at_1b(Y, cxB[1])) / dY : (y1 + X_at_1b(Y, 2 - cxB[1])) / dY;
	cx[2] = (cxB[2] > 0) ? (z1 - X_at_1b(Z, cxB[2])) / dZ : (z1 + X_at_1b(Z, 2 - cxB[2])) / dZ;
#else
	cx[0] = (cxB[0] > 0) ? (x1 - X_at_1b(X, cxB[0])) / dX : (x1 + X_at_1b(X, 1 - cxB[0])) / dX;
	cx[1] = (cxB[1] > 0) ? (y1 - X_at_1b(Y, cxB[1])) / dY : (y1 + X_at_1b(Y, 1 - cxB[1])) / dY;
	cx[2] = (cxB[2] > 0) ? (z1 - X_at_1b(Z, cxB[2])) / dZ : (z1 + X_at_1b(Z, 1 - cxB[2])) / dZ;
#endif
}

__device__ __forceinline__ void global_interp_device(
    const int* ex, const double* X, const double* Y, const double* Z,
    const double* f, double* f_int,
    double x1, double y1, double z1,
    int ORDN, const double* SoA, int symmetry
) {
    if (ORDN > MAX_ORDN) { f_int[0] = NAN; gpu_stop(); return; }

	int cxB[3], cxT[3];
	build_interp_window(ex, X, Y, Z, x1, y1, z1, ORDN, symmetry, cxB, cxT);

	double cx[3];
	compute_local_coords(X, Y, Z, cxB, x1, y1, z1, cx);

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
                double f_val = fetch_f_sym(f, ex, gx, gy, gz, SoA);

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
	if (ii <= 0) { ii = 2 - ii; factor *= SoA[0]; }
	if (jj <= 0) { jj = 2 - jj; factor *= SoA[1]; }
	if (kk <= 0) { kk = 2 - kk; factor *= SoA[2]; }
#else
	if (ii <= 0) { ii = 1 - ii; factor *= SoA[0]; }
	if (jj <= 0) { jj = 1 - jj; factor *= SoA[1]; }
	if (kk <= 0) { kk = 1 - kk; factor *= SoA[2]; }
#endif

	if (ii < 1 || ii > extc[0]) return 0.0;
	if (jj < 1 || jj > extc[1]) return 0.0;
	if (kk < 1 || kk > extc[2]) return 0.0;

	return f_at_1b(func, extc, ii, jj, kk) * factor;
}

__device__ __forceinline__ double warpReduceSum(double val) {
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {
        val += __shfl_down_sync(0xffffffff, val, offset);
    }
    return val;
}

#endif // FMISC_GPU_CUH
