#ifndef ADVECTION_DISSIPATION_SMEM_GPU_CUH
#define ADVECTION_DISSIPATION_SMEM_GPU_CUH

// One-axis lopsided-advection + KO-dissipation weight assembly helpers.
//
// Both BSSN and BSSN-EM advection_dissipation kernels follow the same
// pattern: build per-axis weight arrays Wx/Wy/Wz of length 2*MAX_DIST+1,
// then sweep a fixed list of (field, field_rhs) pairs and accumulate
//
//     f_rhs[idx] += Sum_d W?[WCTR+d] * f[idx + d * stride]
//
// The lopsided stencil shape and the KO stencil shape depend only on GW
// and on the local axis context (v, i, imin, imax, h).  This header
// provides the reusable "assemble one axis" device helpers; the kernel
// itself owns the field list and the final __ldg sweep.
//
// Author: extracted from bssn_advection_dissipation_kernel for reuse by
// BSSN-EM (and later BSSN-EScalar / Z4c).  Numerical formulas unchanged.

#include <cmath>

namespace adv_diss {

// Longest one-sided stencil reach across lopsided + KO for this GW.
//   GW=2: reach 2  GW=3: reach 3  GW=4: reach 5  GW=5: reach 5
// Mirrors `MAX_DIST = (GW <= 3) ? GW : 5` in the BSSN kernel.
template<int GW>
__device__ inline constexpr int max_dist() {
    return (GW <= 3) ? GW : 5;
}

template<int GW>
__device__ inline constexpr int wsize() { return 2 * max_dist<GW>() + 1; }

template<int GW>
__device__ inline constexpr int wctr() { return max_dist<GW>(); }

// KO normalization 2^(2*GW) and sign (-1)^(GW-1).
template<int GW>
__device__ inline constexpr double ko_cof()  { return double(1 << (2 * GW)); }
template<int GW>
__device__ inline constexpr double ko_sign() {
    return ((GW - 1) % 2 == 0) ? 1.0 : -1.0;
}

// Lopsided/upwind advection stencil weights for one axis.
// Accumulates into W in place.
//
//   W:    weight buffer, length wsize<GW>(), centre = wctr<GW>()
//   v:    advection velocity along this axis (beta^i)
//   h:    grid spacing along this axis
//   i:    current logical index along this axis (0-based)
//   imin: lower symmetry-extended bound for this axis
//   imax: upper interior bound for this axis (= ex - 1)
//
// Formulas verbatim from bssn_advection_dissipation_kernel.
template<int GW>
__device__ __forceinline__ void assemble_advect_weights_1d(
    double* W, double v, double h,
    int i, int imin, int imax
) {
    constexpr int WCTR = wctr<GW>();

    if constexpr (GW == 2) {
        const double d2dh = 0.5 / h;
        const double v_dh = v * d2dh;

        if (v > 0.0 && i + 2 <= imax) {
            W[WCTR]   += v_dh * (-3.0);
            W[WCTR+1] += v_dh *  4.0;
            W[WCTR+2] += v_dh * (-1.0);
        } else if (v > 0.0 && i + 1 <= imax) {
            W[WCTR]   += v_dh * (-1.0);
            W[WCTR+1] += v_dh *  1.0;
        }
        if (v < 0.0 && i - 2 >= imin) {
            double nv = -v_dh;
            W[WCTR]   += nv * (-3.0);
            W[WCTR-1] += nv *  4.0;
            W[WCTR-2] += nv * (-1.0);
        } else if (v < 0.0 && i - 1 >= imin) {
            W[WCTR]   += (-v_dh) * (-1.0);
            W[WCTR-1] += (-v_dh) *  1.0;
        }
    } else if constexpr (GW == 3) {
        const double d12dh = 1.0 / (12.0 * h);
        const double v_dh  = v * d12dh;

        int sx = 0;
        if (v > 0.0)      { if (i + 3 <= imax) sx =  3; else if (i + 2 <= imax) sx =  2; else if (i + 1 <= imax) sx =  1; }
        else if (v < 0.0) { if (i - 3 >= imin) sx = -3; else if (i - 2 >= imin) sx = -2; else if (i - 1 >= imin) sx = -1; }

        if (sx == 3 || sx == -1) {
            W[WCTR-1] += v_dh *  -3.0;
            W[WCTR  ] += v_dh * -10.0;
            W[WCTR+1] += v_dh *  18.0;
            W[WCTR+2] += v_dh *  -6.0;
            W[WCTR+3] += v_dh *   1.0;
        } else if (sx == 2 || sx == -2) {
            W[WCTR-2] += v_dh *   1.0;
            W[WCTR-1] += v_dh *  -8.0;
            W[WCTR+1] += v_dh *   8.0;
            W[WCTR+2] += v_dh *  -1.0;
        } else if (sx == 1 || sx == -3) {
            W[WCTR-3] += v_dh *  -1.0;
            W[WCTR-2] += v_dh *   6.0;
            W[WCTR-1] += v_dh * -18.0;
            W[WCTR  ] += v_dh *  10.0;
            W[WCTR+1] += v_dh *   3.0;
        }
    } else if constexpr (GW == 4) {
        const double d60dh = 1.0 / 60.0 / h;
        const double d12dh = 1.0 / 12.0 / h;
        const double d2dh  = 0.5 / h;

        if (v >= 0.0 && i + 4 <= imax && i - 2 >= imin) {
            double c = v * d60dh;
            W[WCTR-2] += c *   2.0; W[WCTR-1] += c * -24.0; W[WCTR  ] += c * -35.0;
            W[WCTR+1] += c *  80.0; W[WCTR+2] += c * -30.0; W[WCTR+3] += c *   8.0; W[WCTR+4] += c *  -1.0;
        } else if (v >= 0.0 && i + 5 <= imax && i - 1 >= imin) {
            double c = v * d60dh;
            W[WCTR-1] += c * -10.0; W[WCTR  ] += c * -77.0; W[WCTR+1] += c * 150.0;
            W[WCTR+2] += c *-100.0; W[WCTR+3] += c *  50.0; W[WCTR+4] += c * -15.0; W[WCTR+5] += c *   2.0;
        } else if (v <= 0.0 && i - 4 >= imin && i + 2 <= imax) {
            double c = -v * d60dh;
            W[WCTR+2] += c *   2.0; W[WCTR+1] += c * -24.0; W[WCTR  ] += c * -35.0;
            W[WCTR-1] += c *  80.0; W[WCTR-2] += c * -30.0; W[WCTR-3] += c *   8.0; W[WCTR-4] += c *  -1.0;
        } else if (v <= 0.0 && i - 5 >= imin && i + 1 <= imax) {
            double c = -v * d60dh;
            W[WCTR+1] += c * -10.0; W[WCTR  ] += c * -77.0; W[WCTR-1] += c * 150.0;
            W[WCTR-2] += c *-100.0; W[WCTR-3] += c *  50.0; W[WCTR-4] += c * -15.0; W[WCTR-5] += c *   2.0;
        } else if (i + 3 <= imax && i - 3 >= imin) {
            double c = v * d60dh;
            W[WCTR-3] += c *  -1.0; W[WCTR-2] += c *   9.0; W[WCTR-1] += c * -45.0;
            W[WCTR+1] += c *  45.0; W[WCTR+2] += c *  -9.0; W[WCTR+3] += c *   1.0;
        } else if (i + 2 <= imax && i - 2 >= imin) {
            double c = v * d12dh;
            W[WCTR-2] += c *   1.0; W[WCTR-1] += c *  -8.0;
            W[WCTR+1] += c *   8.0; W[WCTR+2] += c *  -1.0;
        } else if (i + 1 <= imax && i - 1 >= imin) {
            double c = v * d2dh;
            W[WCTR-1] += c *  -1.0; W[WCTR+1] += c *   1.0;
        }
    } else if constexpr (GW == 5) {
        const double d840dh = 1.0 / 840.0 / h;
        const double d60dh  = 1.0 /  60.0 / h;
        const double d12dh  = 1.0 /  12.0 / h;
        const double d2dh   = 0.5 / h;

        if (v >= 0.0 && i + 5 <= imax && i - 3 >= imin) {
            double c = v * d840dh;
            W[WCTR-3] += c *   -5.0; W[WCTR-2] += c *   60.0; W[WCTR-1] += c * -420.0; W[WCTR  ] += c * -378.0;
            W[WCTR+1] += c * 1050.0; W[WCTR+2] += c * -420.0; W[WCTR+3] += c *  140.0; W[WCTR+4] += c *  -30.0; W[WCTR+5] += c *    3.0;
        } else if (v <= 0.0 && i - 5 >= imin && i + 3 <= imax) {
            double c = -v * d840dh;
            W[WCTR+3] += c *   -5.0; W[WCTR+2] += c *   60.0; W[WCTR+1] += c * -420.0; W[WCTR  ] += c * -378.0;
            W[WCTR-1] += c * 1050.0; W[WCTR-2] += c * -420.0; W[WCTR-3] += c *  140.0; W[WCTR-4] += c *  -30.0; W[WCTR-5] += c *    3.0;
        } else if (i + 4 <= imax && i - 4 >= imin) {
            double c = v * d840dh;
            W[WCTR-4] += c *   3.0; W[WCTR-3] += c *  -32.0; W[WCTR-2] += c *  168.0; W[WCTR-1] += c * -672.0;
            W[WCTR+1] += c * 672.0; W[WCTR+2] += c * -168.0; W[WCTR+3] += c *   32.0; W[WCTR+4] += c *   -3.0;
        } else if (i + 3 <= imax && i - 3 >= imin) {
            double c = v * d60dh;
            W[WCTR-3] += c *  -1.0; W[WCTR-2] += c *   9.0; W[WCTR-1] += c * -45.0;
            W[WCTR+1] += c *  45.0; W[WCTR+2] += c *  -9.0; W[WCTR+3] += c *   1.0;
        } else if (i + 2 <= imax && i - 2 >= imin) {
            double c = v * d12dh;
            W[WCTR-2] += c *   1.0; W[WCTR-1] += c *  -8.0;
            W[WCTR+1] += c *   8.0; W[WCTR+2] += c *  -1.0;
        } else if (i + 1 <= imax && i - 1 >= imin) {
            double c = v * d2dh;
            W[WCTR-1] += c *  -1.0; W[WCTR+1] += c *   1.0;
        }
    }
}

// Kreiss-Oliger dissipation stencil weights for one axis.
// Caller passes edh = ko_sign<GW>() * eps / (ko_cof<GW>() * h).
// Accumulates into W in place.
//
// Coefficients: (D+ D-)^GW, i.e. binomial(2*GW, k) * (-1)^k for k = 0..2*GW.
template<int GW>
__device__ __forceinline__ void assemble_ko_weights_1d(double* W, double edh) {
    constexpr int WCTR = wctr<GW>();
    if constexpr (GW == 2) {
        constexpr double c[5] = {1.0, -4.0, 6.0, -4.0, 1.0};
        #pragma unroll
        for (int d = 0; d < 5; ++d) W[WCTR-2+d] += edh * c[d];
    } else if constexpr (GW == 3) {
        constexpr double c[7] = {1.0, -6.0, 15.0, -20.0, 15.0, -6.0, 1.0};
        #pragma unroll
        for (int d = 0; d < 7; ++d) W[WCTR-3+d] += edh * c[d];
    } else if constexpr (GW == 4) {
        constexpr double c[9] = {1.0, -8.0, 28.0, -56.0, 70.0, -56.0, 28.0, -8.0, 1.0};
        #pragma unroll
        for (int d = 0; d < 9; ++d) W[WCTR-4+d] += edh * c[d];
    } else if constexpr (GW == 5) {
        constexpr double c[11] = {1.0, -10.0, 45.0, -120.0, 210.0, -252.0, 210.0, -120.0, 45.0, -10.0, 1.0};
        #pragma unroll
        for (int d = 0; d < 11; ++d) W[WCTR-5+d] += edh * c[d];
    }
}

// Quick boundary check for safe KO application (stencil fully inside [imin,imax]).
template<int GW>
__device__ __forceinline__ bool ko_fits(int i, int imin, int imax) {
    return (i - GW >= imin) && (i + GW <= imax);
}

} // namespace adv_diss

#endif // ADVECTION_DISSIPATION_SMEM_GPU_CUH
