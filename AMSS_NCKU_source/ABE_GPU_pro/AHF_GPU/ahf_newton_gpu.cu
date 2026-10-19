#include "ahf_newton_gpu.h"
#include <cmath>

namespace {
constexpr int threads = 256;

template<bool Maximum>
__device__ void reduce(double value, int invalid, double* values, int* bad) {
    const int t = threadIdx.x;
    values[t] = value;
    bad[t] = invalid;
    __syncthreads();
    for (int stride = threads / 2; stride; stride /= 2) {
        if (t < stride) {
            if (Maximum)
                values[t] = values[t] > values[t + stride] ? values[t] : values[t + stride];
            else
                values[t] = __dadd_rn(values[t], values[t + stride]);
            bad[t] |= bad[t + stride];
        }
        __syncthreads();
    }
}

__global__ void baseline_kernel(int n, int ng, const int *map,
                                const double *ghost, AhfNewtonState *s) {
    __shared__ double sums[threads];
    __shared__ int bad[threads];
    const int t = threadIdx.x;
    if (t == 0) *s = AhfNewtonState{};
    __syncthreads();
    double sum = 0.0;
    int invalid = 0;
    for (int q = t; q < n; q += threads) {
        const int g = map[q];
        if (g < 0 || g >= ng) {
            invalid = 1;
        } else {
            const double x = ghost[g];
            if (!isfinite(x) || x <= 0.0) invalid = 1;
            else sum = __dadd_rn(sum, x);
        }
    }
    reduce<false>(sum, invalid, sums, bad);
    if (t == 0) {
        if (bad[0]) s->status = 1;
        else {
            s->baseline_mean = __ddiv_rn(sums[0], double(n));
            if (!isfinite(s->baseline_mean)) s->status = 1;
            else s->phase = 1;
        }
    }
}

__global__ void shift_kernel(int n, int ng, const int *map, double *h,
                             double *ghost, AhfNewtonState *s, bool restore) {
    __shared__ int active;
    __shared__ double sums[threads];
    __shared__ int bad[threads];
    const int t = threadIdx.x;
    if (t == 0) {
        active = s->status == 0;
        if (active && s->phase != (restore ? 2 : 1)) {
            s->status = 2;
            active = 0;
        }
    }
    __syncthreads();
    int invalid = 0;
    if (active) {
        for (int g = t; g < ng; g += threads)
            if (!isfinite(ghost[g])) invalid = 1;
        for (int q = t; q < n; q += threads)
            if (map[q] < 0 || map[q] >= ng) invalid = 1;
    }
    const int any_invalid = __syncthreads_or(invalid);
    if (t == 0 && active && any_invalid) {
        s->status = 1;
        active = 0;
    }
    __syncthreads();
    if (active)
        for (int g = t; g < ng; g += threads)
            ghost[g] = __dadd_rn(ghost[g], restore ? -1e-6 : 1e-6);
    __syncthreads();
    double sum = 0.0;
    if (active) {
        for (int q = t; q < n; q += threads) {
            h[q] = ghost[map[q]];
            if (restore) sum = __dadd_rn(sum, h[q]);
        }
    }
    reduce<false>(sum, 0, sums, bad);
    if (t == 0 && active) {
        if (restore) s->restored_mean = __ddiv_rn(sums[0], double(n));
        s->phase = restore ? 3 : 2;
    }
}

__global__ void update_kernel(int n, int ng, const int *map, double *h,
                              double *ghost, const double *delta,
                              AhfNewtonState *s) {
    __shared__ int active;
    __shared__ double maxima[threads];
    __shared__ int bad[threads];
    const int t = threadIdx.x;
    if (t == 0) {
        active = s->status == 0;
        if (active && s->phase != 3) {
            s->status = 2;
            active = 0;
        }
    }
    __syncthreads();
    double max_delta = 0.0;
    int invalid = 0;
    if (active) {
        for (int q = t; q < n; q += threads) {
            const int g = map[q];
            if (g < 0 || g >= ng) invalid = 1;
            else if (!isfinite(ghost[g])) invalid = 1;
            const double d = delta[q];
            if (!isfinite(d)) invalid = 1;
            else {
                const double abs_delta = d > 0.0 ? d : -d;
                max_delta = max_delta > abs_delta ? max_delta : abs_delta;
            }
        }
    }
    reduce<true>(max_delta, invalid, maxima, bad);
    if (t == 0 && active) {
        if (bad[0] || !isfinite(s->baseline_mean) || s->baseline_mean <= 0.0) {
            s->status = 1;
            active = 0;
        } else {
            s->limit = __dmul_rn(0.1, s->baseline_mean);
            s->max_delta = maxima[0];
            s->scale = maxima[0] <= s->limit ? 1.0 : __ddiv_rn(s->limit, maxima[0]);
        }
    }
    __syncthreads();
    invalid = 0;
    if (active)
        for (int q = t; q < n; q += threads)
            if (!isfinite(__dadd_rn(ghost[map[q]], -__dmul_rn(s->scale, delta[q]))))
                invalid = 1;
    const int any_invalid = __syncthreads_or(invalid);
    if (t == 0 && active && any_invalid) {
        s->status = 1;
        active = 0;
    }
    __syncthreads();
    if (active) {
        for (int q = t; q < n; q += threads) {
            const int g = map[q];
            ghost[g] = __dadd_rn(ghost[g], -__dmul_rn(s->scale, delta[q]));
            h[q] = ghost[g];
        }
    }
    __syncthreads();
    if (t == 0 && active) s->phase = 4;
}

bool valid(int n, int ng, const int *map, const double *ghost,
           const AhfNewtonState *s) {
    return n > 0 && ng >= n && map && ghost && s;
}
}
extern "C" cudaError_t ahf_native_baseline_mean(cudaStream_t st, int n, int ng,
                                            const int *map, const double *ghost,
                                            AhfNewtonState *s) {
    if (!valid(n, ng, map, ghost, s)) return cudaErrorInvalidValue;
    baseline_kernel<<<1, threads, 0, st>>>(n, ng, map, ghost, s);
    return cudaGetLastError();
}
extern "C" cudaError_t ahf_native_radial_plus(cudaStream_t st, int n, int ng,
                                          const int *map, double *h,
                                          double *ghost, AhfNewtonState *s) {
    if (!valid(n, ng, map, ghost, s) || !h) return cudaErrorInvalidValue;
    shift_kernel<<<1, threads, 0, st>>>(n, ng, map, h, ghost, s, false);
    return cudaGetLastError();
}
extern "C" cudaError_t ahf_native_radial_restore(cudaStream_t st, int n, int ng,
                                             const int *map, double *h,
                                             double *ghost, AhfNewtonState *s) {
    if (!valid(n, ng, map, ghost, s) || !h) return cudaErrorInvalidValue;
    shift_kernel<<<1, threads, 0, st>>>(n, ng, map, h, ghost, s, true);
    return cudaGetLastError();
}
extern "C" cudaError_t ahf_native_newton_update(cudaStream_t st, int n, int ng,
                                            const int *map, double *h,
                                            double *ghost, const double *delta,
                                            AhfNewtonState *s) {
    if (!valid(n, ng, map, ghost, s) || !h || !delta) return cudaErrorInvalidValue;
    update_kernel<<<1, threads, 0, st>>>(n, ng, map, h, ghost, delta, s);
    return cudaGetLastError();
}
