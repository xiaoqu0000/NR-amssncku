// Task209 original implementation, 2026-09-30. No proprietary code or tables.
// Mathematical provenance: atan(x)=atan(a)+atan((x-a)/(1+a*x));
// atan(t)=sum_{k>=0}(-1)^k*t^(2*k+1)/(2*k+1), |t|<=1/4.
// atan(1/2) and atan(1) hi/lo constants generated with MPFR 4.2.2,
// RN at 256 bits, binary64 rounding independently checked at 512 bits.
// This is a measured high-precision approximation, NOT a universal correctly-
// rounded or universal Intel libimf compatibility guarantee. See host-reference.json.
#ifndef AHF_DEVICE_ATAN_CUH
#define AHF_DEVICE_ATAN_CUH
#include <cuda_runtime.h>
#include <math.h>

namespace ahf_atan_detail {
struct DD { double hi, lo; };

__device__ __forceinline__ DD normalize(double a, double b) {
    double s = __dadd_rn(a, b);
    double v = __dsub_rn(s, a);
    double e = __dadd_rn(__dsub_rn(a, __dsub_rn(s, v)), __dsub_rn(b, v));
    return {s, e};
}
__device__ __forceinline__ DD add(DD a, DD b) {
    DD h = normalize(a.hi, b.hi);
    DD l = normalize(a.lo, b.lo);
    DD m = normalize(h.lo, l.hi);
    DD r = normalize(h.hi, m.hi);
    return normalize(r.hi, __dadd_rn(__dadd_rn(r.lo, m.lo), l.lo));
}
__device__ __forceinline__ DD negate(DD a) { return {-a.hi, -a.lo}; }
__device__ __forceinline__ DD sub(DD a, DD b) { return add(a, negate(b)); }
__device__ __forceinline__ DD mul(DD a, DD b) {
    double p = __dmul_rn(a.hi, b.hi);
    double e = __fma_rn(a.hi, b.hi, -p);
    e = __dadd_rn(e, __dmul_rn(a.hi, b.lo));
    e = __dadd_rn(e, __dmul_rn(a.lo, b.hi));
    e = __dadd_rn(e, __dmul_rn(a.lo, b.lo));
    return normalize(p, e);
}
__device__ __forceinline__ DD div(DD a, DD b) {
    double q1 = __ddiv_rn(a.hi, b.hi);
    DD r = sub(a, mul(b, {q1, 0.0}));
    double q2 = __ddiv_rn(r.hi, b.hi);
    r = sub(r, mul(b, {q2, 0.0}));
    double q3 = __ddiv_rn(__dadd_rn(r.hi, r.lo), b.hi);
    return add(normalize(q1, q2), {q3, 0.0});
}
} // namespace ahf_atan_detail

/** Device-only atan for rounded binary64 shell face ratios in [-1,1].
 * Preserves the sign of zero, including subnormal input. Nonfinite or out-of-
 * domain inputs return NaN (the shell caller must first select a dominant face).
 * No host state, input-dependent lookup, allocation, or runtime initialization.
 * Requires ordinary binary64 CUDA math; do not build with --use_fast_math.
 */
__device__ __forceinline__ double ahf_device_atan(double x) {
    using namespace ahf_atan_detail;
    if (!(fabs(x) <= 1.0)) return NAN;
    // For these inputs the cubic correction is below half a binary64 ulp.
    if (fabs(x) <= 0x1p-27) return x;
    double u = fabs(x);
    DD anchor{0.0, 0.0};
    DD t{u, 0.0};
    if (u > 0.75) {
        anchor = {0x1.921fb54442d18p-1, 0x1.1a62633145c07p-55};
        t = div(sub({u, 0.0}, {1.0, 0.0}), add({u, 0.0}, {1.0, 0.0}));
    } else if (u > 0.25) {
        anchor = {0x1.dac670561bb4fp-2, 0x1.a2b7f222f65e2p-56};
        t = div(sub({u, 0.0}, {0.5, 0.0}), add(mul({u, 0.0}, {0.5, 0.0}), {1.0, 0.0}));
    }
    DD minus_t2 = ahf_atan_detail::negate(mul(t, t));
    DD power = t, sum = t;
    // 28 terms; neglected alternating tail <= (1/4)^57/57 < 2^-119.
    // Fixed mathematical order, independent of fixtures or CPU reference values.
    #pragma unroll 1
    for (int k = 1; k < 28; ++k) {
        power = mul(power, minus_t2);
        sum = add(sum, div(power, {static_cast<double>(2*k+1), 0.0}));
    }
    DD result = add(anchor, sum);
    return copysign(__dadd_rn(result.hi, result.lo), x);
}
#endif
