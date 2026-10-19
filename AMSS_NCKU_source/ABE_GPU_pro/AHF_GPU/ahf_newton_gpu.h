#ifndef AHF_NATIVE_NEWTON_H
#define AHF_NATIVE_NEWTON_H
#include <cuda_runtime_api.h>

// AHF Newton operations with FP64 block reductions and original point arithmetic.
// All numerical pointers below are DEVICE pointers. nominal_ghost maps dynamic
// n nominal points to distinct ghost slots in original pn/rho/sigma order;
// ghost_h contains ng packed ghosted points, including unused corner storage.
// h is the nominal mirror, kept consistent with ghost_h[nominal_ghost[q]].
// Caller owns allocation, ghost synchronization/expansion, and stream ordering.
// No allocations, host numerical work, uploads, downloads or synchronizations.
struct AhfNewtonState {
    double baseline_mean, restored_mean, limit, max_delta, scale;
    int status, phase; // 0 success; status 1 invalid numbers, 2 invalid sequence
};

extern "C" {
// Must execute BEFORE any radial perturbation or expansion_Jacobian work.
cudaError_t ahf_native_baseline_mean(cudaStream_t stream, int n, int ng,
                               const int *nominal_ghost, const double *ghost_h,
                               AhfNewtonState *state);
// Add the original epsilon=1e-6 to ALL ghosted slots and refresh nominal h.
// Caller then performs its normal ghost synchronization and plus expansion.
cudaError_t ahf_native_radial_plus(cudaStream_t stream, int n, int ng,
                             const int *nominal_ghost, double *h,
                             double *ghost_h, AhfNewtonState *state);
// Call only on the original successful radial-expansion path. Subtract 1e-6
// from the CURRENT ghosted h (never saved bitcopy). Do NOT sync afterwards:
// original restore leaves nonnominal ghosts minus epsilon until next expansion.
cudaError_t ahf_native_radial_restore(cudaStream_t stream, int n, int ng,
                                const int *nominal_ghost, double *h,
                                double *ghost_h, AhfNewtonState *state);
// Uses saved baseline_mean; block max(abs(delta)) and parallel point updates.
// Updates nominal slots only, leaving other ghosts unchanged. factor=0.1 is
// fixed to the original scheduler value. No solver/stopping changes.
cudaError_t ahf_native_newton_update(cudaStream_t stream, int n, int ng,
                               const int *nominal_ghost, double *h,
                               double *ghost_h, const double *delta,
                               AhfNewtonState *state);
}
#endif
