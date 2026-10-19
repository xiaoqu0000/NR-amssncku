#include "enforce_algebra.h"

#include <cuda_runtime.h>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "gpu_runtime_coordinator.h"

constexpr int kEnforceGaBlockSize = 256;
constexpr double kUnitValue = 1.0;
constexpr double kOneThird = 1.0 / 3.0;
constexpr double kSymmetricPairFactor = 2.0;

__global__ void project_metric_trace_kernel(
    int cell_count,
    double* conformal_xx, double* conformal_xy, double* conformal_xz,
    double* conformal_yy, double* conformal_yz, double* conformal_zz,
    double* tracefree_xx, double* tracefree_xy, double* tracefree_xz,
    double* tracefree_yy, double* tracefree_yz, double* tracefree_zz
);

inline int flatten_extent(const int extent[3]) {
    return extent[0] * extent[1] * extent[2];
}

inline int launch_blocks_for(int item_count) {
    return (item_count + kEnforceGaBlockSize - 1) / kEnforceGaBlockSize;
}

void trace_enforce_det_summary(
    cudaStream_t stream,
    const char* phase,
    const int ex[3],
    double* d_dxx, double* d_gxy, double* d_gxz,
    double* d_dyy, double* d_gyz, double* d_dzz
) {
    if (!getenv("AMSS_TRAJECTORY_PREFIX_TRACE_ENFORCE_DET"))
        return;

    const int nx = ex[0];
    const int ny = ex[1];
    const int nz = ex[2];
    const int cell_count = flatten_extent(ex);
    std::vector<double> h_dxx(cell_count), h_gxy(cell_count), h_gxz(cell_count);
    std::vector<double> h_dyy(cell_count), h_gyz(cell_count), h_dzz(cell_count);

    CUDA_CHECK(cudaMemcpyAsync(h_dxx.data(), d_dxx, static_cast<size_t>(cell_count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_gxy.data(), d_gxy, static_cast<size_t>(cell_count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_gxz.data(), d_gxz, static_cast<size_t>(cell_count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_dyy.data(), d_dyy, static_cast<size_t>(cell_count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_gyz.data(), d_gyz, static_cast<size_t>(cell_count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_dzz.data(), d_dzz, static_cast<size_t>(cell_count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));

    int nonfinite_metric_count = 0;
    int nonfinite_det_count = 0;
    int nonpositive_det_count = 0;
    int first_i = -1;
    int first_j = -1;
    int first_k = -1;
    double first_det = 0.0;
    double first_metric_xx = 0.0;
    double first_metric_xy = 0.0;
    double first_metric_xz = 0.0;
    double first_metric_yy = 0.0;
    double first_metric_yz = 0.0;
    double first_metric_zz = 0.0;

    for (int idx = 0; idx < cell_count; ++idx) {
        const double metric_xx = h_dxx[idx] + kUnitValue;
        const double metric_xy = h_gxy[idx];
        const double metric_xz = h_gxz[idx];
        const double metric_yy = h_dyy[idx] + kUnitValue;
        const double metric_yz = h_gyz[idx];
        const double metric_zz = h_dzz[idx] + kUnitValue;
        const bool metric_finite = std::isfinite(metric_xx) && std::isfinite(metric_xy) &&
            std::isfinite(metric_xz) && std::isfinite(metric_yy) &&
            std::isfinite(metric_yz) && std::isfinite(metric_zz);
        if (!metric_finite)
            ++nonfinite_metric_count;

        const double det = metric_xx * metric_yy * metric_zz + metric_xy * metric_yz * metric_xz
                         + metric_xz * metric_xy * metric_yz
                         - metric_xz * metric_yy * metric_xz
                         - metric_xy * metric_xy * metric_zz
                         - metric_xx * metric_yz * metric_yz;
        if (!std::isfinite(det))
            ++nonfinite_det_count;
        if (!std::isfinite(det) || det <= 0.0) {
            ++nonpositive_det_count;
            if (first_i < 0) {
                first_i = idx % nx;
                first_j = (idx / nx) % ny;
                first_k = idx / (nx * ny);
                first_det = det;
                first_metric_xx = metric_xx;
                first_metric_xy = metric_xy;
                first_metric_xz = metric_xz;
                first_metric_yy = metric_yy;
                first_metric_yz = metric_yz;
                first_metric_zz = metric_zz;
            }
        }
    }

    std::cout << " enforce_det"
              << " phase=" << phase
              << " shape=(" << nx << "," << ny << "," << nz << ")"
              << " nonfinite_metric_count=" << nonfinite_metric_count
              << " nonfinite_det_count=" << nonfinite_det_count
              << " nonpositive_det_count=" << nonpositive_det_count
              << " first_ijk=(" << first_i << "," << first_j << "," << first_k << ")"
              << " first_det=" << first_det
              << " first_metric=(" << first_metric_xx << "," << first_metric_xy << ","
              << first_metric_xz << "," << first_metric_yy << "," << first_metric_yz
              << "," << first_metric_zz << ")"
              << std::endl;
}

void gpu_enforce_ga_launch(
    cudaStream_t &stream,
    int ex[3],
    double* d_dxx, double* d_gxy, double* d_gxz,
    double* d_dyy, double* d_gyz, double* d_dzz,
    double* d_Axx, double* d_Axy, double* d_Axz,
    double* d_Ayy, double* d_Ayz, double* d_Azz
) {
    const int cell_count = flatten_extent(ex);
    const int grid_span = launch_blocks_for(cell_count);

    trace_enforce_det_summary(
        stream, "before", ex,
        d_dxx, d_gxy, d_gxz,
        d_dyy, d_gyz, d_dzz
    );

    project_metric_trace_kernel<<<grid_span, kEnforceGaBlockSize, 0, stream>>>(
        cell_count,
        d_dxx, d_gxy, d_gxz,
        d_dyy, d_gyz, d_dzz,
        d_Axx, d_Axy, d_Axz, 
        d_Ayy, d_Ayz, d_Azz
    );

    trace_enforce_det_summary(
        stream, "after", ex,
        d_dxx, d_gxy, d_gxz,
        d_dyy, d_gyz, d_dzz
    );
}

__global__ void project_metric_trace_kernel(
    int cell_count,
    double* conformal_xx, double* conformal_xy, double* conformal_xz,
    double* conformal_yy, double* conformal_yz, double* conformal_zz,
    double* tracefree_xx, double* tracefree_xy, double* tracefree_xz,
    double* tracefree_yy, double* tracefree_yz, double* tracefree_zz
) {
    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_id >= cell_count) {
        return;
    }

    double metric_xx = conformal_xx[point_id] + kUnitValue;
    double metric_yy = conformal_yy[point_id] + kUnitValue;
    double metric_zz = conformal_zz[point_id] + kUnitValue;
    double metric_xy = conformal_xy[point_id];
    double metric_xz = conformal_xz[point_id];
    double metric_yz = conformal_yz[point_id];

    const double metric_det = metric_xx * metric_yy * metric_zz + metric_xy * metric_yz * metric_xz
                            + metric_xz * metric_xy * metric_yz
                            - metric_xz * metric_yy * metric_xz
                            - metric_xy * metric_xy * metric_zz
                            - metric_xx * metric_yz * metric_yz;
    const double det_scale = kUnitValue / pow(metric_det, kOneThird);

    metric_xx *= det_scale;
    metric_xy *= det_scale;
    metric_xz *= det_scale;
    metric_yy *= det_scale;
    metric_yz *= det_scale;
    metric_zz *= det_scale;

    conformal_xx[point_id] = metric_xx - kUnitValue;
    conformal_xy[point_id] = metric_xy;
    conformal_xz[point_id] = metric_xz;
    conformal_yy[point_id] = metric_yy - kUnitValue;
    conformal_yz[point_id] = metric_yz;
    conformal_zz[point_id] = metric_zz - kUnitValue;

    const double inverse_xx = metric_yy * metric_zz - metric_yz * metric_yz;
    const double inverse_xy = metric_yz * metric_xz - metric_xy * metric_zz;
    const double inverse_xz = metric_xy * metric_yz - metric_yy * metric_xz;
    const double inverse_yy = metric_xx * metric_zz - metric_xz * metric_xz;
    const double inverse_yz = metric_xy * metric_xz - metric_xx * metric_yz;
    const double inverse_zz = metric_xx * metric_yy - metric_xy * metric_xy;

    const double mixed_trace = inverse_xx * tracefree_xx[point_id]
                             + inverse_yy * tracefree_yy[point_id]
                             + inverse_zz * tracefree_zz[point_id]
                             + kSymmetricPairFactor * (
                                   inverse_xy * tracefree_xy[point_id]
                                 + inverse_xz * tracefree_xz[point_id]
                                 + inverse_yz * tracefree_yz[point_id]
                               );

    tracefree_xx[point_id] -= kOneThird * metric_xx * mixed_trace;
    tracefree_xy[point_id] -= kOneThird * metric_xy * mixed_trace;
    tracefree_xz[point_id] -= kOneThird * metric_xz * mixed_trace;
    tracefree_yy[point_id] -= kOneThird * metric_yy * mixed_trace;
    tracefree_yz[point_id] -= kOneThird * metric_yz * mixed_trace;
    tracefree_zz[point_id] -= kOneThird * metric_zz * mixed_trace;
}
