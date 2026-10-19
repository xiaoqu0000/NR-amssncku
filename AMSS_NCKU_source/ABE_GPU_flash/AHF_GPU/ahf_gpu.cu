#include "ahf_gpu.h"

#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>

namespace {

constexpr int kBlockSize = 256;

__global__ void ahf_radial_residual_newton_kernel(
    int point_count,
    const double* unit_x,
    const double* unit_y,
    const double* unit_z,
    double solver_center_x,
    double solver_center_y,
    double solver_center_z,
    double physical_center_x,
    double physical_center_y,
    double physical_center_z,
    double target_radius,
    double* horizon_radius,
    double* residual_out,
    double* jacobian_diag_out,
    int newton_iterations)
{
    const int point_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_index >= point_count)
        return;

    const double ux = unit_x[point_index];
    const double uy = unit_y[point_index];
    const double uz = unit_z[point_index];

    for (int iteration = 0; iteration < newton_iterations; ++iteration) {
        const double h = horizon_radius[point_index];
        const double px = solver_center_x + h * ux;
        const double py = solver_center_y + h * uy;
        const double pz = solver_center_z + h * uz;
        const double dx = px - physical_center_x;
        const double dy = py - physical_center_y;
        const double dz = pz - physical_center_z;
        const double norm = sqrt(dx * dx + dy * dy + dz * dz);
        const double residual = norm - target_radius;
        const double jacobian_diag = (norm > 0.0) ? (dx * ux + dy * uy + dz * uz) / norm : 1.0;

        residual_out[point_index] = residual;
        jacobian_diag_out[point_index] = jacobian_diag;
        if (fabs(jacobian_diag) > 1.0e-14)
            horizon_radius[point_index] -= residual / jacobian_diag;
    }
}

inline int launch_blocks_for(int point_count)
{
    return (point_count + kBlockSize - 1) / kBlockSize;
}

} // namespace

void gpu_ahf_radial_residual_newton_launch(
    cudaStream_t stream,
    int point_count,
    const double* d_unit_x,
    const double* d_unit_y,
    const double* d_unit_z,
    double solver_center_x,
    double solver_center_y,
    double solver_center_z,
    double physical_center_x,
    double physical_center_y,
    double physical_center_z,
    double target_radius,
    double* d_horizon_radius,
    double* d_residual,
    double* d_jacobian_diag,
    int newton_iterations)
{
    ahf_radial_residual_newton_kernel<<<launch_blocks_for(point_count), kBlockSize, 0, stream>>>(
        point_count,
        d_unit_x, d_unit_y, d_unit_z,
        solver_center_x, solver_center_y, solver_center_z,
        physical_center_x, physical_center_y, physical_center_z,
        target_radius,
        d_horizon_radius,
        d_residual,
        d_jacobian_diag,
        newton_iterations);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "gpu_ahf_radial_residual_newton_launch: %s\n", cudaGetErrorString(err));
    }
}
