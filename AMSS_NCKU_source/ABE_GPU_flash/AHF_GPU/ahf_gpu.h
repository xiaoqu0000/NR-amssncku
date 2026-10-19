#ifndef AHF_GPU_H
#define AHF_GPU_H

#include <cuda_runtime.h>

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
    int newton_iterations);

#endif /* AHF_GPU_H */
