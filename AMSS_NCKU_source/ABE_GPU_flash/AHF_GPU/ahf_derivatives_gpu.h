#ifndef AHF_DERIVATIVES_GPU_H
#define AHF_DERIVATIVES_GPU_H

#include <cuda_runtime.h>

void gpu_ah_prepare_derivatives_launch(
    cudaStream_t& stream,
    int* ex,
    double* X, double* Y, double* Z,
    double* phi,
    double* phi_x, double* phi_y, double* phi_z,
    double* gxx,
    double* gxx_x, double* gxx_y, double* gxx_z,
    double* gxy,
    double* gxy_x, double* gxy_y, double* gxy_z,
    double* gxz,
    double* gxz_x, double* gxz_y, double* gxz_z,
    double* gyy,
    double* gyy_x, double* gyy_y, double* gyy_z,
    double* gyz,
    double* gyz_x, double* gyz_y, double* gyz_z,
    double* gzz,
    double* gzz_x, double* gzz_y, double* gzz_z,
    int symmetry);

#endif /* AHF_DERIVATIVES_GPU_H */
