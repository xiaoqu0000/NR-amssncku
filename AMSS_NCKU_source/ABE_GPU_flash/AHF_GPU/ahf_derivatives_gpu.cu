#include "ahf_derivatives_gpu.h"

#include <cstdlib>
#include <iostream>

#include "diff_new_gpu.cuh"

namespace {

__global__ void ahf_field_derivatives_kernel(
    int nx, int ny, int nz,
    const double* X, const double* Y, const double* Z,
    const double* field,
    double* dx, double* dy, double* dz,
    double sym1, double sym2, double sym3,
    int symmetry)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= nx || j >= ny || k >= nz)
        return;

    const int ex[3] = {nx, ny, nz};
    const int n = (k * ny + j) * nx + i;
    d_fderivs_point(
        ex, field,
        dx + n, dy + n, dz + n,
        X, Y, Z,
        sym1, sym2, sym3,
        symmetry, 0,
        i, j, k);
}

void fail_if_cuda_error(cudaError_t status, const char* op)
{
    if (status != cudaSuccess) {
        std::cerr << op << ": " << cudaGetErrorString(status) << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void launch_field_derivatives(
    cudaStream_t stream,
    int nx, int ny, int nz,
    const double* X, const double* Y, const double* Z,
    const double* field,
    double* dx, double* dy, double* dz,
    double sym1, double sym2, double sym3,
    int symmetry)
{
    const dim3 block(8, 8, 4);
    const dim3 grid(
        (nx + block.x - 1) / block.x,
        (ny + block.y - 1) / block.y,
        (nz + block.z - 1) / block.z);

    ahf_field_derivatives_kernel<<<grid, block, 0, stream>>>(
        nx, ny, nz,
        X, Y, Z,
        field,
        dx, dy, dz,
        sym1, sym2, sym3,
        symmetry);
    fail_if_cuda_error(cudaGetLastError(), "ahf_field_derivatives_kernel");
}

} // namespace

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
    int symmetry)
{
    const int nx = ex[0];
    const int ny = ex[1];
    const int nz = ex[2];

    launch_field_derivatives(stream, nx, ny, nz, X, Y, Z,
                             phi, phi_x, phi_y, phi_z,
                             1.0, 1.0, 1.0, symmetry);
    launch_field_derivatives(stream, nx, ny, nz, X, Y, Z,
                             gxx, gxx_x, gxx_y, gxx_z,
                             1.0, 1.0, 1.0, symmetry);
    launch_field_derivatives(stream, nx, ny, nz, X, Y, Z,
                             gxy, gxy_x, gxy_y, gxy_z,
                             -1.0, -1.0, 1.0, symmetry);
    launch_field_derivatives(stream, nx, ny, nz, X, Y, Z,
                             gxz, gxz_x, gxz_y, gxz_z,
                             -1.0, 1.0, -1.0, symmetry);
    launch_field_derivatives(stream, nx, ny, nz, X, Y, Z,
                             gyy, gyy_x, gyy_y, gyy_z,
                             1.0, 1.0, 1.0, symmetry);
    launch_field_derivatives(stream, nx, ny, nz, X, Y, Z,
                             gyz, gyz_x, gyz_y, gyz_z,
                             1.0, -1.0, -1.0, symmetry);
    launch_field_derivatives(stream, nx, ny, nz, X, Y, Z,
                             gzz, gzz_x, gzz_y, gzz_z,
                             1.0, 1.0, 1.0, symmetry);
}
