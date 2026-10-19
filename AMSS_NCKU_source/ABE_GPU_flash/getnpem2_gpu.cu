#include "getnpem2.h"

#include <cuda_runtime.h>

#include "getnpem2_gpu.cuh"

namespace {

__global__ void getnpem2_kernel(
    int ex0,
    int ex1,
    int ex2,
    const double* x_axis,
    const double* y_axis,
    const double* z_axis,
    const double* chi,
    const double* dxx,
    const double* gxy,
    const double* gxz,
    const double* dyy,
    const double* gyz,
    const double* dzz,
    const double* ex,
    const double* ey,
    const double* ez,
    const double* bx,
    const double* by,
    const double* bz,
    double* rphi2,
    double* iphi2,
    int symmetry
) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) {
        return;
    }

    const int idx = i + j * ex0 + k * ex0 * ex1;
    const amss_gpu_em::NpValue phi2 = amss_gpu_em::compute_phi2_point(
        x_axis[i], y_axis[j], z_axis[k],
        chi[idx],
        dxx[idx], gxy[idx], gxz[idx], dyy[idx], gyz[idx], dzz[idx],
        ex[idx], ey[idx], ez[idx], bx[idx], by[idx], bz[idx]
    );
    rphi2[idx] = phi2.real;
    iphi2[idx] = phi2.imag;
    (void)symmetry;
}

__global__ void getnpem1_kernel(
    int ex0,
    int ex1,
    int ex2,
    const double* x_axis,
    const double* y_axis,
    const double* z_axis,
    const double* chi,
    const double* dxx,
    const double* gxy,
    const double* gxz,
    const double* dyy,
    const double* gyz,
    const double* dzz,
    const double* ex,
    const double* ey,
    const double* ez,
    const double* bx,
    const double* by,
    const double* bz,
    double* rphi1,
    double* iphi1,
    int symmetry
) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) {
        return;
    }

    const int idx = i + j * ex0 + k * ex0 * ex1;
    const amss_gpu_em::NpValue phi1 = amss_gpu_em::compute_phi1_point(
        x_axis[i], y_axis[j], z_axis[k],
        chi[idx],
        dxx[idx], gxy[idx], gxz[idx], dyy[idx], gyz[idx], dzz[idx],
        ex[idx], ey[idx], ez[idx], bx[idx], by[idx], bz[idx]
    );
    rphi1[idx] = phi1.real;
    iphi1[idx] = phi1.imag;
    (void)symmetry;
}

}  // namespace

void gpu_getnpem2_launch(
    cudaStream_t stream,
    int ex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Ex, const double* d_Ey, const double* d_Ez,
    const double* d_Bx, const double* d_By, const double* d_Bz,
    double* d_Rphi2, double* d_Iphi2,
    int symmetry
) {
    const dim3 block_dim(8, 8, 4);
    const dim3 grid_dim(
        (ex[0] + block_dim.x - 1) / block_dim.x,
        (ex[1] + block_dim.y - 1) / block_dim.y,
        (ex[2] + block_dim.z - 1) / block_dim.z
    );

    getnpem2_kernel<<<grid_dim, block_dim, 0, stream>>>(
        ex[0], ex[1], ex[2],
        d_X, d_Y, d_Z,
        d_chi,
        d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
        d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
        d_Rphi2, d_Iphi2,
        symmetry
    );
}

void gpu_getnpem1_launch(
    cudaStream_t stream,
    int ex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Ex, const double* d_Ey, const double* d_Ez,
    const double* d_Bx, const double* d_By, const double* d_Bz,
    double* d_Rphi1, double* d_Iphi1,
    int symmetry
) {
    const dim3 block_dim(8, 8, 4);
    const dim3 grid_dim(
        (ex[0] + block_dim.x - 1) / block_dim.x,
        (ex[1] + block_dim.y - 1) / block_dim.y,
        (ex[2] + block_dim.z - 1) / block_dim.z
    );

    getnpem1_kernel<<<grid_dim, block_dim, 0, stream>>>(
        ex[0], ex[1], ex[2],
        d_X, d_Y, d_Z,
        d_chi,
        d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
        d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
        d_Rphi1, d_Iphi1,
        symmetry
    );
}
