#include "fadmquantites_bssn.h"

#include <cuda_runtime.h>
#include <math.h>

// 假设 fmisc.h 中声明了 d_fderivs_point
#include "fmisc.h" 
#include "gpu_runtime_coordinator.h"
#include "derivatives.h"
#include "diff_new_gpu.cuh"

namespace {

constexpr double kUnitValue = 1.0;
constexpr double kQuarterValue = 0.25;
constexpr double kOneEighthValue = 1.0 / 8.0;
constexpr double kInverseTwoPi = 1.0 / (2.0 * M_PI);
constexpr double kConformalWeightExponent = -1.25;

inline dim3 adm_mass_block_shape() {
    return dim3(8, 8, 8);
}

inline dim3 adm_mass_grid_shape(const int extent[3]) {
    const dim3 block = adm_mass_block_shape();
    return dim3(
        (extent[0] + block.x - 1) / block.x,
        (extent[1] + block.y - 1) / block.y,
        (extent[2] + block.z - 1) / block.z
    );
}

__device__ __forceinline__ int flatten_cell_index(int i0, int j0, int k0, int extent_x, int extent_y) {
    return (k0 * extent_y + j0) * extent_x + i0;
}

__device__ __forceinline__ void fill_inverse_metric(
    double metric_xx,
    double metric_xy,
    double metric_xz,
    double metric_yy,
    double metric_yz,
    double metric_zz,
    double& inverse_xx,
    double& inverse_xy,
    double& inverse_xz,
    double& inverse_yy,
    double& inverse_yz,
    double& inverse_zz
) {
    const double metric_det = metric_xx * metric_yy * metric_zz
                            + 2.0 * metric_xy * metric_xz * metric_yz
                            - metric_yy * metric_xz * metric_xz
                            - metric_zz * metric_xy * metric_xy
                            - metric_xx * metric_yz * metric_yz;

    inverse_xx = (metric_yy * metric_zz - metric_yz * metric_yz) / metric_det;
    inverse_xy = (metric_yz * metric_xz - metric_xy * metric_zz) / metric_det;
    inverse_xz = (metric_xy * metric_yz - metric_yy * metric_xz) / metric_det;
    inverse_yy = (metric_xx * metric_zz - metric_xz * metric_xz) / metric_det;
    inverse_yz = (metric_xy * metric_xz - metric_xx * metric_yz) / metric_det;
    inverse_zz = (metric_xx * metric_yy - metric_xy * metric_xy) / metric_det;
}

}  // namespace

__global__ void assemble_adm_mass_density_kernel(
    int extent_x, int extent_y, int extent_z,
    const double* __restrict__ coord_x, const double* __restrict__ coord_y, const double* __restrict__ coord_z,
    const double* __restrict__ conformal_chi, const double* __restrict__ trace_k,
    const double* __restrict__ metric_xx, const double* __restrict__ metric_xy, const double* __restrict__ metric_xz,
    const double* __restrict__ metric_yy, const double* __restrict__ metric_yz, const double* __restrict__ metric_zz,
    const double* __restrict__ tracefree_xx, const double* __restrict__ tracefree_xy, const double* __restrict__ tracefree_xz,
    const double* __restrict__ tracefree_yy, const double* __restrict__ tracefree_yz, const double* __restrict__ tracefree_zz,
    const double* __restrict__ gamma_x, const double* __restrict__ gamma_y, const double* __restrict__ gamma_z,
    double* __restrict__ mass_density_x, double* __restrict__ mass_density_y, double* __restrict__ mass_density_z,
    int symmetry_mode
);

void gpu_admmass_bssn_launch(
    cudaStream_t stream, const int ext[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* chi, const double* trK,
    const double* dxx, const double* gxy, const double* gxz, 
    const double* dyy, const double* gyz, const double* dzz,
    const double* Axx, const double* Axy, const double* Axz, 
    const double* Ayy, const double* Ayz, const double* Azz, 
    const double* Gamx, const double* Gamy, const double* Gamz,
    double* massx, double* massy, double* massz,
    int symmetry
) {
    const dim3 block = adm_mass_block_shape();
    const dim3 grid = adm_mass_grid_shape(ext);

    assemble_adm_mass_density_kernel<<<grid, block, 0, stream>>>(
        ext[0], ext[1], ext[2], 
        d_X, d_Y, d_Z,
        chi, trK, 
        dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, 
        massx, massy, massz, 
        symmetry
    );
}

__global__ void assemble_adm_mass_density_kernel(
    int extent_x, int extent_y, int extent_z,
    const double* __restrict__ coord_x, const double* __restrict__ coord_y, const double* __restrict__ coord_z,
    const double* __restrict__ conformal_chi, const double* __restrict__ trace_k,
    const double* __restrict__ metric_xx, const double* __restrict__ metric_xy, const double* __restrict__ metric_xz,
    const double* __restrict__ metric_yy, const double* __restrict__ metric_yz, const double* __restrict__ metric_zz,
    const double* __restrict__ tracefree_xx, const double* __restrict__ tracefree_xy, const double* __restrict__ tracefree_xz,
    const double* __restrict__ tracefree_yy, const double* __restrict__ tracefree_yz, const double* __restrict__ tracefree_zz,
    const double* __restrict__ gamma_x, const double* __restrict__ gamma_y, const double* __restrict__ gamma_z,
    double* __restrict__ mass_density_x, double* __restrict__ mass_density_y, double* __restrict__ mass_density_z,
    int symmetry_mode
) {
    const int i0 = blockIdx.x * blockDim.x + threadIdx.x;
    const int j0 = blockIdx.y * blockDim.y + threadIdx.y;
    const int k0 = blockIdx.z * blockDim.z + threadIdx.z;

    if (i0 >= extent_x || j0 >= extent_y || k0 >= extent_z) {
        return;
    }

    const int flat_index = flatten_cell_index(i0, j0, k0, extent_x, extent_y);
    int extent[3] = {extent_x, extent_y, extent_z};

    double chi_gradient_x = 0.0;
    double chi_gradient_y = 0.0;
    double chi_gradient_z = 0.0;

    d_fderivs_point(
        extent,
        conformal_chi,
        &chi_gradient_x,
        &chi_gradient_y,
        &chi_gradient_z,
        coord_x,
        coord_y,
        coord_z,
        kUnitValue,
        kUnitValue,
        kUnitValue,
        symmetry_mode,
        0,
        i0,
        j0,
        k0
    );

    const double physical_metric_xx = metric_xx[flat_index] + kUnitValue;
    const double physical_metric_yy = metric_yy[flat_index] + kUnitValue;
    const double physical_metric_zz = metric_zz[flat_index] + kUnitValue;
    const double physical_metric_xy = metric_xy[flat_index];
    const double physical_metric_xz = metric_xz[flat_index];
    const double physical_metric_yz = metric_yz[flat_index];

    double inverse_metric_xx = 0.0;
    double inverse_metric_xy = 0.0;
    double inverse_metric_xz = 0.0;
    double inverse_metric_yy = 0.0;
    double inverse_metric_yz = 0.0;
    double inverse_metric_zz = 0.0;
    fill_inverse_metric(
        physical_metric_xx,
        physical_metric_xy,
        physical_metric_xz,
        physical_metric_yy,
        physical_metric_yz,
        physical_metric_zz,
        inverse_metric_xx,
        inverse_metric_xy,
        inverse_metric_xz,
        inverse_metric_yy,
        inverse_metric_yz,
        inverse_metric_zz
    );

    const double chi_prefactor = kQuarterValue * pow(conformal_chi[flat_index] + kUnitValue, kConformalWeightExponent);

    (void)trace_k;
    (void)tracefree_xx;
    (void)tracefree_xy;
    (void)tracefree_xz;
    (void)tracefree_yy;
    (void)tracefree_yz;
    (void)tracefree_zz;

    mass_density_x[flat_index] =
        (kOneEighthValue * gamma_x[flat_index]
         + chi_prefactor * (
               inverse_metric_xx * chi_gradient_x
             + inverse_metric_xy * chi_gradient_y
             + inverse_metric_xz * chi_gradient_z)) * kInverseTwoPi;
    mass_density_y[flat_index] =
        (kOneEighthValue * gamma_y[flat_index]
         + chi_prefactor * (
               inverse_metric_xy * chi_gradient_x
             + inverse_metric_yy * chi_gradient_y
             + inverse_metric_yz * chi_gradient_z)) * kInverseTwoPi;
    mass_density_z[flat_index] =
        (kOneEighthValue * gamma_z[flat_index]
         + chi_prefactor * (
               inverse_metric_xz * chi_gradient_x
             + inverse_metric_yz * chi_gradient_y
             + inverse_metric_zz * chi_gradient_z)) * kInverseTwoPi;
}
