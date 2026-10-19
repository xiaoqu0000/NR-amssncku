#include "surface_integral.h"

#include <cuda_runtime.h>
#include <device_functions.h>
#include <device_launch_parameters.h>

#include "data_motion_helper.h"
#include "fadmquantites_bssn.h"
#include "getnpem2_gpu.cuh"
#include "misc_gpu.cuh"

#define PI M_PI

namespace {

constexpr int kLinearBlockSize = 256;
constexpr double kOne = 1.0;
constexpr double kOneEighth = 0.125;

struct RankWorkRange {
    int begin;
    int end;
};

struct DevicePointTriplet {
    double* values[3];
};

int linear_grid_size(int item_count) {
    return (item_count + kLinearBlockSize - 1) / kLinearBlockSize;
}

RankWorkRange split_rank_range(int total_count, int rank_count, int rank_id) {
    const int base_count = total_count / rank_count;
    const int extra_count = total_count - rank_count * base_count;
    if (extra_count > rank_id) {
        const int begin = rank_id * base_count + rank_id;
        return RankWorkRange{begin, begin + base_count};
    }

    const int begin = rank_id * base_count + extra_count;
    return RankWorkRange{begin, begin + base_count - 1};
}

void synchronize_level_streams(cgh* grid_hierarchy, int lev, int myrank) {
    for (MyList<Patch>* patch_node = grid_hierarchy->PatL[lev]; patch_node; patch_node = patch_node->next) {
        for (MyList<Block>* block_node = patch_node->data->blb; block_node; block_node = block_node->next) {
            Block* block = block_node->data;
            if (myrank == block->rank) {
                CUDA_CHECK(cudaStreamSynchronize(block->stream));
            }
            if (block_node == patch_node->data->ble) {
                break;
            }
        }
    }
}

__device__ __forceinline__ void invert_metric_tensor(
    double gxx,
    double gxy,
    double gxz,
    double gyy,
    double gyz,
    double gzz,
    double& gupxx,
    double& gupxy,
    double& gupxz,
    double& gupyy,
    double& gupyz,
    double& gupzz
) {
    const double determinant = gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz
        - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz;
    gupxx = (gyy * gzz - gyz * gyz) / determinant;
    gupxy = -(gxy * gzz - gyz * gxz) / determinant;
    gupxz = (gxy * gyz - gyy * gxz) / determinant;
    gupyy = (gxx * gzz - gxz * gxz) / determinant;
    gupyz = -(gxx * gyz - gxy * gxz) / determinant;
    gupzz = (gxx * gyy - gxy * gxy) / determinant;
}

__device__ __forceinline__ void raise_curvature_tensor(
    double gupxx,
    double gupxy,
    double gupxz,
    double gupyy,
    double gupyz,
    double gupzz,
    double axx,
    double axy,
    double axz,
    double ayy,
    double ayz,
    double azz,
    double& aupxx,
    double& aupxy,
    double& aupxz,
    double& aupyx,
    double& aupyy,
    double& aupyz,
    double& aupzx,
    double& aupzy,
    double& aupzz
) {
    aupxx = gupxx * axx + gupxy * axy + gupxz * axz;
    aupxy = gupxx * axy + gupxy * ayy + gupxz * ayz;
    aupxz = gupxx * axz + gupxy * ayz + gupxz * azz;
    aupyx = gupxy * axx + gupyy * axy + gupyz * axz;
    aupyy = gupxy * axy + gupyy * ayy + gupyz * ayz;
    aupyz = gupxy * axz + gupyy * ayz + gupyz * azz;
    aupzx = gupxz * axx + gupyz * axy + gupzz * axz;
    aupzy = gupxz * axy + gupyz * ayy + gupzz * ayz;
    aupzz = gupxz * axz + gupyz * ayz + gupzz * azz;
}

__device__ __forceinline__ void transform_linear_momentum_terms(
    double chi_scale,
    double trk,
    double gxx,
    double gxy,
    double gxz,
    double gyy,
    double gyz,
    double gzz,
    double& axx,
    double& axy,
    double& axz,
    double& ayy,
    double& ayz,
    double& azz
) {
    axx = chi_scale * (axx + gxx * trk / 3.0) - trk;
    axy = chi_scale * (axy + gxy * trk / 3.0);
    axz = chi_scale * (axz + gxz * trk / 3.0);
    ayy = chi_scale * (ayy + gyy * trk / 3.0) - trk;
    ayz = chi_scale * (ayz + gyz * trk / 3.0);
    azz = chi_scale * (azz + gzz * trk / 3.0) - trk;
}

__device__ __forceinline__ int symmetry_image_count(int symmetry_mode) {
    if (symmetry_mode == 0) {
        return 1;
    }
    if (symmetry_mode == 1) {
        return 2;
    }
    return 8;
}

__global__ void project_surface_points_kernel(
    int point_count,
    double radius,
    const double* nx,
    const double* ny,
    const double* nz,
    double* px,
    double* py,
    double* pz
) {
    const int point_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_index >= point_count) {
        return;
    }

    px[point_index] = radius * nx[point_index];
    py[point_index] = radius * ny[point_index];
    pz[point_index] = radius * nz[point_index];
}

__global__ void accumulate_surface_observables_kernel(
    int begin_index,
    int end_index,
    int phi_count,
    int variable_count,
    int symmetry_mode,
    const double* shell_values,
    const double* px,
    const double* py,
    const double* pz,
    const double* nx,
    const double* ny,
    const double* nz,
    const double* theta_weight,
    double* reductions
) {
    const int point_index = begin_index + blockIdx.x * blockDim.x + threadIdx.x;

    double mass = 0.0;
    double linear_x = 0.0;
    double linear_y = 0.0;
    double linear_z = 0.0;
    double angular_x = 0.0;
    double angular_y = 0.0;
    double angular_z = 0.0;

    if (point_index <= end_index) {
        const int theta_index = point_index / phi_count;
        const double weight = theta_weight[theta_index];

        const double chi = shell_values[variable_count * point_index + 3];
        const double trk = shell_values[variable_count * point_index + 4];
        const double gxx = shell_values[variable_count * point_index + 5] + kOne;
        const double gxy = shell_values[variable_count * point_index + 6];
        const double gxz = shell_values[variable_count * point_index + 7];
        const double gyy = shell_values[variable_count * point_index + 8] + kOne;
        const double gyz = shell_values[variable_count * point_index + 9];
        const double gzz = shell_values[variable_count * point_index + 10] + kOne;
        double axx = shell_values[variable_count * point_index + 11];
        double axy = shell_values[variable_count * point_index + 12];
        double axz = shell_values[variable_count * point_index + 13];
        double ayy = shell_values[variable_count * point_index + 14];
        double ayz = shell_values[variable_count * point_index + 15];
        double azz = shell_values[variable_count * point_index + 16];

        const double chi_scale = 1.0 / (1.0 + chi);
        const double psi = chi_scale * sqrt(chi_scale);
        const double normal_x = nx[point_index];
        const double normal_y = ny[point_index];
        const double normal_z = nz[point_index];

        mass = (
            shell_values[variable_count * point_index] * normal_x +
            shell_values[variable_count * point_index + 1] * normal_y +
            shell_values[variable_count * point_index + 2] * normal_z
        ) * weight;

        double gupxx = 0.0;
        double gupxy = 0.0;
        double gupxz = 0.0;
        double gupyy = 0.0;
        double gupyz = 0.0;
        double gupzz = 0.0;
        invert_metric_tensor(gxx, gxy, gxz, gyy, gyz, gzz, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);

        double aupxx = 0.0;
        double aupxy = 0.0;
        double aupxz = 0.0;
        double aupyx = 0.0;
        double aupyy = 0.0;
        double aupyz = 0.0;
        double aupzx = 0.0;
        double aupzy = 0.0;
        double aupzz = 0.0;
        raise_curvature_tensor(gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, axx, axy, axz, ayy, ayz, azz, aupxx, aupxy, aupxz, aupyx, aupyy, aupyz, aupzx, aupzy, aupzz);

        const double coord_x = px[point_index];
        const double coord_y = py[point_index];
        const double coord_z = pz[point_index];

        if (symmetry_mode == 0) {
            angular_x = kOneEighth * psi * (normal_x * (coord_y * aupxz - coord_z * aupxy) + normal_y * (coord_y * aupyz - coord_z * aupyy) + normal_z * (coord_y * aupzz - coord_z * aupzy)) * weight;
            angular_y = kOneEighth * psi * (normal_x * (coord_z * aupxx - coord_x * aupxz) + normal_y * (coord_z * aupyx - coord_x * aupyz) + normal_z * (coord_z * aupzx - coord_x * aupzz)) * weight;
            angular_z = kOneEighth * psi * (normal_x * (coord_x * aupxy - coord_y * aupxx) + normal_y * (coord_x * aupyy - coord_y * aupyx) + normal_z * (coord_x * aupzy - coord_y * aupzx)) * weight;
        } else if (symmetry_mode == 1) {
            angular_z = kOneEighth * psi * (normal_x * (coord_x * aupxy - coord_y * aupxx) + normal_y * (coord_x * aupyy - coord_y * aupyx) + normal_z * (coord_x * aupzy - coord_y * aupzx)) * weight;
        }

        transform_linear_momentum_terms(chi_scale, trk, gxx, gxy, gxz, gyy, gyz, gzz, axx, axy, axz, ayy, ayz, azz);

        if (symmetry_mode == 0) {
            linear_x = kOneEighth * psi * (normal_x * axx + normal_y * axy + normal_z * axz) * weight;
            linear_y = kOneEighth * psi * (normal_x * axy + normal_y * ayy + normal_z * ayz) * weight;
            linear_z = kOneEighth * psi * (normal_x * axz + normal_y * ayz + normal_z * azz) * weight;
        } else if (symmetry_mode == 1) {
            linear_x = kOneEighth * psi * (normal_x * axx + normal_y * axy + normal_z * axz) * weight;
            linear_y = kOneEighth * psi * (normal_x * axy + normal_y * ayy + normal_z * ayz) * weight;
        }
    }

    __shared__ double block_reduction[7];
    if (threadIdx.x == 0) {
        for (int reduction_slot = 0; reduction_slot < 7; ++reduction_slot) {
            block_reduction[reduction_slot] = 0.0;
        }
    }
    __syncthreads();

    atomicAdd(&block_reduction[0], mass);
    atomicAdd(&block_reduction[1], linear_x);
    atomicAdd(&block_reduction[2], linear_y);
    atomicAdd(&block_reduction[3], linear_z);
    atomicAdd(&block_reduction[4], angular_x);
    atomicAdd(&block_reduction[5], angular_y);
    atomicAdd(&block_reduction[6], angular_z);
    __syncthreads();

    if (threadIdx.x == 0) {
        for (int reduction_slot = 0; reduction_slot < 7; ++reduction_slot) {
            atomicAdd(&reductions[reduction_slot], block_reduction[reduction_slot]);
        }
    }
}

__global__ void accumulate_wave_modes_kernel(
    int begin_index,
    int end_index,
    int phi_count,
    int symmetry_mode,
    int spin_weight,
    int max_l,
    int mode_count,
    const double* shell_values,
    const double* theta_weight,
    const double* cos_theta_axis,
    double dphi,
    double r_soa_x,
    double r_soa_y,
    double r_soa_z,
    double i_soa_x,
    double i_soa_y,
    double i_soa_z,
    double* mode_real,
    double* mode_imag
) {
    const int point_index = begin_index + blockIdx.x * blockDim.x + threadIdx.x;
    if (point_index > end_index) {
        return;
    }

    const int theta_index = point_index / phi_count;
    const int phi_index = point_index - theta_index * phi_count;
    const int image_count = symmetry_image_count(symmetry_mode);
    const double weight = theta_weight[theta_index];
    const double costheta_seed = cos_theta_axis[theta_index];
    const double phi_seed = (phi_index + 0.5) * dphi;
    const double psi4_real_seed = shell_values[2 * point_index];
    const double psi4_imag_seed = shell_values[2 * point_index + 1];

    int mode_slot = 0;
    for (int ell = spin_weight; ell <= max_l; ++ell) {
        for (int emm = -ell; emm <= ell; ++emm) {
            double local_real = 0.0;
            double local_imag = 0.0;

            for (int image_id = 0; image_id < image_count; ++image_id) {
                double costheta = costheta_seed;
                double phi = phi_seed;
                double psi4_real = psi4_real_seed;
                double psi4_imag = psi4_imag_seed;

                switch (image_id) {
                    case 0:
                        break;
                    case 1:
                        costheta = -costheta_seed;
                        psi4_real *= r_soa_z;
                        psi4_imag *= i_soa_z;
                        break;
                    case 2:
                        phi = -phi_seed;
                        psi4_real *= r_soa_y;
                        psi4_imag *= i_soa_y;
                        break;
                    case 3:
                        costheta = -costheta_seed;
                        phi = -phi_seed;
                        psi4_real *= r_soa_z * r_soa_y;
                        psi4_imag *= i_soa_z * i_soa_y;
                        break;
                    case 4:
                        phi = M_PI - phi_seed;
                        psi4_real *= r_soa_x;
                        psi4_imag *= i_soa_x;
                        break;
                    case 5:
                        costheta = -costheta_seed;
                        phi = M_PI - phi_seed;
                        psi4_real *= r_soa_z * r_soa_x;
                        psi4_imag *= i_soa_z * i_soa_x;
                        break;
                    case 6:
                        phi = M_PI + phi_seed;
                        psi4_real *= r_soa_y * r_soa_x;
                        psi4_imag *= i_soa_y * i_soa_x;
                        break;
                    default:
                        costheta = -costheta_seed;
                        phi = M_PI + phi_seed;
                        psi4_real *= r_soa_z * r_soa_y * r_soa_x;
                        psi4_imag *= i_soa_z * i_soa_y * i_soa_x;
                        break;
                }

                const double cos_m_phi = cos(emm * phi);
                const double sin_m_phi = sin(emm * phi);
                const double harmonic_factor = sqrt((2.0 * ell + 1.0) / (4.0 * M_PI)) * misc::wigner_d_device(ell, emm, spin_weight, costheta);

                local_real += harmonic_factor * (psi4_real * cos_m_phi + psi4_imag * sin_m_phi) * weight;
                local_imag += harmonic_factor * (psi4_imag * cos_m_phi - psi4_real * sin_m_phi) * weight;
            }

            atomicAdd(&mode_real[mode_slot], local_real);
            atomicAdd(&mode_imag[mode_slot], local_imag);
            ++mode_slot;
        }
    }

    (void)mode_count;
}

struct EmSurfacePoint {
    double px;
    double py;
    double pz;
    double ex;
    double ey;
    double ez;
    double bx;
    double by;
    double bz;
    double chi;
    double gxx;
    double gxy;
    double gxz;
    double gyy;
    double gyz;
    double gzz;
};

__device__ __forceinline__ EmSurfacePoint load_em_surface_point(
    const double* shell_values,
    int variable_count,
    int point_index,
    const double* px,
    const double* py,
    const double* pz
) {
    EmSurfacePoint point;
    point.px = px[point_index];
    point.py = py[point_index];
    point.pz = pz[point_index];
    point.ex = shell_values[variable_count * point_index];
    point.ey = shell_values[variable_count * point_index + 1];
    point.ez = shell_values[variable_count * point_index + 2];
    point.bx = shell_values[variable_count * point_index + 3];
    point.by = shell_values[variable_count * point_index + 4];
    point.bz = shell_values[variable_count * point_index + 5];
    point.chi = shell_values[variable_count * point_index + 6];
    point.gxx = shell_values[variable_count * point_index + 7];
    point.gxy = shell_values[variable_count * point_index + 8];
    point.gxz = shell_values[variable_count * point_index + 9];
    point.gyy = shell_values[variable_count * point_index + 10];
    point.gyz = shell_values[variable_count * point_index + 11];
    point.gzz = shell_values[variable_count * point_index + 12];
    return point;
}

__device__ __forceinline__ void apply_em_surface_image(
    int image_id,
    double costheta_seed,
    double phi_seed,
    double& costheta,
    double& phi,
    EmSurfacePoint& point
) {
    costheta = costheta_seed;
    phi = phi_seed;

    switch (image_id) {
        case 0:
            break;
        case 1:
            costheta = -costheta_seed;
            point.pz = -point.pz;
            point.ez = -point.ez;
            point.bx = -point.bx;
            point.by = -point.by;
            point.gxz = -point.gxz;
            point.gyz = -point.gyz;
            break;
        case 2:
            phi = -phi_seed;
            point.py = -point.py;
            point.ey = -point.ey;
            point.bx = -point.bx;
            point.bz = -point.bz;
            point.gxy = -point.gxy;
            point.gyz = -point.gyz;
            break;
        case 3:
            costheta = -costheta_seed;
            phi = -phi_seed;
            point.py = -point.py;
            point.pz = -point.pz;
            point.ez = -point.ez;
            point.bz = -point.bz;
            point.gxz = -point.gxz;
            point.ey = -point.ey;
            point.by = -point.by;
            point.gxy = -point.gxy;
            break;
        case 4:
            phi = M_PI - phi_seed;
            point.px = -point.px;
            point.ex = -point.ex;
            point.by = -point.by;
            point.bz = -point.bz;
            point.gxy = -point.gxy;
            point.gxz = -point.gxz;
            break;
        case 5:
            costheta = -costheta_seed;
            phi = M_PI - phi_seed;
            point.pz = -point.pz;
            point.px = -point.px;
            point.ez = -point.ez;
            point.bz = -point.bz;
            point.gyz = -point.gyz;
            point.ex = -point.ex;
            point.bx = -point.bx;
            point.gxy = -point.gxy;
            break;
        case 6:
            phi = M_PI + phi_seed;
            point.px = -point.px;
            point.py = -point.py;
            point.ex = -point.ex;
            point.bx = -point.bx;
            point.gxz = -point.gxz;
            point.ey = -point.ey;
            point.by = -point.by;
            point.gyz = -point.gyz;
            break;
        default:
            costheta = -costheta_seed;
            phi = M_PI + phi_seed;
            point.px = -point.px;
            point.py = -point.py;
            point.pz = -point.pz;
            point.ex = -point.ex;
            point.ey = -point.ey;
            point.ez = -point.ez;
            break;
    }
}

__global__ void accumulate_em_wave_modes_kernel(
    int begin_index,
    int end_index,
    int phi_count,
    int variable_count,
    int symmetry_mode,
    int spin_weight,
    int max_l,
    int mode_count,
    const double* shell_values,
    const double* px,
    const double* py,
    const double* pz,
    const double* theta_weight,
    const double* cos_theta_axis,
    double dphi,
    double* mode_real,
    double* mode_imag
) {
    const int point_index = begin_index + blockIdx.x * blockDim.x + threadIdx.x;
    if (point_index > end_index) {
        return;
    }

    const int theta_index = point_index / phi_count;
    const int phi_index = point_index - theta_index * phi_count;
    const int image_count = symmetry_image_count(symmetry_mode);
    const double weight = theta_weight[theta_index];
    const double costheta_seed = cos_theta_axis[theta_index];
    const double phi_seed = (phi_index + 0.5) * dphi;
    const EmSurfacePoint seed = load_em_surface_point(shell_values, variable_count, point_index, px, py, pz);

    int mode_slot = 0;
    for (int ell = spin_weight; ell <= max_l; ++ell) {
        for (int emm = -ell; emm <= ell; ++emm) {
            double local_real = 0.0;
            double local_imag = 0.0;

            for (int image_id = 0; image_id < image_count; ++image_id) {
                EmSurfacePoint point = seed;
                double costheta = costheta_seed;
                double phi = phi_seed;
                apply_em_surface_image(image_id, costheta_seed, phi_seed, costheta, phi, point);

                const amss_gpu_em::NpValue phi_value = (spin_weight == 1)
                    ? amss_gpu_em::compute_phi2_point(
                        point.px, point.py, point.pz, point.chi,
                        point.gxx, point.gxy, point.gxz, point.gyy, point.gyz, point.gzz,
                        point.ex, point.ey, point.ez, point.bx, point.by, point.bz)
                    : amss_gpu_em::compute_phi1_point(
                        point.px, point.py, point.pz, point.chi,
                        point.gxx, point.gxy, point.gxz, point.gyy, point.gyz, point.gzz,
                        point.ex, point.ey, point.ez, point.bx, point.by, point.bz);

                const double cos_m_phi = cos(emm * phi);
                const double sin_m_phi = sin(emm * phi);
                const double harmonic_factor = sqrt((2.0 * ell + 1.0) / (4.0 * M_PI)) * misc::wigner_d_device(ell, emm, spin_weight, costheta);
                const double chipn1 = point.chi + 1.0;
                const double scaled_real = phi_value.real / (chipn1 * chipn1);
                const double scaled_imag = phi_value.imag / (chipn1 * chipn1);

                local_real += harmonic_factor * (scaled_real * cos_m_phi + scaled_imag * sin_m_phi) * weight;
                local_imag += harmonic_factor * (scaled_imag * cos_m_phi - scaled_real * sin_m_phi) * weight;
            }

            atomicAdd(&mode_real[mode_slot], local_real);
            atomicAdd(&mode_imag[mode_slot], local_imag);
            ++mode_slot;
        }
    }

    (void)mode_count;
}

} 

void surface_integral::ensure_wave_buffer_capacity(int NN) {
    const int required_capacity = 2 * NN;
    if (wave_buffer_capacity >= required_capacity)
        return;

    if (d_wave_buffer)
        GpuRuntimeCoordinator::instance().release_device_buffer(d_wave_buffer, wave_buffer_capacity, stream);
    delete[] h_wave_buffer;

    d_wave_buffer = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(required_capacity, stream);
    h_wave_buffer = new double[required_capacity];
    wave_buffer_capacity = required_capacity;
}

void gpu_scale_normals_launch(
    cudaStream_t stream, int n_tot, double rex,
    const double* d_nx, const double* d_ny, const double* d_nz,
    double* d_pox_0, double* d_pox_1, double* d_pox_2
) {
    project_surface_points_kernel<<<linear_grid_size(n_tot), kLinearBlockSize, 0, stream>>>(
        n_tot,
        rex,
        d_nx,
        d_ny,
        d_nz,
        d_pox_0,
        d_pox_1,
        d_pox_2
    );
}

void gpu_surf_MassPAng_launch(
    cudaStream_t stream,
    int Nmin, int Nmax, int N_phi, int InList, int Symmetry,
    const double* d_shellf,
    const double* d_pox_0, const double* d_pox_1, const double* d_pox_2,
    const double* d_nx, const double* d_ny, const double* d_nz,
    const double* d_wtcostheta, double* d_reductions
) {
    const int active_count = Nmax - Nmin + 1;
    if (active_count <= 0) {
        return;
    }

    accumulate_surface_observables_kernel<<<linear_grid_size(active_count), kLinearBlockSize, 0, stream>>>(
        Nmin,
        Nmax,
        N_phi,
        InList,
        Symmetry,
        d_shellf,
        d_pox_0,
        d_pox_1,
        d_pox_2,
        d_nx,
        d_ny,
        d_nz,
        d_wtcostheta,
        d_reductions
    );
}

void gpu_surf_Wave_launch(
    cudaStream_t stream,
    int Nmin, int Nmax, int N_phi, int Symmetry,
    int spinw, int maxl, int NN,
    const double* d_shellf,
    const double* d_wtcostheta, const double* d_arcostheta,
    double dphi,
    double Rpsi4_SoA_0, double Rpsi4_SoA_1, double Rpsi4_SoA_2,
    double Ipsi4_SoA_0, double Ipsi4_SoA_1, double Ipsi4_SoA_2,
    double* d_RP_out, double* d_IP_out
) {
    const int active_count = Nmax - Nmin + 1;
    if (active_count <= 0) {
        return;
    }

    accumulate_wave_modes_kernel<<<linear_grid_size(active_count), kLinearBlockSize, 0, stream>>>(
        Nmin,
        Nmax,
        N_phi,
        Symmetry,
        spinw,
        maxl,
        NN,
        d_shellf,
        d_wtcostheta,
        d_arcostheta,
        dphi,
        Rpsi4_SoA_0,
        Rpsi4_SoA_1,
        Rpsi4_SoA_2,
        Ipsi4_SoA_0,
        Ipsi4_SoA_1,
        Ipsi4_SoA_2,
        d_RP_out,
        d_IP_out
    );
}

void gpu_surf_Wave_EM_launch(
    cudaStream_t stream,
    int Nmin, int Nmax, int N_phi, int InList, int Symmetry,
    int spinw, int maxl, int NN,
    const double* d_shellf,
    const double* d_pox_0, const double* d_pox_1, const double* d_pox_2,
    const double* d_wtcostheta, const double* d_arcostheta,
    double dphi,
    double* d_RP_out, double* d_IP_out
) {
    const int active_count = Nmax - Nmin + 1;
    if (active_count <= 0) {
        return;
    }

    accumulate_em_wave_modes_kernel<<<linear_grid_size(active_count), kLinearBlockSize, 0, stream>>>(
        Nmin,
        Nmax,
        N_phi,
        InList,
        Symmetry,
        spinw,
        maxl,
        NN,
        d_shellf,
        d_pox_0,
        d_pox_1,
        d_pox_2,
        d_wtcostheta,
        d_arcostheta,
        dphi,
        d_RP_out,
        d_IP_out
    );
}

void surface_integral::gpu_surf_MassPAng(
    double rex, int lev, cgh *GH, var *chi, var *trK,
    var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
    var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
    var *Gmx, var *Gmy, var *Gmz,
    var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs,
    double *Rout, monitor *Monitor
) {
    if (myrank == 0 && GH->grids[lev] != 1) {
        if (Monitor && Monitor->outfile) {
            Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
        } else {
            cout << "WARNING: surface integral on multipatches" << endl;
        }
    }

    const int InList = 17;
    MyList<var> *DG_List = new MyList<var>(Sfx_rhs);
    DG_List->insert(Sfy_rhs);
    DG_List->insert(Sfz_rhs);
    DG_List->insert(chi);
    DG_List->insert(trK);
    DG_List->insert(gxx);
    DG_List->insert(gxy);
    DG_List->insert(gxz);
    DG_List->insert(gyy);
    DG_List->insert(gyz);
    DG_List->insert(gzz);
    DG_List->insert(Axx);
    DG_List->insert(Axy);
    DG_List->insert(Axz);
    DG_List->insert(Ayy);
    DG_List->insert(Ayz);
    DG_List->insert(Azz);

    for (MyList<Patch> *patch_node = GH->PatL[lev]; patch_node; patch_node = patch_node->next) {
        for (MyList<Block> *block_node = patch_node->data->blb; block_node; block_node = block_node->next) {
            Block *block = block_node->data;
            if (myrank == block->rank) {
                gpu_admmass_bssn_launch(
                    block->stream,
                    block->shape, block->d_X[0], block->d_X[1], block->d_X[2],
                    block->d_fgfs[chi->sgfn], block->d_fgfs[trK->sgfn],
                    block->d_fgfs[gxx->sgfn], block->d_fgfs[gxy->sgfn], block->d_fgfs[gxz->sgfn], block->d_fgfs[gyy->sgfn], block->d_fgfs[gyz->sgfn], block->d_fgfs[gzz->sgfn],
                    block->d_fgfs[Axx->sgfn], block->d_fgfs[Axy->sgfn], block->d_fgfs[Axz->sgfn], block->d_fgfs[Ayy->sgfn], block->d_fgfs[Ayz->sgfn], block->d_fgfs[Azz->sgfn],
                    block->d_fgfs[Gmx->sgfn], block->d_fgfs[Gmy->sgfn], block->d_fgfs[Gmz->sgfn],
                    block->d_fgfs[Sfx_rhs->sgfn], block->d_fgfs[Sfy_rhs->sgfn], block->d_fgfs[Sfz_rhs->sgfn],
                    Symmetry
                );
            }
            if (block_node == patch_node->data->ble) {
                break;
            }
        }
    }
    synchronize_level_streams(GH, lev, myrank);

    DevicePointTriplet surface_points{{d_surface_points[0], d_surface_points[1], d_surface_points[2]}};
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, surface_points.values[0], surface_points.values[1], surface_points.values[2]);

    double* d_shellf = d_shellf_buffer;
    GH->PatL[lev]->data->Interp_Points_GPU(this->stream, DG_List, n_tot, surface_points.values, d_shellf, Symmetry);

    const RankWorkRange work_range = split_rank_range(n_tot, cpusize, myrank);
    double* h_reductions = h_reductions_buffer;
    double* d_reductions = d_reductions_buffer;
    CUDA_CHECK(cudaMemsetAsync(d_reductions, 0, 7 * sizeof(double), this->stream));
    gpu_surf_MassPAng_launch(
        this->stream,
        work_range.begin,
        work_range.end,
        N_phi,
        InList,
        Symmetry,
        d_shellf,
        surface_points.values[0],
        surface_points.values[1],
        surface_points.values[2],
        d_nx_g,
        d_ny_g,
        d_nz_g,
        d_wtcostheta,
        d_reductions
    );
    GpuRuntimeCoordinator::instance().copy_to_host(h_reductions, d_reductions, 7, this->stream);
    GpuRuntimeCoordinator::wait_for_stream(this->stream);

    if (cpusize > 1) {
        MPI_Allreduce(MPI_IN_PLACE, h_reductions, 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    double mass = h_reductions[0];
    double px = h_reductions[1];
    double py = h_reductions[2];
    double pz = h_reductions[3];
    double sx = h_reductions[4];
    double sy = h_reductions[5];
    double sz = h_reductions[6];

    mass = mass * rex * rex * dphi * factor;
    px = px * rex * rex * dphi * (1.0 / PI) * factor;
    py = py * rex * rex * dphi * (1.0 / PI) * factor;
    pz = pz * rex * rex * dphi * (1.0 / PI) * factor;
    sx = sx * rex * rex * dphi * (1.0 / PI) * factor;
    sy = sy * rex * rex * dphi * (1.0 / PI) * factor;
    sz = sz * rex * rex * dphi * (1.0 / PI) * factor;

    Rout[0] = mass;
    Rout[1] = px;
    Rout[2] = py;
    Rout[3] = pz;
    Rout[4] = sx;
    Rout[5] = sy;
    Rout[6] = sz;

    DG_List->clearList();
    GpuRuntimeCoordinator::wait_for_stream(this->stream);
}

void surface_integral::gpu_surf_Wave(
    double rex, int lev, cgh *GH, var *Rpsi4, var *Ipsi4,
    int spinw, int maxl, int NN, double *RP, double *IP,
    monitor *Monitor
) {
    if (myrank == 0 && GH->grids[lev] != 1) {
        if (Monitor && Monitor->outfile) {
            Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
        } else {
            cout << "WARNING: surface integral on multipatches" << endl;
        }
    }

    const int InList = 2;
    MyList<var> *DG_List = new MyList<var>(Rpsi4);
    DG_List->insert(Ipsi4);

    DevicePointTriplet surface_points{{d_surface_points[0], d_surface_points[1], d_surface_points[2]}};
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, surface_points.values[0], surface_points.values[1], surface_points.values[2]);

    double* d_shellf = d_shellf_buffer;
    GH->PatL[lev]->data->Interp_Points_GPU(this->stream, DG_List, n_tot, surface_points.values, d_shellf, Symmetry);

    const RankWorkRange work_range = split_rank_range(n_tot, cpusize, myrank);
    ensure_wave_buffer_capacity(NN);
    double *d_wave_out = d_wave_buffer;
    double *d_RP_out = d_wave_out;
    double *d_IP_out = d_wave_out + NN;
    CUDA_CHECK(cudaMemsetAsync(d_wave_out, 0, 2 * NN * sizeof(double), this->stream));

    gpu_surf_Wave_launch(
        this->stream,
        work_range.begin,
        work_range.end,
        N_phi,
        Symmetry,
        spinw,
        maxl,
        NN,
        d_shellf,
        d_wtcostheta,
        d_arcostheta,
        dphi,
        Rpsi4->SoA[0],
        Rpsi4->SoA[1],
        Rpsi4->SoA[2],
        Ipsi4->SoA[0],
        Ipsi4->SoA[1],
        Ipsi4->SoA[2],
        d_RP_out,
        d_IP_out
    );

    double *h_wave_out = h_wave_buffer;
    GpuRuntimeCoordinator::instance().copy_to_host(h_wave_out, d_wave_out, 2 * NN, this->stream);
    GpuRuntimeCoordinator::wait_for_stream(this->stream);

    for (int mode_index = 0; mode_index < NN; ++mode_index) {
        h_wave_out[mode_index] *= rex * dphi;
        h_wave_out[NN + mode_index] *= rex * dphi;
    }

    if (cpusize > 1) {
        MPI_Allreduce(MPI_IN_PLACE, h_wave_out, 2 * NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    for (int mode_index = 0; mode_index < NN; ++mode_index) {
        RP[mode_index] = h_wave_out[mode_index];
        IP[mode_index] = h_wave_out[NN + mode_index];
    }

    DG_List->clearList();
    GpuRuntimeCoordinator::wait_for_stream(this->stream);
}

void surface_integral::gpu_surf_Wave_EM(
    double rex, int lev, cgh *GH,
    var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
    var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
    int spinw, int maxl, int NN, double *RP, double *IP,
    monitor *Monitor
) {
    if (myrank == 0 && GH->grids[lev] != 1) {
        if (Monitor && Monitor->outfile) {
            Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
        } else {
            cout << "WARNING: surface integral on multipatches" << endl;
        }
    }

    const int InList = 13;
    MyList<var> *DG_List = new MyList<var>(Ex);
    DG_List->insert(Ey);
    DG_List->insert(Ez);
    DG_List->insert(Bx);
    DG_List->insert(By);
    DG_List->insert(Bz);
    DG_List->insert(chi);
    DG_List->insert(gxx);
    DG_List->insert(gxy);
    DG_List->insert(gxz);
    DG_List->insert(gyy);
    DG_List->insert(gyz);
    DG_List->insert(gzz);

    for (MyList<Patch> *patch_node = GH->PatL[lev]; patch_node; patch_node = patch_node->next) {
        for (MyList<Block> *block_node = patch_node->data->blb; block_node; block_node = block_node->next) {
            Block *block = block_node->data;
            if (myrank == block->rank) {
                for (MyList<var> *var_node = DG_List; var_node; var_node = var_node->next) {
                    block->require_on_gpu(var_node->data->sgfn);
                }
            }
            if (block_node == patch_node->data->ble) {
                break;
            }
        }
    }

    DevicePointTriplet surface_points{{d_surface_points[0], d_surface_points[1], d_surface_points[2]}};
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, surface_points.values[0], surface_points.values[1], surface_points.values[2]);

    double* d_shellf = d_shellf_buffer;
    GH->PatL[lev]->data->Interp_Points_GPU(this->stream, DG_List, n_tot, surface_points.values, d_shellf, Symmetry);

    const RankWorkRange work_range = split_rank_range(n_tot, cpusize, myrank);
    ensure_wave_buffer_capacity(NN);
    double *d_wave_out = d_wave_buffer;
    double *d_RP_out = d_wave_out;
    double *d_IP_out = d_wave_out + NN;
    CUDA_CHECK(cudaMemsetAsync(d_wave_out, 0, 2 * NN * sizeof(double), this->stream));

    gpu_surf_Wave_EM_launch(
        this->stream,
        work_range.begin,
        work_range.end,
        N_phi,
        InList,
        Symmetry,
        spinw,
        maxl,
        NN,
        d_shellf,
        surface_points.values[0],
        surface_points.values[1],
        surface_points.values[2],
        d_wtcostheta,
        d_arcostheta,
        dphi,
        d_RP_out,
        d_IP_out
    );

    double *h_wave_out = h_wave_buffer;
    GpuRuntimeCoordinator::instance().copy_to_host(h_wave_out, d_wave_out, 2 * NN, this->stream);
    GpuRuntimeCoordinator::wait_for_stream(this->stream);

    for (int mode_index = 0; mode_index < NN; ++mode_index) {
        h_wave_out[mode_index] *= rex * dphi;
        h_wave_out[NN + mode_index] *= rex * dphi;
    }

    if (cpusize > 1) {
        MPI_Allreduce(MPI_IN_PLACE, h_wave_out, 2 * NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    for (int mode_index = 0; mode_index < NN; ++mode_index) {
        RP[mode_index] = h_wave_out[mode_index];
        IP[mode_index] = h_wave_out[NN + mode_index];
    }

    DG_List->clearList();
    GpuRuntimeCoordinator::wait_for_stream(this->stream);
}
