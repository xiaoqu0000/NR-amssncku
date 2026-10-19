#include "sommerfeld_rout.h"

#include <cuda_runtime.h>
#include <math.h>

#include "gpu_runtime_coordinator.h"
#include "fmisc.h"
#include "fmisc_gpu.cuh"

constexpr int kSommerfeldStencilWidth = 6;
constexpr int kNoSymmetry = 0;
constexpr int kEquatorialSymmetry = 1;
constexpr int kOctantSymmetry = 2;
constexpr int kCorrectorStage = 1;

constexpr double kZeroValue = 0.0;
constexpr double kUnitValue = 1.0;
constexpr double kTwoValue = 2.0;

__global__ void apply_sommerfeld_transport_kernel(
    int ex0, int ex1, int ex2,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double* f_rhs,
    const double* f0,
    double velocity,
    const double SYM1, const double SYM2, const double SYM3,
    int symmetry_mode
);

__global__ void apply_sommerfeld_transport_ss_kernel(
    int ex0, int ex1, int ex2,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double* f_rhs,
    const double* f0,
    double velocity,
    const double SYM1, const double SYM2, const double SYM3,
    int symmetry_mode
);

__global__ void apply_sommerfeld_projection_kernel(
    int ex0, int ex1, int ex2,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double dT,
    const double* chi0, const double* Lap0,
    const double* f0, double* f,
    const double SYM1, const double SYM2, const double SYM3,
    int symmetry_mode,
    int predictor_or_corrector
);

inline dim3 sommerfeld_block_shape() {
    return dim3(8, 8, 4);
}

inline dim3 sommerfeld_grid_shape(const int ex[3]) {
    const dim3 block = sommerfeld_block_shape();
    return dim3(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );
}

__device__ __forceinline__ double axis_spacing(const double* axis) {
    return axis[1] - axis[0];
}

__device__ __forceinline__ bool lower_face_is_physical(
    int axis_id,
    double lower_bound,
    double spacing,
    int symmetry_mode
) {
    if (axis_id == 0) {
        return !(symmetry_mode == kOctantSymmetry && fabs(lower_bound) < spacing / kTwoValue);
    }
    if (axis_id == 1) {
        return !(symmetry_mode == kOctantSymmetry && fabs(lower_bound) < spacing / kTwoValue);
    }
    return !(symmetry_mode > kNoSymmetry && fabs(lower_bound) < spacing / kTwoValue);
}

__device__ __forceinline__ bool touches_sommerfeld_surface(
    int i1, int j1, int k1,
    int ex0, int ex1, int ex2,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    int symmetry_mode
) {
    const double spacing_x = axis_spacing(X);
    const double spacing_y = axis_spacing(Y);
    const double spacing_z = axis_spacing(Z);

    const bool on_upper_x = (i1 == ex0 && fabs(X[ex0 - 1] - xmax) < spacing_x);
    const bool on_upper_y = (j1 == ex1 && fabs(Y[ex1 - 1] - ymax) < spacing_y);
    const bool on_upper_z = (k1 == ex2 && fabs(Z[ex2 - 1] - zmax) < spacing_z);

    const bool on_lower_x = (i1 == 1 && fabs(X[0] - xmin) < spacing_x &&
                             lower_face_is_physical(0, xmin, spacing_x, symmetry_mode));
    const bool on_lower_y = (j1 == 1 && fabs(Y[0] - ymin) < spacing_y &&
                             lower_face_is_physical(1, ymin, spacing_y, symmetry_mode));
    const bool on_lower_z = (k1 == 1 && fabs(Z[0] - zmin) < spacing_z &&
                             lower_face_is_physical(2, zmin, spacing_z, symmetry_mode));

    return on_upper_x || on_upper_y || on_upper_z || on_lower_x || on_lower_y || on_lower_z;
}

__device__ __forceinline__ int flat_cell_index(int i0, int j0, int k0, int ex0, int ex1) {
    return k0 * ex1 * ex0 + j0 * ex0 + i0;
}

__device__ __forceinline__ int stencil_window_begin(double offset, int center_1b) {
    return (offset > kZeroValue)
        ? center_1b - static_cast<int>(trunc(offset)) - kSommerfeldStencilWidth / 2
        : center_1b - static_cast<int>(trunc(offset)) - kSommerfeldStencilWidth / 2 + 1;
}

__device__ __forceinline__ double shifted_fraction(double offset) {
    return (offset > kZeroValue)
        ? trunc(offset) - offset + kSommerfeldStencilWidth / 2
        : trunc(offset) - offset + kSommerfeldStencilWidth / 2 - 1;
}

__device__ __forceinline__ void pull_window_inside_domain(
    int lower_corner[3],
    int upper_corner[3],
    double fractional_offset[3],
    const int extent[3],
    int symmetry_mode
) {
    if (symmetry_mode == kNoSymmetry && lower_corner[2] < 1) {
        fractional_offset[2] += (lower_corner[2] - 1);
        upper_corner[2] -= (lower_corner[2] - 1);
        lower_corner[2] = 1;
    }
    if (symmetry_mode < kOctantSymmetry && lower_corner[1] < 1) {
        fractional_offset[1] += (lower_corner[1] - 1);
        upper_corner[1] -= (lower_corner[1] - 1);
        lower_corner[1] = 1;
    }
    if (symmetry_mode < kOctantSymmetry && lower_corner[0] < 1) {
        fractional_offset[0] += (lower_corner[0] - 1);
        upper_corner[0] -= (lower_corner[0] - 1);
        lower_corner[0] = 1;
    }

    for (int axis = 0; axis < 3; ++axis) {
        if (upper_corner[axis] > extent[axis]) {
            fractional_offset[axis] += (upper_corner[axis] - extent[axis]);
            lower_corner[axis] -= (upper_corner[axis] - extent[axis]);
            upper_corner[axis] = extent[axis];
        }
    }
}

__device__ __forceinline__ double symmetry_sample(
    const double* field,
    const int extent[3],
    int i1, int j1, int k1,
    const double symmetry_axes[3]
) {
    return d_symmetry_bd_1b(1, const_cast<int*>(extent), const_cast<double*>(field), i1, j1, k1, symmetry_axes);
}

__device__ __forceinline__ double offset_sample(
    int axis,
    int i1, int j1, int k1,
    int delta,
    const double* field,
    const int extent[3],
    const double symmetry_axes[3]
) {
    const int ii = i1 + (axis == 0 ? delta : 0);
    const int jj = j1 + (axis == 1 ? delta : 0);
    const int kk = k1 + (axis == 2 ? delta : 0);
    return symmetry_sample(field, extent, ii, jj, kk, symmetry_axes);
}

__device__ __forceinline__ double directional_boundary_gradient(
    int axis,
    int i1, int j1, int k1,
    int lower_limit,
    int upper_limit,
    double wave_speed_component,
    double inverse_two_spacing,
    const double* field,
    const int extent[3],
    const double symmetry_axes[3]
) {
    if (wave_speed_component == kZeroValue) {
        return kZeroValue;
    }

    const int center_index = (axis == 0) ? i1 : (axis == 1 ? j1 : k1);
    const int plus_one = center_index + 1;
    const int plus_two = center_index + 2;
    const int minus_one = center_index - 1;
    const int minus_two = center_index - 2;

    if (wave_speed_component > kZeroValue) {
        if (minus_two >= lower_limit) {
            return inverse_two_spacing * (
                3.0 * offset_sample(axis, i1, j1, k1, 0, field, extent, symmetry_axes) -
                4.0 * offset_sample(axis, i1, j1, k1, -1, field, extent, symmetry_axes) +
                offset_sample(axis, i1, j1, k1, -2, field, extent, symmetry_axes)
            );
        }
        if (minus_one >= lower_limit) {
            return inverse_two_spacing * (
                -offset_sample(axis, i1, j1, k1, -1, field, extent, symmetry_axes) +
                offset_sample(axis, i1, j1, k1, 1, field, extent, symmetry_axes)
            );
        }
        return inverse_two_spacing * (
            -offset_sample(axis, i1, j1, k1, 2, field, extent, symmetry_axes) +
            4.0 * offset_sample(axis, i1, j1, k1, 1, field, extent, symmetry_axes) -
            3.0 * offset_sample(axis, i1, j1, k1, 0, field, extent, symmetry_axes)
        );
    }

    if (plus_two <= upper_limit) {
        return inverse_two_spacing * (
            -offset_sample(axis, i1, j1, k1, 2, field, extent, symmetry_axes) +
            4.0 * offset_sample(axis, i1, j1, k1, 1, field, extent, symmetry_axes) -
            3.0 * offset_sample(axis, i1, j1, k1, 0, field, extent, symmetry_axes)
        );
    }
    if (plus_one <= upper_limit) {
        return inverse_two_spacing * (
            -offset_sample(axis, i1, j1, k1, -1, field, extent, symmetry_axes) +
            offset_sample(axis, i1, j1, k1, 1, field, extent, symmetry_axes)
        );
    }
    return inverse_two_spacing * (
        3.0 * offset_sample(axis, i1, j1, k1, 0, field, extent, symmetry_axes) -
        4.0 * offset_sample(axis, i1, j1, k1, -1, field, extent, symmetry_axes) +
        offset_sample(axis, i1, j1, k1, -2, field, extent, symmetry_axes)
    );
}

void gpu_sommerfeld_rout_launch(
    cudaStream_t &stream,
    int ex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double dT, const double* d_chi0, const double* d_Lap0,
    const double* d_f0, double* d_f, const double SoA[3],
    int Symmetry, int precor
) {
    const dim3 block = sommerfeld_block_shape();
    const dim3 grid = sommerfeld_grid_shape(ex);

    apply_sommerfeld_projection_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], d_X, d_Y, d_Z, xmin, ymin, zmin, xmax, ymax, zmax,
        dT, d_chi0, d_Lap0, d_f0, d_f, SoA[0], SoA[1], SoA[2], Symmetry, precor
    );
}

void gpu_sommerfeld_routbam_launch(
    cudaStream_t &stream,
    int ex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double* d_f_rhs, const double* d_f0,
    double velocity, const double SoA[3], int Symmetry
) {
    const dim3 block = sommerfeld_block_shape();
    const dim3 grid = sommerfeld_grid_shape(ex);

    apply_sommerfeld_transport_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], d_X, d_Y, d_Z, xmin, ymin, zmin, xmax, ymax, zmax,
        d_f_rhs, d_f0, velocity, SoA[0], SoA[1], SoA[2], Symmetry
    );
}

void gpu_sommerfeld_routbam_ss_launch(
    cudaStream_t &stream,
    int ex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double* d_f_rhs, const double* d_f0,
    double velocity, const double SoA[3], int Symmetry
) {
    const dim3 block = sommerfeld_block_shape();
    const dim3 grid = sommerfeld_grid_shape(ex);

    apply_sommerfeld_transport_ss_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], d_X, d_Y, d_Z, xmin, ymin, zmin, xmax, ymax, zmax,
        d_f_rhs, d_f0, velocity, SoA[0], SoA[1], SoA[2], Symmetry
    );
}

__global__ void apply_sommerfeld_projection_kernel(
    int ex0, int ex1, int ex2,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double dT,
    const double* chi0, const double* Lap0,
    const double* f0, double* f,
    const double SYM1, const double SYM2, const double SYM3,
    int symmetry_mode,
    int predictor_or_corrector
) {
    const int i0 = blockIdx.x * blockDim.x + threadIdx.x;
    const int j0 = blockIdx.y * blockDim.y + threadIdx.y;
    const int k0 = blockIdx.z * blockDim.z + threadIdx.z;

    if (i0 >= ex0 || j0 >= ex1 || k0 >= ex2) {
        return;
    }

    const int i1 = i0 + 1;
    const int j1 = j0 + 1;
    const int k1 = k0 + 1;

    if (!touches_sommerfeld_surface(
            i1, j1, k1, ex0, ex1, ex2, X, Y, Z,
            xmin, ymin, zmin, xmax, ymax, zmax, symmetry_mode)) {
        return;
    }

    const int linear_index = flat_cell_index(i0, j0, k0, ex0, ex1);
    if (predictor_or_corrector == kCorrectorStage) {
        f[linear_index] = f0[linear_index];
        return;
    }

    const int extent[3] = {ex0, ex1, ex2};
    const double symmetry_axes[3] = {SYM1, SYM2, SYM3};
    const double spacing[3] = {axis_spacing(X), axis_spacing(Y), axis_spacing(Z)};
    const double radial_distance = sqrt(X[i0] * X[i0] + Y[j0] * Y[j0] + Z[k0] * Z[k0]);
    const double transport_scale =
        (Lap0[linear_index] + kUnitValue) * sqrt(kUnitValue + chi0[linear_index]) * dT / radial_distance;
    const double retained_weight = kUnitValue - transport_scale;

    double fractional_offset[3] = {
        transport_scale * X[i0] / spacing[0],
        transport_scale * Y[j0] / spacing[1],
        transport_scale * Z[k0] / spacing[2]
    };

    int lower_corner[3] = {
        stencil_window_begin(fractional_offset[0], i1),
        stencil_window_begin(fractional_offset[1], j1),
        stencil_window_begin(fractional_offset[2], k1)
    };

    for (int axis = 0; axis < 3; ++axis) {
        fractional_offset[axis] = shifted_fraction(fractional_offset[axis]);
    }

    int upper_corner[3];
    for (int axis = 0; axis < 3; ++axis) {
        upper_corner[axis] = lower_corner[axis] + kSommerfeldStencilWidth - 1;
    }

    pull_window_inside_domain(lower_corner, upper_corner, fractional_offset, extent, symmetry_mode);

    double neighborhood[kSommerfeldStencilWidth * kSommerfeldStencilWidth * kSommerfeldStencilWidth];
    d_decide3d(
        const_cast<int*>(extent),
        const_cast<double*>(f0),
        const_cast<double*>(f0),
        lower_corner,
        upper_corner,
        symmetry_axes,
        neighborhood,
        kSommerfeldStencilWidth,
        symmetry_mode
    );

    double interpolation_axis[kSommerfeldStencilWidth];
    for (int idx = 0; idx < kSommerfeldStencilWidth; ++idx) {
        interpolation_axis[idx] = static_cast<double>(idx);
    }

    double ignored_error = 0.0;
    double interpolated_value = 0.0;
    d_polin3_1b(
        interpolation_axis,
        interpolation_axis,
        interpolation_axis,
        neighborhood,
        fractional_offset[0],
        fractional_offset[1],
        fractional_offset[2],
        interpolated_value,
        ignored_error,
        kSommerfeldStencilWidth
    );

    f[linear_index] = interpolated_value * retained_weight;
}

__global__ void apply_sommerfeld_transport_kernel(
    int ex0, int ex1, int ex2,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double* f_rhs,
    const double* f0,
    double velocity,
    const double SYM1, const double SYM2, const double SYM3,
    int symmetry_mode
) {
    const int i0 = blockIdx.x * blockDim.x + threadIdx.x;
    const int j0 = blockIdx.y * blockDim.y + threadIdx.y;
    const int k0 = blockIdx.z * blockDim.z + threadIdx.z;

    if (i0 >= ex0 || j0 >= ex1 || k0 >= ex2) {
        return;
    }

    const int i1 = i0 + 1;
    const int j1 = j0 + 1;
    const int k1 = k0 + 1;

    if (!touches_sommerfeld_surface(
            i1, j1, k1, ex0, ex1, ex2, X, Y, Z,
            xmin, ymin, zmin, xmax, ymax, zmax, symmetry_mode)) {
        return;
    }

    const int linear_index = flat_cell_index(i0, j0, k0, ex0, ex1);
    const int extent[3] = {ex0, ex1, ex2};
    const double symmetry_axes[3] = {SYM1, SYM2, SYM3};

    const double spacing_x = axis_spacing(X);
    const double spacing_y = axis_spacing(Y);
    const double spacing_z = axis_spacing(Z);

    int min_i = 1;
    int min_j = 1;
    int min_k = 1;
    if (symmetry_mode > kNoSymmetry && fabs(Z[0]) < spacing_z) {
        min_k = 0;
    }
    if (symmetry_mode > kEquatorialSymmetry && fabs(X[0]) < spacing_x) {
        min_i = 0;
    }
    if (symmetry_mode > kEquatorialSymmetry && fabs(Y[0]) < spacing_y) {
        min_j = 0;
    }

    const int max_i = ex0;
    const int max_j = ex1;
    const int max_k = ex2;

    const double radius = sqrt(X[i0] * X[i0] + Y[j0] * Y[j0] + Z[k0] * Z[k0]);
    const double velocity_x = velocity * X[i0] / radius;
    const double velocity_y = velocity * Y[j0] / radius;
    const double velocity_z = velocity * Z[k0] / radius;

    const double grad_x = directional_boundary_gradient(
        0, i1, j1, k1, min_i, max_i, velocity_x, kUnitValue / (kTwoValue * spacing_x), f0, extent, symmetry_axes
    );
    const double grad_y = directional_boundary_gradient(
        1, i1, j1, k1, min_j, max_j, velocity_y, kUnitValue / (kTwoValue * spacing_y), f0, extent, symmetry_axes
    );
    const double grad_z = directional_boundary_gradient(
        2, i1, j1, k1, min_k, max_k, velocity_z, kUnitValue / (kTwoValue * spacing_z), f0, extent, symmetry_axes
    );

    f_rhs[linear_index] =
        -velocity * (grad_x * X[i0] + grad_y * Y[j0] + grad_z * Z[k0] + f0[linear_index]) / radius;
}

__global__ void apply_sommerfeld_transport_ss_kernel(
    int ex0, int ex1, int ex2,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    double* f_rhs,
    const double* f0,
    double velocity,
    const double SYM1, const double SYM2, const double SYM3,
    int symmetry_mode
) {
    (void)X;
    (void)Y;
    (void)xmin;
    (void)ymin;
    (void)xmax;
    (void)ymax;
    (void)SYM1;
    (void)SYM2;
    (void)SYM3;
    (void)symmetry_mode;

    const int i0 = blockIdx.x * blockDim.x + threadIdx.x;
    const int j0 = blockIdx.y * blockDim.y + threadIdx.y;
    const int k0 = blockIdx.z * blockDim.z + threadIdx.z;

    if (i0 >= ex0 || j0 >= ex1 || k0 >= ex2) {
        return;
    }

    const double dZ = axis_spacing(Z);
    const bool on_outer_shell = fabs(Z[ex2 - 1] - zmax) < dZ && k0 == ex2 - 1;
    const bool on_inner_shell = fabs(Z[0] - zmin) < dZ && k0 == 0;
    if (!on_outer_shell && !on_inner_shell) {
        return;
    }

    const int linear_index = flat_cell_index(i0, j0, k0, ex0, ex1);
    if (on_inner_shell) {
        f_rhs[linear_index] = kZeroValue;
        return;
    }

    if (velocity == kZeroValue) {
        return;
    }

    const double d2dz = kUnitValue / (kTwoValue * dZ);
    double fz = kZeroValue;

    if (velocity > kZeroValue) {
        if (k0 - 2 >= 0) {
            fz = d2dz * (
                3.0 * f0[linear_index] -
                4.0 * f0[flat_cell_index(i0, j0, k0 - 1, ex0, ex1)] +
                f0[flat_cell_index(i0, j0, k0 - 2, ex0, ex1)]
            );
        } else if (k0 - 1 >= 0 && k0 + 1 < ex2) {
            fz = d2dz * (
                -f0[flat_cell_index(i0, j0, k0 - 1, ex0, ex1)] +
                f0[flat_cell_index(i0, j0, k0 + 1, ex0, ex1)]
            );
        } else if (k0 + 2 < ex2) {
            fz = d2dz * (
                -f0[flat_cell_index(i0, j0, k0 + 2, ex0, ex1)] +
                4.0 * f0[flat_cell_index(i0, j0, k0 + 1, ex0, ex1)] -
                3.0 * f0[linear_index]
            );
        }
    } else {
        if (k0 + 2 < ex2) {
            fz = d2dz * (
                -f0[flat_cell_index(i0, j0, k0 + 2, ex0, ex1)] +
                4.0 * f0[flat_cell_index(i0, j0, k0 + 1, ex0, ex1)] -
                3.0 * f0[linear_index]
            );
        } else if (k0 - 1 >= 0 && k0 + 1 < ex2) {
            fz = d2dz * (
                -f0[flat_cell_index(i0, j0, k0 - 1, ex0, ex1)] +
                f0[flat_cell_index(i0, j0, k0 + 1, ex0, ex1)]
            );
        } else if (k0 - 2 >= 0) {
            fz = d2dz * (
                3.0 * f0[linear_index] -
                4.0 * f0[flat_cell_index(i0, j0, k0 - 1, ex0, ex1)] +
                f0[flat_cell_index(i0, j0, k0 - 2, ex0, ex1)]
            );
        }
    }

    f_rhs[linear_index] = -velocity * (fz + f0[linear_index] / Z[k0]);
}
