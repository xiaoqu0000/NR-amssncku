#include "shell_patch_gpu.h"

#include <cuda_runtime.h>

#include <math.h>

#include "gpu_runtime_coordinator.h"

namespace {

constexpr int kShellTransferBlockSize = 128;

__device__ __forceinline__ int flatten_index_0b(int i, int j, int k, int nx, int ny)
{
    return (k * ny + j) * nx + i;
}

__device__ __forceinline__ int reflect_1b(int index, int extent, double soa, double& factor)
{
    if (index >= 1 && index <= extent)
        return index;

    factor *= soa;
#ifdef Vertex
    if (index < 1)
        index = 2 - index;
    else
        index = 2 * extent - index;
#else
    if (index < 1)
        index = 1 - index;
    else
        index = 2 * extent + 1 - index;
#endif
    if (index < 1 || index > extent)
        return -1;
    return index;
}

__device__ __forceinline__ double fetch_tbd(
    const double* field,
    const int shape[dim],
    int i1b,
    int j1b,
    int k1b,
    const double soa[dim])
{
    double factor = 1.0;
    i1b = reflect_1b(i1b, shape[0], soa[0], factor);
    j1b = reflect_1b(j1b, shape[1], soa[1], factor);
    k1b = reflect_1b(k1b, shape[2], soa[2], factor);
    if (i1b < 1 || j1b < 1 || k1b < 1)
        return NAN;
    return factor * field[flatten_index_0b(i1b - 1, j1b - 1, k1b - 1, shape[0], shape[1])];
}

__device__ __forceinline__ double fetch_sntbd(
    const double* field,
    const int shape[dim],
    int i1b,
    int j1b,
    int k1b,
    int active_axis,
    double soa)
{
    double factor = 1.0;
    if (active_axis == 0)
        i1b = reflect_1b(i1b, shape[0], soa, factor);
    else
        j1b = reflect_1b(j1b, shape[1], soa, factor);
    if (i1b < 1 || i1b > shape[0] || j1b < 1 || j1b > shape[1] || k1b < 1 || k1b > shape[2])
        return NAN;
    return factor * field[flatten_index_0b(i1b - 1, j1b - 1, k1b - 1, shape[0], shape[1])];
}

__device__ __forceinline__ void shell_3d_soa(int sst, double in0, double in1, double in2, double out[dim])
{
    if (sst == -1) {
        out[0] = in0;
        out[1] = in1;
        out[2] = in2;
    } else if (sst == 0 || sst == 1) {
        out[0] = in0;
        out[1] = in1;
        out[2] = 0.0;
    } else if (sst == 2 || sst == 3) {
        out[0] = in1;
        out[1] = in2;
        out[2] = 0.0;
    } else {
        out[0] = in0;
        out[1] = in2;
        out[2] = 0.0;
    }
}

__device__ __forceinline__ double shell_1d_soa(int sst, int dumyd, double in0, double in1, double in2)
{
    if (dumyd == 1) {
        if (sst == 2 || sst == 3)
            return in1;
        return in0;
    }
    if (sst == 0 || sst == 1)
        return in1;
    return in2;
}

__device__ __forceinline__ double interpolate_point_3d(
    const double* field,
    const ShellPatchDevicePoint& point,
    const double* coef,
    double SoA_0,
    double SoA_1,
    double SoA_2)
{
    const int ordn = 2 * ghost_width;
    double soa[dim];
    shell_3d_soa(point.ssst, SoA_0, SoA_1, SoA_2, soa);

    double value = 0.0;
    for (int k = 0; k < ordn; ++k) {
        const int k1b = point.sind[2] + k + 1;
        for (int j = 0; j < ordn; ++j) {
            const int j1b = point.sind[1] + j + 1;
            for (int i = 0; i < ordn; ++i) {
                const int i1b = point.sind[0] + i + 1;
                const double sample = fetch_tbd(field, point.shape, i1b, j1b, k1b, soa);
                if (isnan(sample))
                    return NAN;
                value += coef[i] * coef[ordn + j] * coef[2 * ordn + k] * sample;
            }
        }
    }
    return value;
}

__device__ __forceinline__ double interpolate_point_1d(
    const double* field,
    const ShellPatchDevicePoint& point,
    const double* coef,
    double SoA_0,
    double SoA_1,
    double SoA_2)
{
    const int ordn = 2 * ghost_width;
    const double soa = shell_1d_soa(point.ssst, point.dumyd, SoA_0, SoA_1, SoA_2);
    const int active_axis = (point.dumyd == 1) ? 0 : 1;

    double value = 0.0;
    for (int m = 0; m < ordn; ++m) {
        int i1b = 0;
        int j1b = 0;
        const int k1b = point.sind[2] + 1;
        if (point.dumyd == 1) {
            i1b = point.sind[0] + m + 1;
            j1b = point.sind[1] + 1;
        } else {
            i1b = point.sind[1] + 1;
            j1b = point.sind[0] + m + 1;
        }
        const double sample = fetch_sntbd(field, point.shape, i1b, j1b, k1b, active_axis, soa);
        if (isnan(sample))
            return NAN;
        value += coef[m] * sample;
    }
    return value;
}

__device__ __forceinline__ int pointcopy_index_0b(const ShellPatchDevicePoint& point)
{
    int index[dim];
    for (int axis = 0; axis < dim; ++axis) {
#ifdef Vertex
        const double spacing = (point.bbox[dim + axis] - point.bbox[axis]) /
                               static_cast<double>(point.shape[axis] - 1);
#else
        const double spacing = (point.bbox[dim + axis] - point.bbox[axis]) /
                               static_cast<double>(point.shape[axis]);
#endif
        const int index_1b = static_cast<int>((point.lpox[axis] - point.bbox[axis]) / spacing + 0.4) + 1;
        if (index_1b < 1 || index_1b > point.shape[axis])
            return -1;
        index[axis] = index_1b - 1;
    }
    return flatten_index_0b(index[0], index[1], index[2], point.shape[0], point.shape[1]);
}

__global__ void shell_patch_apply_intertransfer_kernel(
    ShellPatchDeviceTransferList transfer,
    double* const* d_source_fields,
    double* const* d_target_fields,
    int source_block_count,
    int target_block_count,
    double SoA_0,
    double SoA_1,
    double SoA_2,
    int Symmetry,
    int* d_error_flag)
{
    const int point_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_index >= transfer.count)
        return;

    (void)Symmetry;
    const ShellPatchDevicePoint source = transfer.d_src[point_index];
    const ShellPatchDevicePoint target = transfer.d_dst[point_index];
    if (source.block_index < 0 || source.block_index >= source_block_count ||
        target.block_index < 0 || target.block_index >= target_block_count) {
        if (d_error_flag)
            atomicExch(d_error_flag, 1);
        return;
    }

    const double* source_field = d_source_fields[source.block_index];
    double* target_field = d_target_fields[target.block_index];
    if (!source_field || !target_field) {
        if (d_error_flag)
            atomicExch(d_error_flag, 1);
        return;
    }

    const double* coef = transfer.d_coef + source.coef_offset;
    double value = NAN;
    const int ordn = 2 * ghost_width;
    if (source.coef_count == dim * ordn && source.dumyd == -1)
        value = interpolate_point_3d(source_field, source, coef, SoA_0, SoA_1, SoA_2);
    else if (source.coef_count == ordn && (source.dumyd == 0 || source.dumyd == 1))
        value = interpolate_point_1d(source_field, source, coef, SoA_0, SoA_1, SoA_2);
    else if (d_error_flag)
        atomicExch(d_error_flag, 1);

    const int target_index = pointcopy_index_0b(target);
    if (target_index < 0 || isnan(value)) {
        if (d_error_flag)
            atomicExch(d_error_flag, 1);
        return;
    }
    target_field[target_index] = value;
}

int launch_blocks_for(int item_count)
{
    return (item_count + kShellTransferBlockSize - 1) / kShellTransferBlockSize;
}

}  // namespace

void shell_patch_intertransfer_gpu_launch(
    cudaStream_t stream,
    const ShellPatchDeviceTransferList& transfer,
    double* const* d_source_fields,
    double* const* d_target_fields,
    int source_block_count,
    int target_block_count,
    double SoA_0,
    double SoA_1,
    double SoA_2,
    int Symmetry,
    int* d_error_flag)
{
    if (transfer.count <= 0)
        return;
    shell_patch_apply_intertransfer_kernel<<<launch_blocks_for(transfer.count), kShellTransferBlockSize, 0, stream>>>(
        transfer,
        d_source_fields,
        d_target_fields,
        source_block_count,
        target_block_count,
        SoA_0,
        SoA_1,
        SoA_2,
        Symmetry,
        d_error_flag);
    CUDA_CHECK(cudaGetLastError());
}
