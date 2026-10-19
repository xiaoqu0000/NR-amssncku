#include <cuda_runtime.h>
#include <math.h>

#include <iostream>
#include <vector>
#include "gpu_runtime_coordinator.h"
#include "fmisc_gpu.cuh"
#include "fmisc.h"

namespace {

constexpr int kLinearBlockSize = 256;
constexpr double kHalfWeight = 0.5;
constexpr double kThreeQuarterWeight = 0.75;
constexpr double kQuarterWeight = 0.25;
constexpr double kThreeEighthWeight = 3.0 / 8.0;
constexpr double kThreeQuarterBlend = 3.0 / 4.0;
constexpr double kNegativeOneEighthWeight = -1.0 / 8.0;
constexpr double kTwentyOneThirtySecondWeight = 21.0 / 32.0;
constexpr double kSevenSixteenthWeight = 7.0 / 16.0;
constexpr double kNegativeThreeThirtySecondWeight = -3.0 / 32.0;
constexpr double kFiveThirtySecondWeight = 5.0 / 32.0;
constexpr double kFifteenSixteenthWeight = 15.0 / 16.0;

struct WindowFieldScratch {
    cudaStream_t stream;
    const double** device_sources;
    double** device_targets;
    int capacity;
};

struct UnpackWindowScratch {
    cudaStream_t stream;
    GpuUnpackWindowTransfer* device_transfers;
    int capacity;
};

inline int flattened_extent(const int extent[3]) {
    return extent[0] * extent[1] * extent[2];
}

inline int launch_blocks_for(int item_count) {
    return (item_count + kLinearBlockSize - 1) / kLinearBlockSize;
}

UnpackWindowScratch& unpack_window_scratch_for_stream(std::vector<UnpackWindowScratch>& scratches, cudaStream_t stream) {
    for (std::vector<UnpackWindowScratch>::iterator scratch = scratches.begin(); scratch != scratches.end(); ++scratch) {
        if (scratch->stream == stream) {
            return *scratch;
        }
    }

    UnpackWindowScratch scratch = {stream, nullptr, 0};
    scratches.push_back(scratch);
    return scratches.back();
}

std::vector<UnpackWindowScratch>& unpack_window_scratches() {
    static std::vector<UnpackWindowScratch> scratches;
    return scratches;
}

void ensure_unpack_window_capacity(UnpackWindowScratch& scratch, int transfer_count) {
    if (transfer_count <= scratch.capacity) {
        return;
    }

    if (scratch.device_transfers) {
        CUDA_CHECK(cudaFree(scratch.device_transfers));
        scratch.device_transfers = nullptr;
    }

    CUDA_CHECK(cudaMalloc(&scratch.device_transfers, transfer_count * sizeof(GpuUnpackWindowTransfer)));
    scratch.capacity = transfer_count;
}

WindowFieldScratch& window_field_scratch_for_stream(std::vector<WindowFieldScratch>& scratches, cudaStream_t stream) {
    for (std::vector<WindowFieldScratch>::iterator scratch = scratches.begin(); scratch != scratches.end(); ++scratch) {
        if (scratch->stream == stream) {
            return *scratch;
        }
    }

    WindowFieldScratch scratch = {stream, nullptr, nullptr, 0};
    scratches.push_back(scratch);
    return scratches.back();
}

std::vector<WindowFieldScratch>& pack_field_scratches() {
    static std::vector<WindowFieldScratch> scratches;
    return scratches;
}

std::vector<WindowFieldScratch>& unpack_field_scratches() {
    static std::vector<WindowFieldScratch> scratches;
    return scratches;
}

void ensure_window_field_capacity(WindowFieldScratch& scratch, int field_count) {
    if (field_count <= scratch.capacity) {
        return;
    }

    if (scratch.device_sources) {
        CUDA_CHECK(cudaFree(scratch.device_sources));
        scratch.device_sources = nullptr;
    }
    if (scratch.device_targets) {
        CUDA_CHECK(cudaFree(scratch.device_targets));
        scratch.device_targets = nullptr;
    }

    CUDA_CHECK(cudaMalloc(&scratch.device_sources, field_count * sizeof(double*)));
    CUDA_CHECK(cudaMalloc(&scratch.device_targets, field_count * sizeof(double*)));
    scratch.capacity = field_count;
}

__device__ __forceinline__ bool point_is_inside_patch(
    double px, double py, double pz,
    double lower_x, double lower_y, double lower_z,
    double upper_x, double upper_y, double upper_z
) {
    return px >= lower_x && px <= upper_x
        && py >= lower_y && py <= upper_y
        && pz >= lower_z && pz <= upper_z;
}

__device__ __forceinline__ void unpack_linear_offset(
    int linear_index,
    int nx,
    int ny,
    int& i,
    int& j,
    int& k
) {
    i = linear_index % nx;
    j = (linear_index / nx) % ny;
    k = linear_index / (nx * ny);
}

__device__ __forceinline__ int flatten_offset_index(
    int i,
    int j,
    int k,
    int nx,
    int ny
) {
    return k * nx * ny + j * nx + i;
}

__global__ void clamp_lower_bound_kernel(int cell_count, double* chi_field, double floor_value);
__global__ void pack_patch_window_kernel(
    const double* __restrict__ source_field,
    double* __restrict__ packed_buffer,
    int source_nx,
    int source_ny,
    int packed_nx,
    int packed_ny,
    int packed_nz,
    int offset_x,
    int offset_y,
    int offset_z
);
__global__ void unpack_patch_window_kernel(
    const double* __restrict__ packed_buffer,
    double* __restrict__ target_field,
    int target_nx,
    int target_ny,
    int packed_nx,
    int packed_ny,
    int packed_nz,
    int offset_x,
    int offset_y,
    int offset_z
);
__global__ void pack_patch_window_fields_kernel(
    const double* const* __restrict__ source_fields,
    double* const* __restrict__ packed_buffers,
    int field_count,
    int source_nx,
    int source_ny,
    int packed_nx,
    int packed_ny,
    int packed_nz,
    int offset_x,
    int offset_y,
    int offset_z
);
__global__ void unpack_patch_window_fields_kernel(
    const double* const* __restrict__ packed_buffers,
    double* const* __restrict__ target_fields,
    int field_count,
    int target_nx,
    int target_ny,
    int packed_nx,
    int packed_ny,
    int packed_nz,
    int offset_x,
    int offset_y,
    int offset_z
);
__global__ void unpack_patch_window_batch_kernel(
    const GpuUnpackWindowTransfer* __restrict__ transfers,
    int transfer_count,
    int blocks_x
);
__global__ void blend_two_fields_kernel(
    int cell_count,
    const double* first_field,
    const double* second_field,
    double first_weight,
    double second_weight,
    double* blended_field
);
__global__ void blend_three_fields_kernel(
    int cell_count,
    const double* first_field,
    const double* second_field,
    const double* third_field,
    double first_weight,
    double second_weight,
    double third_weight,
    double* blended_field
);
__global__ void accumulate_interp_samples_kernel(
    int sample_count,
    int dimension_count,
    double* sample_x,
    double* sample_y,
    double* sample_z,
    int extent_x,
    int extent_y,
    int extent_z,
    double* axis_x,
    double* axis_y,
    double* axis_z,
    double* field,
    double lower_x,
    double lower_y,
    double lower_z,
    double upper_x,
    double upper_y,
    double upper_z,
    int interpolation_order,
    double symmetry_x,
    double symmetry_y,
    double symmetry_z,
    int symmetry_mode,
    int variable_slot,
    int variable_count,
    double* shell_values,
    int* hit_count
);
__global__ void accumulate_interp_first_match_kernel(
    int sample_count,
    int dimension_count,
    double* sample_x,
    double* sample_y,
    double* sample_z,
    int extent_x,
    int extent_y,
    int extent_z,
    double* axis_x,
    double* axis_y,
    double* axis_z,
    double* field,
    double lower_x,
    double lower_y,
    double lower_z,
    double upper_x,
    double upper_y,
    double upper_z,
    int interpolation_order,
    double symmetry_x,
    double symmetry_y,
    double symmetry_z,
    int symmetry_mode,
    int variable_slot,
    int variable_count,
    double* shell_values,
    int* hit_count,
    int block_tag,
    int* owner
);
__global__ void accumulate_l2_norm_kernel(
    const double* __restrict__ field,
    int imin,
    int imax,
    int jmin,
    int jmax,
    int kmin,
    int kmax,
    int nx,
    int ny,
    int nz,
    double* __restrict__ reduced_sum
);
__global__ void write_unique_interp_samples_kernel(
    int active_count,
    int dimension_count,
    int* active_indices,
    double* sample_x,
    double* sample_y,
    double* sample_z,
    int extent_x,
    int extent_y,
    int extent_z,
    double* axis_x,
    double* axis_y,
    double* axis_z,
    double* field,
    int interpolation_order,
    double symmetry_x,
    double symmetry_y,
    double symmetry_z,
    int symmetry_mode,
    int variable_slot,
    int variable_count,
    double* shell_values
);

}  // namespace

void gpu_global_interp_amr_launch(
    cudaStream_t stream,
    int active_count,
    int DIM,
    int* d_active_indices,
    double* d_XX_0,
    double* d_XX_1,
    double* d_XX_2,
    int shape_0,
    int shape_1,
    int shape_2,
    double* d_X_0,
    double* d_X_1,
    double* d_X_2,
    double* d_field,
    int ordn,
    double SoA_0,
    double SoA_1,
    double SoA_2,
    int Symmetry,
    int var_idx,
    int num_var,
    double* d_shellf
) {
    if (active_count == 0) {
        return;
    }

    write_unique_interp_samples_kernel<<<launch_blocks_for(active_count), kLinearBlockSize, 0, stream>>>(
        active_count,
        DIM,
        d_active_indices,
        d_XX_0,
        d_XX_1,
        d_XX_2,
        shape_0,
        shape_1,
        shape_2,
        d_X_0,
        d_X_1,
        d_X_2,
        d_field,
        ordn,
        SoA_0,
        SoA_1,
        SoA_2,
        Symmetry,
        var_idx,
        num_var,
        d_shellf
    );
}

void gpu_l2normhelper_launch(
    cudaStream_t stream,
    const int* ex,
    const double* X,
    const double* Y,
    const double* Z,
    double xmin,
    double ymin,
    double zmin,
    double xmax,
    double ymax,
    double zmax,
    const double* d_f,
    double& f_out,
    int gw
) {
    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    int imin = gw;
    int jmin = gw;
    int kmin = gw;
    int imax = ex[0] - gw - 1;
    int jmax = ex[1] - gw - 1;
    int kmax = ex[2] - gw - 1;

    if (fabs(X[ex[0] - 1] - xmax) < dX) imax = ex[0] - 1;
    if (fabs(Y[ex[1] - 1] - ymax) < dY) jmax = ex[1] - 1;
    if (fabs(Z[ex[2] - 1] - zmax) < dZ) kmax = ex[2] - 1;
    if (fabs(X[0] - xmin) < dX) imin = 0;
    if (fabs(Y[0] - ymin) < dY) jmin = 0;
    if (fabs(Z[0] - zmin) < dZ) kmin = 0;

    const int nx_proc = imax - imin + 1;
    const int ny_proc = jmax - jmin + 1;
    const int nz_proc = kmax - kmin + 1;
    if (nx_proc <= 0 || ny_proc <= 0 || nz_proc <= 0) {
        f_out = 0.0;
        return;
    }

    double* d_sum = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(1, stream);
    const dim3 block_shape(8, 8, 8);
    const dim3 grid_shape(
        (nx_proc + block_shape.x - 1) / block_shape.x,
        (ny_proc + block_shape.y - 1) / block_shape.y,
        (nz_proc + block_shape.z - 1) / block_shape.z
    );

    accumulate_l2_norm_kernel<<<grid_shape, block_shape, 0, stream>>>(
        d_f,
        imin,
        imax,
        jmin,
        jmax,
        kmin,
        kmax,
        ex[0],
        ex[1],
        ex[2],
        d_sum
    );

    double h_sum = 0.0;
    GpuRuntimeCoordinator::instance().copy_to_host(&h_sum, d_sum, 1, stream);
    GpuRuntimeCoordinator::instance().release_device_buffer(d_sum, 1, stream);
    GpuRuntimeCoordinator::wait_for_stream(stream);
    f_out = h_sum * dX * dY * dZ;
}

void gpu_global_interp_launch(
    cudaStream_t stream,
    int NN,
    int DIM,
    double* d_XX_0,
    double* d_XX_1,
    double* d_XX_2,
    int shape_0,
    int shape_1,
    int shape_2,
    double* d_X_0,
    double* d_X_1,
    double* d_X_2,
    double* d_field,
    double llb_0,
    double llb_1,
    double llb_2,
    double uub_0,
    double uub_1,
    double uub_2,
    int ordn,
    double SoA_0,
    double SoA_1,
    double SoA_2,
    int Symmetry,
    int var_idx,
    int num_var,
    double* d_shellf,
    int* d_weight
) {
    accumulate_interp_samples_kernel<<<launch_blocks_for(NN), kLinearBlockSize, 0, stream>>>(
        NN,
        DIM,
        d_XX_0,
        d_XX_1,
        d_XX_2,
        shape_0,
        shape_1,
        shape_2,
        d_X_0,
        d_X_1,
        d_X_2,
        d_field,
        llb_0,
        llb_1,
        llb_2,
        uub_0,
        uub_1,
        uub_2,
        ordn,
        SoA_0,
        SoA_1,
        SoA_2,
        Symmetry,
        var_idx,
        num_var,
        d_shellf,
        d_weight
    );
}

void gpu_global_interp_first_match_launch(
    cudaStream_t stream,
    int NN,
    int DIM,
    double* d_XX_0,
    double* d_XX_1,
    double* d_XX_2,
    int shape_0,
    int shape_1,
    int shape_2,
    double* d_X_0,
    double* d_X_1,
    double* d_X_2,
    double* d_field,
    double llb_0,
    double llb_1,
    double llb_2,
    double uub_0,
    double uub_1,
    double uub_2,
    int ordn,
    double SoA_0,
    double SoA_1,
    double SoA_2,
    int Symmetry,
    int var_idx,
    int num_var,
    double* d_shellf,
    int* d_weight,
    int block_tag,
    int* d_owner
) {
    accumulate_interp_first_match_kernel<<<launch_blocks_for(NN), kLinearBlockSize, 0, stream>>>(
        NN,
        DIM,
        d_XX_0,
        d_XX_1,
        d_XX_2,
        shape_0,
        shape_1,
        shape_2,
        d_X_0,
        d_X_1,
        d_X_2,
        d_field,
        llb_0,
        llb_1,
        llb_2,
        uub_0,
        uub_1,
        uub_2,
        ordn,
        SoA_0,
        SoA_1,
        SoA_2,
        Symmetry,
        var_idx,
        num_var,
        d_shellf,
        d_weight,
        block_tag,
        d_owner
    );
}

void gpu_average2m_launch(cudaStream_t stream, const int ext[3], const double* d_f1, const double* d_f2, const double* d_f3, double* d_fout) {
    const int cell_count = flattened_extent(ext);
    blend_three_fields_kernel<<<launch_blocks_for(cell_count), kLinearBlockSize, 0, stream>>>(
        cell_count,
        d_f1,
        d_f2,
        d_f3,
        kFiveThirtySecondWeight,
        kFifteenSixteenthWeight,
        kNegativeThreeThirtySecondWeight,
        d_fout
    );
}

void gpu_average2p_launch(cudaStream_t stream, const int ext[3], const double* d_f1, const double* d_f2, const double* d_f3, double* d_fout) {
    const int cell_count = flattened_extent(ext);
    blend_three_fields_kernel<<<launch_blocks_for(cell_count), kLinearBlockSize, 0, stream>>>(
        cell_count,
        d_f1,
        d_f2,
        d_f3,
        kTwentyOneThirtySecondWeight,
        kSevenSixteenthWeight,
        kNegativeThreeThirtySecondWeight,
        d_fout
    );
}

void gpu_average2_launch(cudaStream_t stream, const int ext[3], const double* d_f1, const double* d_f2, const double* d_f3, double* d_fout) {
    const int cell_count = flattened_extent(ext);
    blend_three_fields_kernel<<<launch_blocks_for(cell_count), kLinearBlockSize, 0, stream>>>(
        cell_count,
        d_f1,
        d_f2,
        d_f3,
        kThreeEighthWeight,
        kThreeQuarterBlend,
        kNegativeOneEighthWeight,
        d_fout
    );
}

void gpu_average3_launch(cudaStream_t stream, const int ext[3], const double* d_f1, const double* d_f2, double* d_fout) {
    const int cell_count = flattened_extent(ext);
    blend_two_fields_kernel<<<launch_blocks_for(cell_count), kLinearBlockSize, 0, stream>>>(
        cell_count,
        d_f1,
        d_f2,
        kThreeQuarterWeight,
        kQuarterWeight,
        d_fout
    );
}

void gpu_average_launch(cudaStream_t stream, const int ext[3], const double* d_f1, const double* d_f2, double* d_fout) {
    const int cell_count = flattened_extent(ext);
    blend_two_fields_kernel<<<launch_blocks_for(cell_count), kLinearBlockSize, 0, stream>>>(
        cell_count,
        d_f1,
        d_f2,
        kHalfWeight,
        kHalfWeight,
        d_fout
    );
}

void gpu_unpack_launch(
    cudaStream_t stream,
    const double* d_src_1d,
    double* d_dst_3d,
    int dst_nx,
    int dst_ny,
    int src_nx,
    int src_ny,
    int src_nz,
    int off_x,
    int off_y,
    int off_z
) {
    const int cell_count = src_nx * src_ny * src_nz;
    unpack_patch_window_kernel<<<launch_blocks_for(cell_count), kLinearBlockSize, 0, stream>>>(
        d_src_1d,
        d_dst_3d,
        dst_nx,
        dst_ny,
        src_nx,
        src_ny,
        src_nz,
        off_x,
        off_y,
        off_z
    );
}

void gpu_unpack_fields_launch(
    cudaStream_t stream,
    const double* const* d_src_fields,
    double* const* d_dst_fields,
    int field_count,
    int dst_nx,
    int dst_ny,
    int src_nx,
    int src_ny,
    int src_nz,
    int off_x,
    int off_y,
    int off_z
) {
    if (field_count <= 0) {
        return;
    }

    WindowFieldScratch& scratch = window_field_scratch_for_stream(unpack_field_scratches(), stream);
    ensure_window_field_capacity(scratch, field_count);
    CUDA_CHECK(cudaMemcpyAsync(scratch.device_sources, d_src_fields, field_count * sizeof(double*), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(scratch.device_targets, d_dst_fields, field_count * sizeof(double*), cudaMemcpyHostToDevice, stream));

    const int cell_count = src_nx * src_ny * src_nz;
    const dim3 block(kLinearBlockSize, 1, 1);
    const dim3 grid(launch_blocks_for(cell_count), field_count, 1);
    unpack_patch_window_fields_kernel<<<grid, block, 0, stream>>>(
        scratch.device_sources,
        scratch.device_targets,
        field_count,
        dst_nx,
        dst_ny,
        src_nx,
        src_ny,
        src_nz,
        off_x,
        off_y,
        off_z
    );
}

int gpu_window_transfer_block_count(int cell_count) {
    return launch_blocks_for(cell_count);
}

void gpu_unpack_batch_launch(
    cudaStream_t stream,
    const GpuUnpackWindowTransfer* transfers,
    int transfer_count,
    int max_blocks_x
) {
    if (transfer_count <= 0 || max_blocks_x <= 0) {
        return;
    }

    UnpackWindowScratch& scratch = unpack_window_scratch_for_stream(unpack_window_scratches(), stream);
    ensure_unpack_window_capacity(scratch, transfer_count);
    CUDA_CHECK(cudaMemcpyAsync(
        scratch.device_transfers,
        transfers,
        transfer_count * sizeof(GpuUnpackWindowTransfer),
        cudaMemcpyHostToDevice,
        stream
    ));

    const int max_transfers_per_launch = 65535;
    int offset = 0;
    while (offset < transfer_count) {
        const int batch_count = std::min(max_transfers_per_launch, transfer_count - offset);
        const dim3 grid(max_blocks_x, batch_count, 1);
        unpack_patch_window_batch_kernel<<<grid, kLinearBlockSize, 0, stream>>>(
            scratch.device_transfers + offset,
            batch_count,
            max_blocks_x
        );
        offset += batch_count;
    }
}

void gpu_pack_launch(
    cudaStream_t stream,
    const double* d_src_3d,
    double* d_dst_1d,
    int src_nx,
    int src_ny,
    int dst_nx,
    int dst_ny,
    int dst_nz,
    int off_x,
    int off_y,
    int off_z
) {
    const int cell_count = dst_nx * dst_ny * dst_nz;
    pack_patch_window_kernel<<<launch_blocks_for(cell_count), kLinearBlockSize, 0, stream>>>(
        d_src_3d,
        d_dst_1d,
        src_nx,
        src_ny,
        dst_nx,
        dst_ny,
        dst_nz,
        off_x,
        off_y,
        off_z
    );
}

void gpu_pack_fields_launch(
    cudaStream_t stream,
    const double* const* d_src_fields,
    double* const* d_dst_fields,
    int field_count,
    int src_nx,
    int src_ny,
    int dst_nx,
    int dst_ny,
    int dst_nz,
    int off_x,
    int off_y,
    int off_z
) {
    if (field_count <= 0) {
        return;
    }

    WindowFieldScratch& scratch = window_field_scratch_for_stream(pack_field_scratches(), stream);
    ensure_window_field_capacity(scratch, field_count);
    CUDA_CHECK(cudaMemcpyAsync(scratch.device_sources, d_src_fields, field_count * sizeof(double*), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(scratch.device_targets, d_dst_fields, field_count * sizeof(double*), cudaMemcpyHostToDevice, stream));

    const int cell_count = dst_nx * dst_ny * dst_nz;
    const dim3 block(kLinearBlockSize, 1, 1);
    const dim3 grid(launch_blocks_for(cell_count), field_count, 1);
    pack_patch_window_fields_kernel<<<grid, block, 0, stream>>>(
        scratch.device_sources,
        scratch.device_targets,
        field_count,
        src_nx,
        src_ny,
        dst_nx,
        dst_ny,
        dst_nz,
        off_x,
        off_y,
        off_z
    );
}

void gpu_lowerboundset_launch(cudaStream_t &stream, int ex[3], double* d_chi0, double TINNY) {
    const int cell_count = flattened_extent(ex);
    clamp_lower_bound_kernel<<<launch_blocks_for(cell_count), kLinearBlockSize, 0, stream>>>(cell_count, d_chi0, TINNY);
}

namespace {

__global__ void clamp_lower_bound_kernel(int cell_count, double* chi_field, double floor_value) {
    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_id >= cell_count) {
        return;
    }

    if (chi_field[point_id] < floor_value) {
        chi_field[point_id] = floor_value;
    }
}

__global__ void pack_patch_window_kernel(
    const double* __restrict__ source_field,
    double* __restrict__ packed_buffer,
    int source_nx,
    int source_ny,
    int packed_nx,
    int packed_ny,
    int packed_nz,
    int offset_x,
    int offset_y,
    int offset_z
) {
    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    const int packed_count = packed_nx * packed_ny * packed_nz;
    if (point_id >= packed_count) {
        return;
    }

    int local_i = 0;
    int local_j = 0;
    int local_k = 0;
    unpack_linear_offset(point_id, packed_nx, packed_ny, local_i, local_j, local_k);
    const int source_index = flatten_offset_index(local_i + offset_x, local_j + offset_y, local_k + offset_z, source_nx, source_ny);
    packed_buffer[point_id] = source_field[source_index];
}

__global__ void unpack_patch_window_kernel(
    const double* __restrict__ packed_buffer,
    double* __restrict__ target_field,
    int target_nx,
    int target_ny,
    int packed_nx,
    int packed_ny,
    int packed_nz,
    int offset_x,
    int offset_y,
    int offset_z
) {
    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    const int packed_count = packed_nx * packed_ny * packed_nz;
    if (point_id >= packed_count) {
        return;
    }

    int local_i = 0;
    int local_j = 0;
    int local_k = 0;
    unpack_linear_offset(point_id, packed_nx, packed_ny, local_i, local_j, local_k);
    const int target_index = flatten_offset_index(local_i + offset_x, local_j + offset_y, local_k + offset_z, target_nx, target_ny);
    target_field[target_index] = packed_buffer[point_id];
}

__global__ void pack_patch_window_fields_kernel(
    const double* const* __restrict__ source_fields,
    double* const* __restrict__ packed_buffers,
    int field_count,
    int source_nx,
    int source_ny,
    int packed_nx,
    int packed_ny,
    int packed_nz,
    int offset_x,
    int offset_y,
    int offset_z
) {
    const int field_index = blockIdx.y;
    if (field_index >= field_count) {
        return;
    }

    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    const int packed_count = packed_nx * packed_ny * packed_nz;
    if (point_id >= packed_count) {
        return;
    }

    int local_i = 0;
    int local_j = 0;
    int local_k = 0;
    unpack_linear_offset(point_id, packed_nx, packed_ny, local_i, local_j, local_k);
    const int source_index = flatten_offset_index(local_i + offset_x, local_j + offset_y, local_k + offset_z, source_nx, source_ny);
    packed_buffers[field_index][point_id] = source_fields[field_index][source_index];
}

__global__ void unpack_patch_window_fields_kernel(
    const double* const* __restrict__ packed_buffers,
    double* const* __restrict__ target_fields,
    int field_count,
    int target_nx,
    int target_ny,
    int packed_nx,
    int packed_ny,
    int packed_nz,
    int offset_x,
    int offset_y,
    int offset_z
) {
    const int field_index = blockIdx.y;
    if (field_index >= field_count) {
        return;
    }

    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    const int packed_count = packed_nx * packed_ny * packed_nz;
    if (point_id >= packed_count) {
        return;
    }

    int local_i = 0;
    int local_j = 0;
    int local_k = 0;
    unpack_linear_offset(point_id, packed_nx, packed_ny, local_i, local_j, local_k);
    const int target_index = flatten_offset_index(local_i + offset_x, local_j + offset_y, local_k + offset_z, target_nx, target_ny);
    target_fields[field_index][target_index] = packed_buffers[field_index][point_id];
}

__global__ void unpack_patch_window_batch_kernel(
    const GpuUnpackWindowTransfer* __restrict__ transfers,
    int transfer_count,
    int max_blocks_x
) {
    const int transfer_index = blockIdx.y;
    if (transfer_index >= transfer_count || blockIdx.x >= max_blocks_x) {
        return;
    }

    const GpuUnpackWindowTransfer transfer = transfers[transfer_index];
    if (blockIdx.x >= transfer.blocks_x) {
        return;
    }

    const int packed_count = transfer.packed_nx * transfer.packed_ny * transfer.packed_nz;
    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_id >= packed_count) {
        return;
    }

    int local_i = 0;
    int local_j = 0;
    int local_k = 0;
    unpack_linear_offset(point_id, transfer.packed_nx, transfer.packed_ny, local_i, local_j, local_k);
    const int target_index = flatten_offset_index(
        local_i + transfer.offset_x,
        local_j + transfer.offset_y,
        local_k + transfer.offset_z,
        transfer.target_nx,
        transfer.target_ny
    );
    transfer.target_field[target_index] = transfer.packed_buffer[point_id];
}

__global__ void blend_two_fields_kernel(
    int cell_count,
    const double* first_field,
    const double* second_field,
    double first_weight,
    double second_weight,
    double* blended_field
) {
    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_id >= cell_count) {
        return;
    }

    blended_field[point_id] = first_weight * first_field[point_id] + second_weight * second_field[point_id];
}

__global__ void blend_three_fields_kernel(
    int cell_count,
    const double* first_field,
    const double* second_field,
    const double* third_field,
    double first_weight,
    double second_weight,
    double third_weight,
    double* blended_field
) {
    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_id >= cell_count) {
        return;
    }

    blended_field[point_id] =
        first_weight * first_field[point_id]
      + second_weight * second_field[point_id]
      + third_weight * third_field[point_id];
}

__global__ void accumulate_interp_samples_kernel(
    int sample_count,
    int dimension_count,
    double* sample_x,
    double* sample_y,
    double* sample_z,
    int extent_x,
    int extent_y,
    int extent_z,
    double* axis_x,
    double* axis_y,
    double* axis_z,
    double* field,
    double lower_x,
    double lower_y,
    double lower_z,
    double upper_x,
    double upper_y,
    double upper_z,
    int interpolation_order,
    double symmetry_x,
    double symmetry_y,
    double symmetry_z,
    int symmetry_mode,
    int variable_slot,
    int variable_count,
    double* shell_values,
    int* hit_count
) {
    const int sample_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (sample_index >= sample_count) {
        return;
    }

    const double px = sample_x[sample_index];
    const double py = sample_y[sample_index];
    const double pz = sample_z[sample_index];
    if (!point_is_inside_patch(px, py, pz, lower_x, lower_y, lower_z, upper_x, upper_y, upper_z)) {
        return;
    }

    double* axis_bundle[3] = {axis_x, axis_y, axis_z};
    double symmetry_axes[3] = {symmetry_x, symmetry_y, symmetry_z};
    const int local_extent[3] = {extent_x, extent_y, extent_z};

    double interpolated_value = 0.0;
    global_interp_device(
        local_extent,
        axis_bundle[0],
        axis_bundle[1],
        axis_bundle[2],
        field,
        &interpolated_value,
        px,
        py,
        pz,
        interpolation_order,
        symmetry_axes,
        symmetry_mode
    );

    atomicAdd(&shell_values[sample_index * variable_count + variable_slot], interpolated_value);
    if (variable_slot == 0) {
        atomicAdd(&hit_count[sample_index], 1);
    }

    (void)dimension_count;
}

__global__ void accumulate_l2_norm_kernel(
    const double* __restrict__ field,
    int imin,
    int imax,
    int jmin,
    int jmax,
    int kmin,
    int kmax,
    int nx,
    int ny,
    int nz,
    double* __restrict__ reduced_sum
) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x + imin;
    const int j = blockIdx.y * blockDim.y + threadIdx.y + jmin;
    const int k = blockIdx.z * blockDim.z + threadIdx.z + kmin;

    double local_sum = 0.0;
    if (i <= imax && j <= jmax && k <= kmax) {
        const long long flat_index = (long long)i + (long long)j * nx + (long long)k * nx * ny;
        const double field_value = field[flat_index];
        local_sum = field_value * field_value;
    }

    local_sum = warpReduceSum(local_sum);
    const int lane_id = threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z * blockDim.x * blockDim.y;
    if ((lane_id % warpSize) == 0) {
        atomicAdd(reduced_sum, local_sum);
    }

    (void)nz;
}

__global__ void write_unique_interp_samples_kernel(
    int active_count,
    int dimension_count,
    int* active_indices,
    double* sample_x,
    double* sample_y,
    double* sample_z,
    int extent_x,
    int extent_y,
    int extent_z,
    double* axis_x,
    double* axis_y,
    double* axis_z,
    double* field,
    int interpolation_order,
    double symmetry_x,
    double symmetry_y,
    double symmetry_z,
    int symmetry_mode,
    int variable_slot,
    int variable_count,
    double* shell_values
) {
    const int active_slot = blockIdx.x * blockDim.x + threadIdx.x;
    if (active_slot >= active_count) {
        return;
    }

    const int sample_index = active_indices[active_slot];
    const double px = sample_x[sample_index];
    const double py = sample_y[sample_index];
    const double pz = sample_z[sample_index];

    double* axis_bundle[3] = {axis_x, axis_y, axis_z};
    double symmetry_axes[3] = {symmetry_x, symmetry_y, symmetry_z};
    const int local_extent[3] = {extent_x, extent_y, extent_z};

    double interpolated_value = 0.0;
    global_interp_device(
        local_extent,
        axis_bundle[0],
        axis_bundle[1],
        axis_bundle[2],
        field,
        &interpolated_value,
        px,
        py,
        pz,
        interpolation_order,
        symmetry_axes,
        symmetry_mode
    );

    shell_values[sample_index * variable_count + variable_slot] = interpolated_value;
    (void)dimension_count;
}

__global__ void accumulate_interp_first_match_kernel(
    int sample_count,
    int dimension_count,
    double* sample_x,
    double* sample_y,
    double* sample_z,
    int extent_x,
    int extent_y,
    int extent_z,
    double* axis_x,
    double* axis_y,
    double* axis_z,
    double* field,
    double lower_x,
    double lower_y,
    double lower_z,
    double upper_x,
    double upper_y,
    double upper_z,
    int interpolation_order,
    double symmetry_x,
    double symmetry_y,
    double symmetry_z,
    int symmetry_mode,
    int variable_slot,
    int variable_count,
    double* shell_values,
    int* hit_count,
    int block_tag,
    int* owner
) {
    const int sample_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (sample_index >= sample_count) {
        return;
    }

    const double px = sample_x[sample_index];
    const double py = sample_y[sample_index];
    const double pz = sample_z[sample_index];
    if (!point_is_inside_patch(px, py, pz, lower_x, lower_y, lower_z, upper_x, upper_y, upper_z)) {
        return;
    }

    double* axis_bundle[3] = {axis_x, axis_y, axis_z};
    double symmetry_axes[3] = {symmetry_x, symmetry_y, symmetry_z};
    const int local_extent[3] = {extent_x, extent_y, extent_z};

    double interpolated_value = 0.0;
    global_interp_device(
        local_extent,
        axis_bundle[0],
        axis_bundle[1],
        axis_bundle[2],
        field,
        &interpolated_value,
        px,
        py,
        pz,
        interpolation_order,
        symmetry_axes,
        symmetry_mode
    );

    const int old_owner = atomicCAS(&owner[sample_index], -1, block_tag);
    if (old_owner == -1 || old_owner == block_tag) {
        shell_values[sample_index * variable_count + variable_slot] = interpolated_value;
        if (variable_slot == 0) {
            atomicAdd(&hit_count[sample_index], 1);
        }
    }

    (void)dimension_count;
}

}  // namespace
