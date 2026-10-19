#ifndef SHELLPATCH_GPU_CUH
#define SHELLPATCH_GPU_CUH

#include <cuda_runtime.h>

// ============================================================================
// GPU device functions and kernel declarations for ShellPatch acceleration
// ============================================================================

// Device-side getlocalpox: converts global Cartesian -> local shell coordinates
__device__ void dev_getlocalpox(
    double gx, double gy, double gz,
    int *sst, double *lx, double *ly, double *lz);

// Batched PACK kernel: 1 thread per (point, var) pair
__global__ void batched_interpind_pack_kernel(
    int num_active, int num_vars,
    const int *__restrict__ d_active_indices,
    const double *__restrict__ d_coef,
    const int *__restrict__ d_sind,
    const int *__restrict__ d_dumyd,
    const int *__restrict__ d_ssst,
    const int *__restrict__ d_src_block_idx,
    int max_coef,
    const unsigned long long *__restrict__ d_block_fgfs,
    const int *__restrict__ d_block_shape,
    const int *__restrict__ d_block_shape_pad,
    int ordn, int fngfs_pitch,
    const double *__restrict__ d_SoA,
    const int *__restrict__ d_var_fgfs,
    double *__restrict__ d_out
);

// Batched UNPACK kernel: writes directly from GPU buffer to d_fgfs (no D2H)
__global__ void batched_pointcopy_unpack_kernel(
    int num_active, int num_vars,
    const int *__restrict__ d_active_indices,
    const double *__restrict__ d_dst_lpox,
    const int *__restrict__ d_dst_block_idx,
    const unsigned long long *__restrict__ d_block_fgfs,
    const int *__restrict__ d_block_shape,
    const int *__restrict__ d_block_shape_pad,
    const double *__restrict__ d_block_bbox,
    int fngfs_pitch,
    const int *__restrict__ d_var_fgfs,
    const double *__restrict__ d_src_buf
);

// Phase 2 GPU kernel: Cartesian boundary point search
// 1 thread per Cartesian grid point
__global__ void setup_phase2_gpu_kernel(
    int total_points,
    const int *d_block_offsets,
    const int *d_block_shape,
    const double *d_block_bbox,
    const double *d_block_CDH,
    int cs_w, double Rmin, double Rmax,
    int *d_out_count,
    double *d_out_gpox,
    double *d_out_lpox,
    int *d_out_sst,
    int max_output
);

// Launch wrappers
void gpu_batched_interpind_pack_launch(
    cudaStream_t stream,
    int num_active, int num_vars,
    const int *d_active_indices,
    const double *d_coef, const int *d_sind,
    const int *d_dumyd, const int *d_ssst,
    const int *d_src_block_idx,
    int max_coef,
    const unsigned long long *d_block_fgfs,
    const int *d_block_shape, const int *d_block_shape_pad,
    int ordn, int fngfs_pitch,
    const double *d_SoA, const int *d_var_fgfs,
    double *d_out
);

void gpu_batched_pointcopy_unpack_launch(
    cudaStream_t stream,
    int num_active, int num_vars,
    const int *d_active_indices,
    const double *d_dst_lpox,
    const int *d_dst_block_idx,
    const unsigned long long *d_block_fgfs,
    const int *d_block_shape, const int *d_block_shape_pad,
    const double *d_block_bbox,
    int fngfs_pitch,
    const int *d_var_fgfs,
    const double *d_src_buf
);

#include <cstdint>

void gpu_setup_phase2_launch(
    cudaStream_t stream,
    int total_points,
    const int *d_block_offsets, const int *d_block_shape,
    const double *d_block_bbox, const double *d_block_CDH,
    int cs_w, double Rmin, double Rmax, double DH_radial,
    int *d_out_count, double *d_out_gpox, double *d_out_lpox, int *d_out_sst,
    uint64_t *d_out_sort_key,
    int max_output
);

#endif // SHELLPATCH_GPU_CUH
