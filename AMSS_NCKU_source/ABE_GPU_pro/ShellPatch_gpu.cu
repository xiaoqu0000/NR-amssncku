// ShellPatch_gpu.cu
// GPU kernels for ShellPatch: batched PACK/UNPACK interpolation and
// Phase 2 Cartesian boundary point processing.
// Follows the same pattern as MPatch_gpu.cu.

#include <cuda_runtime.h>
#include <math.h>
#include <stdint.h>

#include "gpu_manager.h"
#include "gpu_common.cuh"
#include "fmisc_gpu.cuh"
#include "ShellPatch_gpu.cuh"
#include "macrodef.h"

// ============================================================================
// Device-side getlocalpox: converts global Cartesian -> local shell coordinates
// ============================================================================

__device__ void dev_getlocalpox(
    double gx, double gy, double gz,
    int *sst, double *lx, double *ly, double *lz)
{
    double R = sqrt(gx*gx + gy*gy + gz*gz);
    *lz = R;
    double ax = fabs(gx), ay = fabs(gy), az = fabs(gz);
    // Use atan (not atan2) to match CPU getlocalpox convention (ShellPatch.C:1235-1280)
    // For negative-denominator faces (sst=1,3,5), atan2 would give quadrant-corrected
    // values in [-pi,pi], while CPU atan(x/z) gives values in [-pi/2,pi/2].
    if (az >= ax && az >= ay) {
        *sst = (gz > 0) ? 0 : 1;
        *lx = atan(gx / gz);
        *ly = atan(gy / gz);
    } else if (ax >= ay && ax >= az) {
        *sst = (gx > 0) ? 2 : 3;
        *lx = atan(gy / gx);
        *ly = atan(gz / gx);
    } else {
        *sst = (gy > 0) ? 4 : 5;
        *lx = atan(gx / gy);
        *ly = atan(gz / gy);
    }
}

// ============================================================================
// Phase 2 GPU kernel: Cartesian boundary point search
// 1 thread per Cartesian grid point — processes all points in parallel.
// Does coordinate computation + getlocalpox + radial filtering.
// CPU post-processing does the fast shell patch block search.
// ============================================================================

__global__ void setup_phase2_gpu_kernel(
    int total_points,
    const int    *d_block_offsets,
    const int    *d_block_shape,
    const double *d_block_bbox,
    const double *d_block_CDH,
    int cs_w, double Rmin, double Rmax, double DH_radial,
    int    *d_out_count,
    double *d_out_gpox,
    double *d_out_lpox,
    int    *d_out_sst,
    uint64_t *d_out_sort_key,
    int max_output
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= total_points) return;

    // 定位当前的 Cartesian Block
    int bid = 0;
    while (bid < 1024 && tid >= d_block_offsets[bid + 1]) bid++;
    if (d_block_offsets[bid + 1] == 0) return;

    int local_tid = tid - d_block_offsets[bid];
    int sh0 = d_block_shape[bid * 3 + 0];
    int sh1 = d_block_shape[bid * 3 + 1];
    int sh2 = d_block_shape[bid * 3 + 2];
    if (sh0 <= 0 || sh1 <= 0 || sh2 <= 0) return;

    // 解析局部三维索引
    int ix = local_tid % sh0;
    int iy = (local_tid / sh0) % sh1;
    int iz = local_tid / (sh0 * sh1);

    // 边界点过滤
    int x_lo = cs_w, x_hi = sh0 - cs_w;
    int y_lo = cs_w, y_hi = sh1 - cs_w;
    int z_lo = cs_w, z_hi = sh2 - cs_w;
    if (x_lo < 0) x_lo = 0; if (x_hi > sh0) x_hi = sh0;
    if (y_lo < 0) y_lo = 0; if (y_hi > sh1) y_hi = sh1;
    if (z_lo < 0) z_lo = 0; if (z_hi > sh2) z_hi = sh2;

    bool boundary = (ix < x_lo || ix >= x_hi || iy < y_lo || iy >= y_hi || iz < z_lo || iz >= z_hi);
    if (!boundary) return;

    double bb0 = d_block_bbox[bid * 6 + 0];
    double bb1 = d_block_bbox[bid * 6 + 1];
    double bb2 = d_block_bbox[bid * 6 + 2];
    double h0 = d_block_CDH[bid * 3 + 0];
    double h1 = d_block_CDH[bid * 3 + 1];
    double h2 = d_block_CDH[bid * 3 + 2];

#ifdef Vertex
    double x = bb0 + ix * h0;
    double y = bb1 + iy * h1;
    double z = bb2 + iz * h2;
#else
    double x = bb0 + (ix + 0.5) * h0;
    double y = bb1 + (iy + 0.5) * h1;
    double z = bb2 + (iz + 0.5) * h2;
#endif

    int sst; double lx, ly, lz;
    dev_getlocalpox(x, y, z, &sst, &lx, &ly, &lz);
    if (lz >= Rmin - 0.0001 * DH_radial && lz <= Rmax + 0.0001 * DH_radial) {
        int idx = atomicAdd(d_out_count, 1);
        if (idx < max_output) {
            d_out_gpox[idx * 3 + 0] = x; d_out_gpox[idx * 3 + 1] = y; d_out_gpox[idx * 3 + 2] = z;
            d_out_lpox[idx * 3 + 0] = lx; d_out_lpox[idx * 3 + 1] = ly; d_out_lpox[idx * 3 + 2] = lz;
            d_out_sst[idx] = sst;

            // [新增] 构造严格的 64 位排序键: bid -> iz -> iy -> ix (复刻 CPU 原生遍历顺序)
            uint64_t key = ((uint64_t)bid << 48) | ((uint64_t)iz << 32) | ((uint64_t)iy << 16) | (uint64_t)ix;
            d_out_sort_key[idx] = key;
        }
    }
}

void gpu_setup_phase2_launch(
    cudaStream_t stream,
    int total_points,
    const int *d_block_offsets, const int *d_block_shape,
    const double *d_block_bbox, const double *d_block_CDH,
    int cs_w, double Rmin, double Rmax, double DH_radial,
    int *d_out_count, double *d_out_gpox, double *d_out_lpox, int *d_out_sst,
    uint64_t *d_out_sort_key,
    int max_output
) {
    if (total_points == 0) return;
    cudaMemsetAsync(d_out_count, 0, sizeof(int), stream);
    int blockSize = 256;
    int gridSize = (total_points + blockSize - 1) / blockSize;
    setup_phase2_gpu_kernel<<<gridSize, blockSize, 0, stream>>>(
        total_points, d_block_offsets, d_block_shape, d_block_bbox, d_block_CDH,
        cs_w, Rmin, Rmax, DH_radial,
        d_out_count, d_out_gpox, d_out_lpox, d_out_sst, d_out_sort_key, max_output
    );
}

// ============================================================================
// Batched PACK Kernel
// Replaces per-pointstru per-var 1-thread kernel launches.
// Each thread processes one (point_index, var_index) pair.
// ============================================================================

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
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total = num_active * num_vars;
    if (tid >= total) return;

    int p = tid / num_vars;
    int v = tid % num_vars;
    int cache_idx = d_active_indices[p];

    int dumyd = d_dumyd[cache_idx];
    int ssst  = d_ssst[cache_idx];
    int blk   = d_src_block_idx[cache_idx];
    int sgfn  = d_var_fgfs[v];

    double *d_f = (double *)(d_block_fgfs[blk * fngfs_pitch + sgfn]);
    if (d_f == nullptr) { d_out[tid] = 0.0; return; }

    int ex[3]  = {d_block_shape[blk * 3 + 0], d_block_shape[blk * 3 + 1], d_block_shape[blk * 3 + 2]};
    int pex[3] = {d_block_shape_pad[blk * 3 + 0], d_block_shape_pad[blk * 3 + 1], d_block_shape_pad[blk * 3 + 2]};
    int inds[3]= {d_sind[cache_idx * 3 + 0], d_sind[cache_idx * 3 + 1], d_sind[cache_idx * 3 + 2]};
    const double *coef = d_coef + cache_idx * max_coef;
    double SoA[3] = {d_SoA[v * 3 + 0], d_SoA[v * 3 + 1], d_SoA[v * 3 + 2]};

    int DIMh = (dumyd == -1) ? dim : 1;
    double val;
    if (DIMh == 3) {
        val = interp3_indexed_device(ex, pex, d_f, inds, coef, ordn, SoA, ssst);
    } else {
        val = interp1_indexed_device(ex, pex, d_f, inds, coef, ordn, SoA, ssst, dumyd);
    }
    d_out[tid] = val;
}

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
) {
    if (num_active == 0 || num_vars == 0) return;
    int total = num_active * num_vars;
    dim3 block(256);
    dim3 grid((total + 255) / 256);
    batched_interpind_pack_kernel<<<grid, block, 0, stream>>>(
        num_active, num_vars, d_active_indices,
        d_coef, d_sind, d_dumyd, d_ssst, d_src_block_idx,
        max_coef, d_block_fgfs, d_block_shape, d_block_shape_pad,
        ordn, fngfs_pitch, d_SoA, d_var_fgfs, d_out
    );
}

// ============================================================================
// Batched UNPACK Kernel
// Writes directly from GPU buffer to d_fgfs — no D2H needed.
// Each thread handles one (point, var) pair.
// ============================================================================

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
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total = num_active * num_vars;
    if (tid >= total) return;

    int p = tid / num_vars;
    int v = tid % num_vars;
    int cache_idx = d_active_indices[p];

    int blk  = d_dst_block_idx[cache_idx];
    int sgfn = d_var_fgfs[v];

    double *d_data = (double *)(d_block_fgfs[blk * fngfs_pitch + sgfn]);
    if (d_data == nullptr) return;

    double lx = d_dst_lpox[cache_idx * 3 + 0];
    double ly = d_dst_lpox[cache_idx * 3 + 1];
    double lz = d_dst_lpox[cache_idx * 3 + 2];

    double llb0 = d_block_bbox[blk * 6 + 0];
    double llb1 = d_block_bbox[blk * 6 + 1];
    double llb2 = d_block_bbox[blk * 6 + 2];
    double uub0 = d_block_bbox[blk * 6 + 3];
    double uub1 = d_block_bbox[blk * 6 + 4];
    double uub2 = d_block_bbox[blk * 6 + 5];

    int ext0 = d_block_shape[blk * 3 + 0];
    int ext1 = d_block_shape[blk * 3 + 1];
    int ext2 = d_block_shape[blk * 3 + 2];
    int pext0 = d_block_shape_pad[blk * 3 + 0];
    int pext1 = d_block_shape_pad[blk * 3 + 1];

#ifdef Vertex
    double h0 = (uub0 - llb0) / (double)(ext0 - 1);
    double h1 = (uub1 - llb1) / (double)(ext1 - 1);
    double h2 = (uub2 - llb2) / (double)(ext2 - 1);
#else
    double h0 = (uub0 - llb0) / (double)ext0;
    double h1 = (uub1 - llb1) / (double)ext1;
    double h2 = (uub2 - llb2) / (double)ext2;
#endif
    int i = (int)((lx - llb0) / h0 + 0.4);
    int j = (int)((ly - llb1) / h1 + 0.4);
    int k = (int)((lz - llb2) / h2 + 0.4);

    if (i < 0 || i >= ext0 || j < 0 || j >= ext1 || k < 0 || k >= ext2) return;

    double val = d_src_buf[tid];
    d_data[(k * pext1 + j) * pext0 + i] = val;
}

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
) {
    if (num_active == 0 || num_vars == 0) return;
    int total = num_active * num_vars;
    dim3 block(256);
    dim3 grid((total + 255) / 256);
    batched_pointcopy_unpack_kernel<<<grid, block, 0, stream>>>(
        num_active, num_vars, d_active_indices,
        d_dst_lpox, d_dst_block_idx, d_block_fgfs,
        d_block_shape, d_block_shape_pad, d_block_bbox,
        fngfs_pitch, d_var_fgfs, d_src_buf
    );
}
