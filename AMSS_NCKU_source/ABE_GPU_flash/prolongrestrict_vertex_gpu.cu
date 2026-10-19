#include "prolongrestrict.h"

#include "macrodef.fh"

#ifdef Vertex

#include "fmisc.h"
#include "fmisc_gpu.cuh"

#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

constexpr int kTransferBlockX = 8;
constexpr int kTransferBlockY = 8;
constexpr int kTransferBlockZ = 4;

// ============================================================
// Vertex prolong3 coefficient tables (symmetric-pair interpolation)
// ============================================================
// Unlike Cell Lagrange interpolation, Vertex uses symmetric pair
// coefficients with NO even/odd reversal.

#if ghost_width == 2
__constant__ double C_PROLONG_VERTEX[2] = {
    -1.0 / 16.0,
     9.0 / 16.0
};
constexpr int kProlongPairs = 2;
#elif ghost_width == 3
__constant__ double C_PROLONG_VERTEX[3] = {
     3.0 / 256.0,
   -25.0 / 256.0,
    75.0 / 128.0
};
constexpr int kProlongPairs = 3;
#elif ghost_width == 4
__constant__ double C_PROLONG_VERTEX[4] = {
    -5.0 / 2048.0,
    49.0 / 2048.0,
  -245.0 / 2048.0,
  1225.0 / 2048.0
};
constexpr int kProlongPairs = 4;
#elif ghost_width == 5
__constant__ double C_PROLONG_VERTEX[5] = {
    35.0 / 65536.0,
  -405.0 / 65536.0,
   567.0 / 16384.0,
 -2205.0 / 16384.0,
 19845.0 / 32768.0
};
constexpr int kProlongPairs = 5;
#else
#error Unsupported ghost_width in prolongrestrict_vertex_gpu.cu; expected ghost_width == 2, 3, 4, or 5.
#endif
constexpr int kVertexStencilWidth = 2 * kProlongPairs;
constexpr int kVertexCenterIndex = ghost_width - 1;

__device__ __forceinline__ double vertex_prolong_source_value(
    const double* __restrict__ source,
    int ext0, int ext1, int ext2,
    int i, int j, int k,
    double SoA0, double SoA1, double SoA2
) {
    return d_symmetry_bd_0b(ghost_width, ext0, ext1, ext2, source,
                            i, j, k, SoA0, SoA1, SoA2);
}

__device__ __forceinline__ int vertex_prolong_coord(bool odd, int cxI, int offset) {
    return odd ? (cxI - 1) : (cxI - ghost_width + offset);
}

__device__ __forceinline__ double vertex_prolong_eval(
    int i, int j, int k,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int extc0, int extc1, int extc2,
    double SoA0, double SoA1, double SoA2,
    const double* __restrict__ source
) {
    int ii = i + lbf0;
    int jj = j + lbf1;
    int kk = k + lbf2;

    int cxI_i = (i + lbf0 + 1) / 2 - lbc0 + 1;
    int cxI_j = (j + lbf1 + 1) / 2 - lbc1 + 1;
    int cxI_k = (k + lbf2 + 1) / 2 - lbc2 + 1;

    bool odd_i = (ii % 2 == 1);
    bool odd_j = (jj % 2 == 1);
    bool odd_k = (kk % 2 == 1);

    double tmp2[kVertexStencilWidth][kVertexStencilWidth];
    for (int cj = 0; cj < kVertexStencilWidth; ++cj) {
        int c_idx_j = vertex_prolong_coord(odd_j, cxI_j, cj);
        for (int ci = 0; ci < kVertexStencilWidth; ++ci) {
            int c_idx_i = vertex_prolong_coord(odd_i, cxI_i, ci);
            if (odd_k) {
                int c_idx_k = cxI_k - 1;
                tmp2[cj][ci] = vertex_prolong_source_value(
                    source, extc0, extc1, extc2,
                    c_idx_i, c_idx_j, c_idx_k, SoA0, SoA1, SoA2);
            } else {
                double val = 0.0;
                #pragma unroll
                for (int p = 0; p < kProlongPairs; ++p) {
                    int left = p;
                    int right = kVertexStencilWidth - 1 - p;
                    int c_left_k = cxI_k - ghost_width + left;
                    int c_right_k = cxI_k - ghost_width + right;
                    double left_val = vertex_prolong_source_value(
                        source, extc0, extc1, extc2,
                        c_idx_i, c_idx_j, c_left_k, SoA0, SoA1, SoA2);
                    double right_val = vertex_prolong_source_value(
                        source, extc0, extc1, extc2,
                        c_idx_i, c_idx_j, c_right_k, SoA0, SoA1, SoA2);
                    val += C_PROLONG_VERTEX[p] * (left_val + right_val);
                }
                tmp2[cj][ci] = val;
            }
        }
    }

    double tmp1[kVertexStencilWidth];
    for (int ci = 0; ci < kVertexStencilWidth; ++ci) {
        if (odd_j) {
            tmp1[ci] = tmp2[kVertexCenterIndex][ci];
        } else {
            double val = 0.0;
            #pragma unroll
            for (int p = 0; p < kProlongPairs; ++p) {
                int left = p;
                int right = kVertexStencilWidth - 1 - p;
                val += C_PROLONG_VERTEX[p] * (tmp2[left][ci] + tmp2[right][ci]);
            }
            tmp1[ci] = val;
        }
    }

    if (odd_i) {
        return tmp1[kVertexCenterIndex];
    }

    double final_val = 0.0;
    #pragma unroll
    for (int p = 0; p < kProlongPairs; ++p) {
        int left = p;
        int right = kVertexStencilWidth - 1 - p;
        final_val += C_PROLONG_VERTEX[p] * (tmp1[left] + tmp1[right]);
    }
    return final_val;
}

// ============================================================
// Vertex restrict3: direct copy (coarse vertex <- coincident fine vertex)
// ============================================================
// In Vertex grids with 2:1 refinement, every coarse vertex coincides
// with a fine vertex at index (2*ic, 2*jc, 2*kc).  Restriction is
// simply:  func_coarse(ic,jc,kc) = funf_fine(2*ic, 2*jc, 2*kc)
//
// This is much simpler than Cell restrict which uses multi-pair
// stencil weighted averaging.

__global__ void restrict3_vertex_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int extc0, int extc1, int extc2,
    int extf0, int extf1, int extf2,
    const double* __restrict__ d_src_f,
    double* __restrict__ d_dst_c
) {
    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= ni || j_local >= nj || k_local >= nk) return;

    int i = i_start + i_local;
    int j = j_start + j_local;
    int k = k_start + k_local;

    int i_coarse = i + lbc0;
    int j_coarse = j + lbc1;
    int k_coarse = k + lbc2;

    int i_fine = 2 * i_coarse - lbf0;
    int j_fine = 2 * j_coarse - lbf1;
    int k_fine = 2 * k_coarse - lbf2;

    int fine_stride_xy = extf0 * extf1;
    int fine_stride_x  = extf0;
    int fine_idx = (k_fine - 1) * fine_stride_xy + (j_fine - 1) * fine_stride_x + (i_fine - 1);
    double val = d_src_f[fine_idx];

    int coarse_stride_xy = extc0 * extc1;
    int coarse_stride_x  = extc0;
    int out_idx = k * coarse_stride_xy + j * coarse_stride_x + i;
    d_dst_c[out_idx] = val;
}

// ============================================================
// Vertex prolong3: symmetric-pair interpolation
// ============================================================
// f = sum Ck * (f_{center-k} + f_{center+k-1}) for k = 0 .. kProlongPairs-1
// Stencil width = 2*kProlongPairs = 2*ghost_width, loaded via d_symmetry_bd_0b.
// No even/odd reversal (unlike Cell Lagrange interpolation).

__global__ void prolong3_vertex_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int extc0, int extc1, int extc2,
    int extf0, int extf1, int extf2,
    double SoA0, double SoA1, double SoA2,
    const double* __restrict__ d_src_c,
    double* __restrict__ d_dst_f
) {
    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= ni || j_local >= nj || k_local >= nk) return;

    int i = i_start + i_local;
    int j = j_start + j_local;
    int k = k_start + k_local;

    double final_val = vertex_prolong_eval(
        i, j, k,
        lbc0, lbc1, lbc2,
        lbf0, lbf1, lbf2,
        extc0, extc1, extc2,
        SoA0, SoA1, SoA2,
        d_src_c);

    int out_idx = k * extf0 * extf1 + j * extf0 + i;
    d_dst_f[out_idx] = final_val;
}

}  // anonymous namespace

// ============================================================
// Host-side: Vertex restrict3 launcher (single-transfer)
// ============================================================
// Uses the Fortran Vertex prolongation spacing: CD = (uubc-llbc)/(extc-1), FD = CD/2.
// and Vertex index bounds: ubf = idint(...)+0.4+1

void gpu_restrict3_launch(
    cudaStream_t stream,
    const double* d_src_f, double* d_dst_c,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
) {
    double FD[3] = {
        (uubf[0] - llbf[0]) / (double)(extf[0] - 1),
        (uubf[1] - llbf[1]) / (double)(extf[1] - 1),
        (uubf[2] - llbf[2]) / (double)(extf[2] - 1)
    };
    double CD[3] = {2.0 * FD[0], 2.0 * FD[1], 2.0 * FD[2]};

    // Base alignment (matching Fortran prolongrestrict_vertex.f90)
    double base[3];
    for (int d = 0; d < 3; ++d) {
        if (llbc[d] <= llbf[d]) {
            base[d] = llbc[d];
        } else {
            int jj = (int)std::trunc((llbc[d] - llbf[d]) / FD[d] + 0.4);
            if ((jj / 2) * 2 == jj) {
                base[d] = llbf[d];
            } else {
                base[d] = llbf[d] - CD[d] / 2.0;
            }
        }
    }

    // Index bounds (Vertex: +1 for upper bounds, matching Fortran)
    int lbf[3], ubf[3], lbc[3], ubc[3];
    for (int d = 0; d < 3; ++d) {
        lbf[d] = (int)std::trunc((llbf[d] - base[d]) / FD[d] + 0.4) + 1;
        ubf[d] = (int)std::trunc((uubf[d] - base[d]) / FD[d] + 0.4) + 1;
        lbc[d] = (int)std::trunc((llbc[d] - base[d]) / CD[d] + 0.4) + 1;
        ubc[d] = (int)std::trunc((uubc[d] - base[d]) / CD[d] + 0.4) + 1;
    }

    int lbr[3], ubr[3], lbrf[3], ubrf[3];
    for (int d = 0; d < 3; ++d) {
        lbr[d] = (int)std::trunc((llbt[d] - base[d]) / CD[d] + 0.4) + 1;
        ubr[d] = (int)std::trunc((uubt[d] - base[d]) / CD[d] + 0.4) + 1;
        lbrf[d] = (int)std::trunc((llbt[d] - base[d]) / FD[d] + 0.4) + 1;
        ubrf[d] = (int)std::trunc((uubt[d] - base[d]) / FD[d] + 0.4) + 1;
    }

    for (int d = 0; d < 3; ++d) {
        if (lbr[d] - lbc[d] < 0 ||
            ubr[d] - lbc[d] >= extc[d] ||
            lbrf[d] - lbf[d] < 0 ||
            ubrf[d] - lbf[d] >= extf[d]) {
            return;
        }
    }

    int i_start = lbr[0] - lbc[0];
    int j_start = lbr[1] - lbc[1];
    int k_start = lbr[2] - lbc[2];
    int i_end   = ubr[0] - lbc[0];
    int j_end   = ubr[1] - lbc[1];
    int k_end   = ubr[2] - lbc[2];

    int ni = i_end - i_start + 1;
    int nj = j_end - j_start + 1;
    int nk = k_end - k_start + 1;

    if (ni <= 0 || nj <= 0 || nk <= 0) return;

    dim3 block(kTransferBlockX, kTransferBlockY, kTransferBlockZ);
    dim3 grid((ni + block.x - 1) / block.x,
              (nj + block.y - 1) / block.y,
              (nk + block.z - 1) / block.z);

    restrict3_vertex_kernel<<<grid, block, 0, stream>>>(
        ni, nj, nk,
        i_start, j_start, k_start,
        lbc[0], lbc[1], lbc[2],
        lbf[0], lbf[1], lbf[2],
        extc[0], extc[1], extc[2],
        extf[0], extf[1], extf[2],
        d_src_f, d_dst_c
    );
}

// ============================================================
// Host-side: Vertex prolong3 launcher (single-transfer)
// ============================================================
// Uses the Fortran Vertex prolongation spacing: CD = (uubc-llbc)/(extc-1), FD = CD/2.
// Vertex transfer segments can span multi-box buffer ranges whose shape-derived spacing
// is not exactly CD/2, so batch builders must not infer FD from the target extent.

void gpu_prolong3_launch(
    cudaStream_t stream,
    const double* d_src_c, double* d_dst_f,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
) {
    double CD[3] = {
        (uubc[0] - llbc[0]) / (double)(extc[0] - 1),
        (uubc[1] - llbc[1]) / (double)(extc[1] - 1),
        (uubc[2] - llbc[2]) / (double)(extc[2] - 1)
    };
    double FD[3] = {CD[0] / 2.0, CD[1] / 2.0, CD[2] / 2.0};

    double base[3];
    for (int d = 0; d < 3; ++d) {
        if (llbc[d] <= llbf[d]) {
            base[d] = llbc[d];
        } else {
            int jj = (int)std::trunc((llbc[d] - llbf[d]) / FD[d] + 0.4);
            if ((jj / 2) * 2 == jj) {
                base[d] = llbf[d];
            } else {
                base[d] = llbf[d] - CD[d] / 2.0;
            }
        }
    }

    int lbf[3], ubf[3], lbc[3], ubc[3];
    for (int d = 0; d < 3; ++d) {
        lbf[d] = (int)std::trunc((llbf[d] - base[d]) / FD[d] + 0.4) + 1;
        ubf[d] = (int)std::trunc((uubf[d] - base[d]) / FD[d] + 0.4) + 1;
        lbc[d] = (int)std::trunc((llbc[d] - base[d]) / CD[d] + 0.4) + 1;
        ubc[d] = (int)std::trunc((uubc[d] - base[d]) / CD[d] + 0.4) + 1;
    }

    int lbp[3], ubp[3], lbpc[3], ubpc[3];
    for (int d = 0; d < 3; ++d) {
        lbp[d] = (int)std::trunc((llbt[d] - base[d]) / FD[d] + 0.4) + 1;
        ubp[d] = (int)std::trunc((uubt[d] - base[d]) / FD[d] + 0.4) + 1;
        lbpc[d] = (int)std::trunc((llbt[d] - base[d]) / CD[d] + 0.4) + 1;
        ubpc[d] = (int)std::trunc((uubt[d] - base[d]) / CD[d] + 0.4) + 1;
    }

    // Match Fortran prolong3 sanity checks: both target and coarse windows must fit.
    constexpr int kProlongUpperMargin = ghost_width - 1;
    for (int d = 0; d < 3; ++d) {
        if (lbp[d] - lbf[d] < 0 ||
            ubp[d] - lbf[d] >= extf[d] ||
            lbpc[d] - lbc[d] < 0 ||
            ubpc[d] - lbc[d] + 1 > extc[d] - kProlongUpperMargin) {
            return;
        }
    }

    int i_start = lbp[0] - lbf[0];
    int j_start = lbp[1] - lbf[1];
    int k_start = lbp[2] - lbf[2];
    int i_end   = ubp[0] - lbf[0];
    int j_end   = ubp[1] - lbf[1];
    int k_end   = ubp[2] - lbf[2];

    int ni = i_end - i_start + 1;
    int nj = j_end - j_start + 1;
    int nk = k_end - k_start + 1;

    if (ni <= 0 || nj <= 0 || nk <= 0) return;

    dim3 block(kTransferBlockX, kTransferBlockY, kTransferBlockZ);
    dim3 grid((ni + block.x - 1) / block.x,
              (nj + block.y - 1) / block.y,
              (nk + block.z - 1) / block.z);

    prolong3_vertex_kernel<<<grid, block, 0, stream>>>(
        ni, nj, nk,
        i_start, j_start, k_start,
        lbc[0], lbc[1], lbc[2],
        lbf[0], lbf[1], lbf[2],
        extc[0], extc[1], extc[2],
        extf[0], extf[1], extf[2],
        SoA[0], SoA[1], SoA[2],
        d_src_c, d_dst_f
    );
}

// ============================================================
// Vertex batch kernels + transfer builders
// ============================================================

__global__ void restrict3_vertex_batch_kernel(
    const GpuStencilTransfer* __restrict__ transfers,
    int transfer_count,
    int blocks_x
) {
    int transfer_index = blockIdx.x / blocks_x;
    if (transfer_index >= transfer_count) return;

    int local_block_x = blockIdx.x - transfer_index * blocks_x;
    const GpuStencilTransfer t = transfers[transfer_index];

    int i_local = local_block_x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= t.ni || j_local >= t.nj || k_local >= t.nk) return;

    int i = t.i_start + i_local;
    int j = t.j_start + j_local;
    int k = t.k_start + k_local;

    int i_coarse = i + t.lbc0;
    int j_coarse = j + t.lbc1;
    int k_coarse = k + t.lbc2;

    int i_fine = 2 * i_coarse - t.lbf0;
    int j_fine = 2 * j_coarse - t.lbf1;
    int k_fine = 2 * k_coarse - t.lbf2;

    int fine_idx = (k_fine - 1) * t.extf0 * t.extf1 + (j_fine - 1) * t.extf0 + (i_fine - 1);
    double val = t.source_field[fine_idx];

    int out_idx = k * t.extc0 * t.extc1 + j * t.extc0 + i;
    t.target_field[out_idx] = val;
}

__global__ void prolong3_vertex_batch_kernel(
    const GpuStencilTransfer* __restrict__ transfers,
    int transfer_count,
    int blocks_x
) {
    int transfer_index = blockIdx.x / blocks_x;
    if (transfer_index >= transfer_count) return;

    int local_block_x = blockIdx.x - transfer_index * blocks_x;
    const GpuStencilTransfer t = transfers[transfer_index];

    int i_local = local_block_x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= t.ni || j_local >= t.nj || k_local >= t.nk) return;

    int i = t.i_start + i_local;
    int j = t.j_start + j_local;
    int k = t.k_start + k_local;

    double final_val = vertex_prolong_eval(
        i, j, k,
        t.lbc0, t.lbc1, t.lbc2,
        t.lbf0, t.lbf1, t.lbf2,
        t.extc0, t.extc1, t.extc2,
        t.SoA0, t.SoA1, t.SoA2,
        t.source_field);

    int out_idx = k * t.extf0 * t.extf1 + j * t.extf0 + i;
    t.target_field[out_idx] = final_val;
}

// ============================================================
// Vertex transfer builders (Vertex spacing formulas)
// ============================================================

bool gpu_build_restrict3_transfer(
    const double* d_src_f, double* d_dst_c,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry,
    GpuStencilTransfer* transfer
) {
    double CD[3], FD[3], base[3];
    int lbc[3], lbf[3], lbr[3], ubr[3], lbrf[3], ubrf[3];

    for (int d = 0; d < 3; d++) {
        FD[d] = (uubf[d] - llbf[d]) / (double)(extf[d] - 1);
        CD[d] = 2.0 * FD[d];
        if (llbc[d] <= llbf[d]) {
            base[d] = llbc[d];
        } else {
            int j_val = (int)std::trunc((llbc[d] - llbf[d]) / FD[d] + 0.4);
            if ((j_val / 2) * 2 == j_val) base[d] = llbf[d];
            else base[d] = llbf[d] - CD[d] / 2.0;
        }

        lbc[d] = (int)std::trunc((llbc[d] - base[d]) / CD[d] + 0.4) + 1;
        lbr[d] = (int)std::trunc((llbt[d] - base[d]) / CD[d] + 0.4) + 1;
        ubr[d] = (int)std::trunc((uubt[d] - base[d]) / CD[d] + 0.4) + 1;
        lbf[d] = (int)std::trunc((llbf[d] - base[d]) / FD[d] + 0.4) + 1;
        lbrf[d] = (int)std::trunc((llbt[d] - base[d]) / FD[d] + 0.4) + 1;
        ubrf[d] = (int)std::trunc((uubt[d] - base[d]) / FD[d] + 0.4) + 1;
    }

    // Match Fortran restrict3 sanity checks for both the coarse output window
    // and the coincident fine-source vertices.
    for (int d = 0; d < 3; ++d) {
        if (lbr[d] - lbc[d] < 0 ||
            ubr[d] - lbc[d] >= extc[d] ||
            lbrf[d] - lbf[d] < 0 ||
            ubrf[d] - lbf[d] >= extf[d]) {
            return false;
        }
    }

    int i_start = lbr[0] - lbc[0];
    int j_start = lbr[1] - lbc[1];
    int k_start = lbr[2] - lbc[2];
    int i_end   = ubr[0] - lbc[0];
    int j_end   = ubr[1] - lbc[1];
    int k_end   = ubr[2] - lbc[2];

    int ni = i_end - i_start + 1;
    int nj = j_end - j_start + 1;
    int nk = k_end - k_start + 1;

    if (ni <= 0 || nj <= 0 || nk <= 0) return false;

    transfer->source_field = d_src_f;
    transfer->target_field = d_dst_c;
    transfer->ni = ni;
    transfer->nj = nj;
    transfer->nk = nk;
    transfer->i_start = i_start;
    transfer->j_start = j_start;
    transfer->k_start = k_start;
    transfer->lbc0 = lbc[0];  transfer->lbc1 = lbc[1];  transfer->lbc2 = lbc[2];
    transfer->lbf0 = lbf[0];  transfer->lbf1 = lbf[1];  transfer->lbf2 = lbf[2];
    transfer->extc0 = extc[0]; transfer->extc1 = extc[1]; transfer->extc2 = extc[2];
    transfer->extf0 = extf[0]; transfer->extf1 = extf[1]; transfer->extf2 = extf[2];
    transfer->SoA0 = SoA[0];  transfer->SoA1 = SoA[1];  transfer->SoA2 = SoA[2];
    transfer->Symmetry = Symmetry;
    return true;
}

bool gpu_build_prolong3_transfer(
    const double* d_src_c, double* d_dst_f,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry,
    GpuStencilTransfer* transfer
) {
    double CD[3], FD[3], base[3];
    int lbc[3], lbf[3], lbp[3], ubp[3], lbpc[3], ubpc[3];

    for (int d = 0; d < 3; d++) {
        CD[d] = (uubc[d] - llbc[d]) / (double)(extc[d] - 1);
        FD[d] = CD[d] / 2.0;
        if (llbc[d] <= llbf[d]) {
            base[d] = llbc[d];
        } else {
            int j_val = (int)std::trunc((llbc[d] - llbf[d]) / FD[d] + 0.4);
            if ((j_val / 2) * 2 == j_val) base[d] = llbf[d];
            else base[d] = llbf[d] - CD[d] / 2.0;
        }

        lbc[d] = (int)std::trunc((llbc[d] - base[d]) / CD[d] + 0.4) + 1;
        lbp[d] = (int)std::trunc((llbt[d] - base[d]) / FD[d] + 0.4) + 1;
        lbf[d] = (int)std::trunc((llbf[d] - base[d]) / FD[d] + 0.4) + 1;
        ubp[d] = (int)std::trunc((uubt[d] - base[d]) / FD[d] + 0.4) + 1;
        lbpc[d] = (int)std::trunc((llbt[d] - base[d]) / CD[d] + 0.4) + 1;
        ubpc[d] = (int)std::trunc((uubt[d] - base[d]) / CD[d] + 0.4) + 1;
    }

    // Match Fortran prolong3 sanity checks: both target and coarse windows must fit.
    constexpr int kProlongUpperMargin = ghost_width - 1;
    for (int d = 0; d < 3; ++d) {
        if (lbp[d] - lbf[d] < 0 ||
            ubp[d] - lbf[d] >= extf[d] ||
            lbpc[d] - lbc[d] < 0 ||
            ubpc[d] - lbc[d] + 1 > extc[d] - kProlongUpperMargin) {
            return false;
        }
    }

    int i_start = lbp[0] - lbf[0];
    int j_start = lbp[1] - lbf[1];
    int k_start = lbp[2] - lbf[2];
    int i_end   = ubp[0] - lbf[0];
    int j_end   = ubp[1] - lbf[1];
    int k_end   = ubp[2] - lbf[2];

    int ni = i_end - i_start + 1;
    int nj = j_end - j_start + 1;
    int nk = k_end - k_start + 1;

    if (ni <= 0 || nj <= 0 || nk <= 0) return false;

    transfer->source_field = d_src_c;
    transfer->target_field = d_dst_f;
    transfer->ni = ni;
    transfer->nj = nj;
    transfer->nk = nk;
    transfer->i_start = i_start;
    transfer->j_start = j_start;
    transfer->k_start = k_start;
    transfer->lbc0 = lbc[0];  transfer->lbc1 = lbc[1];  transfer->lbc2 = lbc[2];
    transfer->lbf0 = lbf[0];  transfer->lbf1 = lbf[1];  transfer->lbf2 = lbf[2];
    transfer->extc0 = extc[0]; transfer->extc1 = extc[1]; transfer->extc2 = extc[2];
    transfer->extf0 = extf[0]; transfer->extf1 = extf[1]; transfer->extf2 = extf[2];
    transfer->SoA0 = SoA[0];  transfer->SoA1 = SoA[1];  transfer->SoA2 = SoA[2];
    transfer->Symmetry = Symmetry;
    return true;
}

struct StencilTransferScratch { cudaStream_t stream; GpuStencilTransfer* device_transfers; int capacity; };
static int block_count_for(int item_count, int block_extent) { return (item_count + block_extent - 1) / block_extent; }
static void max_transfer_grid_blocks(const GpuStencilTransfer* transfers, int transfer_count,
                                     dim3 block, int* blocks_x, int* blocks_y, int* blocks_z) {
    *blocks_x = 1;
    *blocks_y = 1;
    *blocks_z = 1;
    for (int i = 0; i < transfer_count; ++i) {
        *blocks_x = std::max(*blocks_x, block_count_for(transfers[i].ni, block.x));
        *blocks_y = std::max(*blocks_y, block_count_for(transfers[i].nj, block.y));
        *blocks_z = std::max(*blocks_z, block_count_for(transfers[i].nk, block.z));
    }
}
static StencilTransferScratch& transfer_scratch_for_stream(std::vector<StencilTransferScratch>& scratches, cudaStream_t stream) { for (size_t i = 0; i < scratches.size(); ++i) if (scratches[i].stream == stream) return scratches[i]; StencilTransferScratch s = {stream, nullptr, 0}; scratches.push_back(s); return scratches.back(); }
static std::vector<StencilTransferScratch>& restrict_transfer_scratches() { static std::vector<StencilTransferScratch> scratches; return scratches; }
static std::vector<StencilTransferScratch>& prolong_transfer_scratches() { static std::vector<StencilTransferScratch> scratches; return scratches; }
static void ensure_stencil_transfer_capacity(StencilTransferScratch& scratch, int n) { if (n <= scratch.capacity) return; if (scratch.device_transfers) cudaFree(scratch.device_transfers); cudaMalloc(&scratch.device_transfers, n * sizeof(GpuStencilTransfer)); scratch.capacity = n; }

void gpu_restrict3_batch_launch(
    cudaStream_t stream,
    const GpuStencilTransfer* transfers,
    int transfer_count
) {
    if (transfer_count <= 0) return;
    StencilTransferScratch& scratch = transfer_scratch_for_stream(restrict_transfer_scratches(), stream);
    ensure_stencil_transfer_capacity(scratch, transfer_count);
    cudaMemcpyAsync(scratch.device_transfers, transfers, transfer_count * sizeof(GpuStencilTransfer), cudaMemcpyHostToDevice, stream);
    dim3 block(kTransferBlockX, kTransferBlockY, kTransferBlockZ);
    int blocks_x = 1;
    int blocks_y = 1;
    int blocks_z = 1;
    max_transfer_grid_blocks(transfers, transfer_count, block, &blocks_x, &blocks_y, &blocks_z);
    dim3 grid(blocks_x * transfer_count, blocks_y, blocks_z);
    restrict3_vertex_batch_kernel<<<grid, block, 0, stream>>>(scratch.device_transfers, transfer_count, blocks_x);
}

void gpu_prolong3_batch_launch(
    cudaStream_t stream,
    const GpuStencilTransfer* transfers,
    int transfer_count
) {
    if (transfer_count <= 0) return;
    StencilTransferScratch& scratch = transfer_scratch_for_stream(prolong_transfer_scratches(), stream);
    ensure_stencil_transfer_capacity(scratch, transfer_count);
    cudaMemcpyAsync(scratch.device_transfers, transfers, transfer_count * sizeof(GpuStencilTransfer), cudaMemcpyHostToDevice, stream);
    dim3 block(kTransferBlockX, kTransferBlockY, kTransferBlockZ);
    int blocks_x = 1;
    int blocks_y = 1;
    int blocks_z = 1;
    max_transfer_grid_blocks(transfers, transfer_count, block, &blocks_x, &blocks_y, &blocks_z);
    dim3 grid(blocks_x * transfer_count, blocks_y, blocks_z);
    prolong3_vertex_batch_kernel<<<grid, block, 0, stream>>>(scratch.device_transfers, transfer_count, blocks_x);
}

#endif  // Vertex
