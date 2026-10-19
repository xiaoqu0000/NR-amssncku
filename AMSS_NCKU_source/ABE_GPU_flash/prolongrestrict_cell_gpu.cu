#include "prolongrestrict.h"

#include "macrodef.fh"

#include "fmisc.h"
#include "fmisc_gpu.cuh"

#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <vector>

#ifdef Cell

namespace {

constexpr int kTransferBlockX = 8;
constexpr int kTransferBlockY = 8;
#if ghost_width == 5
constexpr int kTransferBlockZ = 1;
#elif ghost_width == 4
constexpr int kTransferBlockZ = 2;
#else
constexpr int kTransferBlockZ = 4;
#endif
constexpr int kTransferThreadCount = kTransferBlockX * kTransferBlockY * kTransferBlockZ;

#if ghost_width == 2
constexpr int kProlongStencilWidth = 4;
constexpr int kProlongSmemX = 8;
constexpr int kProlongSmemY = 8;
constexpr int kProlongSmemZ = 6;
constexpr int kRestrictStencilWidth = 4;
constexpr int kRestrictSmemX = 18;
constexpr int kRestrictSmemY = 18;
constexpr int kRestrictSmemZ = 10;
#elif ghost_width == 3
constexpr int kProlongStencilWidth = 6;
constexpr int kProlongSmemX = 10;
constexpr int kProlongSmemY = 10;
constexpr int kProlongSmemZ = 8;
constexpr int kRestrictStencilWidth = 6;
constexpr int kRestrictSmemX = 20;
constexpr int kRestrictSmemY = 20;
constexpr int kRestrictSmemZ = 12;
#elif ghost_width == 4
constexpr int kProlongStencilWidth = 8;
constexpr int kProlongSmemX = 12;
constexpr int kProlongSmemY = 12;
constexpr int kProlongSmemZ = 10;
constexpr int kRestrictStencilWidth = 8;
constexpr int kRestrictSmemX = 22;
constexpr int kRestrictSmemY = 22;
constexpr int kRestrictSmemZ = 10;
#elif ghost_width == 5
constexpr int kProlongStencilWidth = 10;
constexpr int kProlongSmemX = 14;
constexpr int kProlongSmemY = 14;
constexpr int kProlongSmemZ = 10;
constexpr int kRestrictStencilWidth = 10;
constexpr int kRestrictSmemX = 24;
constexpr int kRestrictSmemY = 24;
constexpr int kRestrictSmemZ = 10;
#else
#error Unsupported ghost_width in prolongrestrict_cell_gpu.cu; expected 2, 3, 4, or 5.
#endif

constexpr int kProlongSmemXY = kProlongSmemX * kProlongSmemY;
constexpr int kProlongSmemTotal = kProlongSmemZ * kProlongSmemXY;
constexpr int kRestrictSmemXY = kRestrictSmemX * kRestrictSmemY;
constexpr int kRestrictSmemTotal = kRestrictSmemZ * kRestrictSmemXY;

struct StencilTransferScratch {
    cudaStream_t stream;
    GpuStencilTransfer* device_transfers;
    int capacity;
};

inline int block_count_for(int item_count, int block_extent) {
    return (item_count + block_extent - 1) / block_extent;
}

StencilTransferScratch& transfer_scratch_for_stream(std::vector<StencilTransferScratch>& scratches, cudaStream_t stream) {
    for (std::vector<StencilTransferScratch>::iterator scratch = scratches.begin(); scratch != scratches.end(); ++scratch) {
        if (scratch->stream == stream) {
            return *scratch;
        }
    }

    StencilTransferScratch scratch = {stream, nullptr, 0};
    scratches.push_back(scratch);
    return scratches.back();
}

std::vector<StencilTransferScratch>& prolong_transfer_scratches() {
    static std::vector<StencilTransferScratch> scratches;
    return scratches;
}

std::vector<StencilTransferScratch>& restrict_transfer_scratches() {
    static std::vector<StencilTransferScratch> scratches;
    return scratches;
}

void ensure_stencil_transfer_capacity(StencilTransferScratch& scratch, int transfer_count) {
    if (transfer_count <= scratch.capacity) {
        return;
    }

    if (scratch.device_transfers) {
        cudaFree(scratch.device_transfers);
    }

    cudaMalloc(&scratch.device_transfers, transfer_count * sizeof(GpuStencilTransfer));
    scratch.capacity = transfer_count;
}

}  // namespace

// Prolongation Coefficients
#if ghost_width == 2
__constant__ double C_PROLONG[10] = {
    -7.0 / 128.0,     // C1
    105.0 / 128.0,    // C2
    35.0 / 128.0,     // C3
    -5.0 / 128.0,     // C4
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
};
#elif ghost_width == 3
__constant__ double C_PROLONG[10] = {
    77.0 / 8192.0,    // C1
    -693.0 / 8192.0,  // C2
    3465.0 / 4096.0,  // C3
    1155.0 / 4096.0,  // C4
    -495.0 / 8192.0,  // C5
    63.0 / 8192.0,    // C6
    0.0,
    0.0,
    0.0,
    0.0
};
#elif ghost_width == 4
__constant__ double C_PROLONG[10] = {
    -495.0 / 262144.0,    // C1
    5005.0 / 262144.0,    // C2
    -27027.0 / 262144.0,  // C3
    225225.0 / 262144.0,  // C4
    75075.0 / 262144.0,   // C5
    -19305.0 / 262144.0,  // C6
    4095.0 / 262144.0,    // C7
    -429.0 / 262144.0,    // C8
    0.0,
    0.0
};
#elif ghost_width == 5
__constant__ double C_PROLONG[10] = {
    13585.0 / 33554432.0,    // C1
    -159885.0 / 33554432.0,  // C2
    230945.0 / 8388608.0,    // C3
    -969969.0 / 8388608.0,   // C4
    14549535.0 / 16777216.0, // C5
    4849845.0 / 16777216.0,  // C6
    -692835.0 / 8388608.0,   // C7
    188955.0 / 8388608.0,    // C8
    -138567.0 / 33554432.0,  // C9
    12155.0 / 33554432.0     // C10
};
#endif

// Restriction Coefficients
#if ghost_width == 2
__constant__ double C_RESTRICT[5] = {
    -1.0 / 16.0,      // C1
    9.0 / 16.0,       // C2
    0.0,
    0.0,
    0.0
};
#elif ghost_width == 3
__constant__ double C_RESTRICT[5] = {
    3.0 / 256.0,      // C1
    -25.0 / 256.0,    // C2
    75.0 / 128.0,     // C3
    0.0,
    0.0
};
#elif ghost_width == 4
__constant__ double C_RESTRICT[5] = {
    -5.0 / 2048.0,    // C1
    49.0 / 2048.0,    // C2
    -245.0 / 2048.0,  // C3
    1225.0 / 2048.0,  // C4
    0.0
};
#elif ghost_width == 5
__constant__ double C_RESTRICT[5] = {
    35.0 / 65536.0,     // C1
    -405.0 / 65536.0,   // C2
    567.0 / 16384.0,    // C3
    -2205.0 / 16384.0,  // C4
    19845.0 / 32768.0   // C5
};
#endif

__global__ void prolong3_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int extc0, int extc1, int extc2,
    int extf0, int extf1, int extf2,
    const double* __restrict__ d_src_c,
    double* __restrict__ d_dst_f,
    double SoA0, double SoA1, double SoA2,
    int Symmetry
) {
    __shared__ double smem[kProlongSmemZ][kProlongSmemY][kProlongSmemX];

    // 1D Block 内部线程 ID，用于协同加载
    int tid = threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;
    
    // 计算当前 Block 在细网格全局空间的起始索引
    int block_i_start = i_start + blockIdx.x * blockDim.x;
    int block_j_start = j_start + blockIdx.y * blockDim.y;
    int block_k_start = k_start + blockIdx.z * blockDim.z;

    constexpr int gw_skip = ghost_width - 1;
    int base_c_i = (block_i_start + lbf0) / 2 - lbc0 - gw_skip;
    int base_c_j = (block_j_start + lbf1) / 2 - lbc1 - gw_skip;
    int base_c_k = (block_k_start + lbf2) / 2 - lbc2 - gw_skip;

    for (int idx = tid; idx < kProlongSmemTotal; idx += kTransferThreadCount) {
        // 1D 索引还原为 Tile 内部的 3D 偏移
        int loc_k = idx / kProlongSmemXY;
        int rem   = idx % kProlongSmemXY;
        int loc_j = rem / kProlongSmemX;
        int loc_i = rem % kProlongSmemX;

        // 调用对称取值函数，直接从 Global Memory 存入 Shared Memory
        smem[loc_k][loc_j][loc_i] = d_symmetry_bd_0b(
            ghost_width, extc0, extc1, extc2, d_src_c,
            base_c_i + loc_i, base_c_j + loc_j, base_c_k + loc_k,
            SoA0, SoA1, SoA2
        );
    }

    __syncthreads();

    // ==========================================
    // 3. 计算阶段 (纯 Shared Memory 读取)
    // ==========================================
    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= ni || j_local >= nj || k_local >= nk) return;

    int i = i_start + i_local;
    int j = j_start + j_local;
    int k = k_start + k_local;

    int ii = i + lbf0;
    int jj = j + lbf1;
    int kk = k + lbf2;

    int cxI_i_0b = (i + lbf0) / 2 - lbc0; 
    int cxI_j_0b = (j + lbf1) / 2 - lbc1;
    int cxI_k_0b = (k + lbf2) / 2 - lbc2;

    bool k_even = ((kk / 2) * 2 == kk);
    bool j_even = ((jj / 2) * 2 == jj);
    bool i_even = ((ii / 2) * 2 == ii);

    // 计算当前线程需要的粗网格点在 Shared Memory 中的起始索引
    int smem_i_base = cxI_i_0b - gw_skip - base_c_i;
    int smem_j_base = cxI_j_0b - gw_skip - base_c_j;
    int smem_k_base = cxI_k_0b - gw_skip - base_c_k;

    double tmp2[kProlongStencilWidth][kProlongStencilWidth];
    double tmp1[kProlongStencilWidth];

    // Z方向插值
    for (int m = 0; m < kProlongStencilWidth; m++) {
        for (int n = 0; n < kProlongStencilWidth; n++) {
            int cur_i_smem = smem_i_base + n;
            int cur_j_smem = smem_j_base + m;
            double val = 0.0;
            
            #pragma unroll
            for (int step = 0; step < kProlongStencilWidth; step++) {
                int c_idx = k_even ? step : (kProlongStencilWidth - 1 - step);
                val += C_PROLONG[c_idx] * smem[smem_k_base + step][cur_j_smem][cur_i_smem];
            }
            tmp2[m][n] = val;
        }
    }

    // Y方向插值
    for (int n = 0; n < kProlongStencilWidth; n++) {
        double val = 0.0;
        #pragma unroll
        for (int step = 0; step < kProlongStencilWidth; step++) {
            int c_idx = j_even ? step : (kProlongStencilWidth - 1 - step);
            val += C_PROLONG[c_idx] * tmp2[step][n];
        }
        tmp1[n] = val;
    }

    // X方向插值
    double final_val = 0.0;
    #pragma unroll
    for (int step = 0; step < kProlongStencilWidth; step++) {
        int c_idx = i_even ? step : (kProlongStencilWidth - 1 - step);
        final_val += C_PROLONG[c_idx] * tmp1[step];
    }

    // 写入细网格显存
    int out_idx = k * (extf0 * extf1) + j * extf0 + i;
    d_dst_f[out_idx] = final_val;
}

__global__ void prolong3_batch_kernel(
    const GpuStencilTransfer* __restrict__ transfers,
    int transfer_count,
    int blocks_x
) {
    const int transfer_index = blockIdx.x / blocks_x;
    if (transfer_index >= transfer_count) return;

    const int local_block_x = blockIdx.x - transfer_index * blocks_x;
    const GpuStencilTransfer transfer = transfers[transfer_index];

    __shared__ double smem[kProlongSmemZ][kProlongSmemY][kProlongSmemX];

    int tid = threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;

    int block_i_start = transfer.i_start + local_block_x * blockDim.x;
    int block_j_start = transfer.j_start + blockIdx.y * blockDim.y;
    int block_k_start = transfer.k_start + blockIdx.z * blockDim.z;

    constexpr int gw_skip = ghost_width - 1;
    int base_c_i = (block_i_start + transfer.lbf0) / 2 - transfer.lbc0 - gw_skip;
    int base_c_j = (block_j_start + transfer.lbf1) / 2 - transfer.lbc1 - gw_skip;
    int base_c_k = (block_k_start + transfer.lbf2) / 2 - transfer.lbc2 - gw_skip;

    for (int idx = tid; idx < kProlongSmemTotal; idx += kTransferThreadCount) {
        int loc_k = idx / kProlongSmemXY;
        int rem   = idx % kProlongSmemXY;
        int loc_j = rem / kProlongSmemX;
        int loc_i = rem % kProlongSmemX;

        smem[loc_k][loc_j][loc_i] = d_symmetry_bd_0b(
            ghost_width, transfer.extc0, transfer.extc1, transfer.extc2, transfer.source_field,
            base_c_i + loc_i, base_c_j + loc_j, base_c_k + loc_k,
            transfer.SoA0, transfer.SoA1, transfer.SoA2
        );
    }

    __syncthreads();

    int i_local = local_block_x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= transfer.ni || j_local >= transfer.nj || k_local >= transfer.nk) return;

    int i = transfer.i_start + i_local;
    int j = transfer.j_start + j_local;
    int k = transfer.k_start + k_local;

    int ii = i + transfer.lbf0;
    int jj = j + transfer.lbf1;
    int kk = k + transfer.lbf2;

    int cxI_i_0b = (i + transfer.lbf0) / 2 - transfer.lbc0;
    int cxI_j_0b = (j + transfer.lbf1) / 2 - transfer.lbc1;
    int cxI_k_0b = (k + transfer.lbf2) / 2 - transfer.lbc2;

    bool k_even = ((kk / 2) * 2 == kk);
    bool j_even = ((jj / 2) * 2 == jj);
    bool i_even = ((ii / 2) * 2 == ii);

    int smem_i_base = cxI_i_0b - gw_skip - base_c_i;
    int smem_j_base = cxI_j_0b - gw_skip - base_c_j;
    int smem_k_base = cxI_k_0b - gw_skip - base_c_k;

    double tmp2[kProlongStencilWidth][kProlongStencilWidth];
    double tmp1[kProlongStencilWidth];

    for (int m = 0; m < kProlongStencilWidth; m++) {
        for (int n = 0; n < kProlongStencilWidth; n++) {
            int cur_i_smem = smem_i_base + n;
            int cur_j_smem = smem_j_base + m;
            double val = 0.0;

            #pragma unroll
            for (int step = 0; step < kProlongStencilWidth; step++) {
                int c_idx = k_even ? step : (kProlongStencilWidth - 1 - step);
                val += C_PROLONG[c_idx] * smem[smem_k_base + step][cur_j_smem][cur_i_smem];
            }
            tmp2[m][n] = val;
        }
    }

    for (int n = 0; n < kProlongStencilWidth; n++) {
        double val = 0.0;
        #pragma unroll
        for (int step = 0; step < kProlongStencilWidth; step++) {
            int c_idx = j_even ? step : (kProlongStencilWidth - 1 - step);
            val += C_PROLONG[c_idx] * tmp2[step][n];
        }
        tmp1[n] = val;
    }

    double final_val = 0.0;
    #pragma unroll
    for (int step = 0; step < kProlongStencilWidth; step++) {
        int c_idx = i_even ? step : (kProlongStencilWidth - 1 - step);
        final_val += C_PROLONG[c_idx] * tmp1[step];
    }

    int out_idx = k * (transfer.extf0 * transfer.extf1) + j * transfer.extf0 + i;
    transfer.target_field[out_idx] = final_val;
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
    int lbc[3], lbf[3], lbp[3], ubp[3];

    for (int d = 0; d < 3; d++) {
        CD[d] = (uubc[d] - llbc[d]) / (double)extc[d];
        FD[d] = (uubf[d] - llbf[d]) / (double)extf[d];
        if (llbc[d] <= llbf[d]) {
            base[d] = llbc[d];
        } else {
            int j_val = (int)std::trunc((llbc[d] - llbf[d]) / FD[d] + 0.4);
            if ((j_val / 2) * 2 == j_val) base[d] = llbf[d];
            else base[d] = llbf[d] - CD[d] / 2.0;
        }

        lbf[d] = (int)std::trunc((llbf[d] - base[d]) / FD[d] + 0.4) + 1;
        lbc[d] = (int)std::trunc((llbc[d] - base[d]) / CD[d] + 0.4) + 1;
        lbp[d] = (int)std::trunc((llbt[d] - base[d]) / FD[d] + 0.4) + 1;
        ubp[d] = (int)std::trunc((uubt[d] - base[d]) / FD[d] + 0.4);
    }

    int i_start = lbp[0] - lbf[0];
    int i_end   = ubp[0] - lbf[0];
    int j_start = lbp[1] - lbf[1];
    int j_end   = ubp[1] - lbf[1];
    int k_start = lbp[2] - lbf[2];
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
    transfer->lbc0 = lbc[0];
    transfer->lbc1 = lbc[1];
    transfer->lbc2 = lbc[2];
    transfer->lbf0 = lbf[0];
    transfer->lbf1 = lbf[1];
    transfer->lbf2 = lbf[2];
    transfer->extc0 = extc[0];
    transfer->extc1 = extc[1];
    transfer->extc2 = extc[2];
    transfer->extf0 = extf[0];
    transfer->extf1 = extf[1];
    transfer->extf2 = extf[2];
    transfer->SoA0 = SoA[0];
    transfer->SoA1 = SoA[1];
    transfer->SoA2 = SoA[2];
    transfer->Symmetry = Symmetry;
    return true;
}

void gpu_prolong3_launch(
    cudaStream_t stream,
    const double* d_src_c, double* d_dst_f,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
) {
    GpuStencilTransfer transfer;
    if (!gpu_build_prolong3_transfer(
            d_src_c,
            d_dst_f,
            llbc,
            uubc,
            extc,
            llbf,
            uubf,
            extf,
            llbt,
            uubt,
            SoA,
            Symmetry,
            &transfer
        )) {
        return;
    }

    dim3 block(kTransferBlockX, kTransferBlockY, kTransferBlockZ);
    dim3 grid(block_count_for(transfer.ni, block.x),
              block_count_for(transfer.nj, block.y),
              block_count_for(transfer.nk, block.z));

    prolong3_kernel<<<grid, block, 0, stream>>>(
        transfer.ni, transfer.nj, transfer.nk,
        transfer.i_start, transfer.j_start, transfer.k_start,
        transfer.lbc0, transfer.lbc1, transfer.lbc2,
        transfer.lbf0, transfer.lbf1, transfer.lbf2,
        transfer.extc0, transfer.extc1, transfer.extc2,
        transfer.extf0, transfer.extf1, transfer.extf2,
        transfer.source_field, transfer.target_field,
        transfer.SoA0, transfer.SoA1, transfer.SoA2,
        transfer.Symmetry
    );
}

void gpu_prolong3_batch_launch(
    cudaStream_t stream,
    const GpuStencilTransfer* transfers,
    int transfer_count
) {
    if (transfer_count <= 0) return;

    StencilTransferScratch& scratch = transfer_scratch_for_stream(prolong_transfer_scratches(), stream);
    ensure_stencil_transfer_capacity(scratch, transfer_count);
    cudaMemcpyAsync(
        scratch.device_transfers,
        transfers,
        transfer_count * sizeof(GpuStencilTransfer),
        cudaMemcpyHostToDevice,
        stream
    );

    dim3 block(kTransferBlockX, kTransferBlockY, kTransferBlockZ);
    const int blocks_x = block_count_for(transfers[0].ni, block.x);
    dim3 grid(blocks_x * transfer_count,
              block_count_for(transfers[0].nj, block.y),
              block_count_for(transfers[0].nk, block.z));
    prolong3_batch_kernel<<<grid, block, 0, stream>>>(scratch.device_transfers, transfer_count, blocks_x);
}

__global__ void restrict3_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int extc0, int extc1, int extc2,
    int extf0, int extf1, int extf2,
    const double* __restrict__ d_src_f,
    double* __restrict__ d_dst_c,
    double SoA0, double SoA1, double SoA2,
    int Symmetry
) {
    __shared__ double smem[kRestrictSmemZ][kRestrictSmemY][kRestrictSmemX];

    int tid = threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;

    // 计算当前 Block 对应的全局粗网格起始索引
    int block_i_start = i_start + blockIdx.x * blockDim.x;
    int block_j_start = j_start + blockIdx.y * blockDim.y;
    int block_k_start = k_start + blockIdx.z * blockDim.z;

    int base_f_i = 2 * (block_i_start + lbc0) - lbf0 - ghost_width;
    int base_f_j = 2 * (block_j_start + lbc1) - lbf1 - ghost_width;
    int base_f_k = 2 * (block_k_start + lbc2) - lbf2 - ghost_width;

    for (int idx = tid; idx < kRestrictSmemTotal; idx += kTransferThreadCount) {
        int loc_k = idx / kRestrictSmemXY;
        int rem   = idx % kRestrictSmemXY;
        int loc_j = rem / kRestrictSmemX;
        int loc_i = rem % kRestrictSmemX;

        smem[loc_k][loc_j][loc_i] = d_symmetry_bd_0b(
            ghost_width - 1, extf0, extf1, extf2, d_src_f,
            base_f_i + loc_i, base_f_j + loc_j, base_f_k + loc_k,
            SoA0, SoA1, SoA2
        );
    }

    __syncthreads();

    // ==========================================
    // 3. 计算阶段
    // ==========================================
    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    // 边界检查必须放在 syncthreads 之后
    if (i_local >= ni || j_local >= nj || k_local >= nk) return;

    int i = i_start + i_local;
    int j = j_start + j_local;
    int k = k_start + k_local;

    int smem_i_start = 2 * threadIdx.x;
    int smem_j_start = 2 * threadIdx.y;
    int smem_k_start = 2 * threadIdx.z;

    double tmp2[kRestrictStencilWidth][kRestrictStencilWidth];
    double tmp1[kRestrictStencilWidth];

    // Z-Direction Restriction
    for (int m = 0; m < kRestrictStencilWidth; m++) {
        for (int n = 0; n < kRestrictStencilWidth; n++) {
            int cur_j_smem = smem_j_start + m;
            int cur_i_smem = smem_i_start + n;
            
            double val = 0.0;
            #pragma unroll
            for (int r = 0; r < kRestrictStencilWidth / 2; r++) {
                val += C_RESTRICT[r] * (
                    smem[smem_k_start + r][cur_j_smem][cur_i_smem] +
                    smem[smem_k_start + kRestrictStencilWidth - 1 - r][cur_j_smem][cur_i_smem]
                );
            }
            
            tmp2[m][n] = val;
        }
    }

    // Y-Direction Restriction
    #pragma unroll
    for (int n = 0; n < kRestrictStencilWidth; n++) {
        double val = 0.0;
        #pragma unroll
        for (int r = 0; r < kRestrictStencilWidth / 2; r++) {
            val += C_RESTRICT[r] * (tmp2[r][n] + tmp2[kRestrictStencilWidth - 1 - r][n]);
        }
        tmp1[n] = val;
    }

    // X-Direction Restriction
    double final_val = 0.0;
    #pragma unroll
    for (int r = 0; r < kRestrictStencilWidth / 2; r++) {
        final_val += C_RESTRICT[r] * (tmp1[r] + tmp1[kRestrictStencilWidth - 1 - r]);
    }

    // 写入目标粗网格显存
    int out_idx = k * (extc0 * extc1) + j * extc0 + i;
    d_dst_c[out_idx] = final_val;
}

__global__ void restrict3_batch_kernel(
    const GpuStencilTransfer* __restrict__ transfers,
    int transfer_count,
    int blocks_x
) {
    const int transfer_index = blockIdx.x / blocks_x;
    if (transfer_index >= transfer_count) return;

    const int local_block_x = blockIdx.x - transfer_index * blocks_x;
    const GpuStencilTransfer transfer = transfers[transfer_index];

    __shared__ double smem[kRestrictSmemZ][kRestrictSmemY][kRestrictSmemX];

    int tid = threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;

    int block_i_start = transfer.i_start + local_block_x * blockDim.x;
    int block_j_start = transfer.j_start + blockIdx.y * blockDim.y;
    int block_k_start = transfer.k_start + blockIdx.z * blockDim.z;

    int base_f_i = 2 * (block_i_start + transfer.lbc0) - transfer.lbf0 - ghost_width;
    int base_f_j = 2 * (block_j_start + transfer.lbc1) - transfer.lbf1 - ghost_width;
    int base_f_k = 2 * (block_k_start + transfer.lbc2) - transfer.lbf2 - ghost_width;

    for (int idx = tid; idx < kRestrictSmemTotal; idx += kTransferThreadCount) {
        int loc_k = idx / kRestrictSmemXY;
        int rem   = idx % kRestrictSmemXY;
        int loc_j = rem / kRestrictSmemX;
        int loc_i = rem % kRestrictSmemX;

        smem[loc_k][loc_j][loc_i] = d_symmetry_bd_0b(
            ghost_width - 1, transfer.extf0, transfer.extf1, transfer.extf2, transfer.source_field,
            base_f_i + loc_i, base_f_j + loc_j, base_f_k + loc_k,
            transfer.SoA0, transfer.SoA1, transfer.SoA2
        );
    }

    __syncthreads();

    int i_local = local_block_x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= transfer.ni || j_local >= transfer.nj || k_local >= transfer.nk) return;

    int i = transfer.i_start + i_local;
    int j = transfer.j_start + j_local;
    int k = transfer.k_start + k_local;

    int smem_i_start = 2 * threadIdx.x;
    int smem_j_start = 2 * threadIdx.y;
    int smem_k_start = 2 * threadIdx.z;

    double tmp2[kRestrictStencilWidth][kRestrictStencilWidth];
    double tmp1[kRestrictStencilWidth];

    for (int m = 0; m < kRestrictStencilWidth; m++) {
        for (int n = 0; n < kRestrictStencilWidth; n++) {
            int cur_j_smem = smem_j_start + m;
            int cur_i_smem = smem_i_start + n;

            double val = 0.0;
            #pragma unroll
            for (int r = 0; r < kRestrictStencilWidth / 2; r++) {
                val += C_RESTRICT[r] * (
                    smem[smem_k_start + r][cur_j_smem][cur_i_smem] +
                    smem[smem_k_start + kRestrictStencilWidth - 1 - r][cur_j_smem][cur_i_smem]
                );
            }

            tmp2[m][n] = val;
        }
    }

    #pragma unroll
    for (int n = 0; n < kRestrictStencilWidth; n++) {
        double val = 0.0;
        #pragma unroll
        for (int r = 0; r < kRestrictStencilWidth / 2; r++) {
            val += C_RESTRICT[r] * (tmp2[r][n] + tmp2[kRestrictStencilWidth - 1 - r][n]);
        }
        tmp1[n] = val;
    }

    double final_val = 0.0;
    #pragma unroll
    for (int r = 0; r < kRestrictStencilWidth / 2; r++) {
        final_val += C_RESTRICT[r] * (tmp1[r] + tmp1[kRestrictStencilWidth - 1 - r]);
    }

    int out_idx = k * (transfer.extc0 * transfer.extc1) + j * transfer.extc0 + i;
    transfer.target_field[out_idx] = final_val;
}

bool gpu_build_restrict3_transfer(
    const double* d_src_f, double* d_dst_c,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry,
    GpuStencilTransfer* transfer
) {
    double CD[3], FD[3], base[3];
    int lbc[3], lbf[3], lbr[3], ubr[3];

    for (int d = 0; d < 3; d++) {
        CD[d] = (uubc[d] - llbc[d]) / (double)extc[d];
        FD[d] = (uubf[d] - llbf[d]) / (double)extf[d];
        if (llbc[d] <= llbf[d]) {
            base[d] = llbc[d];
        } else {
            int j_val = (int)std::trunc((llbc[d] - llbf[d]) / FD[d] + 0.4);
            if ((j_val / 2) * 2 == j_val) base[d] = llbf[d];
            else base[d] = llbf[d] - CD[d] / 2.0;
        }

        lbf[d] = (int)std::trunc((llbf[d] - base[d]) / FD[d] + 0.4) + 1;
        lbc[d] = (int)std::trunc((llbc[d] - base[d]) / CD[d] + 0.4) + 1;
        lbr[d] = (int)std::trunc((llbt[d] - base[d]) / CD[d] + 0.4) + 1;
        ubr[d] = (int)std::trunc((uubt[d] - base[d]) / CD[d] + 0.4);
    }

    int i_start = lbr[0] - lbc[0];
    int i_end   = ubr[0] - lbc[0];
    int j_start = lbr[1] - lbc[1];
    int j_end   = ubr[1] - lbc[1];
    int k_start = lbr[2] - lbc[2];
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
    transfer->lbc0 = lbc[0];
    transfer->lbc1 = lbc[1];
    transfer->lbc2 = lbc[2];
    transfer->lbf0 = lbf[0];
    transfer->lbf1 = lbf[1];
    transfer->lbf2 = lbf[2];
    transfer->extc0 = extc[0];
    transfer->extc1 = extc[1];
    transfer->extc2 = extc[2];
    transfer->extf0 = extf[0];
    transfer->extf1 = extf[1];
    transfer->extf2 = extf[2];
    transfer->SoA0 = SoA[0];
    transfer->SoA1 = SoA[1];
    transfer->SoA2 = SoA[2];
    transfer->Symmetry = Symmetry;
    return true;
}

void gpu_restrict3_launch(
    cudaStream_t stream,
    const double* d_src_f, double* d_dst_c,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
) {
    GpuStencilTransfer transfer;
    if (!gpu_build_restrict3_transfer(
            d_src_f,
            d_dst_c,
            llbc,
            uubc,
            extc,
            llbf,
            uubf,
            extf,
            llbt,
            uubt,
            SoA,
            Symmetry,
            &transfer
        )) {
        return;
    }

    dim3 block(kTransferBlockX, kTransferBlockY, kTransferBlockZ);
    dim3 grid(block_count_for(transfer.ni, block.x),
              block_count_for(transfer.nj, block.y),
              block_count_for(transfer.nk, block.z));

    restrict3_kernel<<<grid, block, 0, stream>>>(
        transfer.ni, transfer.nj, transfer.nk,
        transfer.i_start, transfer.j_start, transfer.k_start,
        transfer.lbc0, transfer.lbc1, transfer.lbc2,
        transfer.lbf0, transfer.lbf1, transfer.lbf2,
        transfer.extc0, transfer.extc1, transfer.extc2,
        transfer.extf0, transfer.extf1, transfer.extf2,
        transfer.source_field, transfer.target_field,
        transfer.SoA0, transfer.SoA1, transfer.SoA2,
        transfer.Symmetry
    );
}

void gpu_restrict3_batch_launch(
    cudaStream_t stream,
    const GpuStencilTransfer* transfers,
    int transfer_count
) {
    if (transfer_count <= 0) return;

    StencilTransferScratch& scratch = transfer_scratch_for_stream(restrict_transfer_scratches(), stream);
    ensure_stencil_transfer_capacity(scratch, transfer_count);
    cudaMemcpyAsync(
        scratch.device_transfers,
        transfers,
        transfer_count * sizeof(GpuStencilTransfer),
        cudaMemcpyHostToDevice,
        stream
    );

    dim3 block(kTransferBlockX, kTransferBlockY, kTransferBlockZ);
    const int blocks_x = block_count_for(transfers[0].ni, block.x);
    dim3 grid(blocks_x * transfer_count,
              block_count_for(transfers[0].nj, block.y),
              block_count_for(transfers[0].nk, block.z));
    restrict3_batch_kernel<<<grid, block, 0, stream>>>(scratch.device_transfers, transfer_count, blocks_x);
}

#endif  // Cell
