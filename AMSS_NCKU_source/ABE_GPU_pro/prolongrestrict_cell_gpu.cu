#include "macrodef.h"
#include "prolongrestrict.h"

#include "fmisc.h"
#include "fmisc_gpu.cuh"
#include "macrodef.fh"

#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>

#ifdef Cell

// ============================================================
//  1. Constant Memory Coefficients for all ghost_width values
// ============================================================

// --- GW=2 ---
__constant__ double C_PROLONG_GW2[4] = {
    -7.0  / 128.0,     // C1
    105.0 / 128.0,     // C2
     35.0 / 128.0,     // C3
    -5.0  / 128.0      // C4
};
__constant__ double C_RESTRICT_GW2[2] = {
    -1.0 / 16.0,       // C1
     9.0 / 16.0        // C2
};

// --- GW=3 ---
__constant__ double C_PROLONG_GW3[6] = {
     77.0  / 8192.0,   // C1
   -693.0  / 8192.0,   // C2
   3465.0  / 4096.0,   // C3
   1155.0  / 4096.0,   // C4
   -495.0  / 8192.0,   // C5
     63.0  / 8192.0    // C6
};
__constant__ double C_RESTRICT_GW3[3] = {
      3.0  / 256.0,    // C1
    -25.0  / 256.0,    // C2
     75.0  / 128.0     // C3
};

// --- GW=4 ---
__constant__ double C_PROLONG_GW4[8] = {
    -495.0    / 262144.0,   // C1
    5005.0    / 262144.0,   // C2
  -27027.0    / 262144.0,   // C3
  225225.0    / 262144.0,   // C4
   75075.0    / 262144.0,   // C5
  -19305.0    / 262144.0,   // C6
    4095.0    / 262144.0,   // C7
    -429.0    / 262144.0    // C8
};
__constant__ double C_RESTRICT_GW4[4] = {
     -5.0  / 2048.0,    // C1
     49.0  / 2048.0,    // C2
   -245.0  / 2048.0,    // C3
   1225.0  / 2048.0     // C4
};

// --- GW=5 ---
__constant__ double C_PROLONG_GW5[10] = {
     13585.0 / 33554432.0,  // C1
   -159885.0 / 33554432.0,  // C2
    230945.0 / 8388608.0,   // C3
    -969969.0 / 8388608.0,  // C4
  14549535.0 / 16777216.0,  // C5
   4849845.0 / 16777216.0,  // C6
    -692835.0 / 8388608.0,  // C7
     188955.0 / 8388608.0,  // C8
    -138567.0 / 33554432.0, // C9
      12155.0 / 33554432.0  // C10
};
__constant__ double C_RESTRICT_GW5[5] = {
      35.0 / 65536.0,    // C1
    -405.0 / 65536.0,    // C2
     567.0 / 16384.0,    // C3
   -2205.0 / 16384.0,    // C4
   19845.0 / 32768.0     // C5
};

// ============================================================
//  2. Coeffs<GW> Trait Struct
//     Uses __device__ functions to return __constant__ memory
//     pointers so they are callable from kernel code.
// ============================================================

template<int GW> struct Coeffs;

template<> struct Coeffs<2> {
    static constexpr int STENCIL = 4;
    static constexpr int HALF    = 2;
    __device__ static const double* prolong()    { return C_PROLONG_GW2; }
    __device__ static const double* restrict_c() { return C_RESTRICT_GW2; }
};

template<> struct Coeffs<3> {
    static constexpr int STENCIL = 6;
    static constexpr int HALF    = 3;
    __device__ static const double* prolong()    { return C_PROLONG_GW3; }
    __device__ static const double* restrict_c() { return C_RESTRICT_GW3; }
};

template<> struct Coeffs<4> {
    static constexpr int STENCIL = 8;
    static constexpr int HALF    = 4;
    __device__ static const double* prolong()    { return C_PROLONG_GW4; }
    __device__ static const double* restrict_c() { return C_RESTRICT_GW4; }
};

template<> struct Coeffs<5> {
    static constexpr int STENCIL = 10;
    static constexpr int HALF    = 5;
    __device__ static const double* prolong()    { return C_PROLONG_GW5; }
    __device__ static const double* restrict_c() { return C_RESTRICT_GW5; }
};

// ============================================================
//  3. Template Kernels
// ============================================================

// -------------------------------------------------------------------
//  prolong_kernel<GW>
//
//  Shared memory dimensions (static):
//    SMEM_Z = 2 + 2*GW   (BZ/2 + 2*GW, BZ=4)
//    SMEM_Y = 4 + 2*GW   (BY/2 + 2*GW, BY=8)
//    SMEM_X = 4 + 2*GW   (BX/2 + 2*GW, BX=8)
//
//  Uses fully-folded 3D interpolation (no tmp2/tmp1 arrays).
// -------------------------------------------------------------------
template<int GW>
__global__ void prolong_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    int extf0, int extf1, int extf2,
    int pextf0, int pextf1, int pextf2,
    const double* __restrict__ d_src_c,
    double* __restrict__ d_dst_f,
    double SoA0, double SoA1, double SoA2,
    int Symmetry
) {
    constexpr int STENCIL = Coeffs<GW>::STENCIL;

    // Shared memory tile: [Z][Y][X]
    constexpr int SMEM_Z = 2 + 2 * GW;  // BZ/2 + 2*GW
    constexpr int SMEM_Y = 4 + 2 * GW;  // BY/2 + 2*GW
    constexpr int SMEM_X = 4 + 2 * GW;  // BX/2 + 2*GW
    constexpr int SMEM_TOTAL = SMEM_Z * SMEM_Y * SMEM_X;

    extern __shared__ double smem_raw[];
    double (*smem)[SMEM_Y][SMEM_X] = reinterpret_cast<double (*)[SMEM_Y][SMEM_X]>(smem_raw);

    int tid = threadIdx.z * blockDim.y * blockDim.x
            + threadIdx.y * blockDim.x
            + threadIdx.x;

    // Block start in fine-grid global space
    int block_i_start = i_start + blockIdx.x * blockDim.x;
    int block_j_start = j_start + blockIdx.y * blockDim.y;
    int block_k_start = k_start + blockIdx.z * blockDim.z;

    // Base coarse-grid index for this block (includes -(GW-1) ghost offset)
    int base_c_i = (block_i_start + lbf0) / 2 - lbc0 - (GW - 1);
    int base_c_j = (block_j_start + lbf1) / 2 - lbc1 - (GW - 1);
    int base_c_k = (block_k_start + lbf2) / 2 - lbc2 - (GW - 1);

    // Cooperative load
    for (int idx = tid; idx < SMEM_TOTAL; idx += 256) {
        int loc_k = idx / (SMEM_Y * SMEM_X);
        int rem   = idx % (SMEM_Y * SMEM_X);
        int loc_j = rem / SMEM_X;
        int loc_i = rem % SMEM_X;

        smem[loc_k][loc_j][loc_i] =
            d_src_c[((base_c_k + loc_k) * pextc1 + (base_c_j + loc_j)) * pextc0
                    + (base_c_i + loc_i)];
    }
    __syncthreads();

    // Compute phase
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

    int cxI_i_0b = ii / 2 - lbc0;
    int cxI_j_0b = jj / 2 - lbc1;
    int cxI_k_0b = kk / 2 - lbc2;

    bool k_even = ((kk / 2) * 2 == kk);
    bool j_even = ((jj / 2) * 2 == jj);
    bool i_even = ((ii / 2) * 2 == ii);

    int smem_i_base = cxI_i_0b - (GW - 1) - base_c_i;
    int smem_j_base = cxI_j_0b - (GW - 1) - base_c_j;
    int smem_k_base = cxI_k_0b - (GW - 1) - base_c_k;

    // Fully-folded 3D interpolation
    double final_val = 0.0;
    const double* C_P = Coeffs<GW>::prolong();

    #pragma unroll
    for (int step_x = 0; step_x < STENCIL; step_x++) {
        int c_idx_x = i_even ? step_x : (STENCIL - 1 - step_x);
        int cur_i = smem_i_base + step_x;
        double val_y = 0.0;
        #pragma unroll
        for (int step_y = 0; step_y < STENCIL; step_y++) {
            int c_idx_y = j_even ? step_y : (STENCIL - 1 - step_y);
            int cur_j = smem_j_base + step_y;
            double val_z = 0.0;
            #pragma unroll
            for (int step_z = 0; step_z < STENCIL; step_z++) {
                int c_idx_z = k_even ? step_z : (STENCIL - 1 - step_z);
                val_z += C_P[c_idx_z] * smem[smem_k_base + step_z][cur_j][cur_i];
            }
            val_y += C_P[c_idx_y] * val_z;
        }
        final_val += C_P[c_idx_x] * val_y;
    }

    // Write output using padded stride
    int out_idx = (k * pextf1 + j) * pextf0 + i;
    d_dst_f[out_idx] = final_val;
}

// -------------------------------------------------------------------
//  restrict_kernel<GW>
//
//  Shared memory dimensions (static):
//    SMEM_Z = 8  + 2*(GW-1)  (2*BZ + 2*(GW-1), BZ=4)
//    SMEM_Y = 16 + 2*(GW-1)  (2*BY + 2*(GW-1), BY=8)
//    SMEM_X = 16 + 2*(GW-1)  (2*BX + 2*(GW-1), BX=8)
//
//  Uses symmetric coefficient folding for restriction stencil.
// -------------------------------------------------------------------
template<int GW>
__global__ void restrict_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    int extf0, int extf1, int extf2,
    int pextf0, int pextf1, int pextf2,
    const double* __restrict__ d_src_f,
    double* __restrict__ d_dst_c,
    double SoA0, double SoA1, double SoA2,
    int Symmetry
) {
    constexpr int STENCIL = Coeffs<GW>::STENCIL;

    // Shared memory tile: [Z][Y][X]
    constexpr int SMEM_Z = 8  + 2 * (GW - 1);  // 2*BZ + 2*(GW-1)
    constexpr int SMEM_Y = 16 + 2 * (GW - 1);  // 2*BY + 2*(GW-1)
    constexpr int SMEM_X = 16 + 2 * (GW - 1);  // 2*BX + 2*(GW-1)
    constexpr int SMEM_TOTAL = SMEM_Z * SMEM_Y * SMEM_X;

    extern __shared__ double smem_raw[];
    double (*smem)[SMEM_Y][SMEM_X] = reinterpret_cast<double (*)[SMEM_Y][SMEM_X]>(smem_raw);

    int tid = threadIdx.z * blockDim.y * blockDim.x
            + threadIdx.y * blockDim.x
            + threadIdx.x;

    int block_i_start = i_start + blockIdx.x * blockDim.x;
    int block_j_start = j_start + blockIdx.y * blockDim.y;
    int block_k_start = k_start + blockIdx.z * blockDim.z;

    // Base fine-grid index (includes -GW ghost offset)
    int base_f_i = 2 * (block_i_start + lbc0) - lbf0 - GW;
    int base_f_j = 2 * (block_j_start + lbc1) - lbf1 - GW;
    int base_f_k = 2 * (block_k_start + lbc2) - lbf2 - GW;

    // Cooperative load
    for (int idx = tid; idx < SMEM_TOTAL; idx += 256) {
        int loc_k = idx / (SMEM_Y * SMEM_X);
        int rem   = idx % (SMEM_Y * SMEM_X);
        int loc_j = rem / SMEM_X;
        int loc_i = rem % SMEM_X;

        smem[loc_k][loc_j][loc_i] =
            d_src_f[((base_f_k + loc_k) * pextf1 + (base_f_j + loc_j)) * pextf0
                    + (base_f_i + loc_i)];
    }
    __syncthreads();

    // Compute phase
    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= ni || j_local >= nj || k_local >= nk) return;

    int i = i_start + i_local;
    int j = j_start + j_local;
    int k = k_start + k_local;

    // Thread-local smem start (independent of GW)
    int smem_i_start = 2 * threadIdx.x;
    int smem_j_start = 2 * threadIdx.y;
    int smem_k_start = 2 * threadIdx.z;

    // Fully-folded 3D restriction with symmetric coefficient folding
    double final_val = 0.0;
    const double* C_R = Coeffs<GW>::restrict_c();

    #pragma unroll
    for (int step_x = 0; step_x < STENCIL; step_x++) {
        int cur_i_smem = smem_i_start + step_x;
        int c_idx_x = (step_x < GW) ? step_x : (STENCIL - 1 - step_x);

        double val_y = 0.0;
        #pragma unroll
        for (int step_y = 0; step_y < STENCIL; step_y++) {
            int cur_j_smem = smem_j_start + step_y;
            int c_idx_y = (step_y < GW) ? step_y : (STENCIL - 1 - step_y);

            double val_z = 0.0;
            #pragma unroll
            for (int step_z = 0; step_z < STENCIL; step_z++) {
                int cur_k_smem = smem_k_start + step_z;
                int c_idx_z = (step_z < GW) ? step_z : (STENCIL - 1 - step_z);

                val_z += C_R[c_idx_z] * smem[cur_k_smem][cur_j_smem][cur_i_smem];
            }
            val_y += C_R[c_idx_y] * val_z;
        }
        final_val += C_R[c_idx_x] * val_y;
    }

    int out_idx = (k * pextc1 + j) * pextc0 + i;
    d_dst_c[out_idx] = final_val;
}

// -------------------------------------------------------------------
//  fused_prolong_kernel<GW>
//
//  Same as prolong_kernel<GW> but loops over vars.num_vars.
//  Loads smem from vars.ptrs[v] each iteration.
//  Writes to d_dst_1d[v * total_per_var + out_idx].
// -------------------------------------------------------------------
template<int GW>
__global__ void fused_prolong_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    int extf0, int extf1, int extf2,
    FusedVars vars,
    double* __restrict__ d_dst_1d,
    int Symmetry
) {
    constexpr int STENCIL = Coeffs<GW>::STENCIL;

    constexpr int SMEM_Z = 2 + 2 * GW;
    constexpr int SMEM_Y = 4 + 2 * GW;
    constexpr int SMEM_X = 4 + 2 * GW;
    constexpr int SMEM_TOTAL = SMEM_Z * SMEM_Y * SMEM_X;

    extern __shared__ double smem_raw[];
    double (*smem)[SMEM_Y][SMEM_X] = reinterpret_cast<double (*)[SMEM_Y][SMEM_X]>(smem_raw);

    int tid = threadIdx.z * blockDim.y * blockDim.x
            + threadIdx.y * blockDim.x
            + threadIdx.x;

    int block_i_start = i_start + blockIdx.x * blockDim.x;
    int block_j_start = j_start + blockIdx.y * blockDim.y;
    int block_k_start = k_start + blockIdx.z * blockDim.z;

    int base_c_i = (block_i_start + lbf0) / 2 - lbc0 - (GW - 1);
    int base_c_j = (block_j_start + lbf1) / 2 - lbc1 - (GW - 1);
    int base_c_k = (block_k_start + lbf2) / 2 - lbc2 - (GW - 1);

    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    int total_per_var = extf0 * extf1 * extf2;

    // Pre-compute per-thread values outside the loop
    bool is_valid = (i_local < ni && j_local < nj && k_local < nk);
    int out_idx;
    int smem_i_base, smem_j_base, smem_k_base;
    bool k_even, j_even, i_even;

    if (is_valid) {
        int i = i_start + i_local;
        int j = j_start + j_local;
        int k = k_start + k_local;

        int ii = i + lbf0;
        int jj = j + lbf1;
        int kk = k + lbf2;

        int cxI_i_0b = ii / 2 - lbc0;
        int cxI_j_0b = jj / 2 - lbc1;
        int cxI_k_0b = kk / 2 - lbc2;

        k_even = ((kk / 2) * 2 == kk);
        j_even = ((jj / 2) * 2 == jj);
        i_even = ((ii / 2) * 2 == ii);

        smem_i_base = cxI_i_0b - (GW - 1) - base_c_i;
        smem_j_base = cxI_j_0b - (GW - 1) - base_c_j;
        smem_k_base = cxI_k_0b - (GW - 1) - base_c_k;

        out_idx = k * (extf0 * extf1) + j * extf0 + i;
    }

    const double* C_P = Coeffs<GW>::prolong();

    for (int v = 0; v < vars.num_vars; ++v) {

        // Cooperative load from vars.ptrs[v]
        for (int idx = tid; idx < SMEM_TOTAL; idx += 256) {
            int loc_k = idx / (SMEM_Y * SMEM_X);
            int rem   = idx % (SMEM_Y * SMEM_X);
            int loc_j = rem / SMEM_X;
            int loc_i = rem % SMEM_X;

            smem[loc_k][loc_j][loc_i] =
                vars.ptrs[v][((base_c_k + loc_k) * pextc1 + (base_c_j + loc_j)) * pextc0
                             + (base_c_i + loc_i)];
        }
        __syncthreads();

        // Fully-folded 3D interpolation
        if (is_valid) {
            double final_val = 0.0;

            #pragma unroll
            for (int step_x = 0; step_x < STENCIL; step_x++) {
                int c_idx_x = i_even ? step_x : (STENCIL - 1 - step_x);
                int cur_i = smem_i_base + step_x;

                double val_y = 0.0;
                #pragma unroll
                for (int step_y = 0; step_y < STENCIL; step_y++) {
                    int c_idx_y = j_even ? step_y : (STENCIL - 1 - step_y);
                    int cur_j = smem_j_base + step_y;

                    double val_z = 0.0;
                    #pragma unroll
                    for (int step_z = 0; step_z < STENCIL; step_z++) {
                        int c_idx_z = k_even ? step_z : (STENCIL - 1 - step_z);
                        val_z += C_P[c_idx_z] * smem[smem_k_base + step_z][cur_j][cur_i];
                    }
                    val_y += C_P[c_idx_y] * val_z;
                }
                final_val += C_P[c_idx_x] * val_y;
            }
            d_dst_1d[v * total_per_var + out_idx] = final_val;
        }
        __syncthreads();
    }
}

// -------------------------------------------------------------------
//  fused_restrict_kernel<GW>
//
//  Same as restrict_kernel<GW> but loops over vars.num_vars.
//  Loads smem from vars.ptrs[v] each iteration.
//  Writes to d_dst_1d[v * total_per_var + out_idx].
// -------------------------------------------------------------------
template<int GW>
__global__ void fused_restrict_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int extc0, int extc1, int extc2,
    int extf0, int extf1, int extf2,
    int pextf0, int pextf1, int pextf2,
    FusedVars vars,
    double* __restrict__ d_dst_1d,
    int Symmetry
) {
    constexpr int STENCIL = Coeffs<GW>::STENCIL;

    constexpr int SMEM_Z = 8  + 2 * (GW - 1);
    constexpr int SMEM_Y = 16 + 2 * (GW - 1);
    constexpr int SMEM_X = 16 + 2 * (GW - 1);
    constexpr int SMEM_TOTAL = SMEM_Z * SMEM_Y * SMEM_X;

    extern __shared__ double smem_raw[];
    double (*smem)[SMEM_Y][SMEM_X] = reinterpret_cast<double (*)[SMEM_Y][SMEM_X]>(smem_raw);

    int tid = threadIdx.z * blockDim.y * blockDim.x
            + threadIdx.y * blockDim.x
            + threadIdx.x;

    int block_i_start = i_start + blockIdx.x * blockDim.x;
    int block_j_start = j_start + blockIdx.y * blockDim.y;
    int block_k_start = k_start + blockIdx.z * blockDim.z;

    int base_f_i = 2 * (block_i_start + lbc0) - lbf0 - GW;
    int base_f_j = 2 * (block_j_start + lbc1) - lbf1 - GW;
    int base_f_k = 2 * (block_k_start + lbc2) - lbf2 - GW;

    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    int total_per_var = extc0 * extc1 * extc2;

    // Pre-compute per-thread values outside the loop
    bool is_valid = (i_local < ni && j_local < nj && k_local < nk);
    int out_idx;
    int smem_i_start, smem_j_start, smem_k_start;

    if (is_valid) {
        int i = i_start + i_local;
        int j = j_start + j_local;
        int k = k_start + k_local;

        smem_i_start = 2 * threadIdx.x;
        smem_j_start = 2 * threadIdx.y;
        smem_k_start = 2 * threadIdx.z;
        out_idx = k * (extc0 * extc1) + j * extc0 + i;
    }

    const double* C_R = Coeffs<GW>::restrict_c();

    for (int v = 0; v < vars.num_vars; ++v) {

        // Cooperative load from vars.ptrs[v]
        for (int idx = tid; idx < SMEM_TOTAL; idx += 256) {
            int loc_k = idx / (SMEM_Y * SMEM_X);
            int rem   = idx % (SMEM_Y * SMEM_X);
            int loc_j = rem / SMEM_X;
            int loc_i = rem % SMEM_X;

            smem[loc_k][loc_j][loc_i] =
                vars.ptrs[v][((base_f_k + loc_k) * pextf1 + (base_f_j + loc_j)) * pextf0
                             + (base_f_i + loc_i)];
        }
        __syncthreads();

        // Fully-folded 3D restriction with symmetric coefficient folding
        if (is_valid) {
            double final_val = 0.0;

            #pragma unroll
            for (int step_x = 0; step_x < STENCIL; step_x++) {
                int cur_i_smem = smem_i_start + step_x;
                int c_idx_x = (step_x < GW) ? step_x : (STENCIL - 1 - step_x);

                double val_y = 0.0;
                #pragma unroll
                for (int step_y = 0; step_y < STENCIL; step_y++) {
                    int cur_j_smem = smem_j_start + step_y;
                    int c_idx_y = (step_y < GW) ? step_y : (STENCIL - 1 - step_y);

                    double val_z = 0.0;
                    #pragma unroll
                    for (int step_z = 0; step_z < STENCIL; step_z++) {
                        int cur_k_smem = smem_k_start + step_z;
                        int c_idx_z = (step_z < GW) ? step_z : (STENCIL - 1 - step_z);

                        val_z += C_R[c_idx_z] * smem[cur_k_smem][cur_j_smem][cur_i_smem];
                    }
                    val_y += C_R[c_idx_y] * val_z;
                }
                final_val += C_R[c_idx_x] * val_y;
            }
            d_dst_1d[v * total_per_var + out_idx] = final_val;
        }
        __syncthreads();
    }
}

// ============================================================
//  4. Template Launch Functions
// ============================================================

template<int GW>
void gpu_prolong_launch_impl(
    cudaStream_t stream,
    const double* d_src_c, double* d_dst_f,
    const double* llbc, const double* uubc, const int* extc, const int* pextc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
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

    if (ni <= 0 || nj <= 0 || nk <= 0) return;

    dim3 block(8, 8, 4);
    dim3 grid((ni + block.x - 1) / block.x,
              (nj + block.y - 1) / block.y,
              (nk + block.z - 1) / block.z);

    constexpr size_t shmem_bytes = (2 + 2*GW) * (4 + 2*GW) * (4 + 2*GW) * sizeof(double);
    prolong_kernel<GW><<<grid, block, shmem_bytes, stream>>>(
        ni, nj, nk,
        i_start, j_start, k_start,
        lbc[0], lbc[1], lbc[2],
        lbf[0], lbf[1], lbf[2],
        extc[0], extc[1], extc[2],
        pextc[0], pextc[1], pextc[2],
        extf[0], extf[1], extf[2],
        pextf[0], pextf[1], pextf[2],
        d_src_c, d_dst_f,
        SoA[0], SoA[1], SoA[2],
        Symmetry
    );
}

template<int GW>
void gpu_restrict_launch_impl(
    cudaStream_t stream,
    const double* d_src_f, double* d_dst_c,
    const double* llbc, const double* uubc, const int* extc, const int* pextc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
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

    if (ni <= 0 || nj <= 0 || nk <= 0) return;

    dim3 block(8, 8, 4);
    dim3 grid((ni + block.x - 1) / block.x,
              (nj + block.y - 1) / block.y,
              (nk + block.z - 1) / block.z);

    constexpr size_t shmem_bytes = (8 + 2*(GW-1)) * (16 + 2*(GW-1)) * (16 + 2*(GW-1)) * sizeof(double);

    if constexpr (GW >= 4) {
        cudaFuncSetAttribute(
            restrict_kernel<GW>,
            cudaFuncAttributeMaxDynamicSharedMemorySize,
            shmem_bytes
        );
        cudaFuncSetAttribute(
            restrict_kernel<GW>,
            cudaFuncAttributePreferredSharedMemoryCarveout,
            cudaSharedmemCarveoutMaxShared
        );
    }

    restrict_kernel<GW><<<grid, block, shmem_bytes, stream>>>(
        ni, nj, nk, i_start, j_start, k_start,
        lbc[0], lbc[1], lbc[2],
        lbf[0], lbf[1], lbf[2],
        extc[0], extc[1], extc[2],
        pextc[0], pextc[1], pextc[2],
        extf[0], extf[1], extf[2],
        pextf[0], pextf[1], pextf[2],
        d_src_f, d_dst_c,
        SoA[0], SoA[1], SoA[2],
        Symmetry
    );
}

template<int GW>
void gpu_fused_prolong_launch_impl(
    cudaStream_t stream,
    FusedVars vars, double* d_dst_1d,
    const double* llbc, const double* uubc, const int* extc, const int* pextc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    int Symmetry
) {
    if (vars.num_vars == 0) return;

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

    if (ni <= 0 || nj <= 0 || nk <= 0) return;

    dim3 block(8, 8, 4);
    dim3 grid((ni + block.x - 1) / block.x,
              (nj + block.y - 1) / block.y,
              (nk + block.z - 1) / block.z);

    constexpr size_t shmem_bytes = (2 + 2*GW) * (4 + 2*GW) * (4 + 2*GW) * sizeof(double);
    fused_prolong_kernel<GW><<<grid, block, shmem_bytes, stream>>>(
        ni, nj, nk,
        i_start, j_start, k_start,
        lbc[0], lbc[1], lbc[2],
        lbf[0], lbf[1], lbf[2],
        extc[0], extc[1], extc[2],
        pextc[0], pextc[1], pextc[2],
        extf[0], extf[1], extf[2],
        vars, d_dst_1d, Symmetry
    );
}

template<int GW>
void gpu_fused_restrict_launch_impl(
    cudaStream_t stream,
    FusedVars vars, double* d_dst_1d,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const double* llbt, const double* uubt,
    int Symmetry
) {
    if (vars.num_vars == 0) return;

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

    if (ni <= 0 || nj <= 0 || nk <= 0) return;

    dim3 block(8, 8, 4);
    dim3 grid((ni + block.x - 1) / block.x,
              (nj + block.y - 1) / block.y,
              (nk + block.z - 1) / block.z);

    constexpr size_t shmem_bytes = (8 + 2*(GW-1)) * (16 + 2*(GW-1)) * (16 + 2*(GW-1)) * sizeof(double);

    if constexpr (GW >= 4) {
        cudaFuncSetAttribute(
            fused_restrict_kernel<GW>,
            cudaFuncAttributeMaxDynamicSharedMemorySize,
            shmem_bytes
        );
        cudaFuncSetAttribute(
            fused_restrict_kernel<GW>,
            cudaFuncAttributePreferredSharedMemoryCarveout,
            cudaSharedmemCarveoutMaxShared
        );
    }

    fused_restrict_kernel<GW><<<grid, block, shmem_bytes, stream>>>(
        ni, nj, nk, i_start, j_start, k_start,
        lbc[0], lbc[1], lbc[2],
        lbf[0], lbf[1], lbf[2],
        extc[0], extc[1], extc[2],
        extf[0], extf[1], extf[2],
        pextf[0], pextf[1], pextf[2],
        vars, d_dst_1d, Symmetry
    );
}

// ============================================================
//  5. API Wrapper Functions
//
//  Compile-time dispatch via ghost_width macro (from macrodef.fh).
//  Signatures match declarations in prolongrestrict.h exactly.
// ============================================================

void gpu_prolong3_launch(
    cudaStream_t stream,
    const double* d_src_c, double* d_dst_f,
    const double* llbc, const double* uubc, const int* extc, const int* pextc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
) {
    gpu_prolong_launch_impl<ghost_width>(stream,
        d_src_c, d_dst_f, llbc, uubc, extc, pextc,
        llbf, uubf, extf, pextf, llbt, uubt, SoA, Symmetry);
}

void gpu_restrict3_launch(
    cudaStream_t stream,
    const double* d_src_f, double* d_dst_c,
    const double* llbc, const double* uubc, const int* extc, const int* pextc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
) {
    gpu_restrict_launch_impl<ghost_width>(stream,
        d_src_f, d_dst_c, llbc, uubc, extc, pextc,
        llbf, uubf, extf, pextf, llbt, uubt, SoA, Symmetry);
}

void gpu_fused_prolong3_launch(
    cudaStream_t stream,
    FusedVars vars, double* d_dst_1d,
    const double* llbc, const double* uubc, const int* extc, const int* pextc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    int Symmetry
) {
    gpu_fused_prolong_launch_impl<ghost_width>(stream,
        vars, d_dst_1d, llbc, uubc, extc, pextc,
        llbf, uubf, extf, llbt, uubt, Symmetry);
}

void gpu_fused_restrict3_launch(
    cudaStream_t stream,
    FusedVars vars, double* d_dst_1d,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const double* llbt, const double* uubt,
    int Symmetry
) {
    gpu_fused_restrict_launch_impl<ghost_width>(stream,
        vars, d_dst_1d, llbc, uubc, extc,
        llbf, uubf, extf, pextf, llbt, uubt, Symmetry);
}

// ============================================================
//  6. Cell prolongcopy3 kernel and launch
// ============================================================

__global__ void prolongcopy_kernel(
    int ni, int nj, int nk,
    int src_i_start, int src_j_start, int src_k_start,
    int dst_i_start, int dst_j_start, int dst_k_start,
    int pextc0, int pextc1,
    int pextf0, int pextf1,
    const double* __restrict__ d_src_c,
    double* __restrict__ d_dst_f
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ni || j >= nj || k >= nk) return;

    int si = src_i_start + i;
    int sj = src_j_start + j;
    int sk = src_k_start + k;

    int di = dst_i_start + i;
    int dj = dst_j_start + j;
    int dk = dst_k_start + k;

    d_dst_f[(dk * pextf1 + dj) * pextf0 + di] =
        d_src_c[(sk * pextc1 + sj) * pextc0 + si];
}

void gpu_prolongcopy_launch_impl(
    cudaStream_t stream,
    const double* d_src_c, double* d_dst_f,
    const double* llbc, const double* uubc, const int* extc, const int* pextc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
) {
    double CD[3], FD[3], base[3];
    int lbc[3], ubc[3], lbf[3], ubf[3], lbp[3], ubp[3], lbpc[3], ubpc[3];

    for (int d = 0; d < 3; d++) {
        CD[d] = (uubc[d] - llbc[d]) / (double)extc[d];
        FD[d] = CD[d] / 2.0;
        if (llbc[d] <= llbf[d]) {
            base[d] = llbc[d];
        } else {
            int j_val = (int)std::trunc((llbc[d] - llbf[d]) / FD[d] + 0.4);
            if ((j_val / 2) * 2 == j_val) base[d] = llbf[d];
            else base[d] = llbf[d] - CD[d] / 2.0;
        }

        lbf[d] = (int)std::trunc((llbf[d] - base[d]) / FD[d] + 0.4) + 1;
        ubf[d] = (int)std::trunc((uubf[d] - base[d]) / FD[d] + 0.4);
        lbc[d] = (int)std::trunc((llbc[d] - base[d]) / CD[d] + 0.4) + 1;
        ubc[d] = (int)std::trunc((uubc[d] - base[d]) / CD[d] + 0.4);
        lbp[d] = (int)std::trunc((llbt[d] - base[d]) / FD[d] + 0.4) + 1;
        lbpc[d] = (int)std::trunc((llbt[d] - base[d]) / CD[d] + 0.4) + 1;
        ubp[d] = (int)std::trunc((uubt[d] - base[d]) / FD[d] + 0.4);
        ubpc[d] = (int)std::trunc((uubt[d] - base[d]) / CD[d] + 0.4);
    }

    int imini = lbpc[0] - lbc[0] + 1 - ghost_width;
    int imaxi = ubpc[0] - lbc[0] + 1 + ghost_width;
    int jmini = lbpc[1] - lbc[1] + 1 - ghost_width;
    int jmaxi = ubpc[1] - lbc[1] + 1 + ghost_width;
    int kmini = lbpc[2] - lbc[2] + 1 - ghost_width;
    int kmaxi = ubpc[2] - lbc[2] + 1 + ghost_width;

    // Convert from 1-based Fortran to 0-based C indices
    int si_start = imini - 1;
    int sj_start = jmini - 1;
    int sk_start = kmini - 1;
    int di_start = 0, dj_start = 0, dk_start = 0;

    int ni = imaxi - imini + 1;
    int nj = jmaxi - jmini + 1;
    int nk = kmaxi - kmini + 1;

    if (ni <= 0 || nj <= 0 || nk <= 0) return;

    // Kernel handles out-of-bounds via the source/dest index mapping
    dim3 block(8, 8, 4);
    dim3 grid((ni + block.x - 1) / block.x,
              (nj + block.y - 1) / block.y,
              (nk + block.z - 1) / block.z);

    prolongcopy_kernel<<<grid, block, 0, stream>>>(
        ni, nj, nk,
        si_start, sj_start, sk_start,
        di_start, dj_start, dk_start,
        pextc[0], pextc[1],
        pextf[0], pextf[1],
        d_src_c, d_dst_f
    );
}

void gpu_prolongcopy3_launch(
    cudaStream_t stream,
    const double* d_src_c, double* d_dst_f,
    const double* llbc, const double* uubc, const int* extc, const int* pextc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
) {
    gpu_prolongcopy_launch_impl(stream,
        d_src_c, d_dst_f,
        llbc, uubc, extc, pextc,
        llbf, uubf, extf, pextf,
        llbt, uubt,
        SoA, Symmetry);
}

// ============================================================
//  7. Cell prolongmix3 kernel and launch
// ============================================================

template<int GW>
__global__ void prolongmix_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int lbpc0, int lbpc1, int lbpc2,
    int ilb0, int ilb1, int ilb2,
    int iub0, int iub1, int iub2,
    int pextf0, int pextf1,
    int pexti0, int pexti1,
    const double* __restrict__ d_src_f,
    const double* __restrict__ d_src_i,
    double* __restrict__ d_dst_f,
    double illb_v0, double illb_v1, double illb_v2,
    double iuub_v0, double iuub_v1, double iuub_v2,
    double llbf_v0, double llbf_v1, double llbf_v2,
    double FD_v0, double FD_v1, double FD_v2,
    int Symmetry
) {
    constexpr int STENCIL = 2 * GW;
    constexpr int NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2;

    double illb_arr[3] = {illb_v0, illb_v1, illb_v2};
    double iuub_arr[3] = {iuub_v0, iuub_v1, iuub_v2};
    double llbf_arr[3] = {llbf_v0, llbf_v1, llbf_v2};
    double FD_arr[3] = {FD_v0, FD_v1, FD_v2};

    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= ni || j_local >= nj || k_local >= nk) return;

    int i = i_start + i_local;
    int j = j_start + j_local;
    int k = k_start + k_local;

    int cxI[3] = {i, j, k};

    // Cell center coordinate (before coarse conversion)
    // i,j,k are 0-based; Fortran uses (cxI_f - 0.5) with 1-based cxI_f = cxI_0 + 1
    double ccp[3];
    for (int d = 0; d < 3; d++)
        ccp[d] = llbf_arr[d] + (cxI[d] + 0.5) * FD_arr[d];

    // fg: which dims have fine data
    int fg[3] = {0, 0, 0};
    for (int d = 0; d < 3; d++)
        if (illb_arr[d] < ccp[d] && iuub_arr[d] > ccp[d]) fg[d] = 1;

    int sum_fg = fg[0] + fg[1] + fg[2];
    if (sum_fg != 2) return;

    // Convert to coarse reference
    int lbf[3] = {lbf0, lbf1, lbf2};
    int lbpc[3] = {lbpc0, lbpc1, lbpc2};
    int ilb[3] = {ilb0, ilb1, ilb2};
    int iub_arr[3] = {iub0, iub1, iub2};

    int cxI_coarse[3];
    for (int d = 0; d < 3; d++)
        cxI_coarse[d] = (cxI[d] + lbf[d]) / 2 - lbpc[d] + 1;

    // Load stencil ya[STENCIL][STENCIL][STENCIL] from funi
    double ya[10][10][10]; // max STENCIL=10
    for (int sz = 0; sz < STENCIL; sz++) {
        for (int sy = 0; sy < STENCIL; sy++) {
            for (int sx = 0; sx < STENCIL; sx++) {
                int gi = cxI_coarse[0] - GW + 1 + sx - 1;
                int gj = cxI_coarse[1] - GW + 1 + sy - 1;
                int gk = cxI_coarse[2] - GW + 1 + sz - 1;
                ya[sz][sy][sx] =
                    d_src_i[((gk) * pexti1 + (gj)) * pexti0 + (gi)];
            }
        }
    }

    int ii = i + lbf0;
    int jj = j + lbf1;
    int kk = k + lbf2;

    bool i_even = ((ii / 2) * 2 == ii);
    bool j_even = ((jj / 2) * 2 == jj);
    bool k_even = ((kk / 2) * 2 == kk);

    const double* C_P = Coeffs<GW>::prolong();

    double X[10], Y[10], Z[10]; // coordinate arrays, max size 10

    // Determine which dimension needs interpolation
    int interp_dim = -1;
    for (int d = 0; d < 3; d++)
        if (fg[d] == 0) interp_dim = d;

    // Build coordinate arrays for boundary dimension
    // and compute cxB, cxT for fine data replacement
    int cxI_fine[3] = {i, j, k};
    int cxB[3], cxT[3];

    for (int d = 0; d < 3; d++) {
        if (d != interp_dim) continue;

        bool even = (d == 0) ? i_even : ((d == 1) ? j_even : k_even);
        double* coord = (d == 0) ? X : ((d == 1) ? Y : Z);

        if (even) {
            for (int n = 1; n <= GW; n++) {
                coord[GW - n]     = -0.5 - (n - 1) * 2;
                coord[GW - 1 + n] =  1.5 + (n - 1) * 2;
            }
        } else {
            for (int n = 1; n <= GW; n++) {
                coord[GW - n]     = -1.5 - (n - 1) * 2;
                coord[GW - 1 + n] =  0.5 + (n - 1) * 2;
            }
        }

        if (cxI_fine[d] > iub_arr[d]) {
            cxB[d] = iub_arr[d] - GW + 1 + (cxI_fine[d] - iub_arr[d] + 1) / 2;
            cxT[d] = iub_arr[d];
        } else if (cxI_fine[d] < ilb[d]) {
            cxB[d] = ilb[d];
            cxT[d] = ilb[d] + GW - 1 - (ilb[d] - cxI_fine[d]) / 2;
        } else {
            cxB[d] = cxI_fine[d];
            cxT[d] = cxI_fine[d];
        }
    }

    // Reduce ya through the two "fine" dimensions and interpolate
    // Dimension mapping: 0=x, 1=y, 2=z
    double result = 0.0;

    if (interp_dim == 0) {
        // Interpolate in x, reduce y and z with prolong coefficients
        double tmp1[10];
        for (int sx = 0; sx < STENCIL; sx++) {
            double tmp2[10];
            for (int sy = 0; sy < STENCIL; sy++) {
                int c_idx_y = j_even ? sy : (STENCIL - 1 - sy);
                double val_z = 0.0;
                for (int sz = 0; sz < STENCIL; sz++) {
                    int c_idx_z = k_even ? sz : (STENCIL - 1 - sz);
                    val_z += C_P[c_idx_z] * ya[sz][sy][sx];
                }
                tmp2[sy] = val_z;
            }
            double val_y = 0.0;
            for (int sy = 0; sy < STENCIL; sy++) {
                int c_idx_y = j_even ? sy : (STENCIL - 1 - sy);
                val_y += C_P[c_idx_y] * tmp2[sy];
            }
            tmp1[sx] = val_y;
        }

        // Replace with fine data
        if (cxI_fine[0] > iub_arr[0]) {
            for (int n = cxB[0]; n <= cxT[0]; n++) {
                int idx = GW - cxI_fine[0] + n;
                if (idx >= 0 && idx < STENCIL) {
                    X[idx] = (double)(n - cxI_fine[0]);
                    tmp1[idx] = d_src_f[(k * pextf1 + j) * pextf0 + n];
                }
            }
        } else if (cxI_fine[0] < ilb[0]) {
            for (int n = cxB[0]; n <= cxT[0]; n++) {
                int idx = GW - cxI_fine[0] + n;
                if (idx >= 0 && idx < STENCIL) {
                    X[idx] = (double)(n - cxI_fine[0]);
                    tmp1[idx] = d_src_f[(k * pextf1 + j) * pextf0 + n];
                }
            }
        }

        double dy;
        result = polint_device(X, tmp1, 0.0, &dy, STENCIL);

    } else if (interp_dim == 1) {
        // Interpolate in y, reduce x and z
        double tmp1[10];
        for (int sy = 0; sy < STENCIL; sy++) {
            double tmp2[10];
            for (int sx = 0; sx < STENCIL; sx++) {
                int c_idx_x = i_even ? sx : (STENCIL - 1 - sx);
                double val_z = 0.0;
                for (int sz = 0; sz < STENCIL; sz++) {
                    int c_idx_z = k_even ? sz : (STENCIL - 1 - sz);
                    val_z += C_P[c_idx_z] * ya[sz][sy][sx];
                }
                tmp2[sx] = val_z;
            }
            double val_x = 0.0;
            for (int sx = 0; sx < STENCIL; sx++) {
                int c_idx_x = i_even ? sx : (STENCIL - 1 - sx);
                val_x += C_P[c_idx_x] * tmp2[sx];
            }
            tmp1[sy] = val_x;
        }

        if (cxI_fine[1] > iub_arr[1]) {
            for (int n = cxB[1]; n <= cxT[1]; n++) {
                int idx = GW - cxI_fine[1] + n;
                if (idx >= 0 && idx < STENCIL) {
                    Y[idx] = (double)(n - cxI_fine[1]);
                    tmp1[idx] = d_src_f[(k * pextf1 + n) * pextf0 + i];
                }
            }
        } else if (cxI_fine[1] < ilb[1]) {
            for (int n = cxB[1]; n <= cxT[1]; n++) {
                int idx = GW - cxI_fine[1] + n;
                if (idx >= 0 && idx < STENCIL) {
                    Y[idx] = (double)(n - cxI_fine[1]);
                    tmp1[idx] = d_src_f[(k * pextf1 + n) * pextf0 + i];
                }
            }
        }

        double dy;
        result = polint_device(Y, tmp1, 0.0, &dy, STENCIL);

    } else { // interp_dim == 2
        // Interpolate in z, reduce x and y
        double tmp1[10];
        for (int sz = 0; sz < STENCIL; sz++) {
            double tmp2[10];
            for (int sx = 0; sx < STENCIL; sx++) {
                int c_idx_x = i_even ? sx : (STENCIL - 1 - sx);
                double val_y = 0.0;
                for (int sy = 0; sy < STENCIL; sy++) {
                    int c_idx_y = j_even ? sy : (STENCIL - 1 - sy);
                    val_y += C_P[c_idx_y] * ya[sz][sy][sx];
                }
                tmp2[sx] = val_y;
            }
            double val_x = 0.0;
            for (int sx = 0; sx < STENCIL; sx++) {
                int c_idx_x = i_even ? sx : (STENCIL - 1 - sx);
                val_x += C_P[c_idx_x] * tmp2[sx];
            }
            tmp1[sz] = val_x;
        }

        if (cxI_fine[2] > iub_arr[2]) {
            for (int n = cxB[2]; n <= cxT[2]; n++) {
                int idx = GW - cxI_fine[2] + n;
                if (idx >= 0 && idx < STENCIL) {
                    Z[idx] = (double)(n - cxI_fine[2]);
                    tmp1[idx] = d_src_f[(n * pextf1 + j) * pextf0 + i];
                }
            }
        } else if (cxI_fine[2] < ilb[2]) {
            for (int n = cxB[2]; n <= cxT[2]; n++) {
                int idx = GW - cxI_fine[2] + n;
                if (idx >= 0 && idx < STENCIL) {
                    Z[idx] = (double)(n - cxI_fine[2]);
                    tmp1[idx] = d_src_f[(n * pextf1 + j) * pextf0 + i];
                }
            }
        }

        double dy;
        result = polint_device(Z, tmp1, 0.0, &dy, STENCIL);
    }

    d_dst_f[(k * pextf1 + j) * pextf0 + i] = result;
}

template<int GW>
void gpu_prolongmix_launch_impl(
    cudaStream_t stream,
    const double* d_src_f, const double* d_src_i, double* d_dst_f,
    const double* llbc, const double* uubc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const int* exti,
    const double* llbp, const double* uubp,
    const double* illb, const double* iuub,
    const double* SoA, int Symmetry
) {
    double FD[3], CD[3], base[3];
    int lbc[3], ubc[3], lbf[3], ubf[3], lbp[3], ubp[3], lbpc[3], ubpc[3];
    int ilb[3], iub_arr[3];

    for (int d = 0; d < 3; d++) {
        FD[d] = (uubf[d] - llbf[d]) / (double)extf[d];
        CD[d] = FD[d] * 2.0;
        if (llbc[d] <= llbf[d]) {
            base[d] = llbc[d];
        } else {
            int j_val = (int)std::trunc((llbc[d] - llbf[d]) / FD[d] + 0.4);
            if ((j_val / 2) * 2 == j_val) base[d] = llbf[d];
            else base[d] = llbf[d] - CD[d] / 2.0;
        }

        lbf[d] = (int)std::trunc((llbf[d] - base[d]) / FD[d] + 0.4) + 1;
        ubf[d] = (int)std::trunc((uubf[d] - base[d]) / FD[d] + 0.4);
        lbc[d] = (int)std::trunc((llbc[d] - base[d]) / CD[d] + 0.4) + 1;
        ubc[d] = (int)std::trunc((uubc[d] - base[d]) / CD[d] + 0.4);
        lbp[d] = (int)std::trunc((llbp[d] - base[d]) / FD[d] + 0.4) + 1;
        lbpc[d] = (int)std::trunc((llbp[d] - base[d]) / CD[d] + 0.4) + 1;
        ubp[d] = (int)std::trunc((uubp[d] - base[d]) / FD[d] + 0.4);
        ubpc[d] = (int)std::trunc((uubp[d] - base[d]) / CD[d] + 0.4);
        ilb[d] = (int)std::trunc((illb[d] - base[d]) / FD[d] + 0.4) + 1;
        iub_arr[d] = (int)std::trunc((iuub[d] - base[d]) / FD[d] + 0.4);
    }

    int imino = lbp[0] - lbf[0] + 1;
    int imaxo = ubp[0] - lbf[0] + 1;
    int jmino = lbp[1] - lbf[1] + 1;
    int jmaxo = ubp[1] - lbf[1] + 1;
    int kmino = lbp[2] - lbf[2] + 1;
    int kmaxo = ubp[2] - lbf[2] + 1;

    // Adjust lbpc/ubpc by ghost_width
    for (int d = 0; d < 3; d++) {
        lbpc[d] -= ghost_width;
        ubpc[d] += ghost_width;
    }

    // Convert ilb/iub to fine-grid relative indices
    int ilb_rel[3], iub_rel[3];
    for (int d = 0; d < 3; d++) {
        ilb_rel[d] = ilb[d] - lbf[d] + 1;
        iub_rel[d] = iub_arr[d] - lbf[d] + 1;
        ilb_rel[d] = std::max(ilb_rel[d], 1);
        iub_rel[d] = std::min(iub_rel[d], extf[d]);
    }

    int ni = imaxo - imino + 1;
    int nj = jmaxo - jmino + 1;
    int nk = kmaxo - kmino + 1;

    if (ni <= 0 || nj <= 0 || nk <= 0) return;

    // Convert from 1-based Fortran to 0-based C indexing
    int i_start = imino - 1;
    int j_start = jmino - 1;
    int k_start = kmino - 1;

    dim3 block(8, 8, 4);
    dim3 grid((ni + block.x - 1) / block.x,
              (nj + block.y - 1) / block.y,
              (nk + block.z - 1) / block.z);

    prolongmix_kernel<GW><<<grid, block, 0, stream>>>(
        ni, nj, nk,
        i_start, j_start, k_start,
        lbc[0], lbc[1], lbc[2],
        lbf[0], lbf[1], lbf[2],
        lbpc[0], lbpc[1], lbpc[2],
        ilb_rel[0] - 1, ilb_rel[1] - 1, ilb_rel[2] - 1,
        iub_rel[0] - 1, iub_rel[1] - 1, iub_rel[2] - 1,
        pextf[0], pextf[1],
        exti[0] + 2 * ghost_width, exti[1] + 2 * ghost_width,
        d_src_f, d_src_i, d_dst_f,
        illb[0], illb[1], illb[2],
        iuub[0], iuub[1], iuub[2],
        llbf[0], llbf[1], llbf[2],
        FD[0], FD[1], FD[2],
        Symmetry
    );
}

void gpu_prolongmix3_launch(
    cudaStream_t stream,
    const double* d_src_f, const double* d_src_i, double* d_dst_f,
    const double* llbc, const double* uubc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const int* exti,
    const double* llbp, const double* uubp,
    const double* illb, const double* iuub,
    const double* SoA, int Symmetry
) {
    gpu_prolongmix_launch_impl<ghost_width>(stream,
        d_src_f, d_src_i, d_dst_f,
        llbc, uubc,
        llbf, uubf, extf, pextf,
        exti,
        llbp, uubp,
        illb, iuub,
        SoA, Symmetry);
}

#endif // Cell
