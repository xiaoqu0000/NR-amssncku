#include "macrodef.h"
// =====================================================================
//  prolongrestrict_vertex_gpu.cu
//
//  GPU implementations of vertex-based prolongation, restriction,
//  prolongcopy, and prolongmix operations.
//
//  KEY DIFFERENCES from cell version:
//    - CD = (uubc-llbc)/(extc-1), FD = CD/2  (vertex spacing)
//    - Vertex prolongation uses restriction coefficients (C_RESTRICT)
//      with parity-based interpolation:
//        even dim -> symmetric coefficient folding over STENCIL points
//        odd  dim -> pick middle value at offset GW-1 (no coeff mult)
//    - Vertex restriction is a direct copy (no interpolation)
//    - Vertex prolongcopy is a block copy with ghost zone handling
//    - Vertex prolongmix uses polint_device for the non-fine dimension
// =====================================================================

#include "prolongrestrict.h"

#include "fmisc.h"
#include "fmisc_gpu.cuh"
#include "macrodef.fh"

#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>

#ifdef Vertex

// ============================================================
//  1. Constant Memory Coefficients
//
//  Vertex prolongation uses what the cell file calls "restriction"
//  coefficients. We define them with CV_PROLONG_GW* names to avoid
//  linker conflicts with the cell .cu file's C_RESTRICT_GW* arrays.
//  The VALUES are identical to the cell file's C_RESTRICT_GW*.
// ============================================================

// --- GW=2 ---
__constant__ double CV_PROLONG_GW2[2] = {
    -1.0 / 16.0,       // C1
     9.0 / 16.0        // C2
};

// --- GW=3 ---
__constant__ double CV_PROLONG_GW3[3] = {
      3.0  / 256.0,    // C1
    -25.0  / 256.0,    // C2
     75.0  / 128.0     // C3
};

// --- GW=4 ---
__constant__ double CV_PROLONG_GW4[4] = {
     -5.0  / 2048.0,   // C1
     49.0  / 2048.0,   // C2
   -245.0  / 2048.0,   // C3
   1225.0  / 2048.0    // C4
};

// --- GW=5 ---
__constant__ double CV_PROLONG_GW5[5] = {
      35.0 / 65536.0,  // C1
    -405.0 / 65536.0,  // C2
     567.0 / 16384.0,  // C3
   -2205.0 / 16384.0,  // C4
   19845.0 / 32768.0   // C5
};

// ============================================================
//  2. VertexCoeffs<GW> Trait Struct
//
//  Returns pointers to the vertex prolongation constant arrays.
//  For vertex, the "prolong" coefficients are the restriction
//  coefficients used with symmetric folding.
// ============================================================

template<int GW> struct VertexCoeffs;

template<> struct VertexCoeffs<2> {
    static constexpr int STENCIL = 4;   // 2 * GW
    static constexpr int HALF    = 2;
    __device__ static const double* prolong() { return CV_PROLONG_GW2; }
};

template<> struct VertexCoeffs<3> {
    static constexpr int STENCIL = 6;
    static constexpr int HALF    = 3;
    __device__ static const double* prolong() { return CV_PROLONG_GW3; }
};

template<> struct VertexCoeffs<4> {
    static constexpr int STENCIL = 8;
    static constexpr int HALF    = 4;
    __device__ static const double* prolong() { return CV_PROLONG_GW4; }
};

template<> struct VertexCoeffs<5> {
    static constexpr int STENCIL = 10;
    static constexpr int HALF    = 5;
    __device__ static const double* prolong() { return CV_PROLONG_GW5; }
};

// ============================================================
//  3. Vertex Common Index Computation (shared by host launchers)
//
//  Computes CD, FD, base, lbc, lbf, lbp, ubp using the vertex
//  formulas: CD = (uubc-llbc)/(extc-1), FD = CD/2.
//  Returns the output iteration bounds [i_start..i_end] etc.
// ============================================================

struct VertexProlongIndex {
    int lbc[3], lbf[3], lbp[3], ubp[3];
    int i_start, i_end, j_start, j_end, k_start, k_end;
    int ni, nj, nk;
};

static inline VertexProlongIndex compute_vertex_prolong_index(
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt
) {
    VertexProlongIndex idx;
    double CD[3], FD[3], base[3];

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

        idx.lbf[d] = (int)std::trunc((llbf[d] - base[d]) / FD[d] + 0.4) + 1;
        idx.lbc[d] = (int)std::trunc((llbc[d] - base[d]) / CD[d] + 0.4) + 1;
        idx.lbp[d] = (int)std::trunc((llbt[d] - base[d]) / FD[d] + 0.4) + 1;
        idx.ubp[d] = (int)std::trunc((uubt[d] - base[d]) / FD[d] + 0.4) + 1;
    }

    idx.i_start = idx.lbp[0] - idx.lbf[0];
    idx.i_end   = idx.ubp[0] - idx.lbf[0];
    idx.j_start = idx.lbp[1] - idx.lbf[1];
    idx.j_end   = idx.ubp[1] - idx.lbf[1];
    idx.k_start = idx.lbp[2] - idx.lbf[2];
    idx.k_end   = idx.ubp[2] - idx.lbf[2];

    idx.ni = idx.i_end - idx.i_start + 1;
    idx.nj = idx.j_end - idx.j_start + 1;
    idx.nk = idx.k_end - idx.k_start + 1;

    return idx;
}

// ============================================================
//  4. Template Kernels
// ============================================================

// -------------------------------------------------------------------
//  prolong_kernel<GW>
//
//  Vertex prolongation kernel with parity-based interpolation.
//
//  Shared memory tile dimensions (same layout as cell):
//    SMEM_Z = 2 + 2*GW
//    SMEM_Y = 4 + 2*GW
//    SMEM_X = 4 + 2*GW
//
//  Algorithm per output point (i,j,k):
//    1. Compute global index ii = i + lbf0 (0-based, same as cell kernel)
//    2. For each dimension d, check parity of ii:
//       - even -> use symmetric coefficient folding (restriction coeffs)
//       - odd  -> pick middle value (offset GW-1 in stencil)
//    3. Fully-folded 3D interpolation across all parities.
// -------------------------------------------------------------------
template<int GW>
__global__ void prolong_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int pextc0, int pextc1,
    int pextf0, int pextf1,
    const double* __restrict__ d_src_c,
    double* __restrict__ d_dst_f,
    int Symmetry
) {
    constexpr int STENCIL = VertexCoeffs<GW>::STENCIL;

    // Shared memory tile: [Z][Y][X]
    constexpr int SMEM_Z = 2 + 2 * GW;
    constexpr int SMEM_Y = 4 + 2 * GW;
    constexpr int SMEM_X = 4 + 2 * GW;
    constexpr int SMEM_TOTAL = SMEM_Z * SMEM_Y * SMEM_X;

    extern __shared__ double smem_raw[];
    double (*smem)[SMEM_Y][SMEM_X] =
        reinterpret_cast<double (*)[SMEM_Y][SMEM_X]>(smem_raw);

    int tid = threadIdx.z * blockDim.y * blockDim.x
            + threadIdx.y * blockDim.x
            + threadIdx.x;

    // Block start in fine-grid global space
    int block_i_start = i_start + blockIdx.x * blockDim.x;
    int block_j_start = j_start + blockIdx.y * blockDim.y;
    int block_k_start = k_start + blockIdx.z * blockDim.z;

    // Base coarse-grid index for this block (includes -(GW-1) ghost offset)
    int base_c_i = (block_i_start + 1 + lbf0) / 2 - lbc0 - (GW - 1);
    int base_c_j = (block_j_start + 1 + lbf1) / 2 - lbc1 - (GW - 1);
    int base_c_k = (block_k_start + 1 + lbf2) / 2 - lbc2 - (GW - 1);

    // Cooperative load into shared memory
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

    // Parity: ii = i + lbf (0-based, same as cell kernel)
    int ii = i + lbf0;
    int jj = j + lbf1;
    int kk = k + lbf2;

    // Coarse index: emulate Fortran's (i_1based + lbf) / 2
    int cxI_i_0b = (i + 1 + lbf0) / 2 - lbc0;
    int cxI_j_0b = (j + 1 + lbf1) / 2 - lbc1;
    int cxI_k_0b = (k + 1 + lbf2) / 2 - lbc2;

    // Parity checks
    bool i_even = ((ii / 2) * 2 == ii);
    bool j_even = ((jj / 2) * 2 == jj);
    bool k_even = ((kk / 2) * 2 == kk);

    // Shared memory base offsets for this thread
    int smem_i_base = cxI_i_0b - (GW - 1) - base_c_i;
    int smem_j_base = cxI_j_0b - (GW - 1) - base_c_j;
    int smem_k_base = cxI_k_0b - (GW - 1) - base_c_k;

    // Fully-folded 3D interpolation with vertex parity logic
    const double* C_V = VertexCoeffs<GW>::prolong();

    double final_val = 0.0;

    #pragma unroll
    for (int step_x = 0; step_x < STENCIL; step_x++) {
        // Vertex parity: even -> symmetric fold, odd -> middle pick
        double coeff_x;
        if (i_even) {
            int c_idx = (step_x < GW) ? step_x : (STENCIL - 1 - step_x);
            coeff_x = C_V[c_idx];
        } else {
            coeff_x = (step_x == GW - 1) ? 1.0 : 0.0;
        }
        int cur_i = smem_i_base + step_x;

        double val_y = 0.0;
        #pragma unroll
        for (int step_y = 0; step_y < STENCIL; step_y++) {
            double coeff_y;
            if (j_even) {
                int c_idx = (step_y < GW) ? step_y : (STENCIL - 1 - step_y);
                coeff_y = C_V[c_idx];
            } else {
                coeff_y = (step_y == GW - 1) ? 1.0 : 0.0;
            }
            int cur_j = smem_j_base + step_y;

            double val_z = 0.0;
            #pragma unroll
            for (int step_z = 0; step_z < STENCIL; step_z++) {
                double coeff_z;
                if (k_even) {
                    int c_idx = (step_z < GW) ? step_z : (STENCIL - 1 - step_z);
                    coeff_z = C_V[c_idx];
                } else {
                    coeff_z = (step_z == GW - 1) ? 1.0 : 0.0;
                }
                val_z += coeff_z * smem[smem_k_base + step_z][cur_j][cur_i];
            }
            val_y += coeff_y * val_z;
        }
        final_val += coeff_x * val_y;
    }

    // Write output using padded stride
    int out_idx = (k * pextf1 + j) * pextf0 + i;
    d_dst_f[out_idx] = final_val;
}

// -------------------------------------------------------------------
//  restrict_kernel
//
//  Vertex restriction is a direct copy: each coarse point maps to
//  a specific fine point. No interpolation or shared memory needed.
//
//  Mapping: fine_C = 2*(coarse_C + lbc) - lbf - 1
// -------------------------------------------------------------------
__global__ void restrict_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int pextc0, int pextc1,
    int pextf0, int pextf1,
    const double* __restrict__ d_src_f,
    double* __restrict__ d_dst_c,
    int Symmetry
) {
    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= ni || j_local >= nj || k_local >= nk) return;

    int i = i_start + i_local;
    int j = j_start + j_local;
    int k = k_start + k_local;

    // Direct copy: fine_C = 2*(coarse_C + lbc) - lbf - 1
    int fine_i = 2 * (i + lbc0) - lbf0 - 1;
    int fine_j = 2 * (j + lbc1) - lbf1 - 1;
    int fine_k = 2 * (k + lbc2) - lbf2 - 1;

    // Read from fine grid, write to coarse grid
    d_dst_c[(k * pextc1 + j) * pextc0 + i] =
        d_src_f[(fine_k * pextf1 + fine_j) * pextf0 + fine_i];
}

// -------------------------------------------------------------------
//  prolongcopy_kernel
//
//  Block copy from coarse grid to fine grid output with ghost zone.
//  Each thread copies one element from src_c[imini+i, jmini+j, kmini+k]
//  to dst_f[i, j, k].
//
//  The host precomputes imini/imaxi etc. using vertex CD/FD formulas.
// -------------------------------------------------------------------
__global__ void prolongcopy_kernel(
    int ni, int nj, int nk,
    int imini, int jmini, int kmini,
    int dst_x_offset, int dst_y_offset, int dst_z_offset,
    int pextc0, int pextc1,
    int pextf0, int pextf1,
    const double* __restrict__ d_src_c,
    double* __restrict__ d_dst_f
) {
    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= ni || j_local >= nj || k_local >= nk) return;

    // Source index in coarse grid
    int src_i = imini + i_local;
    int src_j = jmini + j_local;
    int src_k = kmini + k_local;

    // Destination index in fine output (with ghost zone offset)
    int dst_i = dst_x_offset + i_local;
    int dst_j = dst_y_offset + j_local;
    int dst_k = dst_z_offset + k_local;

    double val = d_src_c[(src_k * pextc1 + src_j) * pextc0 + src_i];
    d_dst_f[(dst_k * pextf1 + dst_j) * pextf0 + dst_i] = val;
}

// -------------------------------------------------------------------
//  prolongmix_kernel<GW>
//
//  Mixed prolongation for points where exactly 2 of 3 dimensions
//  have fine-grid data available. The remaining dimension uses
//  polint_device (polynomial interpolation).
//
//  Per-thread algorithm:
//    1. Load local stencil ya[STENCIL][STENCIL][STENCIL] from funi
//    2. Compute fg (which dims have fine data in [cxB, cxT])
//    3. For the two dims with fine data:
//       - Apply restriction coefficient folding (even/odd parity)
//       - Replace with actual fine values at appropriate positions
//    4. Call polint_device for the remaining dimension at x=0
// -------------------------------------------------------------------
template<int GW>
__global__ void prolongmix_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int lbpc0, int lbpc1, int lbpc2,
    int ilb0, int ilb1, int ilb2,
    int iub0, int iub1, int iub2,
    int pextc0, int pextc1,
    int pextf0, int pextf1,
    int pexti0, int pexti1,
    const double* __restrict__ d_src_f,   // fine grid data (funf)
    const double* __restrict__ d_src_i,   // coarse stencil data (funi)
    double* __restrict__ d_dst_f,         // fine grid output (funf, in-place)
    int Symmetry
) {
    constexpr int STENCIL = 2 * GW;

    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= ni || j_local >= nj || k_local >= nk) return;

    int i = i_start + i_local;
    int j = j_start + j_local;
    int k = k_start + k_local;

    // Parity (0-based: ii = i + lbf)
    int ii = i + lbf0;
    int jj = j + lbf1;
    int kk = k + lbf2;

    // Coarse-level index: emulate Fortran's (i_1based + lbf) / 2
    int cxI_i = (i + 1 + lbf0) / 2;
    int cxI_j = (j + 1 + lbf1) / 2;
    int cxI_k = (k + 1 + lbf2) / 2;

    // Convert to funi array index: cxI = cxI - lbpc + 1
    cxI_i = cxI_i - lbpc0 + 1;
    cxI_j = cxI_j - lbpc1 + 1;
    cxI_k = cxI_k - lbpc2 + 1;

    // Load local stencil from funi (same pattern as cell kernel)
    double ya[STENCIL][STENCIL][STENCIL];
    for (int sz = 0; sz < STENCIL; sz++) {
        for (int sy = 0; sy < STENCIL; sy++) {
            for (int sx = 0; sx < STENCIL; sx++) {
                int src_i = cxI_i - GW + sx;
                int src_j = cxI_j - GW + sy;
                int src_k = cxI_k - GW + sz;
                ya[sz][sy][sx] =
                    d_src_i[(src_k * pexti1 + src_j) * pexti0 + src_i];
            }
        }
    }

    // Determine fine data bounds for this point
    int cxB[3], cxT[3];
    cxB[0] = (i - GW > ilb0) ? i - GW : ilb0;
    cxB[1] = (j - GW > ilb1) ? j - GW : ilb1;
    cxB[2] = (k - GW > ilb2) ? k - GW : ilb2;
    cxT[0] = (i + GW - 1 < iub0) ? i + GW - 1 : iub0;
    cxT[1] = (j + GW - 1 < iub1) ? j + GW - 1 : iub1;
    cxT[2] = (k + GW - 1 < iub2) ? k + GW - 1 : iub2;

    // fg[d] = 1 if dimension d has fine data covering this point
    int fg[3];
    fg[0] = (cxT[0] >= i && cxB[0] <= i) ? 1 : 0;
    fg[1] = (cxT[1] >= j && cxB[1] <= j) ? 1 : 0;
    fg[2] = (cxT[2] >= k && cxB[2] <= k) ? 1 : 0;

    int sum_fg = fg[0] + fg[1] + fg[2];

    // Only process points with exactly 2 fine-grid dimensions
    // and valid fine data range
    if (sum_fg != 2 || cxT[0] < cxB[0] || cxT[1] < cxB[1] || cxT[2] < cxB[2]) {
        return;
    }

    // Build X/Y/Z arrays for polint_device
    double X[STENCIL], Y[STENCIL], Z[STENCIL];
    // Build parity-based coordinate arrays
    for (int n = 0; n < GW; n++) {
        if (ii / 2 * 2 == ii) {
            X[GW - 1 - n] = -1.0 - (double)(n) * 2.0;
            X[GW + n]     =  1.0 + (double)(n) * 2.0;
        } else {
            X[GW - 1 - n] = -(double)(n) * 2.0;
            X[GW + n]     =  (double)(n + 1) * 2.0;
        }
    }
    for (int n = 0; n < GW; n++) {
        if (jj / 2 * 2 == jj) {
            Y[GW - 1 - n] = -1.0 - (double)(n) * 2.0;
            Y[GW + n]     =  1.0 + (double)(n) * 2.0;
        } else {
            Y[GW - 1 - n] = -(double)(n) * 2.0;
            Y[GW + n]     =  (double)(n + 1) * 2.0;
        }
    }
    for (int n = 0; n < GW; n++) {
        if (kk / 2 * 2 == kk) {
            Z[GW - 1 - n] = -1.0 - (double)(n) * 2.0;
            Z[GW + n]     =  1.0 + (double)(n) * 2.0;
        } else {
            Z[GW - 1 - n] = -(double)(n) * 2.0;
            Z[GW + n]     =  (double)(n + 1) * 2.0;
        }
    }

    // Replace coordinates for dimensions that have fine data
    // Fortran: X(n) = n - GW (1-based n); GPU 0-based: X[n] = n - GW + 1
    for (int n = cxB[0] + GW - i; n <= cxT[0] + GW - i; n++) {
        X[n] = (double)(n - GW + 1);
    }
    for (int n = cxB[1] + GW - j; n <= cxT[1] + GW - j; n++) {
        Y[n] = (double)(n - GW + 1);
    }
    for (int n = cxB[2] + GW - k; n <= cxT[2] + GW - k; n++) {
        Z[n] = (double)(n - GW + 1);
    }

    const double* C_V = VertexCoeffs<GW>::prolong();
    bool k_even = ((kk / 2) * 2 == kk);
    bool j_even = ((jj / 2) * 2 == jj);
    bool i_even = ((ii / 2) * 2 == ii);

    double result;
    double ddy;  // error estimate from polint_device (unused)

    // Case 1: x-dimension has no fine data (fg[0] == 0)
    // Collapse ya over z and y using restriction coefficients,
    // then polint in x.
    if (fg[0] == 0) {
        // Step 1: Reduce ya[STENCIL][STENCIL][STENCIL] -> tmp2[STENCIL][STENCIL]
        //         along z dimension
        double tmp2[STENCIL][STENCIL];
        #pragma unroll
        for (int sy = 0; sy < STENCIL; sy++) {
            #pragma unroll
            for (int sx = 0; sx < STENCIL; sx++) {
                double val = 0.0;
                if (k_even) {
                    #pragma unroll
                    for (int sz = 0; sz < STENCIL; sz++) {
                        int c_idx = (sz < GW) ? sz : (STENCIL - 1 - sz);
                        val += C_V[c_idx] * ya[sz][sy][sx];
                    }
                } else {
                    val = ya[GW - 1][sy][sx];
                }
                tmp2[sy][sx] = val;
            }
        }

        // Step 2: Reduce tmp2[STENCIL][STENCIL] -> tmp1[STENCIL]
        //         along y dimension
        double tmp1[STENCIL];
        #pragma unroll
        for (int sx = 0; sx < STENCIL; sx++) {
            double val = 0.0;
            if (j_even) {
                #pragma unroll
                for (int sy = 0; sy < STENCIL; sy++) {
                    int c_idx = (sy < GW) ? sy : (STENCIL - 1 - sy);
                    val += C_V[c_idx] * tmp2[sy][sx];
                }
            } else {
                val = tmp2[GW - 1][sx];
            }
            tmp1[sx] = val;
        }

        // Step 3: Replace with actual fine data at positions where fg[0] exists
        // n runs over fine data positions in the stencil: n = cxB[0]+GW-i .. cxT[0]+GW-i
        // The corresponding fine grid index is i + (n - GW)
        for (int n = cxB[0] + GW - i; n <= cxT[0] + GW - i; n++) {
            int fine_i = i + (n - GW);
            tmp1[n] = d_src_f[(k * pextf1 + j) * pextf0 + fine_i];
        }

        // Step 4: polint in x at x=0
        result = polint_device(X, tmp1, 0.0, &ddy, STENCIL);
    }
    // Case 2: y-dimension has no fine data (fg[1] == 0)
    else if (fg[1] == 0) {
        // Reduce ya over z, then x
        double tmp2[STENCIL][STENCIL];
        #pragma unroll
        for (int sy = 0; sy < STENCIL; sy++) {
            #pragma unroll
            for (int sx = 0; sx < STENCIL; sx++) {
                double val = 0.0;
                if (k_even) {
                    #pragma unroll
                    for (int sz = 0; sz < STENCIL; sz++) {
                        int c_idx = (sz < GW) ? sz : (STENCIL - 1 - sz);
                        val += C_V[c_idx] * ya[sz][sy][sx];
                    }
                } else {
                    val = ya[GW - 1][sy][sx];
                }
                tmp2[sy][sx] = val;
            }
        }

        double tmp1[STENCIL];
        #pragma unroll
        for (int sy = 0; sy < STENCIL; sy++) {
            double val = 0.0;
            if (i_even) {
                #pragma unroll
                for (int sx = 0; sx < STENCIL; sx++) {
                    int c_idx = (sx < GW) ? sx : (STENCIL - 1 - sx);
                    val += C_V[c_idx] * tmp2[sy][sx];
                }
            } else {
                val = tmp2[sy][GW - 1];
            }
            tmp1[sy] = val;
        }

        // Replace with actual fine data at positions where fg[1] exists
        for (int n = cxB[1] + GW - j; n <= cxT[1] + GW - j; n++) {
            int fine_j = j + (n - GW);
            tmp1[n] = d_src_f[(k * pextf1 + fine_j) * pextf0 + i];
        }

        result = polint_device(Y, tmp1, 0.0, &ddy, STENCIL);
    }
    // Case 3: z-dimension has no fine data (fg[2] == 0)
    else {
        // Reduce ya over x, then y
        double tmp2[STENCIL][STENCIL];
        #pragma unroll
        for (int sz = 0; sz < STENCIL; sz++) {
            #pragma unroll
            for (int sy = 0; sy < STENCIL; sy++) {
                double val = 0.0;
                if (i_even) {
                    #pragma unroll
                    for (int sx = 0; sx < STENCIL; sx++) {
                        int c_idx = (sx < GW) ? sx : (STENCIL - 1 - sx);
                        val += C_V[c_idx] * ya[sz][sy][sx];
                    }
                } else {
                    val = ya[sz][sy][GW - 1];
                }
                tmp2[sz][sy] = val;
            }
        }

        double tmp1[STENCIL];
        #pragma unroll
        for (int sz = 0; sz < STENCIL; sz++) {
            double val = 0.0;
            if (j_even) {
                #pragma unroll
                for (int sy = 0; sy < STENCIL; sy++) {
                    int c_idx = (sy < GW) ? sy : (STENCIL - 1 - sy);
                    val += C_V[c_idx] * tmp2[sz][sy];
                }
            } else {
                val = tmp2[sz][GW - 1];
            }
            tmp1[sz] = val;
        }

        // Replace with actual fine data at positions where fg[2] exists
        for (int n = cxB[2] + GW - k; n <= cxT[2] + GW - k; n++) {
            int fine_k = k + (n - GW);
            tmp1[n] = d_src_f[(fine_k * pextf1 + j) * pextf0 + i];
        }

        result = polint_device(Z, tmp1, 0.0, &ddy, STENCIL);
    }

    d_dst_f[(k * pextf1 + j) * pextf0 + i] = result;
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
    int pextc0, int pextc1,
    int extf0, int extf1, int extf2,
    FusedVars vars,
    double* __restrict__ d_dst_1d,
    int Symmetry
) {
    constexpr int STENCIL = VertexCoeffs<GW>::STENCIL;

    constexpr int SMEM_Z = 2 + 2 * GW;
    constexpr int SMEM_Y = 4 + 2 * GW;
    constexpr int SMEM_X = 4 + 2 * GW;
    constexpr int SMEM_TOTAL = SMEM_Z * SMEM_Y * SMEM_X;

    extern __shared__ double smem_raw[];
    double (*smem)[SMEM_Y][SMEM_X] =
        reinterpret_cast<double (*)[SMEM_Y][SMEM_X]>(smem_raw);

    int tid = threadIdx.z * blockDim.y * blockDim.x
            + threadIdx.y * blockDim.x
            + threadIdx.x;

    int block_i_start = i_start + blockIdx.x * blockDim.x;
    int block_j_start = j_start + blockIdx.y * blockDim.y;
    int block_k_start = k_start + blockIdx.z * blockDim.z;

    int base_c_i = (block_i_start + 1 + lbf0) / 2 - lbc0 - (GW - 1);
    int base_c_j = (block_j_start + 1 + lbf1) / 2 - lbc1 - (GW - 1);
    int base_c_k = (block_k_start + 1 + lbf2) / 2 - lbc2 - (GW - 1);

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

        int cxI_i_0b = (i + 1 + lbf0) / 2 - lbc0;
        int cxI_j_0b = (j + 1 + lbf1) / 2 - lbc1;
        int cxI_k_0b = (k + 1 + lbf2) / 2 - lbc2;

        k_even = ((kk / 2) * 2 == kk);
        j_even = ((jj / 2) * 2 == jj);
        i_even = ((ii / 2) * 2 == ii);

        smem_i_base = cxI_i_0b - (GW - 1) - base_c_i;
        smem_j_base = cxI_j_0b - (GW - 1) - base_c_j;
        smem_k_base = cxI_k_0b - (GW - 1) - base_c_k;

        out_idx = (k * extf1 + j) * extf0 + i;
    }

    const double* C_V = VertexCoeffs<GW>::prolong();

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

        // Fully-folded 3D interpolation with vertex parity logic
        if (is_valid) {
            double final_val = 0.0;

            #pragma unroll
            for (int step_x = 0; step_x < STENCIL; step_x++) {
                double coeff_x;
                if (i_even) {
                    int c_idx = (step_x < GW) ? step_x : (STENCIL - 1 - step_x);
                    coeff_x = C_V[c_idx];
                } else {
                    coeff_x = (step_x == GW - 1) ? 1.0 : 0.0;
                }
                int cur_i = smem_i_base + step_x;

                double val_y = 0.0;
                #pragma unroll
                for (int step_y = 0; step_y < STENCIL; step_y++) {
                    double coeff_y;
                    if (j_even) {
                        int c_idx = (step_y < GW) ? step_y : (STENCIL - 1 - step_y);
                        coeff_y = C_V[c_idx];
                    } else {
                        coeff_y = (step_y == GW - 1) ? 1.0 : 0.0;
                    }
                    int cur_j = smem_j_base + step_y;

                    double val_z = 0.0;
                    #pragma unroll
                    for (int step_z = 0; step_z < STENCIL; step_z++) {
                        double coeff_z;
                        if (k_even) {
                            int c_idx = (step_z < GW) ? step_z : (STENCIL - 1 - step_z);
                            coeff_z = C_V[c_idx];
                        } else {
                            coeff_z = (step_z == GW - 1) ? 1.0 : 0.0;
                        }
                        val_z += coeff_z * smem[smem_k_base + step_z][cur_j][cur_i];
                    }
                    val_y += coeff_y * val_z;
                }
                final_val += coeff_x * val_y;
            }
            d_dst_1d[v * total_per_var + out_idx] = final_val;
        }
        __syncthreads();
    }
}

// -------------------------------------------------------------------
//  fused_restrict_kernel
//
//  Same as restrict_kernel but loops over vars.num_vars.
//  Direct copy per variable, no interpolation needed.
// -------------------------------------------------------------------
__global__ void fused_restrict_kernel(
    int ni, int nj, int nk,
    int i_start, int j_start, int k_start,
    int lbc0, int lbc1, int lbc2,
    int lbf0, int lbf1, int lbf2,
    int pextf0, int pextf1,
    int extc0, int extc1, int extc2,
    FusedVars vars,
    double* __restrict__ d_dst_1d,
    int Symmetry
) {
    int i_local = blockIdx.x * blockDim.x + threadIdx.x;
    int j_local = blockIdx.y * blockDim.y + threadIdx.y;
    int k_local = blockIdx.z * blockDim.z + threadIdx.z;

    if (i_local >= ni || j_local >= nj || k_local >= nk) return;

    int i = i_start + i_local;
    int j = j_start + j_local;
    int k = k_start + k_local;

    int fine_i = 2 * (i + lbc0) - lbf0 - 1;
    int fine_j = 2 * (j + lbc1) - lbf1 - 1;
    int fine_k = 2 * (k + lbc2) - lbf2 - 1;

    int total_per_var = extc0 * extc1 * extc2;
    int out_idx = (k * extc1 + j) * extc0 + i;

    for (int v = 0; v < vars.num_vars; ++v) {
        d_dst_1d[v * total_per_var + out_idx] =
            vars.ptrs[v][(fine_k * pextf1 + fine_j) * pextf0 + fine_i];
    }
}

// ============================================================
//  5. Template Launch Functions (host side)
// ============================================================

// -------------------------------------------------------------------
//  gpu_prolong_launch_impl<GW>
//
//  Host launch function for vertex prolongation.
//  Uses vertex CD/FD: CD = (uubc-llbc)/(extc-1), FD = CD/2.
// -------------------------------------------------------------------
template<int GW>
void gpu_prolong_launch_impl(
    cudaStream_t stream,
    const double* d_src_c, double* d_dst_f,
    const double* llbc, const double* uubc, const int* extc, const int* pextc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
) {
    VertexProlongIndex idx = compute_vertex_prolong_index(
        llbc, uubc, extc, llbf, uubf, extf, llbt, uubt);

    if (idx.ni <= 0 || idx.nj <= 0 || idx.nk <= 0) return;

    dim3 block(8, 8, 4);
    dim3 grid((idx.ni + block.x - 1) / block.x,
              (idx.nj + block.y - 1) / block.y,
              (idx.nk + block.z - 1) / block.z);

    constexpr size_t shmem_bytes =
        (2 + 2 * GW) * (4 + 2 * GW) * (4 + 2 * GW) * sizeof(double);

    prolong_kernel<GW><<<grid, block, shmem_bytes, stream>>>(
        idx.ni, idx.nj, idx.nk,
        idx.i_start, idx.j_start, idx.k_start,
        idx.lbc[0], idx.lbc[1], idx.lbc[2],
        idx.lbf[0], idx.lbf[1], idx.lbf[2],
        pextc[0], pextc[1],
        pextf[0], pextf[1],
        d_src_c, d_dst_f,
        Symmetry
    );
}

// -------------------------------------------------------------------
//  gpu_restrict_launch_impl
//
//  Host launch function for vertex restriction (direct copy).
//  Uses vertex CD/FD: FD = (uubf-llbf)/(extf-1), CD = 2*FD.
// -------------------------------------------------------------------
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
        FD[d] = (uubf[d] - llbf[d]) / (double)(extf[d] - 1);
        CD[d] = 2.0 * FD[d];

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
        ubr[d] = (int)std::trunc((uubt[d] - base[d]) / CD[d] + 0.4) + 1;
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

    restrict_kernel<<<grid, block, 0, stream>>>(
        ni, nj, nk,
        i_start, j_start, k_start,
        lbc[0], lbc[1], lbc[2],
        lbf[0], lbf[1], lbf[2],
        pextc[0], pextc[1],
        pextf[0], pextf[1],
        d_src_f, d_dst_c,
        Symmetry
    );
}

// -------------------------------------------------------------------
//  gpu_prolongcopy_launch_impl
//
//  Host launch for vertex prolongcopy (block copy with ghost zones).
//  Uses vertex CD/FD: CD = (uubc-llbc)/(extc-1), FD = CD/2.
// -------------------------------------------------------------------
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
        CD[d] = (uubc[d] - llbc[d]) / (double)(extc[d] - 1);
        FD[d] = CD[d] / 2.0;

        if (llbc[d] <= llbf[d]) {
            base[d] = llbc[d];
        } else {
            int j_val = (int)std::trunc((llbc[d] - llbf[d]) / FD[d] + 0.4);
            if ((j_val / 2) * 2 == j_val) base[d] = llbf[d];
            else base[d] = llbf[d] - CD[d] / 2.0;
        }

        lbf[d]  = (int)std::trunc((llbf[d] - base[d]) / FD[d] + 0.4) + 1;
        ubf[d]  = (int)std::trunc((uubf[d] - base[d]) / FD[d] + 0.4) + 1;
        lbc[d]  = (int)std::trunc((llbc[d] - base[d]) / CD[d] + 0.4) + 1;
        ubc[d]  = (int)std::trunc((uubc[d] - base[d]) / CD[d] + 0.4) + 1;
        lbp[d]  = (int)std::trunc((llbt[d] - base[d]) / FD[d] + 0.4) + 1;
        lbpc[d] = (int)std::trunc((llbt[d] - base[d]) / CD[d] + 0.4) + 1;
        ubp[d]  = (int)std::trunc((uubt[d] - base[d]) / FD[d] + 0.4) + 1;
        ubpc[d] = (int)std::trunc((uubt[d] - base[d]) / CD[d] + 0.4) + 1;
    }

    // Source region in coarse grid (with ghost zones)
    int imini = lbpc[0] - lbc[0] + 1 - ghost_width;
    int imaxi = ubpc[0] - lbc[0] + 1 + ghost_width;
    int jmini = lbpc[1] - lbc[1] + 1 - ghost_width;
    int jmaxi = ubpc[1] - lbc[1] + 1 + ghost_width;
    int kmini = lbpc[2] - lbc[2] + 1 - ghost_width;
    int kmaxi = ubpc[2] - lbc[2] + 1 + ghost_width;

    // Clamp to valid range
    if (imini < 1) imini = 1;
    if (jmini < 1) jmini = 1;
    if (kmini < 1) kmini = 1;
    if (imaxi > extc[0]) imaxi = extc[0];
    if (jmaxi > extc[1]) jmaxi = extc[1];
    if (kmaxi > extc[2]) kmaxi = extc[2];

    // Destination region in output (fine) grid
    // Convert: imini in coarse -> (imini + lbc - 1) in coarse coords -> fine index
    // dst offset = 1-based from 1 in fine output
    int dst_x = 1;  // destination starts at 1
    int dst_y = 1;
    int dst_z = 1;

    int ni = imaxi - imini + 1;
    int nj = jmaxi - jmini + 1;
    int nk = kmaxi - kmini + 1;

    if (ni <= 0 || nj <= 0 || nk <= 0) return;

    dim3 block(8, 8, 4);
    dim3 grid((ni + block.x - 1) / block.x,
              (nj + block.y - 1) / block.y,
              (nk + block.z - 1) / block.z);

    // Convert from 1-based Fortran to 0-based C indices
    // imini, jmini, kmini are 1-based Fortran indices into extc array
    // dst_x, dst_y, dst_z are 1-based Fortran indices into output array
    // Kernel expects 0-based indices
    prolongcopy_kernel<<<grid, block, 0, stream>>>(
        ni, nj, nk,
        imini - 1, jmini - 1, kmini - 1,  // 0-based src offset
        dst_x - 1, dst_y - 1, dst_z - 1,  // 0-based dst offset
        pextc[0], pextc[1],
        pextf[0], pextf[1],
        d_src_c, d_dst_f
    );
}

// -------------------------------------------------------------------
//  gpu_prolongmix_launch_impl<GW>
//
//  Host launch for vertex prolongmix (mixed interpolation).
//  Uses vertex CD/FD: FD = (uubf-llbf)/(extf-1), CD = 2*FD.
// -------------------------------------------------------------------
template<int GW>
void gpu_prolongmix_launch_impl(
    cudaStream_t stream,
    const double* d_src_f,   // funf (fine grid, in/out)
    const double* d_src_i,   // funi (coarse stencil data with ghost zones)
    double* d_dst_f,         // funf output (same as d_src_f for in-place)
    const double* llbc, const double* uubc,
    const double* llbf, const double* uubf, const int* extf,
    const int* exti,         // exti = interior coarse extent (kernel adds 2*GW for padding)
    const double* llbp, const double* uubp,
    const double* illb, const double* iuub,  // fine grid inner bounds
    const double* SoA, int Symmetry
) {
    double CD[3], FD[3], base[3];
    int lbc[3], ubc[3], lbf[3], ubf[3];
    int lbp[3], ubp[3], lbpc[3], ubpc[3];
    int ilb[3], iub[3];

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

        lbf[d]  = (int)std::trunc((llbf[d] - base[d]) / FD[d] + 0.4) + 1;
        ubf[d]  = (int)std::trunc((uubf[d] - base[d]) / FD[d] + 0.4) + 1;
        lbc[d]  = (int)std::trunc((llbc[d] - base[d]) / CD[d] + 0.4) + 1;
        ubc[d]  = (int)std::trunc((uubc[d] - base[d]) / CD[d] + 0.4) + 1;
        lbp[d]  = (int)std::trunc((llbp[d] - base[d]) / FD[d] + 0.4) + 1;
        lbpc[d] = (int)std::trunc((llbp[d] - base[d]) / CD[d] + 0.4) + 1;
        ubp[d]  = (int)std::trunc((uubp[d] - base[d]) / FD[d] + 0.4) + 1;
        ubpc[d] = (int)std::trunc((uubp[d] - base[d]) / CD[d] + 0.4) + 1;

        ilb[d] = (int)std::trunc((illb[d] - base[d]) / FD[d] + 0.4) + 1;
        iub[d] = (int)std::trunc((iuub[d] - base[d]) / FD[d] + 0.4) + 1;
    }

    // Output iteration range
    int imino = lbp[0] - lbf[0] + 1;
    int imaxo = ubp[0] - lbf[0] + 1;
    int jmino = lbp[1] - lbf[1] + 1;
    int jmaxo = ubp[1] - lbf[1] + 1;
    int kmino = lbp[2] - lbf[2] + 1;
    int kmaxo = ubp[2] - lbf[2] + 1;

    // Ghost zone for coarse array reference
    lbpc[0] -= ghost_width; lbpc[1] -= ghost_width; lbpc[2] -= ghost_width;
    ubpc[0] += ghost_width; ubpc[1] += ghost_width; ubpc[2] += ghost_width;

    // Convert ilb/iub to fine-grid array indices
    ilb[0] = ilb[0] - lbf[0] + 1;
    ilb[1] = ilb[1] - lbf[1] + 1;
    ilb[2] = ilb[2] - lbf[2] + 1;
    iub[0] = iub[0] - lbf[0] + 1;
    iub[1] = iub[1] - lbf[1] + 1;
    iub[2] = iub[2] - lbf[2] + 1;

    // Clamp to valid range (domain division by parallelization)
    ilb[0] = (ilb[0] > 1) ? ilb[0] : 1;
    ilb[1] = (ilb[1] > 1) ? ilb[1] : 1;
    ilb[2] = (ilb[2] > 1) ? ilb[2] : 1;
    iub[0] = (iub[0] < extf[0]) ? iub[0] : extf[0];
    iub[1] = (iub[1] < extf[1]) ? iub[1] : extf[1];
    iub[2] = (iub[2] < extf[2]) ? iub[2] : extf[2];

    // Convert from 1-based Fortran to 0-based C indexing
    // Output range
    int i_start = imino - 1;  // 0-based
    int j_start = jmino - 1;
    int k_start = kmino - 1;
    int ni = imaxo - imino + 1;
    int nj = jmaxo - jmino + 1;
    int nk = kmaxo - kmino + 1;

    if (ni <= 0 || nj <= 0 || nk <= 0) return;

    dim3 block(8, 8, 4);
    dim3 grid((ni + block.x - 1) / block.x,
              (nj + block.y - 1) / block.y,
              (nk + block.z - 1) / block.z);

    prolongmix_kernel<GW><<<grid, block, 0, stream>>>(
        ni, nj, nk,
        i_start, j_start, k_start,
        lbc[0], lbc[1], lbc[2],       // lbc (1-based Fortran)
        lbf[0], lbf[1], lbf[2],       // lbf (1-based Fortran)
        lbpc[0], lbpc[1], lbpc[2],             // lbpc (adjusted, same as cell)
        ilb[0] - 1, ilb[1] - 1, ilb[2] - 1,     // ilb (0-based C)
        iub[0] - 1, iub[1] - 1, iub[2] - 1,     // iub (0-based C)
        exti[0] + 2*GW, exti[1] + 2*GW,     // pextc (padded coarse dims)
        extf[0] + 2*GW, extf[1] + 2*GW,     // pextf (padded fine dims)
        exti[0] + 2*GW, exti[1] + 2*GW,     // pexti (padded funi strides)
        d_src_f, d_src_i, d_dst_f,
        Symmetry
    );
}

// -------------------------------------------------------------------
//  gpu_fused_prolong_launch_impl<GW>
// -------------------------------------------------------------------
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

    VertexProlongIndex idx = compute_vertex_prolong_index(
        llbc, uubc, extc, llbf, uubf, extf, llbt, uubt);

    if (idx.ni <= 0 || idx.nj <= 0 || idx.nk <= 0) return;

    dim3 block(8, 8, 4);
    dim3 grid((idx.ni + block.x - 1) / block.x,
              (idx.nj + block.y - 1) / block.y,
              (idx.nk + block.z - 1) / block.z);

    constexpr size_t shmem_bytes =
        (2 + 2 * GW) * (4 + 2 * GW) * (4 + 2 * GW) * sizeof(double);

    fused_prolong_kernel<GW><<<grid, block, shmem_bytes, stream>>>(
        idx.ni, idx.nj, idx.nk,
        idx.i_start, idx.j_start, idx.k_start,
        idx.lbc[0], idx.lbc[1], idx.lbc[2],
        idx.lbf[0], idx.lbf[1], idx.lbf[2],
        pextc[0], pextc[1],
        extf[0], extf[1], extf[2],
        vars, d_dst_1d, Symmetry
    );
}

// -------------------------------------------------------------------
//  gpu_fused_restrict_launch_impl
// -------------------------------------------------------------------
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
        FD[d] = (uubf[d] - llbf[d]) / (double)(extf[d] - 1);
        CD[d] = 2.0 * FD[d];

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
        ubr[d] = (int)std::trunc((uubt[d] - base[d]) / CD[d] + 0.4) + 1;
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

    fused_restrict_kernel<<<grid, block, 0, stream>>>(
        ni, nj, nk,
        i_start, j_start, k_start,
        lbc[0], lbc[1], lbc[2],
        lbf[0], lbf[1], lbf[2],
        pextf[0], pextf[1],
        extc[0], extc[1], extc[2],
        vars, d_dst_1d, Symmetry
    );
}

// ============================================================
//  6. API Wrapper Functions
//
//  Compile-time dispatch via ghost_width macro (from macrodef.fh).
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

void gpu_prolongmix3_launch(
    cudaStream_t stream,
    const double* d_src_f,   // funf
    const double* d_src_i,   // funi
    double* d_dst_f,         // funf output
    const double* llbc, const double* uubc,
    const double* llbf, const double* uubf, const int* extf, const int* pextf,
    const int* exti,
    const double* llbp, const double* uubp,
    const double* illb, const double* iuub,
    const double* SoA, int Symmetry
) {
    gpu_prolongmix_launch_impl<ghost_width>(stream,
        d_src_f, d_src_i, d_dst_f, llbc, uubc, llbf, uubf, extf,
        exti, llbp, uubp, illb, iuub, SoA, Symmetry);
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

#endif // Vertex
