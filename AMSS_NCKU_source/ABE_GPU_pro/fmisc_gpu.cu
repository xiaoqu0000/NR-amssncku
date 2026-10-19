#include "macrodef.h"
#include <cuda_runtime.h>
#include <math.h>

#include <iostream>
#include "gpu_manager.h"
#include "gpu_common.cuh"
#include "fmisc_gpu.cuh"
#include "gpu_fusion.h"
#include "fmisc.h"
#include "gpu_info.h"
#include "MyList.h"
#include "var.h"
#include "MPatch.h"

__global__ void lowerboundset_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    double* chi0, double TINNY
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = i + j * pex0 + k * pex0 * pex1;
    if (chi0[idx] < TINNY) chi0[idx] = TINNY;
}

void gpu_lowerboundset_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    double* d_chi0, double TINNY
) {
    int total = ex[0] * ex[1] * ex[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);

    lowerboundset_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], d_chi0, TINNY
    );
}

__global__ void gpu_pack_kernel(
    const double* __restrict__ src_3d, double* __restrict__ dst_1d,
    int src_pnx, int src_pny,
    int dst_nx, int dst_ny, int dst_nz,
    int off_x, int off_y, int off_z
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= dst_nx * dst_ny * dst_nz) return;
    int i, j, k;
    linear_to_ijk(tid, dst_nx, dst_ny, i, j, k);

    int src_idx = (i + off_x)
                   + src_pnx * ((j + off_y)
                   + src_pny * (k + off_z));

    int dst_idx = i
                   + dst_nx * (j + dst_ny * k);

    dst_1d[dst_idx] = src_3d[src_idx];
}

__global__ void gpu_unpack_kernel(
    const double* __restrict__ src_1d, double* __restrict__ dst_3d,
    int dst_pnx, int dst_pny,
    int src_nx, int src_ny, int src_nz,
    int off_x, int off_y, int off_z
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= src_nx * src_ny * src_nz) return;
    int i, j, k;
    linear_to_ijk(tid, src_nx, src_ny, i, j, k);

    int src_idx = i
                   + src_nx * (j + src_ny * k);

    int dst_idx = (i + off_x)
                   + dst_pnx * ((j + off_y)
                   + dst_pny * (k + off_z));

    dst_3d[dst_idx] = src_1d[src_idx];
}

void gpu_pack_launch(
    cudaStream_t stream, const double* d_src_3d, double* d_dst_1d,
    int src_pnx, int src_pny, int dst_nx, int dst_ny, int dst_nz,
    int off_x, int off_y, int off_z
) {
    int total = dst_nx * dst_ny * dst_nz;
    dim3 block(256);
    dim3 grid((total + 255) / 256);

    gpu_pack_kernel<<<grid, block, 0, stream>>>(
        d_src_3d, d_dst_1d, src_pnx, src_pny, dst_nx, dst_ny, dst_nz, off_x, off_y, off_z
    );
}

void gpu_unpack_launch(
    cudaStream_t stream, const double* d_src_1d, double* d_dst_3d,
    int dst_pnx, int dst_pny, int src_nx, int src_ny, int src_nz,
    int off_x, int off_y, int off_z
) {
    int total = src_nx * src_ny * src_nz;
    dim3 block(256);
    dim3 grid((total + 255) / 256);

    gpu_unpack_kernel<<<grid, block, 0, stream>>>(
        d_src_1d, d_dst_3d, dst_pnx, dst_pny, src_nx, src_ny, src_nz, off_x, off_y, off_z
    );
}

__global__ void set_value_kernel(int ex0, int ex1, int ex2, int pex0, int pex1, double* data, double value) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    data[(k * pex1 + j) * pex0 + i] = value;
}

__global__ void add_value_kernel(int ex0, int ex1, int ex2, int pex0, int pex1, double* data, double value) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    data[(k * pex1 + j) * pex0 + i] += value;
}

__global__ void array_copy_kernel(int ex0, int ex1, int ex2, int pex0, int pex1, double* dst, const double* src) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = (k * pex1 + j) * pex0 + i;
    dst[idx] = src[idx];
}

__global__ void array_add_kernel(int ex0, int ex1, int ex2, int pex0, int pex1, double* dst, const double* src) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = (k * pex1 + j) * pex0 + i;
    dst[idx] += src[idx];
}

__global__ void array_subtract_kernel(int ex0, int ex1, int ex2, int pex0, int pex1, double* dst, const double* src) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = (k * pex1 + j) * pex0 + i;
    dst[idx] -= src[idx];
}

void gpu_set_value_launch(cudaStream_t stream, const int ex[3], const int pex[3], double* d_data, double value) {
    int total = ex[0] * ex[1] * ex[2];
    if (total <= 0) return;
    dim3 block(256), grid((total + 255) / 256);
    set_value_kernel<<<grid, block, 0, stream>>>(ex[0], ex[1], ex[2], pex[0], pex[1], d_data, value);
}

void gpu_add_value_launch(cudaStream_t stream, const int ex[3], const int pex[3], double* d_data, double value) {
    int total = ex[0] * ex[1] * ex[2];
    if (total <= 0) return;
    dim3 block(256), grid((total + 255) / 256);
    add_value_kernel<<<grid, block, 0, stream>>>(ex[0], ex[1], ex[2], pex[0], pex[1], d_data, value);
}

void gpu_array_copy_launch(cudaStream_t stream, const int ex[3], const int pex[3], double* d_dst, const double* d_src) {
    int total = ex[0] * ex[1] * ex[2];
    if (total <= 0) return;
    dim3 block(256), grid((total + 255) / 256);
    array_copy_kernel<<<grid, block, 0, stream>>>(ex[0], ex[1], ex[2], pex[0], pex[1], d_dst, d_src);
}

void gpu_array_add_launch(cudaStream_t stream, const int ex[3], const int pex[3], double* d_dst, const double* d_src) {
    int total = ex[0] * ex[1] * ex[2];
    if (total <= 0) return;
    dim3 block(256), grid((total + 255) / 256);
    array_add_kernel<<<grid, block, 0, stream>>>(ex[0], ex[1], ex[2], pex[0], pex[1], d_dst, d_src);
}

void gpu_array_subtract_launch(cudaStream_t stream, const int ex[3], const int pex[3], double* d_dst, const double* d_src) {
    int total = ex[0] * ex[1] * ex[2];
    if (total <= 0) return;
    dim3 block(256), grid((total + 255) / 256);
    array_subtract_kernel<<<grid, block, 0, stream>>>(ex[0], ex[1], ex[2], pex[0], pex[1], d_dst, d_src);
}

__global__ void scan_nan_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* X, const double* Y, const double* Z,
    const double* f, int* count, double* first_xyzf
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = (k * pex1 + j) * pex0 + i;
    double v = f[idx];
    if (isnan(v)) {
        int old = atomicAdd(count, 1);
        if (old == 0) {
            first_xyzf[0] = X[i];
            first_xyzf[1] = Y[j];
            first_xyzf[2] = Z[k];
            first_xyzf[3] = v;
        }
    }
}

void gpu_scan_nan_launch(
    cudaStream_t stream, const int ex[3], const int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_f, int* d_count, double* d_first_xyzf
) {
    int total = ex[0] * ex[1] * ex[2];
    cudaMemsetAsync(d_count, 0, sizeof(int), stream);
    cudaMemsetAsync(d_first_xyzf, 0, 4 * sizeof(double), stream);
    if (total <= 0) return;
    dim3 block(256), grid((total + 255) / 256);
    scan_nan_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], d_X, d_Y, d_Z, d_f, d_count, d_first_xyzf
    );
}

__global__ void pointcopy_cell_kernel(
    double* __restrict__ data_out,
    int ext0, int ext1, int ext2,
    int pext0, int pext1,
    double llb0, double llb1, double llb2,
    double uub0, double uub1, double uub2,
    double xx, double yy, double zz, double dv
) {
#ifdef Vertex
    double h0 = (uub0 - llb0) / (double)(ext0 - 1);
    double h1 = (uub1 - llb1) / (double)(ext1 - 1);
    double h2 = (uub2 - llb2) / (double)(ext2 - 1);
#else
    double h0 = (uub0 - llb0) / (double)ext0;
    double h1 = (uub1 - llb1) / (double)ext1;
    double h2 = (uub2 - llb2) / (double)ext2;
#endif
    int i = (int)((xx - llb0) / h0 + 0.4);
    int j = (int)((yy - llb1) / h1 + 0.4);
    int k = (int)((zz - llb2) / h2 + 0.4);
    if (i < 0 || i >= ext0 || j < 0 || j >= ext1 || k < 0 || k >= ext2) return;
    data_out[(k * pext1 + j) * pext0 + i] = dv;
}

void gpu_pointcopy_cell_launch(
    cudaStream_t stream,
    const int ext_out[3], const int pext_out[3],
    double* d_data_out,
    const double llbout[3], const double uubout[3],
    double xx, double yy, double zz, double dv
) {
    pointcopy_cell_kernel<<<1, 1, 0, stream>>>(
        d_data_out,
        ext_out[0], ext_out[1], ext_out[2],
        pext_out[0], pext_out[1],
        llbout[0], llbout[1], llbout[2],
        uubout[0], uubout[1], uubout[2],
        xx, yy, zz, dv
    );
}

__global__ void copy_cell_kernel(
    double* __restrict__ data_out,
    const double* __restrict__ data_in,
    int out0, int out1, int out2,
    int pout0, int pout1,
    int in0, int in1, int in2,
    int pin0, int pin1,
    double llbo0, double llbo1, double llbo2,
    double uubo0, double uubo1, double uubo2,
    double llbi0, double llbi1, double llbi2,
    double uubi0, double uubi1, double uubi2,
    double lcopy0, double lcopy1, double lcopy2,
    double ucopy0, double ucopy1, double ucopy2
) {
#ifdef Vertex
    double ho0 = (uubo0 - llbo0) / (double)(out0 - 1);
    double ho1 = (uubo1 - llbo1) / (double)(out1 - 1);
    double ho2 = (uubo2 - llbo2) / (double)(out2 - 1);
    double hi0 = (uubi0 - llbi0) / (double)(in0 - 1);
    double hi1 = (uubi1 - llbi1) / (double)(in1 - 1);
    double hi2 = (uubi2 - llbi2) / (double)(in2 - 1);
#else
    double ho0 = (uubo0 - llbo0) / (double)out0;
    double ho1 = (uubo1 - llbo1) / (double)out1;
    double ho2 = (uubo2 - llbo2) / (double)out2;
    double hi0 = (uubi0 - llbi0) / (double)in0;
    double hi1 = (uubi1 - llbi1) / (double)in1;
    double hi2 = (uubi2 - llbi2) / (double)in2;
#endif

    int illo0 = (int)((lcopy0 - llbo0) / ho0 + 0.4);
    int illo1 = (int)((lcopy1 - llbo1) / ho1 + 0.4);
    int illo2 = (int)((lcopy2 - llbo2) / ho2 + 0.4);
    int iuuo0 = out0 - 1 - (int)((uubo0 - ucopy0) / ho0 + 0.4);
    int iuuo1 = out1 - 1 - (int)((uubo1 - ucopy1) / ho1 + 0.4);
    int iuuo2 = out2 - 1 - (int)((uubo2 - ucopy2) / ho2 + 0.4);
    int illi0 = (int)((lcopy0 - llbi0) / hi0 + 0.4);
    int illi1 = (int)((lcopy1 - llbi1) / hi1 + 0.4);
    int illi2 = (int)((lcopy2 - llbi2) / hi2 + 0.4);

    int nx = iuuo0 - illo0 + 1;
    int ny = iuuo1 - illo1 + 1;
    int nz = iuuo2 - illo2 + 1;
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (nx <= 0 || ny <= 0 || nz <= 0 || tid >= nx * ny * nz) return;

    int i, j, k;
    linear_to_ijk(tid, nx, ny, i, j, k);
    int oi = illo0 + i;
    int oj = illo1 + j;
    int ok = illo2 + k;
    int ii = illi0 + i;
    int ij = illi1 + j;
    int ik = illi2 + k;
    if (oi < 0 || oi >= out0 || oj < 0 || oj >= out1 || ok < 0 || ok >= out2) return;
    if (ii < 0 || ii >= in0 || ij < 0 || ij >= in1 || ik < 0 || ik >= in2) return;
    data_out[(ok * pout1 + oj) * pout0 + oi] = data_in[(ik * pin1 + ij) * pin0 + ii];
}

void gpu_copy_cell_launch(
    cudaStream_t stream,
    const int ext_out[3], const int pext_out[3], double* d_data_out,
    const double llbout[3], const double uubout[3],
    const int ext_in[3], const int pext_in[3], const double* d_data_in,
    const double llbin[3], const double uubin[3],
    const double lcopy[3], const double ucopy[3]
) {
#ifdef Vertex
    int n0 = (int)((ucopy[0] - lcopy[0]) / ((uubout[0] - llbout[0]) / (double)(ext_out[0] - 1)) + 1.4);
    int n1 = (int)((ucopy[1] - lcopy[1]) / ((uubout[1] - llbout[1]) / (double)(ext_out[1] - 1)) + 1.4);
    int n2 = (int)((ucopy[2] - lcopy[2]) / ((uubout[2] - llbout[2]) / (double)(ext_out[2] - 1)) + 1.4);
#else
    int n0 = (int)((ucopy[0] - lcopy[0]) / ((uubout[0] - llbout[0]) / (double)ext_out[0]) + 1.4);
    int n1 = (int)((ucopy[1] - lcopy[1]) / ((uubout[1] - llbout[1]) / (double)ext_out[1]) + 1.4);
    int n2 = (int)((ucopy[2] - lcopy[2]) / ((uubout[2] - llbout[2]) / (double)ext_out[2]) + 1.4);
#endif
    
    int total = n0 * n1 * n2;
    if (total <= 0) return;
    dim3 block(256);
    dim3 grid((total + 255) / 256);
    copy_cell_kernel<<<grid, block, 0, stream>>>(
        d_data_out, d_data_in,
        ext_out[0], ext_out[1], ext_out[2], pext_out[0], pext_out[1],
        ext_in[0], ext_in[1], ext_in[2], pext_in[0], pext_in[1],
        llbout[0], llbout[1], llbout[2], uubout[0], uubout[1], uubout[2],
        llbin[0], llbin[1], llbin[2], uubin[0], uubin[1], uubin[2],
        lcopy[0], lcopy[1], lcopy[2], ucopy[0], ucopy[1], ucopy[2]
    );
}

// =====================================================================
// Time Level Interpolation Kernels
// =====================================================================

__global__ void average_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* f1, const double* f2, double* fout
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = i + j * pex0 + k * pex0 * pex1;
    fout[idx] = 0.5 * (f1[idx] + f2[idx]);
}

__global__ void average3_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* f1, const double* f2, double* fout
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = i + j * pex0 + k * pex0 * pex1;
    fout[idx] = 0.75 * f1[idx] + 0.25 * f2[idx];
}

__global__ void average2_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* f1, const double* f2, const double* f3, double* fout
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = i + j * pex0 + k * pex0 * pex1;
    fout[idx] = (3.0 / 8.0) * f1[idx] + (3.0 / 4.0) * f2[idx] - (1.0 / 8.0) * f3[idx];
}

__global__ void average2p_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* f1, const double* f2, const double* f3, double* fout
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = i + j * pex0 + k * pex0 * pex1;
    fout[idx] = (5.0 / 32.0) * f1[idx] + (15.0 / 16.0) * f2[idx] - (3.0 / 32.0) * f3[idx];
}

__global__ void average2m_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* f1, const double* f2, const double* f3, double* fout
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = i + j * pex0 + k * pex0 * pex1;
    fout[idx] = (5.0 / 32.0) * f1[idx] + (15.0 / 16.0) * f2[idx] - (3.0 / 32.0) * f3[idx];
}

void gpu_average_launch(cudaStream_t stream, const int ex[3], const int pex[3], const double* d_f1, const double* d_f2, double* d_fout) {
    int total = ex[0] * ex[1] * ex[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);
    average_kernel<<<grid, block, 0, stream>>>(ex[0], ex[1], ex[2], pex[0], pex[1], d_f1, d_f2, d_fout);
}

void gpu_average3_launch(cudaStream_t stream, const int ex[3], const int pex[3], const double* d_f1, const double* d_f2, double* d_fout) {
    int total = ex[0] * ex[1] * ex[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);
    average3_kernel<<<grid, block, 0, stream>>>(ex[0], ex[1], ex[2], pex[0], pex[1], d_f1, d_f2, d_fout);
}

void gpu_average2_launch(cudaStream_t stream, const int ex[3], const int pex[3], const double* d_f1, const double* d_f2, const double* d_f3, double* d_fout) {
    int total = ex[0] * ex[1] * ex[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);
    average2_kernel<<<grid, block, 0, stream>>>(ex[0], ex[1], ex[2], pex[0], pex[1], d_f1, d_f2, d_f3, d_fout);
}

void gpu_average2p_launch(cudaStream_t stream, const int ex[3], const int pex[3], const double* d_f1, const double* d_f2, const double* d_f3, double* d_fout) {
    int total = ex[0] * ex[1] * ex[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);
    average2p_kernel<<<grid, block, 0, stream>>>(ex[0], ex[1], ex[2], pex[0], pex[1], d_f1, d_f2, d_f3, d_fout);
}

void gpu_average2m_launch(cudaStream_t stream, const int ex[3], const int pex[3], const double* d_f1, const double* d_f2, const double* d_f3, double* d_fout) {
    int total = ex[0] * ex[1] * ex[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);
    average2m_kernel<<<grid, block, 0, stream>>>(ex[0], ex[1], ex[2], pex[0], pex[1], d_f1, d_f2, d_f3, d_fout);
}

__global__ void interp_2_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ f, double* __restrict__ f_int,
    int il, int jl, int kl, double Dx, double Dy, double Dz,
    int ordn, double SoA0, double SoA1, double SoA2
) {
    if (threadIdx.x || blockIdx.x) return;
    if (ordn <= 0 || ordn > MAX_ORDN) {
        *f_int = NAN;
        gpu_stop();
        return;
    }
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const double SoA[3] = {SoA0, SoA1, SoA2};
    double wx[MAX_ORDN], wy[MAX_ORDN], wz[MAX_ORDN];
    compute_lagrange_weights(Dx, ordn, wx);
    compute_lagrange_weights(Dy, ordn, wy);
    compute_lagrange_weights(Dz, ordn, wz);

    double val = 0.0;
    for (int k = 0; k < ordn; ++k) {
        for (int j = 0; j < ordn; ++j) {
            for (int i = 0; i < ordn; ++i) {
                val += wx[i] * wy[j] * wz[k] *
                       fetch_f_sym_0b_pad(f, ex, pex, il + i, jl + j, kl + k, SoA);
            }
        }
    }
    *f_int = val;
}

void gpu_interp_2_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_f, double* d_f_int,
    int il, int jl, int kl,
    double Dx, double Dy, double Dz,
    int ordn, const double SoA[3]
) {
    interp_2_kernel<<<1, 1, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_f, d_f_int, il, jl, kl, Dx, Dy, Dz, ordn, SoA[0], SoA[1], SoA[2]
    );
}

__global__ void polint_kernel(const double* xa, const double* ya, double x, double* y, double* dy, int ordn) {
    if (threadIdx.x || blockIdx.x) return;
    *y = polint_device(xa, ya, x, dy, ordn);
}

void gpu_polint_launch(cudaStream_t stream, const double* d_xa, const double* d_ya, double x, double* d_y, double* d_dy, int ordn) {
    polint_kernel<<<1, 1, 0, stream>>>(d_xa, d_ya, x, d_y, d_dy, ordn);
}

__global__ void polin2_kernel(
    const double* x1a, const double* x2a, const double* ya,
    double x1, double x2, double* y, double* dy, int ordn
) {
    if (threadIdx.x || blockIdx.x) return;
    *y = polin2_device(x1a, x2a, ya, x1, x2, dy, ordn);
}

void gpu_polin2_launch(
    cudaStream_t stream,
    const double* d_x1a, const double* d_x2a, const double* d_ya,
    double x1, double x2, double* d_y, double* d_dy, int ordn
) {
    polin2_kernel<<<1, 1, 0, stream>>>(d_x1a, d_x2a, d_ya, x1, x2, d_y, d_dy, ordn);
}

__global__ void polin3_kernel(
    const double* x1a, const double* x2a, const double* x3a, const double* ya,
    double x1, double x2, double x3, double* y, double* dy, int ordn
) {
    if (threadIdx.x || blockIdx.x) return;
    *y = polin3_device(x1a, x2a, x3a, ya, x1, x2, x3, dy, ordn);
}

void gpu_polin3_launch(
    cudaStream_t stream,
    const double* d_x1a, const double* d_x2a, const double* d_x3a, const double* d_ya,
    double x1, double x2, double x3, double* d_y, double* d_dy, int ordn
) {
    polin3_kernel<<<1, 1, 0, stream>>>(d_x1a, d_x2a, d_x3a, d_ya, x1, x2, x3, d_y, d_dy, ordn);
}

__global__ void d2dump_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* __restrict__ data_in, double* __restrict__ data_out,
    double llb2, double uub2, int gord, double SoA2
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1) return;
    int i = tid % ex0;
    int j = tid / ex0;
#ifdef Vertex
    double dZ = (uub2 - llb2) / (double)(ex2 - 1);
#else
    double dZ = (uub2 - llb2) / (double)ex2;
#endif
    int k = (int)((0.0 - llb2) / dZ + 0.4);
    if (k < 0) return;

    auto at = [&](int kk) {
        int rk = kk;
        double factor = 1.0;
#ifdef Vertex
        if (kk < 0) {
            rk = -kk;
            factor = SoA2;
        }
        if (kk >= ex2) {
            rk = 2 * ex2 - kk - 2;
            factor = SoA2;
        }
#else
        if (kk < 0) {
            rk = -kk - 1; // Cell 边界穿透映射
            factor = SoA2;
        }
        if (kk >= ex2) {
            rk = 2 * ex2 - kk - 1;
            factor = SoA2;
        }
#endif
        if (rk < 0 || rk >= ex2) return 0.0;
        return factor * data_in[(rk * pex1 + j) * pex0 + i];
    };

    double out = 0.0;
#ifdef Vertex
    // Vertex 模式下，网格点直接落在边界上，无需权重插值
    out = at(k);
#else
    // Cell 模式下，通过周围点插值到物理边界
    if (gord == 2) {
        out = 0.5625 * (at(k) + at(k - 1)) - 0.0625 * (at(k + 1) + at(k - 2));
    } else if (gord == 3) {
        out = 0.5859375 * (at(k) + at(k - 1))
            - 0.09765625 * (at(k + 1) + at(k - 2))
            + 0.01171875 * (at(k + 2) + at(k - 3));
    } else if (gord == 4) {
        out = 0.5981445312 * (at(k) + at(k - 1))
            - 0.1196289063 * (at(k + 1) + at(k - 2))
            + 0.02392578125 * (at(k + 2) + at(k - 3))
            - 0.00244140625 * (at(k + 3) + at(k - 4));
    } else if (gord == 5) {
        out = 0.6056213378 * (at(k) + at(k - 1))
            - 0.1345825196 * (at(k + 1) + at(k - 2))
            + 0.03460693359 * (at(k + 2) + at(k - 3))
            - 0.006179809571 * (at(k + 3) + at(k - 4))
            + 0.0005340576171 * (at(k + 4) + at(k - 5));
    } else {
        out = at(k);
    }
#endif
    data_out[j * ex0 + i] = out;
}

void gpu_d2dump_launch(
    cudaStream_t stream, const int ex[3], const int pex[3],
    const double* d_data_in, double* d_data_out,
    const double llb[3], const double uub[3],
    int gord, const double SoA[3]
) {
    int total = ex[0] * ex[1];
    if (total <= 0) return;
    dim3 block(256), grid((total + 255) / 256);
    d2dump_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1],
        d_data_in, d_data_out, llb[2], uub[2], gord, SoA[2]
    );
}

__global__ void global_interp_kernel(
    int NN, int DIM,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    double* d_field,
    double llb_0, double llb_1, double llb_2,
    double uub_0, double uub_1, double uub_2,
    double DH_0, double DH_1, double DH_2,
    int ordn, double SoA_0, double SoA_1, double SoA_2, int Symmetry,
    int var_idx, int num_var,
    double* d_shellf, int* d_weight
) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= NN) return;

    double px = d_XX_0[j], py = d_XX_1[j], pz = d_XX_2[j];

    double tol_0 = DH_0 / 2.0;
    double tol_1 = DH_1 / 2.0;
    double tol_2 = DH_2 / 2.0;

    if (px - llb_0 < -tol_0 || px - uub_0 > tol_0) return;
    if (DIM > 1 && (py - llb_1 < -tol_1 || py - uub_1 > tol_1)) return;
    if (DIM > 2 && (pz - llb_2 < -tol_2 || pz - uub_2 > tol_2)) return;

    double SoA_arr[3] = {SoA_0, SoA_1, SoA_2};
    const int ex[3]  = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    double val = 0.0;
    global_interp_device(
        ex, pex, d_X_0, d_X_1, d_X_2,
        d_field, &val,
        px, py, pz,
        ordn, SoA_arr, Symmetry
    );

    atomicAdd(&d_shellf[j * num_var + var_idx], val);
    if (var_idx == 0) atomicAdd(&d_weight[j], 1);
}

void gpu_global_interp_launch(
    cudaStream_t stream,
    int NN, int DIM,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    double* d_field,
    double llb_0, double llb_1, double llb_2,
    double uub_0, double uub_1, double uub_2,
    double DH_0, double DH_1, double DH_2,
    int ordn, double SoA_0, double SoA_1, double SoA_2,
    int Symmetry, int var_idx, int num_var,
    double* d_shellf, int* d_weight
) {
    int blockSize = 256;
    int gridSize = (NN + blockSize - 1) / blockSize;

    global_interp_kernel<<<gridSize, blockSize, 0, stream>>>(
        NN, DIM,
        d_XX_0, d_XX_1, d_XX_2,
        shape_0, shape_1, shape_2,
        pshape_0, pshape_1, pshape_2,
        d_X_0, d_X_1, d_X_2, d_field,
        llb_0, llb_1, llb_2, uub_0, uub_1, uub_2,
        DH_0, DH_1, DH_2,
        ordn, SoA_0, SoA_1, SoA_2, Symmetry,
        var_idx, num_var, d_shellf, d_weight
    );
}

// ============================================================================
// batch_interp_ss_kernel + gpu_batch_interp_ss_launch (ShellPatch Interp_Points_GPU)
// N-point batch version for ShellPatch::Interp_Points_GPU
// ============================================================================

__global__ void batch_interp_ss_kernel(
    int NN, int DIM, int sst,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    double* d_field,
    double llb_0, double llb_1, double llb_2,
    double uub_0, double uub_1, double uub_2,
    double DH_0, double DH_1, double DH_2,
    int ordn, double SoA_0, double SoA_1, double SoA_2, int Symmetry,
    int var_idx, int num_var,
    double* d_shellf, int* d_weight,
    int* d_point_sst
) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= NN) return;

    // Only interpolate if this point belongs to the current sector
    if (d_point_sst[j] != sst) return;

    double px = d_XX_0[j], py = d_XX_1[j], pz = d_XX_2[j];
    if (px - llb_0 < -DH_0 / 2.0 || px - uub_0 > DH_0 / 2.0) return;
    if (DIM > 1 && (py - llb_1 < -DH_1 / 2.0 || py - uub_1 > DH_1 / 2.0)) return;
    if (DIM > 2 && (pz - llb_2 < -DH_2 / 2.0 || pz - uub_2 > DH_2 / 2.0)) return;

    double SoA_arr[3] = {SoA_0, SoA_1, SoA_2};
    const int ex[3]  = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    double val = 0.0;
    global_interp_ss_device(
        ex, pex, d_X_0, d_X_1, d_X_2,
        d_field, &val,
        px, py, pz,
        ordn, SoA_arr, Symmetry, sst
    );

    atomicAdd(&d_shellf[j * num_var + var_idx], val);
    if (var_idx == 0) atomicAdd(&d_weight[j], 1);
}

void gpu_batch_interp_ss_launch(
    cudaStream_t stream,
    int NN, int DIM, int sst,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    double* d_field,
    double llb_0, double llb_1, double llb_2,
    double uub_0, double uub_1, double uub_2,
    double DH_0, double DH_1, double DH_2,
    int ordn, double SoA_0, double SoA_1, double SoA_2,
    int Symmetry, int var_idx, int num_var,
    double* d_shellf, int* d_weight,
    int* d_point_sst
) {
    int blockSize = 256;
    int gridSize = (NN + blockSize - 1) / blockSize;

    batch_interp_ss_kernel<<<gridSize, blockSize, 0, stream>>>(
        NN, DIM, sst,
        d_XX_0, d_XX_1, d_XX_2,
        shape_0, shape_1, shape_2,
        pshape_0, pshape_1, pshape_2,
        d_X_0, d_X_1, d_X_2, d_field,
        llb_0, llb_1, llb_2, uub_0, uub_1, uub_2,
        DH_0, DH_1, DH_2,
        ordn, SoA_0, SoA_1, SoA_2, Symmetry,
        var_idx, num_var, d_shellf, d_weight, d_point_sst
    );
}

// ============================================================================
// batch_interp_ss_2d_kernel + gpu_batch_interp_ss_2d_launch (NullShellPatch Interp_Points_2D_GPU)
// N-point batch 2D version for NullShellPatch
// ============================================================================

__global__ void batch_interp_ss_2d_kernel(
    int NN, int sst, int indZ,
    double* d_XX_0, double* d_XX_1,
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    double* d_X_0, double* d_X_1,
    double* d_field,
    double llb_0, double llb_1,
    double uub_0, double uub_1,
    int ordn, double SoA_0, double SoA_1, double SoA_2, int Symmetry,
    int var_idx, int num_var,
    double* d_shellf, int* d_weight
) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= NN) return;

    double px = d_XX_0[j], py = d_XX_1[j];
    if (px < llb_0 || px > uub_0) return;
    if (py < llb_1 || py > uub_1) return;

    double SoA_arr[3] = {SoA_0, SoA_1, SoA_2};
    const int ex[3]  = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    double val = 0.0;
    global_interp_ss_2d_device(
        ex, pex, d_X_0, d_X_1,
        d_field, &val,
        px, py,
        ordn, SoA_arr, Symmetry, sst, indZ
    );

    atomicAdd(&d_shellf[j * num_var + var_idx], val);
    if (var_idx == 0) atomicAdd(&d_weight[j], 1);
}

void gpu_batch_interp_ss_2d_launch(
    cudaStream_t stream,
    int NN, int sst, int indZ,
    double* d_XX_0, double* d_XX_1,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    double* d_X_0, double* d_X_1,
    double* d_field,
    double llb_0, double llb_1,
    double uub_0, double uub_1,
    int ordn, double SoA_0, double SoA_1, double SoA_2,
    int Symmetry, int var_idx, int num_var,
    double* d_shellf, int* d_weight
) {
    int blockSize = 256;
    int gridSize = (NN + blockSize - 1) / blockSize;

    batch_interp_ss_2d_kernel<<<gridSize, blockSize, 0, stream>>>(
        NN, sst, indZ,
        d_XX_0, d_XX_1,
        shape_0, shape_1, shape_2,
        pshape_0, pshape_1, pshape_2,
        d_X_0, d_X_1, d_field,
        llb_0, llb_1, uub_0, uub_1,
        ordn, SoA_0, SoA_1, SoA_2, Symmetry,
        var_idx, num_var, d_shellf, d_weight
    );
}

__global__ void global_interp_scalar_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* X, const double* Y, const double* Z,
    const double* __restrict__ f, double* __restrict__ f_int,
    double x1, double y1, double z1,
    int ordn, double SoA0, double SoA1, double SoA2, int symmetry
) {
    if (threadIdx.x || blockIdx.x) return;
    if (ordn <= 0 || ordn > MAX_ORDN) {
        *f_int = NAN;
        gpu_stop();
        return;
    }
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const double SoA[3] = {SoA0, SoA1, SoA2};
    global_interp_device(ex, pex, X, Y, Z, f, f_int, x1, y1, z1, ordn, SoA, symmetry);
}

void gpu_global_interp_scalar_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_f, double* d_f_int,
    double x1, double y1, double z1,
    int ordn, const double SoA[3], int symmetry
) {
    global_interp_scalar_kernel<<<1, 1, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_X, d_Y, d_Z, d_f, d_f_int,
        x1, y1, z1, ordn, SoA[0], SoA[1], SoA[2], symmetry
    );
}

__global__ void global_interp_ss_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* X, const double* Y, const double* Z,
    const double* __restrict__ f, double* __restrict__ f_int,
    double x1, double y1, double z1,
    int ordn, double SoA0, double SoA1, double SoA2, int sst
) {
    if (threadIdx.x || blockIdx.x) return;
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const double SoA[3] = {SoA0, SoA1, SoA2};
    // 直接调用已完美修复的双模 Device 函数
    global_interp_ss_device(ex, pex, X, Y, Z, f, f_int, x1, y1, z1, ordn, SoA, 0, sst);
}

void gpu_global_interp_ss_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_f, double* d_f_int,
    double x1, double y1, double z1,
    int ordn, const double SoA[3], int sst
) {
    global_interp_ss_kernel<<<1, 1, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_X, d_Y, d_Z, d_f, d_f_int,
        x1, y1, z1, ordn, SoA[0], SoA[1], SoA[2], sst
    );
}

__global__ void global_interp_ss_2d_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* X, const double* Y,
    const double* __restrict__ f, double* __restrict__ f_int,
    double x1, double y1, int indZ,
    int ordn, double SoA0, double SoA1, double SoA2, int sst
) {
    if (threadIdx.x || blockIdx.x) return;
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const double SoA[3] = {SoA0, SoA1, SoA2};
    // 直接调用已完美修复的双模 Device 函数
    global_interp_ss_2d_device(ex, pex, X, Y, f, f_int, x1, y1, ordn, SoA, 0, sst, indZ);
}

void gpu_global_interp_ss_2d_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_X, const double* d_Y,
    const double* d_f, double* d_f_int,
    double x1, double y1, int indZ,
    int ordn, const double SoA[3], int sst
) {
    global_interp_ss_2d_kernel<<<1, 1, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_X, d_Y, d_f, d_f_int,
        x1, y1, indZ, ordn, SoA[0], SoA[1], SoA[2], sst
    );
}

__global__ void l2normhelper_kernel(
    const double* __restrict__ f,
    int imin, int imax,
    int jmin, int jmax,
    int kmin, int kmax,
    int pnx, int pny,
    double dv,
    double* __restrict__ d_out
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + imin;
    int j = blockIdx.y * blockDim.y + threadIdx.y + jmin;
    int k = blockIdx.z * blockDim.z + threadIdx.z + kmin;

    double my_val = 0.0;

    if (i <= imax && j <= jmax && k <= kmax) {
        long long idx = (long long)i + (long long)j * pnx + (long long)k * pnx * pny;
        double val = f[idx];
        my_val = val * val * dv;
    }

    my_val = warpReduceSum(my_val);

    int tid = threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z * blockDim.x * blockDim.y;

    if ((tid % warpSize) == 0) {
        atomicAdd(d_out, my_val);
    }
}

void gpu_l2normhelper_launch(
    cudaStream_t stream,
    const int* ex, const int* pex,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    const double* d_f, double* d_sum, int gw
) {
    double dX = X[1] - X[0];
    double dY = Y[1] - Y[0];
    double dZ = Z[1] - Z[0];
    double dv = dX * dY * dZ;

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

    int nx_proc = imax - imin + 1;
    int ny_proc = jmax - jmin + 1;
    int nz_proc = kmax - kmin + 1;

    if (nx_proc <= 0 || ny_proc <= 0 || nz_proc <= 0) {
        return;
    }

    dim3 blockDim(8, 8, 8);
    dim3 gridDim((nx_proc + blockDim.x - 1) / blockDim.x,
                 (ny_proc + blockDim.y - 1) / blockDim.y,
                 (nz_proc + blockDim.z - 1) / blockDim.z);

    l2normhelper_kernel<<<gridDim, blockDim, 0, stream>>>(
        d_f,
        imin, imax, jmin, jmax, kmin, kmax,
        pex[0], pex[1],
        dv, d_sum
    );
}

void gpu_l2normhelper_sh_launch(
    cudaStream_t stream,
    const int* ex, const int* pex,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    const double* d_f, double* d_sum, int gw, int ogw, int Symmetry
) {
    double dX = X[1] - X[0];
    double dY = Y[1] - Y[0];
    double dZ = Z[1] - Z[0];
    double dv = dX * dY * dZ;
    double PIo4 = acos(-1.0) / 4.0;

    int imin = gw;
    int jmin = gw;
    int kmin = gw;
    int imax = ex[0] - gw - 1;
    int jmax = ex[1] - gw - 1;
    int kmax = ex[2] - gw - 1;

    if (fabs(X[ex[0] - 1] - xmax) < dX) imax = (X[ex[0] - 1] - PIo4 > dX) ? ex[0] - ogw - 1 : ex[0] - 1;
    if (fabs(Y[ex[1] - 1] - ymax) < dY) jmax = (Y[ex[1] - 1] - PIo4 > dY) ? ex[1] - ogw - 1 : ex[1] - 1;
    if (fabs(Z[ex[2] - 1] - zmax) < dZ) kmax = ex[2] - 1;

    if (fabs(X[0] - xmin) < dX) imin = (X[0] + PIo4 < dX) ? ogw : 0;
    if (fabs(Y[0] - ymin) < dY) jmin = (Y[0] + PIo4 < dY) ? ogw : 0;
    if (fabs(Z[0] - zmin) < dZ) kmin = 0;

    if (Symmetry == 1) {
        if (fabs(ymin + gw * dY) < dY && Y[0] < 0.0) jmin = gw;
        if (fabs(ymax - gw * dY) < dY && Y[ex[1] - 1] > 0.0) jmax = ex[1] - gw - 1;
    }
    if (Symmetry == 2) {
        if (fabs(xmin + gw * dX) < dX && X[0] < 0.0) imin = gw;
        if (fabs(ymin + gw * dY) < dY && Y[0] < 0.0) jmin = gw;
    }

    int nx_proc = imax - imin + 1;
    int ny_proc = jmax - jmin + 1;
    int nz_proc = kmax - kmin + 1;
    if (nx_proc <= 0 || ny_proc <= 0 || nz_proc <= 0) return;

    dim3 blockDim(8, 8, 8);
    dim3 gridDim((nx_proc + blockDim.x - 1) / blockDim.x,
                 (ny_proc + blockDim.y - 1) / blockDim.y,
                 (nz_proc + blockDim.z - 1) / blockDim.z);
    l2normhelper_kernel<<<gridDim, blockDim, 0, stream>>>(
        d_f, imin, imax, jmin, jmax, kmin, kmax, pex[0], pex[1], dv, d_sum
    );
}

__global__ void l2normhelper_rms_kernel(
    const double* __restrict__ f,
    int imin, int imax, int jmin, int jmax, int kmin, int kmax,
    int pnx, int pny,
    double* __restrict__ d_out, int* __restrict__ d_count
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + imin;
    int j = blockIdx.y * blockDim.y + threadIdx.y + jmin;
    int k = blockIdx.z * blockDim.z + threadIdx.z + kmin;

    double my_val = 0.0;
    int my_count = 0;
    if (i <= imax && j <= jmax && k <= kmax) {
        long long idx = (long long)i + (long long)j * pnx + (long long)k * pnx * pny;
        double val = f[idx];
        my_val = val * val;
        my_count = 1;
    }
    my_val = warpReduceSum(my_val);
    my_count = warpReduceSumInt(my_count);
    int tid = threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z * blockDim.x * blockDim.y;
    if ((tid % warpSize) == 0) {
        atomicAdd(d_out, my_val);
        atomicAdd(d_count, my_count);
    }
}

void gpu_l2normhelper_sh_rms_launch(
    cudaStream_t stream,
    const int* ex, const int* pex,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    const double* d_f, double* d_sum, int gw, int ogw, int Symmetry, int* d_count
) {
    double dX = X[1] - X[0];
    double dY = Y[1] - Y[0];
    double dZ = Z[1] - Z[0];
    double PIo4 = acos(-1.0) / 4.0;

    int imin = gw, jmin = gw, kmin = gw;
    int imax = ex[0] - gw - 1, jmax = ex[1] - gw - 1, kmax = ex[2] - gw - 1;
    if (fabs(X[ex[0] - 1] - xmax) < dX) imax = (X[ex[0] - 1] - PIo4 > dX) ? ex[0] - ogw - 1 : ex[0] - 1;
    if (fabs(Y[ex[1] - 1] - ymax) < dY) jmax = (Y[ex[1] - 1] - PIo4 > dY) ? ex[1] - ogw - 1 : ex[1] - 1;
    if (fabs(Z[ex[2] - 1] - zmax) < dZ) kmax = ex[2] - 1;
    if (fabs(X[0] - xmin) < dX) imin = (X[0] + PIo4 < dX) ? ogw : 0;
    if (fabs(Y[0] - ymin) < dY) jmin = (Y[0] + PIo4 < dY) ? ogw : 0;
    if (fabs(Z[0] - zmin) < dZ) kmin = 0;
    if (Symmetry == 1) {
        if (fabs(ymin + gw * dY) < dY && Y[0] < 0.0) jmin = gw;
        if (fabs(ymax - gw * dY) < dY && Y[ex[1] - 1] > 0.0) jmax = ex[1] - gw - 1;
    }
    if (Symmetry == 2) {
        if (fabs(xmin + gw * dX) < dX && X[0] < 0.0) imin = gw;
        if (fabs(ymin + gw * dY) < dY && Y[0] < 0.0) jmin = gw;
    }

    int nx_proc = imax - imin + 1;
    int ny_proc = jmax - jmin + 1;
    int nz_proc = kmax - kmin + 1;
    if (nx_proc <= 0 || ny_proc <= 0 || nz_proc <= 0) return;
    dim3 blockDim(8, 8, 8);
    dim3 gridDim((nx_proc + blockDim.x - 1) / blockDim.x,
                 (ny_proc + blockDim.y - 1) / blockDim.y,
                 (nz_proc + blockDim.z - 1) / blockDim.z);
    l2normhelper_rms_kernel<<<gridDim, blockDim, 0, stream>>>(
        d_f, imin, imax, jmin, jmax, kmin, kmax, pex[0], pex[1], d_sum, d_count
    );
}

__device__ __forceinline__ unsigned long long atomicMaxDoublePositive(double* address, double val) {
    unsigned long long* ull = reinterpret_cast<unsigned long long*>(address);
    return atomicMax(ull, __double_as_longlong(val));
}

__global__ void find_maximum_init_kernel(double* max_val, int* max_idx) {
    *max_val = 0.0;
    max_idx[0] = -1;
    max_idx[1] = -1;
    max_idx[2] = -1;
}

__global__ void find_maximum_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* fun,
    int llb0, int llb1, int llb2,
    int uub0, int uub1, int uub2,
    double* max_val
) {
    int imin = llb0;
    int jmin = llb1;
    int kmin = llb2;
    int imax = ex0 - uub0 - 1;
    int jmax = ex1 - uub1 - 1;
    int kmax = ex2 - uub2 - 1;
    int nx = imax - imin + 1;
    int ny = jmax - jmin + 1;
    int nz = kmax - kmin + 1;
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (nx <= 0 || ny <= 0 || nz <= 0 || tid >= nx * ny * nz) return;

    int ii, jj, kk;
    linear_to_ijk(tid, nx, ny, ii, jj, kk);
    int i = imin + ii;
    int j = jmin + jj;
    int k = kmin + kk;
    int idx = (k * pex1 + j) * pex0 + i;
    double v = fabs(fun[idx]);
    if (!isfinite(v)) return;
    atomicMaxDoublePositive(max_val, v);
}

__global__ void find_maximum_idx_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* fun,
    int llb0, int llb1, int llb2,
    int uub0, int uub1, int uub2,
    const double* max_val, int* max_idx
) {
    int imin = llb0;
    int jmin = llb1;
    int kmin = llb2;
    int imax = ex0 - uub0 - 1;
    int jmax = ex1 - uub1 - 1;
    int kmax = ex2 - uub2 - 1;
    int nx = imax - imin + 1;
    int ny = jmax - jmin + 1;
    int nz = kmax - kmin + 1;
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (nx <= 0 || ny <= 0 || nz <= 0 || tid >= nx * ny * nz) return;

    int ii, jj, kk;
    linear_to_ijk(tid, nx, ny, ii, jj, kk);
    int i = imin + ii;
    int j = jmin + jj;
    int k = kmin + kk;
    int idx = (k * pex1 + j) * pex0 + i;
    double v = fabs(fun[idx]);
    if (!isfinite(v) || v != *max_val) return;
    if (atomicCAS(&max_idx[0], -1, i) == -1) {
        max_idx[1] = j;
        max_idx[2] = k;
    }
}

__global__ void find_maximum_pos_kernel(
    const double* X, const double* Y, const double* Z,
    const int* max_idx, double* pos
) {
    if (max_idx[0] < 0) {
        pos[0] = NAN;
        pos[1] = NAN;
        pos[2] = NAN;
        return;
    }
    pos[0] = X[max_idx[0]];
    pos[1] = Y[max_idx[1]];
    pos[2] = Z[max_idx[2]];
}

void gpu_find_maximum_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_fun,
    double* d_val, double* d_pos,
    const int llb[3], const int uub[3],
    int* d_tmp_idx
) {
    int nx = ex[0] - llb[0] - uub[0];
    int ny = ex[1] - llb[1] - uub[1];
    int nz = ex[2] - llb[2] - uub[2];
    int total = nx * ny * nz;
    if (total <= 0) return;
    dim3 block(256), grid((total + 255) / 256);
    find_maximum_init_kernel<<<1, 1, 0, stream>>>(d_val, d_tmp_idx);
    find_maximum_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1],
        d_fun,
        llb[0], llb[1], llb[2], uub[0], uub[1], uub[2],
        d_val
    );
    find_maximum_idx_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1],
        d_fun,
        llb[0], llb[1], llb[2], uub[0], uub[1], uub[2],
        d_val, d_tmp_idx
    );
    find_maximum_pos_kernel<<<1, 1, 0, stream>>>(d_X, d_Y, d_Z, d_tmp_idx, d_pos);
}

__global__ void wigner_d_kernel(int l, int m, int s, double costheta, double* out) {
    if (threadIdx.x || blockIdx.x) return;
    *out = wigner_d_device(l, m, s, costheta);
}

void gpu_wigner_d_launch(cudaStream_t stream, int l, int m, int s, double costheta, double* d_out) {
    wigner_d_kernel<<<1, 1, 0, stream>>>(l, m, s, costheta, d_out);
}

__global__ void ffact_kernel(int n, double* out) {
    if (threadIdx.x || blockIdx.x) return;
    *out = ffact_device(n);
}

void gpu_ffact_launch(cudaStream_t stream, int n, double* d_out) {
    ffact_kernel<<<1, 1, 0, stream>>>(n, d_out);
}

__global__ void yslm_kernel(int s, int l, int m, double theta, double phi, double* real_out, double* imag_out) {
    if (threadIdx.x || blockIdx.x) return;
    yslm_device(s, l, m, theta, phi, real_out, imag_out);
}

void gpu_yslm_launch(
    cudaStream_t stream, int s, int l, int m, double theta, double phi,
    double* d_real_out, double* d_imag_out
) {
    yslm_kernel<<<1, 1, 0, stream>>>(s, l, m, theta, phi, d_real_out, d_imag_out);
}

__global__ void global_interpind_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ f, double* __restrict__ f_int,
    int ordn, double SoA0, double SoA1, double SoA2,
    const int* __restrict__ inds, const double* __restrict__ coef, int sst
) {
    if (threadIdx.x || blockIdx.x) return;
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const int indv[3] = {inds[0], inds[1], inds[2]};
    const double SoA[3] = {SoA0, SoA1, SoA2};
    *f_int = interp3_indexed_device(ex, pex, f, indv, coef, ordn, SoA, sst);
}

void gpu_global_interpind_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_f, double* d_f_int,
    int ordn, const double SoA[3],
    const int* d_inds, const double* d_coef, int sst
) {
    global_interpind_kernel<<<1, 1, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_f, d_f_int, ordn, SoA[0], SoA[1], SoA[2], d_inds, d_coef, sst
    );
}

__global__ void global_interpind2d_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ f, double* __restrict__ f_int,
    int ordn, double SoA0, double SoA1, double SoA2,
    const int* __restrict__ inds, const double* __restrict__ coef, int sst
) {
    if (threadIdx.x || blockIdx.x) return;
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const int indv[3] = {inds[0], inds[1], inds[2]};
    const double SoA[3] = {SoA0, SoA1, SoA2};
    *f_int = interp2_indexed_device(ex, pex, f, indv, coef, ordn, SoA, sst);
}

void gpu_global_interpind2d_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_f, double* d_f_int,
    int ordn, const double SoA[3],
    const int* d_inds, const double* d_coef, int sst
) {
    global_interpind2d_kernel<<<1, 1, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_f, d_f_int, ordn, SoA[0], SoA[1], SoA[2], d_inds, d_coef, sst
    );
}

__global__ void global_interpind1d_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    const double* __restrict__ f, double* __restrict__ f_int,
    int ordn, double SoA0, double SoA1, double SoA2,
    const int* __restrict__ inds, const double* __restrict__ coef, int sst, int dumyd
) {
    if (threadIdx.x || blockIdx.x) return;
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const int indv[3] = {inds[0], inds[1], inds[2]};
    const double SoA[3] = {SoA0, SoA1, SoA2};
    *f_int = interp1_indexed_device(ex, pex, f, indv, coef, ordn, SoA, sst, dumyd);
}

void gpu_global_interpind1d_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_f, double* d_f_int,
    int ordn, const double SoA[3],
    const int* d_inds, const double* d_coef, int sst, int dumyd
) {
    global_interpind1d_kernel<<<1, 1, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_f, d_f_int, ordn, SoA[0], SoA[1], SoA[2], d_inds, d_coef, sst, dumyd
    );
}

__global__ void global_interp_amr_kernel(
    int active_count, int DIM,
    int* d_active_indices,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    double* d_field,
    int ordn, double SoA_0, double SoA_1, double SoA_2, int Symmetry,
    int var_idx, int num_var,
    double* d_shellf
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= active_count) return;

    int j = d_active_indices[idx];
    double px = d_XX_0[j], py = d_XX_1[j], pz = d_XX_2[j];

    double SoA_arr[3] = {SoA_0, SoA_1, SoA_2};
    const int ex[3]  = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    double val = 0.0;
    global_interp_device(
        ex, pex, d_X_0, d_X_1, d_X_2,
        d_field, &val,
        px, py, pz,
        ordn, SoA_arr, Symmetry
    );

    d_shellf[j * num_var + var_idx] = val;
}

void gpu_global_interp_amr_launch(
    cudaStream_t stream,
    int active_count, int DIM,
    int* d_active_indices,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    double* d_field,
    int ordn, double SoA_0, double SoA_1, double SoA_2,
    int Symmetry, int var_idx, int num_var,
    double* d_shellf
) {
    if (active_count == 0) return;
    int blockSize = 256;
    int gridSize = (active_count + blockSize - 1) / blockSize;

    global_interp_amr_kernel<<<gridSize, blockSize, 0, stream>>>(
        active_count, DIM, d_active_indices,
        d_XX_0, d_XX_1, d_XX_2,
        shape_0, shape_1, shape_2,
        pshape_0, pshape_1, pshape_2,
        d_X_0, d_X_1, d_X_2, d_field,
        ordn, SoA_0, SoA_1, SoA_2, Symmetry,
        var_idx, num_var, d_shellf
    );
}

__global__ void apply_symmetry_bd_kernel(
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    double* __restrict__ f,
    double SYM0, double SYM1, double SYM2
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= pextc0 * pextc1 * pextc2) return;
    int i_pad, j_pad, k_pad;
    linear_to_ijk(tid, pextc0, pextc1, i_pad, j_pad, k_pad);

    int i = i_pad - GPU_PAD_X_LEFT;
    int j = j_pad - ghost_width;
    int k = k_pad - ghost_width;

    if (i < -ghost_width) return;
    if (i >= 0 && j >= 0 && k >= 0) return;
    if (i >= extc0 || j >= extc1 || k >= extc2) return; // 上/下 ghost 混合角点，跳过

    int out_idx = (k * pextc1 + j) * pextc0 + i;

    bool i_neg = (i < 0);
    bool j_neg = (j < 0);
    bool k_neg = (k < 0);
#ifdef Vertex
    // Vertex: symmetry plane is ON a grid point (i=0 at x=0).
    // Ghost point i=-1 reflects to i=1, i=-2 to i=2, etc.
    int ri = i_neg ? (-i) : i;
    int rj = j_neg ? (-j) : j;
    int rk = k_neg ? (-k) : k;
#else
    // Cell: symmetry plane is BETWEEN grid points.
    // Ghost point i=-1 reflects to i=0, i=-2 to i=1, etc.
    int ri = i_neg ? (-i - 1) : i;
    int rj = j_neg ? (-j - 1) : j;
    int rk = k_neg ? (-k - 1) : k;
#endif
    int src_idx = (rk * pextc1 + rj) * pextc0 + ri;

    double factor = (i_neg ? SYM0 : 1.0)
                  * (j_neg ? SYM1 : 1.0)
                  * (k_neg ? SYM2 : 1.0);
    f[out_idx] = f[src_idx] * factor;
}

// sym_type: 0 = none, 1 = equatorial (z only), 2 = octant (x,y,z)
__global__ void apply_symmetry_bd_fused_kernel(
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    int sym_type,
    double* __restrict__ chi, double* __restrict__ trK,
    double* __restrict__ dxx, double* __restrict__ gxy, double* __restrict__ gxz,
    double* __restrict__ dyy, double* __restrict__ gyz, double* __restrict__ dzz,
    double* __restrict__ Axx, double* __restrict__ Axy, double* __restrict__ Axz,
    double* __restrict__ Ayy, double* __restrict__ Ayz, double* __restrict__ Azz,
    double* __restrict__ Gamx, double* __restrict__ Gamy, double* __restrict__ Gamz,
    double* __restrict__ Lap,
    double* __restrict__ betax, double* __restrict__ betay, double* __restrict__ betaz,
    double* __restrict__ dtSfx, double* __restrict__ dtSfy, double* __restrict__ dtSfz
) {
    if (sym_type == 0) return;
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= pextc0 * pextc1 * pextc2) return;
    int i_pad, j_pad, k_pad;
    linear_to_ijk(tid, pextc0, pextc1, i_pad, j_pad, k_pad);

    int i = i_pad - GPU_PAD_X_LEFT;
    int j = j_pad - ghost_width;
    int k = k_pad - ghost_width;

    if (i < -ghost_width) return;
    if (i >= 0 && j >= 0 && k >= 0) return;
    if (i >= extc0 || j >= extc1 || k >= extc2) return;

    // equatorial: only k < 0 ghost zones need symmetry; skip i/j ghost zones
    if (sym_type == 1 && k >= 0) return;

    int out_idx = (k * pextc1 + j) * pextc0 + i;

    bool i_neg = (i < 0);
    bool j_neg = (j < 0);
    bool k_neg = (k < 0);
#ifdef Vertex
    int ri = i_neg ? (-i) : i;
    int rj = j_neg ? (-j) : j;
    int rk = k_neg ? (-k) : k;
#else
    int ri = i_neg ? (-i - 1) : i;
    int rj = j_neg ? (-j - 1) : j;
    int rk = k_neg ? (-k - 1) : k;
#endif
    int src_idx = (rk * pextc1 + rj) * pextc0 + ri;

    double fi = i_neg ? -1.0 : 1.0;
    double fj = j_neg ? -1.0 : 1.0;
    double fk = k_neg ? -1.0 : 1.0;

    double f_sss = 1.0;
    double f_aas = fi * fj;
    double f_asa = fi * fk;
    double f_saa = fj * fk;
    double f_ass = fi;
    double f_sas = fj;
    double f_ssa = fk;

    chi[out_idx]   = chi[src_idx]   * f_sss;
    trK[out_idx]   = trK[src_idx]   * f_sss;
    dxx[out_idx]   = dxx[src_idx]   * f_sss;
    dyy[out_idx]   = dyy[src_idx]   * f_sss;
    dzz[out_idx]   = dzz[src_idx]   * f_sss;

    gxy[out_idx]   = gxy[src_idx]   * f_aas;
    gxz[out_idx]   = gxz[src_idx]   * f_asa;
    gyz[out_idx]   = gyz[src_idx]   * f_saa;

    Axx[out_idx]   = Axx[src_idx]   * f_sss;
    Ayy[out_idx]   = Ayy[src_idx]   * f_sss;
    Azz[out_idx]   = Azz[src_idx]   * f_sss;
    Axy[out_idx]   = Axy[src_idx]   * f_aas;
    Axz[out_idx]   = Axz[src_idx]   * f_asa;
    Ayz[out_idx]   = Ayz[src_idx]   * f_saa;

    Gamx[out_idx]  = Gamx[src_idx]  * f_ass;
    Gamy[out_idx]  = Gamy[src_idx]  * f_sas;
    Gamz[out_idx]  = Gamz[src_idx]  * f_ssa;

    Lap[out_idx]   = Lap[src_idx]   * f_sss;

    betax[out_idx] = betax[src_idx] * f_ass;
    betay[out_idx] = betay[src_idx] * f_sas;
    betaz[out_idx] = betaz[src_idx] * f_ssa;

    dtSfx[out_idx] = dtSfx[src_idx] * f_ass;
    dtSfy[out_idx] = dtSfy[src_idx] * f_sas;
    dtSfz[out_idx] = dtSfz[src_idx] * f_ssa;
}

void gpu_apply_symmetry_bd(
    cudaStream_t stream,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    double* __restrict__ func,
    double SoA0, double SoA1, double SoA2
) {
    int total = pextc0 * pextc1 * pextc2;
    dim3 block(256);
    dim3 grid((total + 255) / 256);

    apply_symmetry_bd_kernel<<<grid, block, 0, stream>>>(
        extc0, extc1, extc2, pextc0, pextc1, pextc2,
        func, SoA0, SoA1, SoA2
    );
}

void gpu_apply_symmetry_bd_fused_launch(
    cudaStream_t stream,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    int sym_type,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz
) {
    if (sym_type == 0) return;
    int total = pextc0 * pextc1 * pextc2;
    dim3 block(256);
    dim3 grid((total + 255) / 256);

    apply_symmetry_bd_fused_kernel<<<grid, block, 0, stream>>>(
        extc0, extc1, extc2, pextc0, pextc1, pextc2, sym_type,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap,
        betax, betay, betaz, dtSfx, dtSfy, dtSfz
    );
}

void gpu_apply_scalar_symmetry_launch(
    cudaStream_t stream,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    double* Sphi, double* Spi
) {
    // Sphi and Spi are scalars with SoA = (1, 1, 1): even parity under all reflections.
    // Reuse the generic single-variable symmetry kernel with SYM=1.0 on all axes.
    apply_symmetry_bd_kernel<<<(pextc0*pextc1*pextc2+255)/256, 256, 0, stream>>>(
        extc0, extc1, extc2, pextc0, pextc1, pextc2,
        Sphi, 1.0, 1.0, 1.0
    );
    apply_symmetry_bd_kernel<<<(pextc0*pextc1*pextc2+255)/256, 256, 0, stream>>>(
        extc0, extc1, extc2, pextc0, pextc1, pextc2,
        Spi, 1.0, 1.0, 1.0
    );
}

// ==========================================
// Fused Pack Kernel (Corrected Padding)
// ==========================================
__global__ void gpu_fused_pack_kernel(
    FusedVars vars,
    double* __restrict__ dst_1d,
    int src_pnx, int src_pny,  // 必须传入物理维度的跨度
    int dst_nx, int dst_ny, int dst_nz,
    int off_x, int off_y, int off_z
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_per_var = dst_nx * dst_ny * dst_nz;
    if (idx >= total_per_var) return;

    int i = idx % dst_nx;
    int j = (idx / dst_nx) % dst_ny;
    int k = idx / (dst_nx * dst_ny);

    // 使用 pnx 和 pny 计算源地址偏移
    int src_idx = (i + off_x) + src_pnx * ((j + off_y) + src_pny * (k + off_z));

    for (int v = 0; v < vars.num_vars; ++v) {
        dst_1d[v * total_per_var + idx] = vars.ptrs[v][src_idx];
    }
}

void gpu_fused_pack_launch(
    cudaStream_t stream,
    FusedVars vars, double* d_dst_1d,
    int src_pnx, int src_pny, // 修复 Launch 签名
    int dst_nx, int dst_ny, int dst_nz,
    int off_x, int off_y, int off_z
) {
    if (vars.num_vars == 0) return;
    int n = dst_nx * dst_ny * dst_nz;
    int block = 256;
    int grid = (n + block - 1) / block;
    gpu_fused_pack_kernel<<<grid, block, 0, stream>>>(
        vars, d_dst_1d,
        src_pnx, src_pny, dst_nx, dst_ny, dst_nz,
        off_x, off_y, off_z
    );
}

// ==========================================
// Fused Unpack Kernel (Corrected Padding)
// ==========================================
__global__ void gpu_fused_unpack_kernel(
    const double* __restrict__ src_1d,
    FusedVars vars,
    int dst_pnx, int dst_pny, // 必须传入物理维度的跨度
    int src_nx, int src_ny, int src_nz,
    int off_x, int off_y, int off_z
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_per_var = src_nx * src_ny * src_nz;
    if (idx >= total_per_var) return;

    int i = idx % src_nx;
    int j = (idx / src_nx) % src_ny;
    int k = idx / (src_nx * src_ny);

    // 使用 pnx 和 pny 计算目标地址偏移
    int dst_idx = (i + off_x) + dst_pnx * ((j + off_y) + dst_pny * (k + off_z));

    for (int v = 0; v < vars.num_vars; ++v) {
        vars.ptrs[v][dst_idx] = src_1d[v * total_per_var + idx];
    }
}

void gpu_fused_unpack_launch(
    cudaStream_t stream,
    const double* d_src_1d, FusedVars vars,
    int dst_pnx, int dst_pny, // 修复 Launch 签名
    int src_nx, int src_ny, int src_nz,
    int off_x, int off_y, int off_z
) {
    if (vars.num_vars == 0) return;
    int n = src_nx * src_ny * src_nz;
    int block = 256;
    int grid = (n + block - 1) / block;
    gpu_fused_unpack_kernel<<<grid, block, 0, stream>>>(
        d_src_1d, vars,
        dst_pnx, dst_pny, src_nx, src_ny, src_nz,
        off_x, off_y, off_z
    );
}

// ==========================================
// Fused Time Level Interpolation Kernels (from main_merge)
// ==========================================

__global__ void gpu_fused_time_interp_3_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, FusedTimeVars3 vars, int tindex
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = i + pex0 * (j + pex1 * k);
    for (int v = 0; v < vars.num_vars; ++v) {
        double f1 = vars.ptrs1[v][idx];
        double f2 = vars.ptrs2[v][idx];
        if (tindex == 0) {
            vars.ptrs3[v][idx] = 0.5 * (f1 + f2);
        } else if (tindex == 1) {
            vars.ptrs3[v][idx] = 0.75 * f1 + 0.25 * f2;
        } else if (tindex == -1) {
            vars.ptrs3[v][idx] = 0.75 * f2 + 0.25 * f1;
        }
    }
}

__global__ void gpu_fused_time_interp_4_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, FusedTimeVars4 vars, int tindex
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0 * ex1 * ex2) return;
    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    int idx = i + pex0 * (j + pex1 * k);
    for (int v = 0; v < vars.num_vars; ++v) {
        double f1 = vars.ptrs1[v][idx];
        double f2 = vars.ptrs2[v][idx];
        double f3 = vars.ptrs3[v][idx];
        if (tindex == 0) {
            vars.ptrs4[v][idx] = (3.0 / 8.0) * f1 + (3.0 / 4.0) * f2 - (1.0 / 8.0) * f3;
        } else if (tindex == 1) {
            vars.ptrs4[v][idx] = (5.0 / 32.0) * f1 + (15.0 / 16.0) * f2 - (3.0 / 32.0) * f3;
        } else if (tindex == -1) {
            vars.ptrs4[v][idx] = (5.0 / 32.0) * f1 + (15.0 / 16.0) * f2 - (3.0 / 32.0) * f3;
        }
    }
}

void gpu_fused_time_interp_3_launch(cudaStream_t stream, const int ext[3], const int pext[3], FusedTimeVars3 vars, int tindex) {
    if (vars.num_vars == 0) return;
    int total = ext[0] * ext[1] * ext[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);
    gpu_fused_time_interp_3_kernel<<<grid, block, 0, stream>>>(
        ext[0], ext[1], ext[2], pext[0], pext[1], vars, tindex
    );
}

void gpu_fused_time_interp_4_launch(cudaStream_t stream, const int ext[3], const int pext[3], FusedTimeVars4 vars, int tindex) {
    if (vars.num_vars == 0) return;
    int total = ext[0] * ext[1] * ext[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);
    gpu_fused_time_interp_4_kernel<<<grid, block, 0, stream>>>(
        ext[0], ext[1], ext[2], pext[0], pext[1], vars, tindex
    );
}

// =====================================================================
// Fused Global Interpolation Kernels (from main_merge)
// =====================================================================

__global__ void fused_global_interp_kernel(
    int NN, int DIM,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    FusedVars vars,
    double llb_0, double llb_1, double llb_2,
    double uub_0, double uub_1, double uub_2,
    double DH_0, double DH_1, double DH_2,
    int ordn, int Symmetry,
    int base_var_idx, int total_num_var,
    double* d_shellf, int* d_weight
) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= NN) return;

    double px = d_XX_0[j];
    double py = d_XX_1[j];
    double pz = d_XX_2[j];

    if (px - llb_0 < -DH_0 / 2.0 || px - uub_0 > DH_0 / 2.0) return;
    if (DIM > 1 && (py - llb_1 < -DH_1 / 2.0 || py - uub_1 > DH_1 / 2.0)) return;
    if (DIM > 2 && (pz - llb_2 < -DH_2 / 2.0 || pz - uub_2 > DH_2 / 2.0)) return;

    double* d_X_arr[3] = {d_X_0, d_X_1, d_X_2};
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    for (int v = 0; v < vars.num_vars; ++v) {
        double SoA_arr[3] = {vars.SoA0[v], vars.SoA1[v], vars.SoA2[v]};
        double val = 0.0;

        global_interp_device(
            ex, pex, d_X_arr[0], d_X_arr[1], d_X_arr[2],
            vars.ptrs[v], &val,
            px, py, pz,
            ordn, SoA_arr, Symmetry
        );

        atomicAdd(&d_shellf[j * total_num_var + (base_var_idx + v)], val);
    }

    if (base_var_idx == 0) {
        atomicAdd(&d_weight[j], 1);
    }
}

void gpu_fused_global_interp_launch(
    cudaStream_t stream,
    int NN, int DIM,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    FusedVars vars,
    double llb_0, double llb_1, double llb_2,
    double uub_0, double uub_1, double uub_2,
    double DH_0, double DH_1, double DH_2,
    int ordn, int Symmetry,
    int base_var_idx, int total_num_var,
    double* d_shellf, int* d_weight
) {
    if (vars.num_vars == 0) return;
    int blockSize = 256;
    int gridSize = (NN + blockSize - 1) / blockSize;

    fused_global_interp_kernel<<<gridSize, blockSize, 0, stream>>>(
        NN, DIM,
        d_XX_0, d_XX_1, d_XX_2,
        shape_0, shape_1, shape_2,
        pshape_0, pshape_1, pshape_2,
        d_X_0, d_X_1, d_X_2,
        vars,
        llb_0, llb_1, llb_2, uub_0, uub_1, uub_2,
        DH_0, DH_1, DH_2,
        ordn, Symmetry,
        base_var_idx, total_num_var, d_shellf, d_weight
    );
}

__global__ void fused_global_interp_amr_kernel(
    int active_count, int DIM,
    int* d_active_indices,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    FusedVars vars,
    int ordn, int Symmetry,
    int base_var_idx, int total_num_var,
    double* d_shellf
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= active_count) return;

    int j = d_active_indices[idx];

    double px = d_XX_0[j];
    double py = d_XX_1[j];
    double pz = d_XX_2[j];

    double* d_X_arr[3] = {d_X_0, d_X_1, d_X_2};
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};

    for (int v = 0; v < vars.num_vars; ++v) {
        double SoA_arr[3] = {vars.SoA0[v], vars.SoA1[v], vars.SoA2[v]};
        double val = 0.0;

        global_interp_device(
            ex, pex, d_X_arr[0], d_X_arr[1], d_X_arr[2],
            vars.ptrs[v], &val,
            px, py, pz,
            ordn, SoA_arr, Symmetry
        );

        d_shellf[j * total_num_var + (base_var_idx + v)] = val;
    }
}

void gpu_fused_global_interp_amr_launch(
    cudaStream_t stream,
    int active_count, int DIM,
    int* d_active_indices,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    FusedVars vars,
    int ordn, int Symmetry,
    int base_var_idx, int total_num_var,
    double* d_shellf
) {
    if (active_count == 0 || vars.num_vars == 0) return;
    int blockSize = 256;
    int gridSize = (active_count + blockSize - 1) / blockSize;

    fused_global_interp_amr_kernel<<<gridSize, blockSize, 0, stream>>>(
        active_count, DIM, d_active_indices,
        d_XX_0, d_XX_1, d_XX_2,
        shape_0, shape_1, shape_2,
        pshape_0, pshape_1, pshape_2,
        d_X_0, d_X_1, d_X_2,
        vars,
        ordn, Symmetry,
        base_var_idx, total_num_var, d_shellf
    );
}

// =====================================================================
// Fused Apply Symmetry BD Kernel (task-queue style, general FusedVars)
// =====================================================================

__global__ void apply_symmetry_bd_fused_vars_kernel(
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    int sym_type,
    FusedVars vars
) {
    if (sym_type == 0) return;
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= pextc0 * pextc1 * pextc2) return;
    int i_pad, j_pad, k_pad;
    linear_to_ijk(tid, pextc0, pextc1, i_pad, j_pad, k_pad);

    int i = i_pad - GPU_PAD_X_LEFT;
    int j = j_pad - ghost_width;
    int k = k_pad - ghost_width;

    if (i < -ghost_width) return;
    if (i >= 0 && j >= 0 && k >= 0) return;

    if (sym_type == 1 && k >= 0) return;

    int out_idx = (k * pextc1 + j) * pextc0 + i;

    if (i >= extc0 || j >= extc1 || k >= extc2) return;

    bool i_neg = (i < 0);
    bool j_neg = (j < 0);
    bool k_neg = (k < 0);
#ifdef Vertex
    // Vertex: symmetry plane is ON a grid point (i=0 at x=0).
    // Ghost point i=-1 reflects to i=1, i=-2 to i=2, etc.
    int ri = i_neg ? (-i) : i;
    int rj = j_neg ? (-j) : j;
    int rk = k_neg ? (-k) : k;
#else
    // Cell: symmetry plane is BETWEEN grid points.
    // Ghost point i=-1 reflects to i=0, i=-2 to i=1, etc.
    int ri = i_neg ? (-i - 1) : i;
    int rj = j_neg ? (-j - 1) : j;
    int rk = k_neg ? (-k - 1) : k;
#endif
    int src_idx = (rk * pextc1 + rj) * pextc0 + ri;

    for (int v = 0; v < vars.num_vars; ++v) {
        double factor = (i_neg ? vars.SoA0[v] : 1.0)
                      * (j_neg ? vars.SoA1[v] : 1.0)
                      * (k_neg ? vars.SoA2[v] : 1.0);
        vars.ptrs[v][out_idx] = vars.ptrs[v][src_idx] * factor;
    }
}

void gpu_apply_symmetry_bd_vars_launch(
    cudaStream_t stream,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    int sym_type,
    FusedVars vars
) {
    if (sym_type == 0 || vars.num_vars == 0) return;
    int total = pextc0 * pextc1 * pextc2;
    dim3 block(256);
    dim3 grid((total + 255) / 256);

    apply_symmetry_bd_fused_vars_kernel<<<grid, block, 0, stream>>>(
        extc0, extc1, extc2, pextc0, pextc1, pextc2, sym_type, vars
    );
}

