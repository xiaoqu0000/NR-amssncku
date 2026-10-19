#include "macrodef.h"
#include "bssn_rhs.h"
#include "rhs_smem_common.cuh"
#include <cstdio>

#include "fmisc_gpu.cuh"
#include "diff_new_gpu.cuh"
#include "kodiss_gpu.cuh"
#include "lopsidediff_gpu.cuh"
#include "gpu_manager.h"

#include <cuda_runtime.h>
#include <mpi.h>
#include <math.h>
#include <device_launch_parameters.h>
#include <iostream>
#include <cuda/pipeline>

#define IDX3D(i, j, k, nx, ny, nz) ((i) + (nx) * ((j) + (ny) * (k)))
#define IDX3D_PAD(i, j, k, pnx, pny) ((i) + (pnx) * ((j) + (pny) * (k)))

// BH data globals (defined in transpbh.C)
extern int BHN;
extern double Mass[3];
extern double PBH[9];


// ================= 预计算常量 =================
constexpr double ONE = 1.0;
constexpr double TWO = 2.0;
constexpr double FOUR_PI  = 4.0 * M_PI;
// ===============================================

constexpr double F1o3 = 1.0 / 3.0;
constexpr double F2o3 = 2.0 / 3.0;
constexpr double F3o2 = 1.5;
constexpr double HALF = 0.5;
constexpr double FF = 0.75;
constexpr double eta = 2.0;

template<int GW>
__launch_bounds__(256, 2)
__global__ void bssn_ricci_kernel(
    const int ex0, const int ex1, const int ex2, 
    const int pex0, const int pex1, const int pex2,
    const int symmetry, const int lev,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    double* __restrict__ Gamxxx_out, double* __restrict__ Gamxxy_out, double* __restrict__ Gamxxz_out,
    double* __restrict__ Gamxyy_out, double* __restrict__ Gamxyz_out, double* __restrict__ Gamxzz_out,
    double* __restrict__ Gamyxx_out, double* __restrict__ Gamyxy_out, double* __restrict__ Gamyxz_out,
    double* __restrict__ Gamyyy_out, double* __restrict__ Gamyyz_out, double* __restrict__ Gamyzz_out,
    double* __restrict__ Gamzxx_out, double* __restrict__ Gamzxy_out, double* __restrict__ Gamzxz_out,
    double* __restrict__ Gamzyy_out, double* __restrict__ Gamzyz_out, double* __restrict__ Gamzzz_out,
    double* __restrict__ Rxx_out, double* __restrict__ Rxy_out, double* __restrict__ Rxz_out,
    double* __restrict__ Ryy_out, double* __restrict__ Ryz_out, double* __restrict__ Rzz_out
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int pdims[3] = {pex0, pex1, pex2};
    double dX = X[1] - X[0]; double dY = Y[1] - Y[0]; double dZ = Z[1] - Z[0];

    int imax = ex0 - 1; int jmax = ex1 - 1; int kmax = ex2 - 1;
    int imin = 0, jmin = 0, kmin = 0;
    constexpr int ORD = GW - 1;
    if (symmetry > 0 && fabs(Z[0]) < dZ) kmin = -ORD;
    if (symmetry > 1 && fabs(X[0]) < dX) imin = -ORD;
    if (symmetry > 1 && fabs(Y[0]) < dY) jmin = -ORD;

    int block_i = blockIdx.x * blockDim.x;
    int block_j = blockIdx.y * blockDim.y;
    int block_k = blockIdx.z * blockDim.z;

    bool is_valid = (i < ex0 && j < ex1 && k < ex2);
    int idx = is_valid ? IDX3D_PAD(i, j, k, pex0, pex1) : 0;

    const double HALF = 0.5; const double TWO = 2.0; const double ONE = 1.0;

    // =========================================================================
    // 0. 前置计算逆度规
    // =========================================================================
    double l_gxx = 0, l_gxy = 0, l_gxz = 0, l_gyy = 0, l_gyz = 0, l_gzz = 0;
    double gupxx = 0, gupxy = 0, gupxz = 0, gupyy = 0, gupyz = 0, gupzz = 0;

    if (is_valid) {
        l_gxx = dxx[idx] + ONE; l_gxy = gxy[idx];       l_gxz = gxz[idx];
        l_gyy = dyy[idx] + ONE; l_gyz = gyz[idx];       l_gzz = dzz[idx] + ONE;

        double det = l_gxx * l_gyy * l_gzz + TWO * l_gxy * l_gyz * l_gxz - l_gxz * l_gyy * l_gxz - l_gxx * l_gyz * l_gyz - l_gzz * l_gxy * l_gxy;
        double inv_det = 1.0 / det;

        gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) * inv_det;
        gupxy = (l_gyz * l_gxz - l_gxy * l_gzz) * inv_det;
        gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) * inv_det;
        gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) * inv_det;
        gupyz = (l_gxy * l_gxz - l_gxx * l_gyz) * inv_det;
        gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) * inv_det;
    }

    // Dynamic shared memory: 3 tiles of [Z][Y][X] doubles.
    // For GW=5 this is 3*3072*8 = 73728 bytes, exceeding the 48KB static smem limit on sm_80.
    extern __shared__ double smem_raw[];
    double (*smem3)[TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X] =
        reinterpret_cast<double(*)[TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X]>(smem_raw);

    // =========================================================================
    // 1. 批量 Load，串行 Compute
    // =========================================================================
    
    // --- Batch 1: dxx, gxy, gxz ---
    load_3_fields_to_smem<GW>(smem3, dxx, gxy, gxz, pdims, block_i, block_j, block_k);
    __syncthreads(); // 仅需一次同步等待三个场加载完成

    double dgxx_x, dgxx_y, dgxx_z, contract_xx = 0;
    {
        double d2_xx, d2_xy, d2_xz, d2_yy, d2_yz, d2_zz;
        compute_all_derivs_smem<GW>(smem3[0], &dgxx_x, &dgxx_y, &dgxx_z, &d2_xx, &d2_xy, &d2_xz, &d2_yy, &d2_yz, &d2_zz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
        if (is_valid) contract_xx = gupxx * d2_xx + gupyy * d2_yy + gupzz * d2_zz + (gupxy * d2_xy + gupxz * d2_xz + gupyz * d2_yz) * TWO;
    } // d2 寄存器在此处安全死亡

    double dgxy_x, dgxy_y, dgxy_z, contract_xy = 0;
    {
        double d2_xx, d2_xy, d2_xz, d2_yy, d2_yz, d2_zz;
        compute_all_derivs_smem<GW>(smem3[1], &dgxy_x, &dgxy_y, &dgxy_z, &d2_xx, &d2_xy, &d2_xz, &d2_yy, &d2_yz, &d2_zz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
        if (is_valid) contract_xy = gupxx * d2_xx + gupyy * d2_yy + gupzz * d2_zz + (gupxy * d2_xy + gupxz * d2_xz + gupyz * d2_yz) * TWO;
    }

    double dgxz_x, dgxz_y, dgxz_z, contract_xz = 0;
    {
        double d2_xx, d2_xy, d2_xz, d2_yy, d2_yz, d2_zz;
        compute_all_derivs_smem<GW>(smem3[2], &dgxz_x, &dgxz_y, &dgxz_z, &d2_xx, &d2_xy, &d2_xz, &d2_yy, &d2_yz, &d2_zz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
        if (is_valid) contract_xz = gupxx * d2_xx + gupyy * d2_yy + gupzz * d2_zz + (gupxy * d2_xy + gupxz * d2_xz + gupyz * d2_yz) * TWO;
    }
    __syncthreads(); // 等待所有线程算完 Batch 1，才能覆写 smem3

    // --- Batch 2: dyy, gyz, dzz ---
    load_3_fields_to_smem<GW>(smem3, dyy, gyz, dzz, pdims, block_i, block_j, block_k);
    __syncthreads();

    double dgyy_x, dgyy_y, dgyy_z, contract_yy = 0;
    {
        double d2_xx, d2_xy, d2_xz, d2_yy, d2_yz, d2_zz;
        compute_all_derivs_smem<GW>(smem3[0], &dgyy_x, &dgyy_y, &dgyy_z, &d2_xx, &d2_xy, &d2_xz, &d2_yy, &d2_yz, &d2_zz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
        if (is_valid) contract_yy = gupxx * d2_xx + gupyy * d2_yy + gupzz * d2_zz + (gupxy * d2_xy + gupxz * d2_xz + gupyz * d2_yz) * TWO;
    }

    double dgyz_x, dgyz_y, dgyz_z, contract_yz = 0;
    {
        double d2_xx, d2_xy, d2_xz, d2_yy, d2_yz, d2_zz;
        compute_all_derivs_smem<GW>(smem3[1], &dgyz_x, &dgyz_y, &dgyz_z, &d2_xx, &d2_xy, &d2_xz, &d2_yy, &d2_yz, &d2_zz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
        if (is_valid) contract_yz = gupxx * d2_xx + gupyy * d2_yy + gupzz * d2_zz + (gupxy * d2_xy + gupxz * d2_xz + gupyz * d2_yz) * TWO;
    }

    double dgzz_x, dgzz_y, dgzz_z, contract_zz = 0;
    {
        double d2_xx, d2_xy, d2_xz, d2_yy, d2_yz, d2_zz;
        compute_all_derivs_smem<GW>(smem3[2], &dgzz_x, &dgzz_y, &dgzz_z, &d2_xx, &d2_xy, &d2_xz, &d2_yy, &d2_yz, &d2_zz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
        if (is_valid) contract_zz = gupxx * d2_xx + gupyy * d2_yy + gupzz * d2_zz + (gupxy * d2_xy + gupxz * d2_xz + gupyz * d2_yz) * TWO;
    }
    __syncthreads();

    // --- Batch 3: Gamx, Gamy, Gamz ---
    load_3_fields_to_smem<GW>(smem3, Gamx, Gamy, Gamz, pdims, block_i, block_j, block_k);
    __syncthreads();

    double dGamxx, dGamxy, dGamxz;
    compute_first_derivs_smem<GW>(smem3[0], &dGamxx, &dGamxy, &dGamxz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double dGamyx, dGamyy, dGamyz;
    compute_first_derivs_smem<GW>(smem3[1], &dGamyx, &dGamyy, &dGamyz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double dGamzx, dGamzy, dGamzz;
    compute_first_derivs_smem<GW>(smem3[2], &dGamzx, &dGamzy, &dGamzz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    // 此处不需要 __syncthreads()，因为 smem 已经完成了它的历史使命

    // =========================================================================
    // 2. 组装 Christoffel 与 Ricci 张量 (与原逻辑完全一致，无任何变更)
    // =========================================================================
    double G1_xxx = HALF * dgxx_x;
    double G1_xxy = HALF * dgxx_y;
    double G1_xxz = HALF * dgxx_z;
    double G1_xyy = HALF * (TWO * dgxy_y - dgyy_x);
    double G1_xyz = HALF * (dgxy_z + dgxz_y - dgyz_x);
    double G1_xzz = HALF * (TWO * dgxz_z - dgzz_x);

    double G1_yxx = HALF * (TWO * dgxy_x - dgxx_y);
    double G1_yxy = HALF * dgyy_x;
    double G1_yxz = HALF * (dgxy_z + dgyz_x - dgxz_y);
    double G1_yyy = HALF * dgyy_y;
    double G1_yyz = HALF * dgyy_z;
    double G1_yzz = HALF * (TWO * dgyz_z - dgzz_y);

    double G1_zxx = HALF * (TWO * dgxz_x - dgxx_z);
    double G1_zxy = HALF * (dgxz_y + dgyz_x - dgxy_z);
    double G1_zxz = HALF * dgzz_x;
    double G1_zyy = HALF * (TWO * dgyz_y - dgyy_z);
    double G1_zyz = HALF * dgzz_y;
    double G1_zzz = HALF * dgzz_z;

    double G2_xxx = gupxx * G1_xxx + gupxy * G1_yxx + gupxz * G1_zxx;
    double G2_xxy = gupxx * G1_xxy + gupxy * G1_yxy + gupxz * G1_zxy;
    double G2_xxz = gupxx * G1_xxz + gupxy * G1_yxz + gupxz * G1_zxz;
    double G2_xyy = gupxx * G1_xyy + gupxy * G1_yyy + gupxz * G1_zyy;
    double G2_xyz = gupxx * G1_xyz + gupxy * G1_yyz + gupxz * G1_zyz;
    double G2_xzz = gupxx * G1_xzz + gupxy * G1_yzz + gupxz * G1_zzz;

    double G2_yxx = gupxy * G1_xxx + gupyy * G1_yxx + gupyz * G1_zxx;
    double G2_yxy = gupxy * G1_xxy + gupyy * G1_yxy + gupyz * G1_zxy;
    double G2_yxz = gupxy * G1_xxz + gupyy * G1_yxz + gupyz * G1_zxz;
    double G2_yyy = gupxy * G1_xyy + gupyy * G1_yyy + gupyz * G1_zyy;
    double G2_yyz = gupxy * G1_xyz + gupyy * G1_yyz + gupyz * G1_zyz;
    double G2_yzz = gupxy * G1_xzz + gupyy * G1_yzz + gupyz * G1_zzz;

    double G2_zxx = gupxz * G1_xxx + gupyz * G1_yxx + gupzz * G1_zxx;
    double G2_zxy = gupxz * G1_xxy + gupyz * G1_yxy + gupzz * G1_zxy;
    double G2_zxz = gupxz * G1_xxz + gupyz * G1_yxz + gupzz * G1_zxz;
    double G2_zyy = gupxz * G1_xyy + gupyz * G1_yyy + gupzz * G1_zyy;
    double G2_zyz = gupxz * G1_xyz + gupyz * G1_yyz + gupzz * G1_zyz;
    double G2_zzz = gupxz * G1_xzz + gupyz * G1_yzz + gupzz * G1_zzz;

    double Gamxa = gupxx * G2_xxx + gupyy * G2_xyy + gupzz * G2_xzz + TWO * (gupxy * G2_xxy + gupxz * G2_xxz + gupyz * G2_xyz);
    double Gamya = gupxx * G2_yxx + gupyy * G2_yyy + gupzz * G2_yzz + TWO * (gupxy * G2_yxy + gupxz * G2_yxz + gupyz * G2_yyz);
    double Gamza = gupxx * G2_zxx + gupyy * G2_zyy + gupzz * G2_zzz + TWO * (gupxy * G2_zxy + gupxz * G2_zxz + gupyz * G2_zyz);

    if (is_valid) {
        Gamxxx_out[idx] = G2_xxx; Gamxxy_out[idx] = G2_xxy; Gamxxz_out[idx] = G2_xxz;
        Gamxyy_out[idx] = G2_xyy; Gamxyz_out[idx] = G2_xyz; Gamxzz_out[idx] = G2_xzz;
        Gamyxx_out[idx] = G2_yxx; Gamyxy_out[idx] = G2_yxy; Gamyxz_out[idx] = G2_yxz;
        Gamyyy_out[idx] = G2_yyy; Gamyyz_out[idx] = G2_yyz; Gamyzz_out[idx] = G2_yzz;
        Gamzxx_out[idx] = G2_zxx; Gamzxy_out[idx] = G2_zxy; Gamzxz_out[idx] = G2_zxz;
        Gamzyy_out[idx] = G2_zyy; Gamzyz_out[idx] = G2_zyz; Gamzzz_out[idx] = G2_zzz;

        double l_Rxx = -HALF * contract_xx + l_gxx * dGamxx + l_gxy * dGamyx + l_gxz * dGamzx + Gamxa * G1_xxx + Gamya * G1_xxy + Gamza * G1_xxz + 
              gupxx * (TWO*(G2_xxx*G1_xxx + G2_yxx*G1_xxy + G2_zxx*G1_xxz) + G2_xxx*G1_xxx + G2_yxx*G1_yxx + G2_zxx*G1_zxx) +
              gupxy * (TWO*(G2_xxx*G1_xxy + G2_yxx*G1_xyy + G2_zxx*G1_xyz + G2_xxy*G1_xxx + G2_yxy*G1_xxy + G2_zxy*G1_xxz) + G2_xxy*G1_xxx + G2_yxy*G1_yxx + G2_zxy*G1_zxx + G2_xxx*G1_xxy + G2_yxx*G1_yxy + G2_zxx*G1_zxy) + 
              gupxz * (TWO*(G2_xxx*G1_xxz + G2_yxx*G1_xyz + G2_zxx*G1_xzz + G2_xxz*G1_xxx + G2_yxz*G1_xxy + G2_zxz*G1_xxz) + G2_xxz*G1_xxx + G2_yxz*G1_yxx + G2_zxz*G1_zxx + G2_xxx*G1_xxz + G2_yxx*G1_yxz + G2_zxx*G1_zxz) + 
              gupyy * (TWO*(G2_xxy*G1_xxy + G2_yxy*G1_xyy + G2_zxy*G1_xyz) + G2_xxy*G1_xxy + G2_yxy*G1_yxy + G2_zxy*G1_zxy) + 
              gupyz * (TWO*(G2_xxy*G1_xxz + G2_yxy*G1_xyz + G2_zxy*G1_xzz + G2_xxz*G1_xxy + G2_yxz*G1_xyy + G2_zxz*G1_xyz) + G2_xxz*G1_xxy + G2_yxz*G1_yxy + G2_zxz*G1_zxy + G2_xxy*G1_xxz + G2_yxy*G1_yxz + G2_zxy*G1_zxz) + 
              gupzz * (TWO*(G2_xxz*G1_xxz + G2_yxz*G1_xyz + G2_zxz*G1_xzz) + G2_xxz*G1_xxz + G2_yxz*G1_yxz + G2_zxz*G1_zxz);

        double l_Ryy = -HALF * contract_yy + l_gxy * dGamxy + l_gyy * dGamyy + l_gyz * dGamzy + Gamxa * G1_yxy + Gamya * G1_yyy + Gamza * G1_yyz + 
              gupxx * (TWO*(G2_xxy*G1_yxx + G2_yxy*G1_yxy + G2_zxy*G1_yxz) + G2_xxy*G1_xxy + G2_yxy*G1_yxy + G2_zxy*G1_zxy) + 
              gupxy * (TWO*(G2_xxy*G1_yxy + G2_yxy*G1_yyy + G2_zxy*G1_yyz + G2_xyy*G1_yxx + G2_yyy*G1_yxy + G2_zyy*G1_yxz) + G2_xyy*G1_xxy + G2_yyy*G1_yxy + G2_zyy*G1_zxy + G2_xxy*G1_xyy + G2_yxy*G1_yyy + G2_zxy*G1_zyy) + 
              gupxz * (TWO*(G2_xxy*G1_yxz + G2_yxy*G1_yyz + G2_zxy*G1_yzz + G2_xyz*G1_yxx + G2_yyz*G1_yxy + G2_zyz*G1_yxz) + G2_xyz*G1_xxy + G2_yyz*G1_yxy + G2_zyz*G1_zxy + G2_xxy*G1_xyz + G2_yxy*G1_yyz + G2_zxy*G1_zyz) + 
              gupyy * (TWO*(G2_xyy*G1_yxy + G2_yyy*G1_yyy + G2_zyy*G1_yyz) + G2_xyy*G1_xyy + G2_yyy*G1_yyy + G2_zyy*G1_zyy) + 
              gupyz * (TWO*(G2_xyy*G1_yxz + G2_yyy*G1_yyz + G2_zyy*G1_yzz + G2_xyz*G1_yxy + G2_yyz*G1_yyy + G2_zyz*G1_yyz) + G2_xyz*G1_xyy + G2_yyz*G1_yyy + G2_zyz*G1_zyy + G2_xyy*G1_xyz + G2_yyy*G1_yyz + G2_zyy*G1_zyz) + 
              gupzz * (TWO*(G2_xyz*G1_yxz + G2_yyz*G1_yyz + G2_zyz*G1_yzz) + G2_xyz*G1_xyz + G2_yyz*G1_yyz + G2_zyz*G1_zyz);

        double l_Rzz = -HALF * contract_zz + l_gxz * dGamxz + l_gyz * dGamyz + l_gzz * dGamzz + Gamxa * G1_zxz + Gamya * G1_zyz + Gamza * G1_zzz + 
              gupxx * (TWO*(G2_xxz*G1_zxx + G2_yxz*G1_zxy + G2_zxz*G1_zxz) + G2_xxz*G1_xxz + G2_yxz*G1_yxz + G2_zxz*G1_zxz) + 
              gupxy * (TWO*(G2_xxz*G1_zxy + G2_yxz*G1_zyy + G2_zxz*G1_zyz + G2_xyz*G1_zxx + G2_yyz*G1_zxy + G2_zyz*G1_zxz) + G2_xyz*G1_xxz + G2_yyz*G1_yxz + G2_zyz*G1_zxz + G2_xxz*G1_xyz + G2_yxz*G1_yyz + G2_zxz*G1_zyz) + 
              gupxz * (TWO*(G2_xxz*G1_zxz + G2_yxz*G1_zyz + G2_zxz*G1_zzz + G2_xzz*G1_zxx + G2_yzz*G1_zxy + G2_zzz*G1_zxz) + G2_xzz*G1_xxz + G2_yzz*G1_yxz + G2_zzz*G1_zxz + G2_xxz*G1_xzz + G2_yxz*G1_yzz + G2_zxz*G1_zzz) + 
              gupyy * (TWO*(G2_xyz*G1_zxy + G2_yyz*G1_zyy + G2_zyz*G1_zyz) + G2_xyz*G1_xyz + G2_yyz*G1_yyz + G2_zyz*G1_zyz) + 
              gupyz * (TWO*(G2_xyz*G1_zxz + G2_yyz*G1_zyz + G2_zyz*G1_zzz + G2_xzz*G1_zxy + G2_yzz*G1_zyy + G2_zzz*G1_zyz) + G2_xzz*G1_xyz + G2_yzz*G1_yyz + G2_zzz*G1_zyz + G2_xyz*G1_xzz + G2_yyz*G1_yzz + G2_zyz*G1_zzz) + 
              gupzz * (TWO*(G2_xzz*G1_zxz + G2_yzz*G1_zyz + G2_zzz*G1_zzz) + G2_xzz*G1_xzz + G2_yzz*G1_yzz + G2_zzz*G1_zzz);

        double l_Rxy = HALF * ( - contract_xy + l_gxx * dGamxy + l_gxy * dGamyy + l_gxz * dGamzy + l_gxy * dGamxx + l_gyy * dGamyx + l_gyz * dGamzx + Gamxa * G1_xxy + Gamya * G1_xyy + Gamza * G1_xyz + Gamxa * G1_yxx + Gamya * G1_yxy + Gamza * G1_yxz) + 
              gupxx * (G2_xxx*G1_yxx + G2_yxx*G1_yxy + G2_zxx*G1_yxz + G2_xxy*G1_xxx + G2_yxy*G1_xxy + G2_zxy*G1_xxz + G2_xxx*G1_xxy + G2_yxx*G1_yxy + G2_zxx*G1_zxy) + 
              gupxy * (G2_xxx*G1_yxy + G2_yxx*G1_yyy + G2_zxx*G1_yyz + G2_xxy*G1_xxy + G2_yxy*G1_xyy + G2_zxy*G1_xyz + G2_xxy*G1_xxy + G2_yxy*G1_yxy + G2_zxy*G1_zxy + G2_xxy*G1_yxx + G2_yxy*G1_yxy + G2_zxy*G1_yxz + G2_xyy*G1_xxx + G2_yyy*G1_xxy + G2_zyy*G1_xxz + G2_xxx*G1_xyy + G2_yxx*G1_yyy + G2_zxx*G1_zyy) + 
              gupxz * (G2_xxx*G1_yxz + G2_yxx*G1_yyz + G2_zxx*G1_yzz + G2_xxy*G1_xxz + G2_yxy*G1_xyz + G2_zxy*G1_xzz + G2_xxz*G1_xxy + G2_yxz*G1_yxy + G2_zxz*G1_zxy + G2_xxz*G1_yxx + G2_yxz*G1_yxy + G2_zxz*G1_yxz + G2_xyz*G1_xxx + G2_yyz*G1_xxy + G2_zyz*G1_xxz + G2_xxx*G1_xyz + G2_yxx*G1_yyz + G2_zxx*G1_zyz) + 
              gupyy * (G2_xxy*G1_yxy + G2_yxy*G1_yyy + G2_zxy*G1_yyz + G2_xyy*G1_xxy + G2_yyy*G1_xyy + G2_zyy*G1_xyz + G2_xxy*G1_xyy + G2_yxy*G1_yyy + G2_zxy*G1_zyy) + 
              gupyz * (G2_xxy*G1_yxz + G2_yxy*G1_yyz + G2_zxy*G1_yzz + G2_xyy*G1_xxz + G2_yyy*G1_xyz + G2_zyy*G1_xzz + G2_xxz*G1_xyy + G2_yxz*G1_yyy + G2_zxz*G1_zyy + G2_xxz*G1_yxy + G2_yxz*G1_yyy + G2_zxz*G1_yyz + G2_xyz*G1_xxy + G2_yyz*G1_xyy + G2_zyz*G1_xyz + G2_xxy*G1_xyz + G2_yxy*G1_yyz + G2_zxy*G1_zyz) + 
              gupzz * (G2_xxz*G1_yxz + G2_yxz*G1_yyz + G2_zxz*G1_yzz + G2_xyz*G1_xxz + G2_yyz*G1_xyz + G2_zyz*G1_xzz + G2_xxz*G1_xyz + G2_yxz*G1_yyz + G2_zxz*G1_zyz);

        double l_Rxz = HALF * ( - contract_xz + l_gxx * dGamxz + l_gxy * dGamyz + l_gxz * dGamzz + l_gxz * dGamxx + l_gyz * dGamyx + l_gzz * dGamzx + Gamxa * G1_xxz + Gamya * G1_xyz + Gamza * G1_xzz + Gamxa * G1_zxx + Gamya * G1_zxy + Gamza * G1_zxz) + 
              gupxx * (G2_xxx*G1_zxx + G2_yxx*G1_zxy + G2_zxx*G1_zxz + G2_xxz*G1_xxx + G2_yxz*G1_xxy + G2_zxz*G1_xxz + G2_xxx*G1_xxz + G2_yxx*G1_yxz + G2_zxx*G1_zxz) + 
              gupxy * (G2_xxx*G1_zxy + G2_yxx*G1_zyy + G2_zxx*G1_zyz + G2_xxz*G1_xxy + G2_yxz*G1_xyy + G2_zxz*G1_xyz + G2_xxy*G1_xxz + G2_yxy*G1_yxz + G2_zxy*G1_zxz + G2_xxy*G1_zxx + G2_yxy*G1_zxy + G2_zxy*G1_zxz + G2_xyz*G1_xxx + G2_yyz*G1_xxy + G2_zyz*G1_xxz + G2_xxx*G1_xyz + G2_yxx*G1_yyz + G2_zxx*G1_zyz) + 
              gupxz * (G2_xxx*G1_zxz + G2_yxx*G1_zyz + G2_zxx*G1_zzz + G2_xxz*G1_xxz + G2_yxz*G1_xyz + G2_zxz*G1_xzz + G2_xxz*G1_xxz + G2_yxz*G1_yxz + G2_zxz*G1_zxz + G2_xxz*G1_zxx + G2_yxz*G1_zxy + G2_zxz*G1_zxz + G2_xzz*G1_xxx + G2_yzz*G1_xxy + G2_zzz*G1_xxz + G2_xxx*G1_xzz + G2_yxx*G1_yzz + G2_zxx*G1_zzz) + 
              gupyy * (G2_xxy*G1_zxy + G2_yxy*G1_zyy + G2_zxy*G1_zyz + G2_xyz*G1_xxy + G2_yyz*G1_xyy + G2_zyz*G1_xyz + G2_xxy*G1_xyz + G2_yxy*G1_yyz + G2_zxy*G1_zyz) + 
              gupyz * (G2_xxy*G1_zxz + G2_yxy*G1_zyz + G2_zxy*G1_zzz + G2_xyz*G1_xxz + G2_yyz*G1_xyz + G2_zyz*G1_xzz + G2_xxz*G1_xyz + G2_yxz*G1_yyz + G2_zxz*G1_zyz + G2_xxz*G1_zxy + G2_yxz*G1_zyy + G2_zxz*G1_zyz + G2_xzz*G1_xxy + G2_yzz*G1_xyy + G2_zzz*G1_xyz + G2_xxy*G1_xzz + G2_yxy*G1_yzz + G2_zxy*G1_zzz) + 
              gupzz * (G2_xxz*G1_zxz + G2_yxz*G1_zyz + G2_zxz*G1_zzz + G2_xzz*G1_xxz + G2_yzz*G1_xyz + G2_zzz*G1_xzz + G2_xxz*G1_xzz + G2_yxz*G1_yzz + G2_zxz*G1_zzz);

        double l_Ryz = HALF * ( - contract_yz + l_gxy * dGamxz + l_gyy * dGamyz + l_gyz * dGamzz + l_gxz * dGamxy + l_gyz * dGamyy + l_gzz * dGamzy + Gamxa * G1_yxz + Gamya * G1_yyz + Gamza * G1_yzz + Gamxa * G1_zxy + Gamya * G1_zyy + Gamza * G1_zyz) + 
              gupxx * (G2_xxy*G1_zxx + G2_yxy*G1_zxy + G2_zxy*G1_zxz + G2_xxz*G1_yxx + G2_yxz*G1_yxy + G2_zxz*G1_yxz + G2_xxy*G1_xxz + G2_yxy*G1_yxz + G2_zxy*G1_zxz) + 
              gupxy * (G2_xxy*G1_zxy + G2_yxy*G1_zyy + G2_zxy*G1_zyz + G2_xxz*G1_yxy + G2_yxz*G1_yyy + G2_zxz*G1_yyz + G2_xyy*G1_xxz + G2_yyy*G1_yxz + G2_zyy*G1_zxz + G2_xyy*G1_zxx + G2_yyy*G1_zxy + G2_zyy*G1_zxz + G2_xyz*G1_yxx + G2_yyz*G1_yxy + G2_zyz*G1_yxz + G2_xxy*G1_xyz + G2_yxy*G1_yyz + G2_zxy*G1_zyz) + 
              gupxz * (G2_xxy*G1_zxz + G2_yxy*G1_zyz + G2_zxy*G1_zzz + G2_xxz*G1_yxz + G2_yxz*G1_yyz + G2_zxz*G1_yzz + G2_xyz*G1_xxz + G2_yyz*G1_yxz + G2_zyz*G1_zxz + G2_xyz*G1_zxx + G2_yyz*G1_zxy + G2_zyz*G1_zxz + G2_xzz*G1_yxx + G2_yzz*G1_yxy + G2_zzz*G1_yxz + G2_xxy*G1_xzz + G2_yxy*G1_yzz + G2_zxy*G1_zzz) + 
              gupyy * (G2_xyy*G1_zxy + G2_yyy*G1_zyy + G2_zyy*G1_zyz + G2_xyz*G1_yxy + G2_yyz*G1_yyy + G2_zyz*G1_yyz + G2_xyy*G1_xyz + G2_yyy*G1_yyz + G2_zyy*G1_zyz) + 
              gupyz * (G2_xyy*G1_zxz + G2_yyy*G1_zyz + G2_zyy*G1_zzz + G2_xyz*G1_yxz + G2_yyz*G1_yyz + G2_zyz*G1_yzz + G2_xyz*G1_xyz + G2_yyz*G1_yyz + G2_zyz*G1_zyz + G2_xyz*G1_zxy + G2_yyz*G1_zyy + G2_zyz*G1_zyz + G2_xzz*G1_yxy + G2_yzz*G1_yyy + G2_zzz*G1_yyz + G2_xyy*G1_xzz + G2_yyy*G1_yzz + G2_zyy*G1_zzz) + 
              gupzz * (G2_xyz*G1_zxz + G2_yyz*G1_zyz + G2_zyz*G1_zzz + G2_xzz*G1_yxz + G2_yzz*G1_yyz + G2_zzz*G1_yzz + G2_xyz*G1_xzz + G2_yyz*G1_yzz + G2_zyz*G1_zzz);

        Rxx_out[idx] = l_Rxx; Rxy_out[idx] = l_Rxy; Rxz_out[idx] = l_Rxz;
        Ryy_out[idx] = l_Ryy; Ryz_out[idx] = l_Ryz; Rzz_out[idx] = l_Rzz;
    }
}

template<int GW, int GAUGE_PARAM>
__launch_bounds__(256, 1)
__global__ void bssn_rhs_eval_kernel(
    const int ex0, const int ex1, const int ex2,
    const int pex0, const int pex1, const int pex2,
    const int symmetry, const int lev,
    const int BHN,
    const double Porg1, const double Porg2, const double Porg3,
    const double Porg4, const double Porg5, const double Porg6,
    const double Mass1, const double Mass2,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Axx, const double* __restrict__ Axy, const double* __restrict__ Axz,
    const double* __restrict__ Ayy, const double* __restrict__ Ayz, const double* __restrict__ Azz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    const double* __restrict__ Lap, 
    const double* __restrict__ betax, const double* __restrict__ betay, const double* __restrict__ betaz,
    const double* __restrict__ dtSfx, const double* __restrict__ dtSfy, const double* __restrict__ dtSfz,
    const double* __restrict__ rho, const double* __restrict__ Sx, const double* __restrict__ Sy, const double* __restrict__ Sz,
    const double* __restrict__ Sxx, const double* __restrict__ Sxy, const double* __restrict__ Sxz, 
    const double* __restrict__ Syy, const double* __restrict__ Syz, const double* __restrict__ Szz,
    const double* __restrict__ Gamxxx_in, const double* __restrict__ Gamxxy_in, const double* __restrict__ Gamxxz_in,
    const double* __restrict__ Gamxyy_in, const double* __restrict__ Gamxyz_in, const double* __restrict__ Gamxzz_in,
    const double* __restrict__ Gamyxx_in, const double* __restrict__ Gamyxy_in, const double* __restrict__ Gamyxz_in,
    const double* __restrict__ Gamyyy_in, const double* __restrict__ Gamyyz_in, const double* __restrict__ Gamyzz_in,
    const double* __restrict__ Gamzxx_in, const double* __restrict__ Gamzxy_in, const double* __restrict__ Gamzxz_in,
    const double* __restrict__ Gamzyy_in, const double* __restrict__ Gamzyz_in, const double* __restrict__ Gamzzz_in,
    const double* __restrict__ Rxx_in, const double* __restrict__ Rxy_in, const double* __restrict__ Rxz_in,
    const double* __restrict__ Ryy_in, const double* __restrict__ Ryz_in, const double* __restrict__ Rzz_in,
    double* __restrict__ chi_rhs, double* __restrict__ trK_rhs,
    double* __restrict__ gxx_rhs, double* __restrict__ gxy_rhs, double* __restrict__ gxz_rhs,
    double* __restrict__ gyy_rhs, double* __restrict__ gyz_rhs, double* __restrict__ gzz_rhs,
    double* __restrict__ Axx_rhs, double* __restrict__ Axy_rhs, double* __restrict__ Axz_rhs,
    double* __restrict__ Ayy_rhs, double* __restrict__ Ayz_rhs, double* __restrict__ Azz_rhs,
    double* __restrict__ Gamx_rhs, double* __restrict__ Gamy_rhs, double* __restrict__ Gamz_rhs,
    double* __restrict__ Lap_rhs, 
    double* __restrict__ betax_rhs, double* __restrict__ betay_rhs, double* __restrict__ betaz_rhs,
    double* __restrict__ dtSfx_rhs, double* __restrict__ dtSfy_rhs, double* __restrict__ dtSfz_rhs
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int dims[3]  = {ex0, ex1, ex2};
    int pdims[3] = {pex0, pex1, pex2};
    double dX = X[1] - X[0]; double dY = Y[1] - Y[0]; double dZ = Z[1] - Z[0];

    int imax = ex0 - 1; int jmax = ex1 - 1; int kmax = ex2 - 1;
    int imin = 0, jmin = 0, kmin = 0;
    constexpr int ORD = GW - 1;
    if (symmetry > 0 && fabs(Z[0]) < dZ) kmin = -ORD;
    if (symmetry > 1 && fabs(X[0]) < dX) imin = -ORD;
    if (symmetry > 1 && fabs(Y[0]) < dY) jmin = -ORD;

    int block_i = blockIdx.x * blockDim.x;
    int block_j = blockIdx.y * blockDim.y;
    int block_k = blockIdx.z * blockDim.z;

    // Dynamic shared memory: 3 tiles of [Z][Y][X] doubles.
    // For GW=5 this is 3*3072*8 = 73728 bytes, exceeding the 48KB static smem limit on sm_80.
    extern __shared__ double smem_raw[];
    double (*smem3)[TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X] =
        reinterpret_cast<double(*)[TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X]>(smem_raw);

    // =========================================================================
    // 1. Batch 1: \beta^x, \beta^y, \beta^z
    // =========================================================================
    load_3_fields_to_smem<GW>(smem3, betax, betay, betaz, pdims, block_i, block_j, block_k);
    __syncthreads();

    double betaxx, betaxy, betaxz, bx_gxxx, bx_gxyx, bx_gxzx, bx_gyyx, bx_gyzx, bx_gzzx;
    compute_all_derivs_smem<GW>(smem3[0], &betaxx, &betaxy, &betaxz, &bx_gxxx, &bx_gxyx, &bx_gxzx, &bx_gyyx, &bx_gyzx, &bx_gzzx, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double betayx, betayy, betayz, by_gxxy, by_gxyy, by_gxzy, by_gyyy, by_gyzy, by_gzzy;
    compute_all_derivs_smem<GW>(smem3[1], &betayx, &betayy, &betayz, &by_gxxy, &by_gxyy, &by_gxzy, &by_gyyy, &by_gyzy, &by_gzzy, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double betazx, betazy, betazz, bz_gxxz, bz_gxyz, bz_gxzz, bz_gyyz, bz_gyzz, bz_gzzz;
    compute_all_derivs_smem<GW>(smem3[2], &betazx, &betazy, &betazz, &bz_gxxz, &bz_gxyz, &bz_gxzz, &bz_gyyz, &bz_gyzz, &bz_gzzz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    
    __syncthreads(); // 等待所有线程算完 Batch 1，才能覆写 smem3

    // =========================================================================
    // 2. Batch 2: \alpha (Lap), \chi, K (trK)
    // =========================================================================
    load_3_fields_to_smem<GW>(smem3, Lap, chi, trK, pdims, block_i, block_j, block_k);
    __syncthreads();

    double Lapx, Lapy, Lapz, fxx_Lap, fxy_Lap, fxz_Lap, fyy_Lap, fyz_Lap, fzz_Lap;
    compute_all_derivs_smem<GW>(smem3[0], &Lapx, &Lapy, &Lapz, &fxx_Lap, &fxy_Lap, &fxz_Lap, &fyy_Lap, &fyz_Lap, &fzz_Lap, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double chix, chiy, chiz, chixx, chixy, chixz, chiyy, chiyz, chizz;
    compute_all_derivs_smem<GW>(smem3[1], &chix, &chiy, &chiz, &chixx, &chixy, &chixz, &chiyy, &chiyz, &chizz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double Kx, Ky, Kz;
    compute_first_derivs_smem<GW>(smem3[2], &Kx, &Ky, &Kz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    // 此处不需要 __syncthreads()，因为 smem 的使命已经结束

    // =========================================================================
    // 3. 物理组装 (完全保留原有数学逻辑)
    // =========================================================================
    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = IDX3D_PAD(i, j, k, pex0, pex1);

    const double ONE = 1.0; const double TWO = 2.0; const double F1o3 = 1.0 / 3.0; const double F2o3 = 2.0 / 3.0;
    const double F3o2 = 1.5; const double HALF = 0.5; const double EIGHT = 8.0; const double FOUR = 4.0; const double F16 = 16.0; const double PI = M_PI;

    double val_Lap = Lap[idx]; double val_chi = chi[idx]; double val_trK = trK[idx];
    double alpn1 = val_Lap + ONE; double chin1 = val_chi + ONE;
    double inv_chin1 = 1.0 / chin1;

    double l_gxx = dxx[idx] + ONE; double l_gxy = gxy[idx]; double l_gxz = gxz[idx];
    double l_gyy = dyy[idx] + ONE; double l_gyz = gyz[idx]; double l_gzz = dzz[idx] + ONE;
    double l_Axx = Axx[idx]; double l_Axy = Axy[idx]; double l_Axz = Axz[idx];
    double l_Ayy = Ayy[idx]; double l_Ayz = Ayz[idx]; double l_Azz = Azz[idx];

    double div_beta = betaxx + betayy + betazz;

    chi_rhs[idx] = F2o3 * chin1 * (alpn1 * val_trK - div_beta);

    gxx_rhs[idx] = -TWO * alpn1 * l_Axx - F2o3 * l_gxx * div_beta + TWO * (l_gxx * betaxx + l_gxy * betayx + l_gxz * betazx);
    gyy_rhs[idx] = -TWO * alpn1 * l_Ayy - F2o3 * l_gyy * div_beta + TWO * (l_gxy * betaxy + l_gyy * betayy + l_gyz * betazy);
    gzz_rhs[idx] = -TWO * alpn1 * l_Azz - F2o3 * l_gzz * div_beta + TWO * (l_gxz * betaxz + l_gyz * betayz + l_gzz * betazz);
    gxy_rhs[idx] = -TWO * alpn1 * l_Axy + F1o3 * l_gxy * div_beta + l_gxx * betaxy + l_gxz * betazy + l_gyy * betayx + l_gyz * betazx - l_gxy * betazz;
    gyz_rhs[idx] = -TWO * alpn1 * l_Ayz + F1o3 * l_gyz * div_beta + l_gxy * betaxz + l_gyy * betayz + l_gxz * betaxy + l_gzz * betazy - l_gyz * betaxx;
    gxz_rhs[idx] = -TWO * alpn1 * l_Axz + F1o3 * l_gxz * div_beta + l_gxx * betaxz + l_gxy * betayz + l_gyz * betayx + l_gzz * betazx - l_gxz * betayy;

    double det = l_gxx * l_gyy * l_gzz + TWO * l_gxy * l_gyz * l_gxz - l_gxz * l_gyy * l_gxz - l_gxx * l_gyz * l_gyz - l_gzz * l_gxy * l_gxy;
    double inv_det = 1.0 / det;

    double gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) * inv_det;
    double gupxy = (l_gyz * l_gxz - l_gxy * l_gzz) * inv_det;
    double gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) * inv_det;
    double gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) * inv_det;
    double gupyz = (l_gxy * l_gxz - l_gxx * l_gyz) * inv_det;
    double gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) * inv_det;

    // 提前消费 18 个 β 二阶导 + 3 个 fXX_beta，在 l_GamXXX 加载前使其死亡
    // 节省约 21 个寄存器的 peak liveness
    double fxx_beta = bx_gxxx + by_gxyy + bz_gxzz;
    double fxy_beta = bx_gxyx + by_gyyy + bz_gyzz;
    double fxz_beta = bx_gxzx + by_gyzy + bz_gzzz;
    double partial_Gamx = F1o3 * (gupxx * fxx_beta + gupxy * fxy_beta + gupxz * fxz_beta)
        + gupxx * bx_gxxx + gupyy * bx_gyyx + gupzz * bx_gzzx
        + TWO * (gupxy * bx_gxyx + gupxz * bx_gxzx + gupyz * bx_gyzx);
    double partial_Gamy = F1o3 * (gupxy * fxx_beta + gupyy * fxy_beta + gupyz * fxz_beta)
        + gupxx * by_gxxy + gupyy * by_gyyy + gupzz * by_gzzy
        + TWO * (gupxy * by_gxyy + gupxz * by_gxzy + gupyz * by_gyzy);
    double partial_Gamz = F1o3 * (gupxz * fxx_beta + gupyz * fxy_beta + gupzz * fxz_beta)
        + gupxx * bz_gxxz + gupyy * bz_gyyz + gupzz * bz_gzzz
        + TWO * (gupxy * bz_gxyz + gupxz * bz_gxzz + gupyz * bz_gyzz);

    double Au_d_xx = gupxx * l_Axx + gupxy * l_Axy + gupxz * l_Axz;
    double Au_d_xy = gupxx * l_Axy + gupxy * l_Ayy + gupxz * l_Ayz;
    double Au_d_xz = gupxx * l_Axz + gupxy * l_Ayz + gupxz * l_Azz;
    double Au_d_yx = gupxy * l_Axx + gupyy * l_Axy + gupyz * l_Axz;
    double Au_d_yy = gupxy * l_Axy + gupyy * l_Ayy + gupyz * l_Ayz;
    double Au_d_yz = gupxy * l_Axz + gupyy * l_Ayz + gupyz * l_Azz;
    double Au_d_zx = gupxz * l_Axx + gupyz * l_Axy + gupzz * l_Axz;
    double Au_d_zy = gupxz * l_Axy + gupyz * l_Ayy + gupzz * l_Ayz;
    double Au_d_zz = gupxz * l_Axz + gupyz * l_Ayz + gupzz * l_Azz;

    double Aup_xx = Au_d_xx * gupxx + Au_d_xy * gupxy + Au_d_xz * gupxz;
    double Aup_xy = Au_d_xx * gupxy + Au_d_xy * gupyy + Au_d_xz * gupyz;
    double Aup_xz = Au_d_xx * gupxz + Au_d_xy * gupyz + Au_d_xz * gupzz;
    double Aup_yy = Au_d_yx * gupxy + Au_d_yy * gupyy + Au_d_yz * gupyz;
    double Aup_yz = Au_d_yx * gupxz + Au_d_yy * gupyz + Au_d_yz * gupzz;
    double Aup_zz = Au_d_zx * gupxz + Au_d_zy * gupyz + Au_d_zz * gupzz;

    // 提前消费 Au_d_XX (9 个)，节省净 3 个寄存器
    double term_xx = l_Axx * Au_d_xx + l_Axy * Au_d_yx + l_Axz * Au_d_zx;
    double term_xy = l_Axx * Au_d_xy + l_Axy * Au_d_yy + l_Axz * Au_d_zy;
    double term_xz = l_Axx * Au_d_xz + l_Axy * Au_d_yz + l_Axz * Au_d_zz;
    double term_yy = l_Axy * Au_d_xy + l_Ayy * Au_d_yy + l_Ayz * Au_d_zy;
    double term_yz = l_Axy * Au_d_xz + l_Ayy * Au_d_yz + l_Ayz * Au_d_zz;
    double term_zz = l_Axz * Au_d_xz + l_Ayz * Au_d_yz + l_Azz * Au_d_zz;

    const double l_Gamxxx = Gamxxx_in[idx]; const double l_Gamxxy = Gamxxy_in[idx]; const double l_Gamxxz = Gamxxz_in[idx];
    const double l_Gamxyy = Gamxyy_in[idx]; const double l_Gamxyz = Gamxyz_in[idx]; const double l_Gamxzz = Gamxzz_in[idx];
    const double l_Gamyxx = Gamyxx_in[idx]; const double l_Gamyxy = Gamyxy_in[idx]; const double l_Gamyxz = Gamyxz_in[idx];
    const double l_Gamyyy = Gamyyy_in[idx]; const double l_Gamyyz = Gamyyz_in[idx]; const double l_Gamyzz = Gamyzz_in[idx];
    const double l_Gamzxx = Gamzxx_in[idx]; const double l_Gamzxy = Gamzxy_in[idx]; const double l_Gamzxz = Gamzxz_in[idx];
    const double l_Gamzyy = Gamzyy_in[idx]; const double l_Gamzyz = Gamzyz_in[idx]; const double l_Gamzzz = Gamzzz_in[idx];

    double l_Rxx = Rxx_in[idx]; double l_Rxy = Rxy_in[idx]; double l_Rxz = Rxz_in[idx];
    double l_Ryy = Ryy_in[idx]; double l_Ryz = Ryz_in[idx]; double l_Rzz = Rzz_in[idx];

    double val_Sx = Sx[idx]; double val_Sy = Sy[idx]; double val_Sz = Sz[idx];

    double val_Gamx_rhs = - TWO * (Lapx * Aup_xx + Lapy * Aup_xy + Lapz * Aup_xz) + 
        TWO * alpn1 * (-F3o2 * inv_chin1 * (chix * Aup_xx + chiy * Aup_xy + chiz * Aup_xz) - gupxx * (F2o3 * Kx + EIGHT * PI * val_Sx) - gupxy * (F2o3 * Ky + EIGHT * PI * val_Sy) - gupxz * (F2o3 * Kz + EIGHT * PI * val_Sz) + 
        l_Gamxxx * Aup_xx + l_Gamxyy * Aup_yy + l_Gamxzz * Aup_zz + TWO * (l_Gamxxy * Aup_xy + l_Gamxxz * Aup_xz + l_Gamxyz * Aup_yz));

    double val_Gamy_rhs = - TWO * (Lapx * Aup_xy + Lapy * Aup_yy + Lapz * Aup_yz) + 
        TWO * alpn1 * (-F3o2 * inv_chin1 * (chix * Aup_xy + chiy * Aup_yy + chiz * Aup_yz) - gupxy * (F2o3 * Kx + EIGHT * PI * val_Sx) - gupyy * (F2o3 * Ky + EIGHT * PI * val_Sy) - gupyz * (F2o3 * Kz + EIGHT * PI * val_Sz) + 
        l_Gamyxx * Aup_xx + l_Gamyyy * Aup_yy + l_Gamyzz * Aup_zz + TWO * (l_Gamyxy * Aup_xy + l_Gamyxz * Aup_xz + l_Gamyyz * Aup_yz));

    double val_Gamz_rhs = - TWO * (Lapx * Aup_xz + Lapy * Aup_yz + Lapz * Aup_zz) + 
        TWO * alpn1 * (-F3o2 * inv_chin1 * (chix * Aup_xz + chiy * Aup_yz + chiz * Aup_zz) - gupxz * (F2o3 * Kx + EIGHT * PI * val_Sx) - gupyz * (F2o3 * Ky + EIGHT * PI * val_Sy) - gupzz * (F2o3 * Kz + EIGHT * PI * val_Sz) + 
        l_Gamzxx * Aup_xx + l_Gamzyy * Aup_yy + l_Gamzzz * Aup_zz + TWO * (l_Gamzxy * Aup_xy + l_Gamzxz * Aup_xz + l_Gamzyz * Aup_yz));

    double Gamxa = gupxx * l_Gamxxx + gupyy * l_Gamxyy + gupzz * l_Gamxzz + TWO * (gupxy * l_Gamxxy + gupxz * l_Gamxxz + gupyz * l_Gamxyz);
    double Gamya = gupxx * l_Gamyxx + gupyy * l_Gamyyy + gupzz * l_Gamyzz + TWO * (gupxy * l_Gamyxy + gupxz * l_Gamyxz + gupyz * l_Gamyyz);
    double Gamza = gupxx * l_Gamzxx + gupyy * l_Gamzyy + gupzz * l_Gamzzz + TWO * (gupxy * l_Gamzxy + gupxz * l_Gamzxz + gupyz * l_Gamzyz);

    val_Gamx_rhs += partial_Gamx + F2o3 * Gamxa * div_beta - (Gamxa * betaxx + Gamya * betaxy + Gamza * betaxz);
    val_Gamy_rhs += partial_Gamy + F2o3 * Gamya * div_beta - (Gamxa * betayx + Gamya * betayy + Gamza * betayz);
    val_Gamz_rhs += partial_Gamz + F2o3 * Gamza * div_beta - (Gamxa * betazx + Gamya * betazy + Gamza * betazz);
    
    double fxx = chixx - l_Gamxxx * chix - l_Gamyxx * chiy - l_Gamzxx * chiz;
    double fxy = chixy - l_Gamxxy * chix - l_Gamyxy * chiy - l_Gamzxy * chiz;
    double fxz = chixz - l_Gamxxz * chix - l_Gamyxz * chiy - l_Gamzxz * chiz;
    double fyy = chiyy - l_Gamxyy * chix - l_Gamyyy * chiy - l_Gamzyy * chiz;
    double fyz = chiyz - l_Gamxyz * chix - l_Gamyyz * chiy - l_Gamzyz * chiz;
    double fzz = chizz - l_Gamxzz * chix - l_Gamyzz * chiy - l_Gamzzz * chiz;

    double f_scalar = gupxx * (fxx - F3o2*inv_chin1 * chix * chix) + gupyy * (fyy - F3o2*inv_chin1 * chiy * chiy) + gupzz * (fzz - F3o2*inv_chin1 * chiz * chiz) + 
                      TWO * (gupxy * (fxy - F3o2*inv_chin1 * chix * chiy) + gupxz * (fxz - F3o2*inv_chin1 * chix * chiz) + gupyz * (fyz - F3o2*inv_chin1 * chiy * chiz));
    
    l_Rxx += (fxx - chix*chix*inv_chin1*HALF + l_gxx * f_scalar)*inv_chin1*HALF;
    l_Ryy += (fyy - chiy*chiy*inv_chin1*HALF + l_gyy * f_scalar)*inv_chin1*HALF;
    l_Rzz += (fzz - chiz*chiz*inv_chin1*HALF + l_gzz * f_scalar)*inv_chin1*HALF;
    l_Rxy += (fxy - chix*chiy*inv_chin1*HALF + l_gxy * f_scalar)*inv_chin1*HALF;
    l_Rxz += (fxz - chix*chiz*inv_chin1*HALF + l_gxz * f_scalar)*inv_chin1*HALF;
    l_Ryz += (fyz - chiy*chiz*inv_chin1*HALF + l_gyz * f_scalar)*inv_chin1*HALF;

    double gx_phy = (gupxx * chix + gupxy * chiy + gupxz * chiz) * inv_chin1;
    double gy_phy = (gupxy * chix + gupyy * chiy + gupyz * chiz) * inv_chin1;
    double gz_phy = (gupxz * chix + gupyz * chiy + gupzz * chiz) * inv_chin1;
    double chi_dot_Lap = gx_phy * Lapx + gy_phy * Lapy + gz_phy * Lapz;

    double D2_Lap_xx = fxx_Lap - l_Gamxxx*Lapx - l_Gamyxx*Lapy - l_Gamzxx*Lapz;
    double D2_Lap_yy = fyy_Lap - l_Gamxyy*Lapx - l_Gamyyy*Lapy - l_Gamzyy*Lapz;
    double D2_Lap_zz = fzz_Lap - l_Gamxzz*Lapx - l_Gamyzz*Lapy - l_Gamzzz*Lapz;
    double D2_Lap_xy = fxy_Lap - l_Gamxxy*Lapx - l_Gamyxy*Lapy - l_Gamzxy*Lapz;
    double D2_Lap_xz = fxz_Lap - l_Gamxxz*Lapx - l_Gamyxz*Lapy - l_Gamzxz*Lapz;
    double D2_Lap_yz = fyz_Lap - l_Gamxyz*Lapx - l_Gamyyz*Lapy - l_Gamzyz*Lapz;

    double fxx_Lap_phy = D2_Lap_xx + (chix * Lapx) * inv_chin1 - HALF * l_gxx * chi_dot_Lap;
    double fyy_Lap_phy = D2_Lap_yy + (chiy * Lapy) * inv_chin1 - HALF * l_gyy * chi_dot_Lap;
    double fzz_Lap_phy = D2_Lap_zz + (chiz * Lapz) * inv_chin1 - HALF * l_gzz * chi_dot_Lap;
    double fxy_Lap_phy = D2_Lap_xy + HALF * (chix * Lapy + chiy * Lapx) * inv_chin1 - HALF * l_gxy * chi_dot_Lap;
    double fxz_Lap_phy = D2_Lap_xz + HALF * (chix * Lapz + chiz * Lapx) * inv_chin1 - HALF * l_gxz * chi_dot_Lap;
    double fyz_Lap_phy = D2_Lap_yz + HALF * (chiy * Lapz + chiz * Lapy) * inv_chin1 - HALF * l_gyz * chi_dot_Lap;

    double conf_trK_rhs = gupxx * D2_Lap_xx + gupyy * D2_Lap_yy + gupzz * D2_Lap_zz + TWO * (gupxy * D2_Lap_xy + gupxz * D2_Lap_xz + gupyz * D2_Lap_yz);
    double trK_rhs_val = conf_trK_rhs - HALF * chi_dot_Lap; 

    double l_Sxx = Sxx[idx], l_Sxy = Sxy[idx], l_Sxz = Sxz[idx], l_Syy = Syy[idx], l_Syz = Syz[idx], l_Szz = Szz[idx];
    double S = chin1 * (gupxx * l_Sxx + gupyy * l_Syy + gupzz * l_Szz + TWO * (gupxy * l_Sxy + gupxz * l_Sxz + gupyz * l_Syz));

    double trA2 = l_Axx * Aup_xx + l_Ayy * Aup_yy + l_Azz * Aup_zz + TWO * (l_Axy * Aup_xy + l_Axz * Aup_xz + l_Ayz * Aup_yz);

    double f = F2o3 * val_trK * val_trK - trA2 - F16*PI*rho[idx] + EIGHT*PI*S;
    double f_trace = -F1o3 * (trK_rhs_val + alpn1*inv_chin1 * f);

    double src_xx = alpn1 * (l_Rxx - EIGHT*PI*l_Sxx) - fxx_Lap_phy;
    double src_yy = alpn1 * (l_Ryy - EIGHT*PI*l_Syy) - fyy_Lap_phy;
    double src_zz = alpn1 * (l_Rzz - EIGHT*PI*l_Szz) - fzz_Lap_phy;
    double src_xy = alpn1 * (l_Rxy - EIGHT*PI*l_Sxy) - fxy_Lap_phy;
    double src_xz = alpn1 * (l_Rxz - EIGHT*PI*l_Sxz) - fxz_Lap_phy;
    double src_yz = alpn1 * (l_Ryz - EIGHT*PI*l_Syz) - fyz_Lap_phy;

    Axx_rhs[idx] = chin1 * (src_xx - l_gxx * f_trace) + alpn1 * (val_trK * l_Axx - TWO * term_xx) + TWO * (l_Axx * betaxx + l_Axy * betayx + l_Axz * betazx) - F2o3 * l_Axx * div_beta;
    Ayy_rhs[idx] = chin1 * (src_yy - l_gyy * f_trace) + alpn1 * (val_trK * l_Ayy - TWO * term_yy) + TWO * (l_Axy * betaxy + l_Ayy * betayy + l_Ayz * betazy) - F2o3 * l_Ayy * div_beta;
    Azz_rhs[idx] = chin1 * (src_zz - l_gzz * f_trace) + alpn1 * (val_trK * l_Azz - TWO * term_zz) + TWO * (l_Axz * betaxz + l_Ayz * betayz + l_Azz * betazz) - F2o3 * l_Azz * div_beta;
    
    Axy_rhs[idx] = chin1 * (src_xy - l_gxy * f_trace) + alpn1 * (val_trK * l_Axy - TWO * term_xy) + l_Axx * betaxy + l_Axz * betazy + l_Ayy * betayx + l_Ayz * betazx - l_Axy * betazz + F1o3 * l_Axy * div_beta;
    Ayz_rhs[idx] = chin1 * (src_yz - l_gyz * f_trace) + alpn1 * (val_trK * l_Ayz - TWO * term_yz) + l_Axy * betaxz + l_Ayy * betayz + l_Axz * betaxy + l_Azz * betazy - l_Ayz * betaxx + F1o3 * l_Ayz * div_beta;
    Axz_rhs[idx] = chin1 * (src_xz - l_gxz * f_trace) + alpn1 * (val_trK * l_Axz - TWO * term_xz) + l_Axx * betaxz + l_Axy * betayz + l_Ayz * betayx + l_Azz * betazx - l_Axz * betayy + F1o3 * l_Axz * div_beta;

    trK_rhs[idx] = -chin1 * trK_rhs_val + alpn1 * (F1o3 * val_trK * val_trK + trA2 + FOUR_PI * (rho[idx] + S));

    Lap_rhs[idx] = -TWO * alpn1 * val_trK;

    if constexpr (GAUGE_PARAM == 0) {
        betax_rhs[idx] = FF * dtSfx[idx];
        betay_rhs[idx] = FF * dtSfy[idx];
        betaz_rhs[idx] = FF * dtSfz[idx];
        dtSfx_rhs[idx] = val_Gamx_rhs - eta * dtSfx[idx];
        dtSfy_rhs[idx] = val_Gamy_rhs - eta * dtSfy[idx];
        dtSfz_rhs[idx] = val_Gamz_rhs - eta * dtSfz[idx];
    } else if constexpr (GAUGE_PARAM == 1) {
        betax_rhs[idx] = Gamx[idx] - eta * betax[idx];
        betay_rhs[idx] = Gamy[idx] - eta * betay[idx];
        betaz_rhs[idx] = Gamz[idx] - eta * betaz[idx];
        dtSfx_rhs[idx] = 0.0;
        dtSfy_rhs[idx] = 0.0;
        dtSfz_rhs[idx] = 0.0;
    } else if constexpr (GAUGE_PARAM == 2 || GAUGE_PARAM == 3) {
        betax_rhs[idx] = FF * dtSfx[idx];
        betay_rhs[idx] = FF * dtSfy[idx];
        betaz_rhs[idx] = FF * dtSfz[idx];
        double reta = gupxx * chix*chix + gupyy * chiy*chiy + gupzz * chiz*chiz
                     + TWO * (gupxy * chix*chiy + gupxz * chix*chiz + gupyz * chiy*chiz);
        if constexpr (GAUGE_PARAM == 2) {
            double sq = sqrt(chin1);
            reta = 1.31 / 2.0 * sqrt(reta / chin1) / ((1.0 - sq) * (1.0 - sq));
        } else {
            reta = 1.31 / 2.0 * sqrt(reta / chin1) / ((1.0 - chin1) * (1.0 - chin1));
        }
        dtSfx_rhs[idx] = val_Gamx_rhs - reta * dtSfx[idx];
        dtSfy_rhs[idx] = val_Gamy_rhs - reta * dtSfy[idx];
        dtSfz_rhs[idx] = val_Gamz_rhs - reta * dtSfz[idx];
    } else if constexpr (GAUGE_PARAM == 4 || GAUGE_PARAM == 5) {
        double reta = gupxx * chix*chix + gupyy * chiy*chiy + gupzz * chiz*chiz
                     + TWO * (gupxy * chix*chiy + gupxz * chix*chiz + gupyz * chiy*chiz);
        if constexpr (GAUGE_PARAM == 4) {
            double sq = sqrt(chin1);
            reta = 1.31 / 2.0 * sqrt(reta / chin1) / ((1.0 - sq) * (1.0 - sq));
        } else {
            reta = 1.31 / 2.0 * sqrt(reta / chin1) / ((1.0 - chin1) * (1.0 - chin1));
        }
        betax_rhs[idx] = FF * Gamx[idx] - reta * betax[idx];
        betay_rhs[idx] = FF * Gamy[idx] - reta * betay[idx];
        betaz_rhs[idx] = FF * Gamz[idx] - reta * betaz[idx];
        dtSfx_rhs[idx] = 0.0;
        dtSfy_rhs[idx] = 0.0;
        dtSfz_rhs[idx] = 0.0;
    } else if constexpr (GAUGE_PARAM == 6 || GAUGE_PARAM == 7) {
        betax_rhs[idx] = FF * dtSfx[idx];
        betay_rhs[idx] = FF * dtSfy[idx];
        betaz_rhs[idx] = FF * dtSfz[idx];
        double reta = eta;
        if (BHN == 2) {
            double M = Mass1 + Mass2;
            double A = 2.0 / M;
            double w1 = 12.0, w2 = 12.0;
            double C1 = 1.0 / Mass1 - A;
            double C2 = 1.0 / Mass2 - A;
            double sep2 = (Porg1-Porg4)*(Porg1-Porg4) + (Porg2-Porg5)*(Porg2-Porg5) + (Porg3-Porg6)*(Porg3-Porg6);
            double r1 = ((Porg1-X[i])*(Porg1-X[i]) + (Porg2-Y[j])*(Porg2-Y[j]) + (Porg3-Z[k])*(Porg3-Z[k])) / sep2;
            double r2 = ((Porg4-X[i])*(Porg4-X[i]) + (Porg5-Y[j])*(Porg5-Y[j]) + (Porg6-Z[k])*(Porg6-Z[k])) / sep2;
            if constexpr (GAUGE_PARAM == 6) {
                reta = A + C1 / (1.0 + w1 * r1) + C2 / (1.0 + w2 * r2);
            } else {
                reta = A + C1 * exp(-w1 * r1) + C2 * exp(-w2 * r2);
            }
        }
        dtSfx_rhs[idx] = val_Gamx_rhs - reta * dtSfx[idx];
        dtSfy_rhs[idx] = val_Gamy_rhs - reta * dtSfy[idx];
        dtSfz_rhs[idx] = val_Gamz_rhs - reta * dtSfz[idx];
    }

    Gamx_rhs[idx] = val_Gamx_rhs; Gamy_rhs[idx] = val_Gamy_rhs; Gamz_rhs[idx] = val_Gamz_rhs;
}

// 将 24 个变量完全平铺的指针分配宏 (保持不变)
#define GET_SINGLE_STEP_POINTERS(step_idx) \
    f = nullptr; f_rhs = nullptr; \
    switch (step_idx) { \
        case 0:  f = dxx; f_rhs = gxx_rhs; break; \
        case 1:  f = gxy; f_rhs = gxy_rhs; break; \
        case 2:  f = gxz; f_rhs = gxz_rhs; break; \
        case 3:  f = dyy; f_rhs = gyy_rhs; break; \
        case 4:  f = gyz; f_rhs = gyz_rhs; break; \
        case 5:  f = dzz; f_rhs = gzz_rhs; break; \
        case 6:  f = Axx; f_rhs = Axx_rhs; break; \
        case 7:  f = Axy; f_rhs = Axy_rhs; break; \
        case 8:  f = Axz; f_rhs = Axz_rhs; break; \
        case 9:  f = Ayy; f_rhs = Ayy_rhs; break; \
        case 10: f = Ayz; f_rhs = Ayz_rhs; break; \
        case 11: f = Azz; f_rhs = Azz_rhs; break; \
        case 12: f = chi; f_rhs = chi_rhs; break; \
        case 13: f = trK; f_rhs = trK_rhs; break; \
        case 14: f = Gamx; f_rhs = Gamx_rhs; break; \
        case 15: f = Gamy; f_rhs = Gamy_rhs; break; \
        case 16: f = Gamz; f_rhs = Gamz_rhs; break; \
        case 17: f = Lap;  f_rhs = Lap_rhs;  break; \
        case 18: f = betax; f_rhs = betax_rhs; break; \
        case 19: f = betay; f_rhs = betay_rhs; break; \
        case 20: f = betaz; f_rhs = betaz_rhs; break; \
        case 21: f = dtSfx; f_rhs = dtSfx_rhs; break; \
        case 22: f = dtSfy; f_rhs = dtSfy_rhs; break; \
        case 23: f = dtSfz; f_rhs = dtSfz_rhs; break; \
    }

// K3: Advection + KO Dissipation unified fused path for all GW
template<int GW>
__launch_bounds__(256, 2)
__global__ void bssn_advection_dissipation_kernel(
    const int ex0, const int ex1, const int ex2,
    const int pex0, const int pex1, const int pex2,
    const int symmetry, const int lev, const double eps,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ betax, const double* __restrict__ betay, const double* __restrict__ betaz,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Axx, const double* __restrict__ Axy, const double* __restrict__ Axz,
    const double* __restrict__ Ayy, const double* __restrict__ Ayz, const double* __restrict__ Azz,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    const double* __restrict__ Lap,
    const double* __restrict__ dtSfx, const double* __restrict__ dtSfy, const double* __restrict__ dtSfz,
    double* __restrict__ gxx_rhs, double* __restrict__ gxy_rhs, double* __restrict__ gxz_rhs,
    double* __restrict__ gyy_rhs, double* __restrict__ gyz_rhs, double* __restrict__ gzz_rhs,
    double* __restrict__ Axx_rhs, double* __restrict__ Axy_rhs, double* __restrict__ Axz_rhs,
    double* __restrict__ Ayy_rhs, double* __restrict__ Ayz_rhs, double* __restrict__ Azz_rhs,
    double* __restrict__ chi_rhs, double* __restrict__ trK_rhs,
    double* __restrict__ Gamx_rhs, double* __restrict__ Gamy_rhs, double* __restrict__ Gamz_rhs,
    double* __restrict__ Lap_rhs,
    double* __restrict__ betax_rhs, double* __restrict__ betay_rhs, double* __restrict__ betaz_rhs,
    double* __restrict__ dtSfx_rhs, double* __restrict__ dtSfy_rhs, double* __restrict__ dtSfz_rhs
) {
    // 核心修复 1：根据实际差分逻辑，GW=4和5的高阶迎风最多会越过5个格点
    constexpr int MAX_DIST = (GW <= 3) ? GW : 5;
    constexpr int WSIZE = 2 * MAX_DIST + 1;
    constexpr int WCTR = MAX_DIST;

    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;

    const int imax = ex0 - 1;
    const int jmax = ex1 - 1;
    const int kmax = ex2 - 1;

    if (i > imax || j > jmax || k > kmax) return;

    const int idx = i + pex0 * (j + pex1 * k);

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    int imin = 0, jmin = 0, kmin = 0;
    if (symmetry > 0 && fabs(Z[0]) < dZ) kmin = -GW;
    if (symmetry > 1 && fabs(X[0]) < dX) imin = -GW;
    if (symmetry > 1 && fabs(Y[0]) < dY) jmin = -GW;

    const double vx = __ldg(&betax[idx]);
    const double vy = __ldg(&betay[idx]);
    const double vz = __ldg(&betaz[idx]);

    // 此时 Wx, Wy, Wz 绝对安全，不会再被 GW=4/5 的长 Stencil 越界覆写
    double Wx[WSIZE] = {0};
    double Wy[WSIZE] = {0};
    double Wz[WSIZE] = {0};

    // =========================================================================
    // Stencil selection + weight assembly
    // =========================================================================
    if constexpr (GW == 3) {
        const double d12dx = 1.0 / (12.0 * dX);
        const double d12dy = 1.0 / (12.0 * dY);
        const double d12dz = 1.0 / (12.0 * dZ);

        const double v_dx = vx * d12dx;
        const double v_dy = vy * d12dy;
        const double v_dz = vz * d12dz;

        int sx = 0, sy = 0, sz = 0;
        if (vx > 0.0) { if (i + 3 <= imax) sx = 3; else if (i + 2 <= imax) sx = 2; else if (i + 1 <= imax) sx = 1; }
        else if (vx < 0.0) { if (i - 3 >= imin) sx = -3; else if (i - 2 >= imin) sx = -2; else if (i - 1 >= imin) sx = -1; }
        
        if (vy > 0.0) { if (j + 3 <= jmax) sy = 3; else if (j + 2 <= jmax) sy = 2; else if (j + 1 <= jmax) sy = 1; }
        else if (vy < 0.0) { if (j - 3 >= jmin) sy = -3; else if (j - 2 >= jmin) sy = -2; else if (j - 1 >= jmin) sy = -1; }
        
        if (vz > 0.0) { if (k + 3 <= kmax) sz = 3; else if (k + 2 <= kmax) sz = 2; else if (k + 1 <= kmax) sz = 1; }
        else if (vz < 0.0) { if (k - 3 >= kmin) sz = -3; else if (k - 2 >= kmin) sz = -2; else if (k - 1 >= kmin) sz = -1; }

        if (sx == 3 || sx == -1) {
            Wx[WCTR-1] += v_dx * -3.0; Wx[WCTR] += v_dx * -10.0; Wx[WCTR+1] += v_dx * 18.0; Wx[WCTR+2] += v_dx * -6.0; Wx[WCTR+3] += v_dx * 1.0;
        } else if (sx == 2 || sx == -2) {
            Wx[WCTR-2] += v_dx * 1.0;  Wx[WCTR-1] += v_dx * -8.0;  Wx[WCTR+1] += v_dx * 8.0;  Wx[WCTR+2] += v_dx * -1.0;
        } else if (sx == 1 || sx == -3) {
            Wx[WCTR-3] += v_dx * -1.0; Wx[WCTR-2] += v_dx * 6.0;   Wx[WCTR-1] += v_dx * -18.0; Wx[WCTR] += v_dx * 10.0; Wx[WCTR+1] += v_dx * 3.0;
        }
        if (sy == 3 || sy == -1) {
            Wy[WCTR-1] += v_dy * -3.0; Wy[WCTR] += v_dy * -10.0; Wy[WCTR+1] += v_dy * 18.0; Wy[WCTR+2] += v_dy * -6.0; Wy[WCTR+3] += v_dy * 1.0;
        } else if (sy == 2 || sy == -2) {
            Wy[WCTR-2] += v_dy * 1.0;  Wy[WCTR-1] += v_dy * -8.0;  Wy[WCTR+1] += v_dy * 8.0;  Wy[WCTR+2] += v_dy * -1.0;
        } else if (sy == 1 || sy == -3) {
            Wy[WCTR-3] += v_dy * -1.0; Wy[WCTR-2] += v_dy * 6.0;   Wy[WCTR-1] += v_dy * -18.0; Wy[WCTR] += v_dy * 10.0; Wy[WCTR+1] += v_dy * 3.0;
        }
        if (sz == 3 || sz == -1) {
            Wz[WCTR-1] += v_dz * -3.0; Wz[WCTR] += v_dz * -10.0; Wz[WCTR+1] += v_dz * 18.0; Wz[WCTR+2] += v_dz * -6.0; Wz[WCTR+3] += v_dz * 1.0;
        } else if (sz == 2 || sz == -2) {
            Wz[WCTR-2] += v_dz * 1.0;  Wz[WCTR-1] += v_dz * -8.0;  Wz[WCTR+1] += v_dz * 8.0;  Wz[WCTR+2] += v_dz * -1.0;
        } else if (sz == 1 || sz == -3) {
            Wz[WCTR-3] += v_dz * -1.0; Wz[WCTR-2] += v_dz * 6.0;   Wz[WCTR-1] += v_dz * -18.0; Wz[WCTR] += v_dz * 10.0; Wz[WCTR+1] += v_dz * 3.0;
        }

    } else if constexpr (GW == 2) {
        const double d2dx = 0.5 / dX;
        const double d2dy = 0.5 / dY;
        const double d2dz = 0.5 / dZ;

        const double v_dx = vx * d2dx;
        const double v_dy = vy * d2dy;
        const double v_dz = vz * d2dz;

        if (vx > 0.0 && i + 2 <= imax) {
            Wx[WCTR] += v_dx * (-3.0); Wx[WCTR+1] += v_dx * 4.0; Wx[WCTR+2] += v_dx * (-1.0);
        } else if (vx > 0.0 && i + 1 <= imax) {
            Wx[WCTR] += v_dx * (-1.0); Wx[WCTR+1] += v_dx * 1.0;
        }
        if (vx < 0.0 && i - 2 >= imin) {
            double nv_dx = -v_dx;
            Wx[WCTR] += nv_dx * (-3.0); Wx[WCTR-1] += nv_dx * 4.0; Wx[WCTR-2] += nv_dx * (-1.0);
        } else if (vx < 0.0 && i - 1 >= imin) {
            Wx[WCTR] += (-v_dx) * (-1.0); Wx[WCTR-1] += (-v_dx) * 1.0;
        }
        
        if (vy > 0.0 && j + 2 <= jmax) {
            Wy[WCTR] += v_dy * (-3.0); Wy[WCTR+1] += v_dy * 4.0; Wy[WCTR+2] += v_dy * (-1.0);
        } else if (vy > 0.0 && j + 1 <= jmax) {
            Wy[WCTR] += v_dy * (-1.0); Wy[WCTR+1] += v_dy * 1.0;
        }
        if (vy < 0.0 && j - 2 >= jmin) {
            double nv_dy = -v_dy;
            Wy[WCTR] += nv_dy * (-3.0); Wy[WCTR-1] += nv_dy * 4.0; Wy[WCTR-2] += nv_dy * (-1.0);
        } else if (vy < 0.0 && j - 1 >= jmin) {
            Wy[WCTR] += (-v_dy) * (-1.0); Wy[WCTR-1] += (-v_dy) * 1.0;
        }

        if (vz > 0.0 && k + 2 <= kmax) {
            Wz[WCTR] += v_dz * (-3.0); Wz[WCTR+1] += v_dz * 4.0; Wz[WCTR+2] += v_dz * (-1.0);
        } else if (vz > 0.0 && k + 1 <= kmax) {
            Wz[WCTR] += v_dz * (-1.0); Wz[WCTR+1] += v_dz * 1.0;
        }
        if (vz < 0.0 && k - 2 >= kmin) {
            double nv_dz = -v_dz;
            Wz[WCTR] += nv_dz * (-3.0); Wz[WCTR-1] += nv_dz * 4.0; Wz[WCTR-2] += nv_dz * (-1.0);
        } else if (vz < 0.0 && k - 1 >= kmin) {
            Wz[WCTR] += (-v_dz) * (-1.0); Wz[WCTR-1] += (-v_dz) * 1.0;
        }

    } else if constexpr (GW == 4) {
        const double d60dx = 1.0 / 60.0 / dX; const double d60dy = 1.0 / 60.0 / dY; const double d60dz = 1.0 / 60.0 / dZ;
        const double d12dx_local = 1.0 / 12.0 / dX; const double d12dy_local = 1.0 / 12.0 / dY; const double d12dz_local = 1.0 / 12.0 / dZ;
        const double d2dx_local = 0.5 / dX; const double d2dy_local = 0.5 / dY; const double d2dz_local = 0.5 / dZ;

        if (vx >= 0.0 && i + 4 <= imax && i - 2 >= imin) {
            double c = vx * d60dx;
            Wx[WCTR-2] += c * 2.0; Wx[WCTR-1] += c * -24.0; Wx[WCTR] += c * -35.0;
            Wx[WCTR+1] += c * 80.0; Wx[WCTR+2] += c * -30.0; Wx[WCTR+3] += c * 8.0; Wx[WCTR+4] += c * -1.0;
        } else if (vx >= 0.0 && i + 5 <= imax && i - 1 >= imin) {
            double c = vx * d60dx;
            Wx[WCTR-1] += c * -10.0; Wx[WCTR] += c * -77.0; Wx[WCTR+1] += c * 150.0;
            Wx[WCTR+2] += c * -100.0; Wx[WCTR+3] += c * 50.0; Wx[WCTR+4] += c * -15.0; Wx[WCTR+5] += c * 2.0;
        } else if (vx <= 0.0 && i - 4 >= imin && i + 2 <= imax) {
            double c = -vx * d60dx;
            Wx[WCTR+2] += c * 2.0; Wx[WCTR+1] += c * -24.0; Wx[WCTR] += c * -35.0;
            Wx[WCTR-1] += c * 80.0; Wx[WCTR-2] += c * -30.0; Wx[WCTR-3] += c * 8.0; Wx[WCTR-4] += c * -1.0;
        } else if (vx <= 0.0 && i - 5 >= imin && i + 1 <= imax) {
            double c = -vx * d60dx;
            Wx[WCTR+1] += c * -10.0; Wx[WCTR] += c * -77.0; Wx[WCTR-1] += c * 150.0;
            Wx[WCTR-2] += c * -100.0; Wx[WCTR-3] += c * 50.0; Wx[WCTR-4] += c * -15.0; Wx[WCTR-5] += c * 2.0;
        } else if (i + 3 <= imax && i - 3 >= imin) {
            double c = vx * d60dx;
            Wx[WCTR-3] += c * -1.0; Wx[WCTR-2] += c * 9.0; Wx[WCTR-1] += c * -45.0;
            Wx[WCTR+1] += c * 45.0; Wx[WCTR+2] += c * -9.0; Wx[WCTR+3] += c * 1.0;
        } else if (i + 2 <= imax && i - 2 >= imin) {
            double c = vx * d12dx_local;
            Wx[WCTR-2] += c * 1.0; Wx[WCTR-1] += c * -8.0; Wx[WCTR+1] += c * 8.0; Wx[WCTR+2] += c * -1.0;
        } else if (i + 1 <= imax && i - 1 >= imin) {
            double c = vx * d2dx_local; Wx[WCTR-1] += c * -1.0; Wx[WCTR+1] += c * 1.0;
        }

        if (vy >= 0.0 && j + 4 <= jmax && j - 2 >= jmin) {
            double c = vy * d60dy;
            Wy[WCTR-2] += c * 2.0; Wy[WCTR-1] += c * -24.0; Wy[WCTR] += c * -35.0;
            Wy[WCTR+1] += c * 80.0; Wy[WCTR+2] += c * -30.0; Wy[WCTR+3] += c * 8.0; Wy[WCTR+4] += c * -1.0;
        } else if (vy >= 0.0 && j + 5 <= jmax && j - 1 >= jmin) {
            double c = vy * d60dy;
            Wy[WCTR-1] += c * -10.0; Wy[WCTR] += c * -77.0; Wy[WCTR+1] += c * 150.0;
            Wy[WCTR+2] += c * -100.0; Wy[WCTR+3] += c * 50.0; Wy[WCTR+4] += c * -15.0; Wy[WCTR+5] += c * 2.0;
        } else if (vy <= 0.0 && j - 4 >= jmin && j + 2 <= jmax) {
            double c = -vy * d60dy;
            Wy[WCTR+2] += c * 2.0; Wy[WCTR+1] += c * -24.0; Wy[WCTR] += c * -35.0;
            Wy[WCTR-1] += c * 80.0; Wy[WCTR-2] += c * -30.0; Wy[WCTR-3] += c * 8.0; Wy[WCTR-4] += c * -1.0;
        } else if (vy <= 0.0 && j - 5 >= jmin && j + 1 <= jmax) {
            double c = -vy * d60dy;
            Wy[WCTR+1] += c * -10.0; Wy[WCTR] += c * -77.0; Wy[WCTR-1] += c * 150.0;
            Wy[WCTR-2] += c * -100.0; Wy[WCTR-3] += c * 50.0; Wy[WCTR-4] += c * -15.0; Wy[WCTR-5] += c * 2.0;
        } else if (j + 3 <= jmax && j - 3 >= jmin) {
            double c = vy * d60dy;
            Wy[WCTR-3] += c * -1.0; Wy[WCTR-2] += c * 9.0; Wy[WCTR-1] += c * -45.0;
            Wy[WCTR+1] += c * 45.0; Wy[WCTR+2] += c * -9.0; Wy[WCTR+3] += c * 1.0;
        } else if (j + 2 <= jmax && j - 2 >= jmin) {
            double c = vy * d12dy_local;
            Wy[WCTR-2] += c * 1.0; Wy[WCTR-1] += c * -8.0; Wy[WCTR+1] += c * 8.0; Wy[WCTR+2] += c * -1.0;
        } else if (j + 1 <= jmax && j - 1 >= jmin) {
            double c = vy * d2dy_local; Wy[WCTR-1] += c * -1.0; Wy[WCTR+1] += c * 1.0;
        }

        if (vz >= 0.0 && k + 4 <= kmax && k - 2 >= kmin) {
            double c = vz * d60dz;
            Wz[WCTR-2] += c * 2.0; Wz[WCTR-1] += c * -24.0; Wz[WCTR] += c * -35.0;
            Wz[WCTR+1] += c * 80.0; Wz[WCTR+2] += c * -30.0; Wz[WCTR+3] += c * 8.0; Wz[WCTR+4] += c * -1.0;
        } else if (vz >= 0.0 && k + 5 <= kmax && k - 1 >= kmin) {
            double c = vz * d60dz;
            Wz[WCTR-1] += c * -10.0; Wz[WCTR] += c * -77.0; Wz[WCTR+1] += c * 150.0;
            Wz[WCTR+2] += c * -100.0; Wz[WCTR+3] += c * 50.0; Wz[WCTR+4] += c * -15.0; Wz[WCTR+5] += c * 2.0;
        } else if (vz <= 0.0 && k - 4 >= kmin && k + 2 <= kmax) {
            double c = -vz * d60dz;
            Wz[WCTR+2] += c * 2.0; Wz[WCTR+1] += c * -24.0; Wz[WCTR] += c * -35.0;
            Wz[WCTR-1] += c * 80.0; Wz[WCTR-2] += c * -30.0; Wz[WCTR-3] += c * 8.0; Wz[WCTR-4] += c * -1.0;
        } else if (vz <= 0.0 && k - 5 >= kmin && k + 1 <= kmax) {
            double c = -vz * d60dz;
            Wz[WCTR+1] += c * -10.0; Wz[WCTR] += c * -77.0; Wz[WCTR-1] += c * 150.0;
            Wz[WCTR-2] += c * -100.0; Wz[WCTR-3] += c * 50.0; Wz[WCTR-4] += c * -15.0; Wz[WCTR-5] += c * 2.0;
        } else if (k + 3 <= kmax && k - 3 >= kmin) {
            double c = vz * d60dz;
            Wz[WCTR-3] += c * -1.0; Wz[WCTR-2] += c * 9.0; Wz[WCTR-1] += c * -45.0;
            Wz[WCTR+1] += c * 45.0; Wz[WCTR+2] += c * -9.0; Wz[WCTR+3] += c * 1.0;
        } else if (k + 2 <= kmax && k - 2 >= kmin) {
            double c = vz * d12dz_local;
            Wz[WCTR-2] += c * 1.0; Wz[WCTR-1] += c * -8.0; Wz[WCTR+1] += c * 8.0; Wz[WCTR+2] += c * -1.0;
        } else if (k + 1 <= kmax && k - 1 >= kmin) {
            double c = vz * d2dz_local; Wz[WCTR-1] += c * -1.0; Wz[WCTR+1] += c * 1.0;
        }

    } else if constexpr (GW == 5) {
        const double d840dx = 1.0 / 840.0 / dX; const double d840dy = 1.0 / 840.0 / dY; const double d840dz = 1.0 / 840.0 / dZ;
        const double d60dx_local = 1.0 / 60.0 / dX; const double d60dy_local = 1.0 / 60.0 / dY; const double d60dz_local = 1.0 / 60.0 / dZ;
        const double d12dx_local = 1.0 / 12.0 / dX; const double d12dy_local = 1.0 / 12.0 / dY; const double d12dz_local = 1.0 / 12.0 / dZ;
        const double d2dx_local = 0.5 / dX; const double d2dy_local = 0.5 / dY; const double d2dz_local = 0.5 / dZ;

        if (vx >= 0.0 && i + 5 <= imax && i - 3 >= imin) {
            double c = vx * d840dx;
            Wx[WCTR-3] += c * -5.0; Wx[WCTR-2] += c * 60.0; Wx[WCTR-1] += c * -420.0; Wx[WCTR] += c * -378.0;
            Wx[WCTR+1] += c * 1050.0; Wx[WCTR+2] += c * -420.0; Wx[WCTR+3] += c * 140.0; Wx[WCTR+4] += c * -30.0; Wx[WCTR+5] += c * 3.0;
        } else if (vx <= 0.0 && i - 5 >= imin && i + 3 <= imax) {
            double c = -vx * d840dx;
            Wx[WCTR+3] += c * -5.0; Wx[WCTR+2] += c * 60.0; Wx[WCTR+1] += c * -420.0; Wx[WCTR] += c * -378.0;
            Wx[WCTR-1] += c * 1050.0; Wx[WCTR-2] += c * -420.0; Wx[WCTR-3] += c * 140.0; Wx[WCTR-4] += c * -30.0; Wx[WCTR-5] += c * 3.0;
        } else if (i + 4 <= imax && i - 4 >= imin) {
            double c = vx * d840dx;
            Wx[WCTR-4] += c * 3.0; Wx[WCTR-3] += c * -32.0; Wx[WCTR-2] += c * 168.0; Wx[WCTR-1] += c * -672.0;
            Wx[WCTR+1] += c * 672.0; Wx[WCTR+2] += c * -168.0; Wx[WCTR+3] += c * 32.0; Wx[WCTR+4] += c * -3.0;
        } else if (i + 3 <= imax && i - 3 >= imin) {
            double c = vx * d60dx_local;
            Wx[WCTR-3] += c * -1.0; Wx[WCTR-2] += c * 9.0; Wx[WCTR-1] += c * -45.0; Wx[WCTR+1] += c * 45.0; Wx[WCTR+2] += c * -9.0; Wx[WCTR+3] += c * 1.0;
        } else if (i + 2 <= imax && i - 2 >= imin) {
            double c = vx * d12dx_local; Wx[WCTR-2] += c * 1.0; Wx[WCTR-1] += c * -8.0; Wx[WCTR+1] += c * 8.0; Wx[WCTR+2] += c * -1.0;
        } else if (i + 1 <= imax && i - 1 >= imin) {
            double c = vx * d2dx_local; Wx[WCTR-1] += c * -1.0; Wx[WCTR+1] += c * 1.0;
        }

        if (vy >= 0.0 && j + 5 <= jmax && j - 3 >= jmin) {
            double c = vy * d840dy;
            Wy[WCTR-3] += c * -5.0; Wy[WCTR-2] += c * 60.0; Wy[WCTR-1] += c * -420.0; Wy[WCTR] += c * -378.0;
            Wy[WCTR+1] += c * 1050.0; Wy[WCTR+2] += c * -420.0; Wy[WCTR+3] += c * 140.0; Wy[WCTR+4] += c * -30.0; Wy[WCTR+5] += c * 3.0;
        } else if (vy <= 0.0 && j - 5 >= jmin && j + 3 <= jmax) {
            double c = -vy * d840dy;
            Wy[WCTR+3] += c * -5.0; Wy[WCTR+2] += c * 60.0; Wy[WCTR+1] += c * -420.0; Wy[WCTR] += c * -378.0;
            Wy[WCTR-1] += c * 1050.0; Wy[WCTR-2] += c * -420.0; Wy[WCTR-3] += c * 140.0; Wy[WCTR-4] += c * -30.0; Wy[WCTR-5] += c * 3.0;
        } else if (j + 4 <= jmax && j - 4 >= jmin) {
            double c = vy * d840dy;
            Wy[WCTR-4] += c * 3.0; Wy[WCTR-3] += c * -32.0; Wy[WCTR-2] += c * 168.0; Wy[WCTR-1] += c * -672.0;
            Wy[WCTR+1] += c * 672.0; Wy[WCTR+2] += c * -168.0; Wy[WCTR+3] += c * 32.0; Wy[WCTR+4] += c * -3.0;
        } else if (j + 3 <= jmax && j - 3 >= jmin) {
            double c = vy * d60dy_local;
            Wy[WCTR-3] += c * -1.0; Wy[WCTR-2] += c * 9.0; Wy[WCTR-1] += c * -45.0; Wy[WCTR+1] += c * 45.0; Wy[WCTR+2] += c * -9.0; Wy[WCTR+3] += c * 1.0;
        } else if (j + 2 <= jmax && j - 2 >= jmin) {
            double c = vy * d12dy_local; Wy[WCTR-2] += c * 1.0; Wy[WCTR-1] += c * -8.0; Wy[WCTR+1] += c * 8.0; Wy[WCTR+2] += c * -1.0;
        } else if (j + 1 <= jmax && j - 1 >= jmin) {
            double c = vy * d2dy_local; Wy[WCTR-1] += c * -1.0; Wy[WCTR+1] += c * 1.0;
        }

        if (vz >= 0.0 && k + 5 <= kmax && k - 3 >= kmin) {
            double c = vz * d840dz;
            Wz[WCTR-3] += c * -5.0; Wz[WCTR-2] += c * 60.0; Wz[WCTR-1] += c * -420.0; Wz[WCTR] += c * -378.0;
            Wz[WCTR+1] += c * 1050.0; Wz[WCTR+2] += c * -420.0; Wz[WCTR+3] += c * 140.0; Wz[WCTR+4] += c * -30.0; Wz[WCTR+5] += c * 3.0;
        } else if (vz <= 0.0 && k - 5 >= kmin && k + 3 <= kmax) {
            double c = -vz * d840dz;
            Wz[WCTR+3] += c * -5.0; Wz[WCTR+2] += c * 60.0; Wz[WCTR+1] += c * -420.0; Wz[WCTR] += c * -378.0;
            Wz[WCTR-1] += c * 1050.0; Wz[WCTR-2] += c * -420.0; Wz[WCTR-3] += c * 140.0; Wz[WCTR-4] += c * -30.0; Wz[WCTR-5] += c * 3.0;
        } else if (k + 4 <= kmax && k - 4 >= kmin) {
            double c = vz * d840dz;
            Wz[WCTR-4] += c * 3.0; Wz[WCTR-3] += c * -32.0; Wz[WCTR-2] += c * 168.0; Wz[WCTR-1] += c * -672.0;
            Wz[WCTR+1] += c * 672.0; Wz[WCTR+2] += c * -168.0; Wz[WCTR+3] += c * 32.0; Wz[WCTR+4] += c * -3.0;
        } else if (k + 3 <= kmax && k - 3 >= kmin) {
            double c = vz * d60dz_local;
            Wz[WCTR-3] += c * -1.0; Wz[WCTR-2] += c * 9.0; Wz[WCTR-1] += c * -45.0; Wz[WCTR+1] += c * 45.0; Wz[WCTR+2] += c * -9.0; Wz[WCTR+3] += c * 1.0;
        } else if (k + 2 <= kmax && k - 2 >= kmin) {
            double c = vz * d12dz_local; Wz[WCTR-2] += c * 1.0; Wz[WCTR-1] += c * -8.0; Wz[WCTR+1] += c * 8.0; Wz[WCTR+2] += c * -1.0;
        } else if (k + 1 <= kmax && k - 1 >= kmin) {
            double c = vz * d2dz_local; Wz[WCTR-1] += c * -1.0; Wz[WCTR+1] += c * 1.0;
        }
    }

    // ---- KO dissipation fusion ----
    const bool do_dissipation = (eps > 0.0 &&
        i - GW >= imin && i + GW <= imax &&
        j - GW >= jmin && j + GW <= jmax &&
        k - GW >= kmin && k + GW <= kmax);

    if (do_dissipation) {
        constexpr double cof = double(1 << (2 * GW));
        constexpr double sign = ((GW - 1) % 2 == 0) ? 1.0 : -1.0;
        const double edx = sign * eps / (cof * dX);
        const double edy = sign * eps / (cof * dY);
        const double edz = sign * eps / (cof * dZ);

        if constexpr (GW == 2) {
            constexpr double c[5] = {1, -4, 6, -4, 1};
            #pragma unroll
            for (int d = 0; d < 5; d++) { Wx[WCTR-2+d] += edx * c[d]; Wy[WCTR-2+d] += edy * c[d]; Wz[WCTR-2+d] += edz * c[d]; }
        } else if constexpr (GW == 3) {
            constexpr double c[7] = {1, -6, 15, -20, 15, -6, 1};
            #pragma unroll
            for (int d = 0; d < 7; d++) { Wx[WCTR-3+d] += edx * c[d]; Wy[WCTR-3+d] += edy * c[d]; Wz[WCTR-3+d] += edz * c[d]; }
        } else if constexpr (GW == 4) {
            constexpr double c[9] = {1, -8, 28, -56, 70, -56, 28, -8, 1};
            #pragma unroll
            for (int d = 0; d < 9; d++) { Wx[WCTR-4+d] += edx * c[d]; Wy[WCTR-4+d] += edy * c[d]; Wz[WCTR-4+d] += edz * c[d]; }
        } else if constexpr (GW == 5) {
            constexpr double c[11] = {1, -10, 45, -120, 210, -252, 210, -120, 45, -10, 1};
            #pragma unroll
            for (int d = 0; d < 11; d++) { Wx[WCTR-5+d] += edx * c[d]; Wy[WCTR-5+d] += edy * c[d]; Wz[WCTR-5+d] += edz * c[d]; }
        }
    }

    // =========================================================================
    // 24-field batched __ldg loop
    // =========================================================================
    const int dx = 1;
    const int dy = pex0;
    const int dz = pex0 * pex1;

    const double center_weight = Wx[WCTR] + Wy[WCTR] + Wz[WCTR];

    const double *f;
    double *f_rhs;

    // Stencil read loop. The if-elseif chain above only sets a Wx/Wy/Wz
    // weight slot when the corresponding stencil branch's boundary check
    // passes, so weights for stencils that would reach beyond the data
    // are 0. The loop reads f[idx ± d*stride] unconditionally for
    // d = 1..MAX_DIST; safe iff every (idx ± d*stride) is mapped.
    //
    // x/y negative reads stay inside d_fgfs_base because their stride is
    // small (1 and pex0); the multi-dim flat-index offset to interior
    // origin is large enough to absorb the shift.
    //
    // z negative is different: stride pex0*pex1 is huge, and for GW=4
    // (this codebase's "6th-order") MAX_DIST=5 > GW so d=5 at k=0 would
    // read pex0*pex1 ≈ 6000 doubles *before* d_fgfs_base — an unmapped
    // page, kernel faults. Clamp d_zn = min(d, k + GW) only for that
    // case; the corresponding Wz[WCTR-d] is 0 (the +5/-5 upwind branches
    // require k - d >= kmin = -GW, i.e., d ≤ k + GW), so the numerical
    // result is unchanged. For GW != 4 the if-constexpr drops out so
    // there is zero runtime cost.
    #pragma unroll 1
    for (int step = 0; step < 24; ++step) {
#if (GAUGE == 1 || GAUGE == 4 || GAUGE == 5)
        if (step >= 21) continue;
#endif
        GET_SINGLE_STEP_POINTERS(step);

        double center_val = __ldg(&f[idx]);
        double res = center_weight * center_val;

        #pragma unroll
        for (int d = 1; d <= MAX_DIST; ++d) {
            res += Wx[WCTR-d] * __ldg(&f[idx - d * dx]) + Wx[WCTR+d] * __ldg(&f[idx + d * dx]);
            res += Wy[WCTR-d] * __ldg(&f[idx - d * dy]) + Wy[WCTR+d] * __ldg(&f[idx + d * dy]);
            int d_zn = d;
            if constexpr (MAX_DIST > GW) {
                d_zn = min(d, k + GW);
            }
            res += Wz[WCTR-d] * __ldg(&f[idx - d_zn * dz]) + Wz[WCTR+d] * __ldg(&f[idx + d * dz]);
        }

        f_rhs[idx] += res;
    }
}


template<int GW>
__launch_bounds__(256, 2)
__global__ void bssn_constraints_kernel(
    const int ex0, const int ex1, const int ex2,
    const int pex0, const int pex1, const int pex2,
    const int symmetry, const int lev,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Axx, const double* __restrict__ Axy, const double* __restrict__ Axz,
    const double* __restrict__ Ayy, const double* __restrict__ Ayz, const double* __restrict__ Azz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    const double* __restrict__ rho, 
    const double* __restrict__ Sx, const double* __restrict__ Sy, const double* __restrict__ Sz,
    const double* __restrict__ Gamxxx_in, const double* __restrict__ Gamxxy_in, const double* __restrict__ Gamxxz_in,
    const double* __restrict__ Gamxyy_in, const double* __restrict__ Gamxyz_in, const double* __restrict__ Gamxzz_in,
    const double* __restrict__ Gamyxx_in, const double* __restrict__ Gamyxy_in, const double* __restrict__ Gamyxz_in,
    const double* __restrict__ Gamyyy_in, const double* __restrict__ Gamyyz_in, const double* __restrict__ Gamyzz_in,
    const double* __restrict__ Gamzxx_in, const double* __restrict__ Gamzxy_in, const double* __restrict__ Gamzxz_in,
    const double* __restrict__ Gamzyy_in, const double* __restrict__ Gamzyz_in, const double* __restrict__ Gamzzz_in,
    const double* __restrict__ Rxx_in, const double* __restrict__ Rxy_in, const double* __restrict__ Rxz_in,
    const double* __restrict__ Ryy_in, const double* __restrict__ Ryz_in, const double* __restrict__ Rzz_in,
    double* __restrict__ ham_Res, 
    double* __restrict__ movx_Res, double* __restrict__ movy_Res, double* __restrict__ movz_Res,
    double* __restrict__ Gmx_Res, double* __restrict__ Gmy_Res, double* __restrict__ Gmz_Res
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int pdims[3] = {pex0, pex1, pex2};
    double dX = X[1] - X[0]; double dY = Y[1] - Y[0]; double dZ = Z[1] - Z[0];

    int imax = ex0 - 1; int jmax = ex1 - 1; int kmax = ex2 - 1;
    int imin = 0, jmin = 0, kmin = 0;
    constexpr int ORD = GW - 1;
    if (symmetry > 0 && fabs(Z[0]) < dZ) kmin = -ORD;
    if (symmetry > 1 && fabs(X[0]) < dX) imin = -ORD;
    if (symmetry > 1 && fabs(Y[0]) < dY) jmin = -ORD;

    int block_i = blockIdx.x * blockDim.x;
    int block_j = blockIdx.y * blockDim.y;
    int block_k = blockIdx.z * blockDim.z;

    // Dynamic shared memory: 3 tiles of [Z][Y][X] doubles.
    // For GW=5 this is 3*3072*8 = 73728 bytes, exceeding the 48KB static smem limit on sm_80.
    extern __shared__ double smem_raw[];
    double (*smem3)[TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X] =
        reinterpret_cast<double(*)[TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X]>(smem_raw);

    // =========================================================================
    // 1. Batch 1: Axx, Axy, Axz
    // =========================================================================
    load_3_fields_to_smem<GW>(smem3, Axx, Axy, Axz, pdims, block_i, block_j, block_k);
    __syncthreads();

    double d_Axx_x, d_Axx_y, d_Axx_z;
    compute_first_derivs_smem<GW>(smem3[0], &d_Axx_x, &d_Axx_y, &d_Axx_z, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double d_Axy_x, d_Axy_y, d_Axy_z;
    compute_first_derivs_smem<GW>(smem3[1], &d_Axy_x, &d_Axy_y, &d_Axy_z, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double d_Axz_x, d_Axz_y, d_Axz_z;
    compute_first_derivs_smem<GW>(smem3[2], &d_Axz_x, &d_Axz_y, &d_Axz_z, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    
    __syncthreads(); // 等待所有线程算完 Batch 1

    // =========================================================================
    // 2. Batch 2: Ayy, Ayz, Azz
    // =========================================================================
    load_3_fields_to_smem<GW>(smem3, Ayy, Ayz, Azz, pdims, block_i, block_j, block_k);
    __syncthreads();

    double d_Ayy_x, d_Ayy_y, d_Ayy_z;
    compute_first_derivs_smem<GW>(smem3[0], &d_Ayy_x, &d_Ayy_y, &d_Ayy_z, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double d_Ayz_x, d_Ayz_y, d_Ayz_z;
    compute_first_derivs_smem<GW>(smem3[1], &d_Ayz_x, &d_Ayz_y, &d_Ayz_z, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double d_Azz_x, d_Azz_y, d_Azz_z;
    compute_first_derivs_smem<GW>(smem3[2], &d_Azz_x, &d_Azz_y, &d_Azz_z, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    __syncthreads();

    // =========================================================================
    // 3. Batch 3: chi, trK (借用 load_3_fields，第三个参数传 trK 占位)
    // =========================================================================
    load_3_fields_to_smem<GW>(smem3, chi, trK, trK, pdims, block_i, block_j, block_k);
    __syncthreads();

    double chix, chiy, chiz, chixx, chixy, chixz, chiyy, chiyz, chizz;
    compute_all_derivs_smem<GW>(smem3[0], &chix, &chiy, &chiz, &chixx, &chixy, &chixz, &chiyy, &chiyz, &chizz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double Kx, Ky, Kz;
    compute_first_derivs_smem<GW>(smem3[1], &Kx, &Ky, &Kz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    // 此处无需 __syncthreads()，数据提取完毕

    // =========================================================================
    // 4. 数学降维与物理组装 (完全保留原有数学逻辑)
    // =========================================================================
    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = IDX3D_PAD(i, j, k, pex0, pex1);

    const double HALF = 0.5; const double TWO = 2.0; const double F2o3 = 2.0 / 3.0; 
    const double EIGHT = 8.0; const double F16 = 16.0; const double PI = M_PI;

    double chin1 = chi[idx] + 1.0; double val_trK = trK[idx];
    double inv_chin1 = 1.0 / chin1;

    double l_gxx = dxx[idx] + 1.0; double l_gxy = gxy[idx]; double l_gxz = gxz[idx];
    double l_gyy = dyy[idx] + 1.0; double l_gyz = gyz[idx]; double l_gzz = dzz[idx] + 1.0;
    double l_Axx = Axx[idx]; double l_Axy = Axy[idx]; double l_Axz = Axz[idx];
    double l_Ayy = Ayy[idx]; double l_Ayz = Ayz[idx]; double l_Azz = Azz[idx];

    double det = l_gxx * l_gyy * l_gzz + TWO * l_gxy * l_gyz * l_gxz - l_gxz * l_gyy * l_gxz - l_gxx * l_gyz * l_gyz - l_gzz * l_gxy * l_gxy;
    double inv_det = 1.0 / det;

    double gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) * inv_det;
    double gupxy = (l_gyz * l_gxz - l_gxy * l_gzz) * inv_det;
    double gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) * inv_det;
    double gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) * inv_det;
    double gupyz = (l_gxy * l_gxz - l_gxx * l_gyz) * inv_det;
    double gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) * inv_det;

    const double l_Gamxxx = Gamxxx_in[idx]; const double l_Gamxxy = Gamxxy_in[idx]; const double l_Gamxxz = Gamxxz_in[idx];
    const double l_Gamxyy = Gamxyy_in[idx]; const double l_Gamxyz = Gamxyz_in[idx]; const double l_Gamxzz = Gamxzz_in[idx];
    const double l_Gamyxx = Gamyxx_in[idx]; const double l_Gamyxy = Gamyxy_in[idx]; const double l_Gamyxz = Gamyxz_in[idx];
    const double l_Gamyyy = Gamyyy_in[idx]; const double l_Gamyyz = Gamyyz_in[idx]; const double l_Gamyzz = Gamyzz_in[idx];
    const double l_Gamzxx = Gamzxx_in[idx]; const double l_Gamzxy = Gamzxy_in[idx]; const double l_Gamzxz = Gamzxz_in[idx];
    const double l_Gamzyy = Gamzyy_in[idx]; const double l_Gamzyz = Gamzyz_in[idx]; const double l_Gamzzz = Gamzzz_in[idx];

    double Gamxa = gupxx * l_Gamxxx + gupyy * l_Gamxyy + gupzz * l_Gamxzz + TWO * (gupxy * l_Gamxxy + gupxz * l_Gamxxz + gupyz * l_Gamxyz);
    double Gamya = gupxx * l_Gamyxx + gupyy * l_Gamyyy + gupzz * l_Gamyzz + TWO * (gupxy * l_Gamyxy + gupxz * l_Gamyxz + gupyz * l_Gamyyz);
    double Gamza = gupxx * l_Gamzxx + gupyy * l_Gamzyy + gupzz * l_Gamzzz + TWO * (gupxy * l_Gamzxy + gupxz * l_Gamzxz + gupyz * l_Gamzyz);

    Gmx_Res[idx] = Gamx[idx] - Gamxa; Gmy_Res[idx] = Gamy[idx] - Gamya; Gmz_Res[idx] = Gamz[idx] - Gamza;

    double Au_d_xx = gupxx * l_Axx + gupxy * l_Axy + gupxz * l_Axz;
    double Au_d_xy = gupxx * l_Axy + gupxy * l_Ayy + gupxz * l_Ayz;
    double Au_d_xz = gupxx * l_Axz + gupxy * l_Ayz + gupxz * l_Azz;
    double Au_d_yx = gupxy * l_Axx + gupyy * l_Axy + gupyz * l_Axz;
    double Au_d_yy = gupxy * l_Axy + gupyy * l_Ayy + gupyz * l_Ayz;
    double Au_d_yz = gupxy * l_Axz + gupyy * l_Ayz + gupyz * l_Azz;
    double Au_d_zx = gupxz * l_Axx + gupyz * l_Axy + gupzz * l_Axz;
    double Au_d_zy = gupxz * l_Axy + gupyz * l_Ayy + gupzz * l_Ayz;
    double Au_d_zz = gupxz * l_Axz + gupyz * l_Ayz + gupzz * l_Azz;

    double trA2 = Au_d_xx*Au_d_xx + Au_d_yy*Au_d_yy + Au_d_zz*Au_d_zz + TWO*(Au_d_xy*Au_d_yx + Au_d_xz*Au_d_zx + Au_d_yz*Au_d_zy);

    double conf_R = gupxx * Rxx_in[idx] + gupyy * Ryy_in[idx] + gupzz * Rzz_in[idx] + TWO * (gupxy * Rxy_in[idx] + gupxz * Rxz_in[idx] + gupyz * Ryz_in[idx]);
    
    double fxx = chixx - l_Gamxxx * chix - l_Gamyxx * chiy - l_Gamzxx * chiz;
    double fyy = chiyy - l_Gamxyy * chix - l_Gamyyy * chiy - l_Gamzyy * chiz;
    double fzz = chizz - l_Gamxzz * chix - l_Gamyzz * chiy - l_Gamzzz * chiz;
    double fxy = chixy - l_Gamxxy * chix - l_Gamyxy * chiy - l_Gamzxy * chiz;
    double fxz = chixz - l_Gamxxz * chix - l_Gamyxz * chiy - l_Gamzxz * chiz;
    double fyz = chiyz - l_Gamxyz * chix - l_Gamyyz * chiy - l_Gamzyz * chiz;

    double D2_chi = gupxx * fxx + gupyy * fyy + gupzz * fzz + TWO * (gupxy * fxy + gupxz * fxz + gupyz * fyz);
    double grad_chi_sq = gupxx * chix*chix + gupyy * chiy*chiy + gupzz * chiz*chiz + TWO * (gupxy * chix*chiy + gupxz * chix*chiz + gupyz * chiy*chiz);
    
    double phys_R = chin1 * conf_R + TWO * D2_chi - 2.5 * inv_chin1 * grad_chi_sq;
    ham_Res[idx] = phys_R + F2o3 * val_trK * val_trK - trA2 - F16 * PI * rho[idx];

    double div_A_x = gupxx * d_Axx_x + gupyy * d_Axy_y + gupzz * d_Axz_z + gupxy * (d_Axx_y + d_Axy_x) + gupxz * (d_Axx_z + d_Axz_x) + gupyz * (d_Axy_z + d_Axz_y);
    double div_A_y = gupxx * d_Axy_x + gupyy * d_Ayy_y + gupzz * d_Ayz_z + gupxy * (d_Axy_y + d_Ayy_x) + gupxz * (d_Axy_z + d_Ayz_x) + gupyz * (d_Ayy_z + d_Ayz_y);
    double div_A_z = gupxx * d_Axz_x + gupyy * d_Ayz_y + gupzz * d_Azz_z + gupxy * (d_Axz_y + d_Ayz_x) + gupxz * (d_Axz_z + d_Azz_x) + gupyz * (d_Ayz_z + d_Azz_y);

    double Gam_A_x = Gamxa * l_Axx + Gamya * l_Axy + Gamza * l_Axz;
    double Gam_A_y = Gamxa * l_Axy + Gamya * l_Ayy + Gamza * l_Ayz;
    double Gam_A_z = Gamxa * l_Axz + Gamya * l_Ayz + Gamza * l_Azz;

    double A_Gam_x = Au_d_xx * l_Gamxxx + Au_d_yx * l_Gamxxy + Au_d_zx * l_Gamxxz + Au_d_xy * l_Gamyxx + Au_d_yy * l_Gamyxy + Au_d_zy * l_Gamyxz + Au_d_xz * l_Gamzxx + Au_d_yz * l_Gamzxy + Au_d_zz * l_Gamzxz;
    double A_Gam_y = Au_d_xx * l_Gamxxy + Au_d_yx * l_Gamxyy + Au_d_zx * l_Gamxyz + Au_d_xy * l_Gamyxy + Au_d_yy * l_Gamyyy + Au_d_zy * l_Gamyyz + Au_d_xz * l_Gamzxy + Au_d_yz * l_Gamzyy + Au_d_zz * l_Gamzyz;
    double A_Gam_z = Au_d_xx * l_Gamxxz + Au_d_yx * l_Gamxyz + Au_d_zx * l_Gamxzz + Au_d_xy * l_Gamyxz + Au_d_yy * l_Gamyyz + Au_d_zy * l_Gamyzz + Au_d_xz * l_Gamzxz + Au_d_yz * l_Gamzyz + Au_d_zz * l_Gamzzz;

    double chi_A_x = 1.5 * inv_chin1 * (Au_d_xx * chix + Au_d_yx * chiy + Au_d_zx * chiz);
    double chi_A_y = 1.5 * inv_chin1 * (Au_d_xy * chix + Au_d_yy * chiy + Au_d_zy * chiz);
    double chi_A_z = 1.5 * inv_chin1 * (Au_d_xz * chix + Au_d_yz * chiy + Au_d_zz * chiz);

    movx_Res[idx] = div_A_x - Gam_A_x - A_Gam_x - chi_A_x - F2o3 * Kx - EIGHT * PI * Sx[idx];
    movy_Res[idx] = div_A_y - Gam_A_y - A_Gam_y - chi_A_y - F2o3 * Ky - EIGHT * PI * Sy[idx];
    movz_Res[idx] = div_A_z - Gam_A_z - A_Gam_z - chi_A_z - F2o3 * Kz - EIGHT * PI * Sz[idx];
}

void gpu_compute_rhs_bssn_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
) {
    // Legacy Block: 给需要用到 SMEM 的 Ricci 和 Eval (依然是 8, 8, 4)
    dim3 block_legacy(DIMX, DIMY, DIMZ);
    dim3 grid_legacy(
        (ex[0] + block_legacy.x - 1) / block_legacy.x,
        (ex[1] + block_legacy.y - 1) / block_legacy.y,
        (ex[2] + block_legacy.z - 1) / block_legacy.z
    );

    // Dynamic shared memory size for K1, K2, K4: 3 tiles of doubles
    constexpr int smem_size = 3 * TileParams<ghost_width>::SIZE * sizeof(double);
    static bool smem_attr_set = false;
    if (!smem_attr_set) {
        cudaFuncSetAttribute(bssn_ricci_kernel<ghost_width>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        cudaFuncSetAttribute(bssn_rhs_eval_kernel<ghost_width, GAUGE>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        cudaFuncSetAttribute(bssn_constraints_kernel<ghost_width>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        smem_attr_set = true;
    }

    // 1. 生成共形几何
    bssn_ricci_kernel<ghost_width><<<grid_legacy, block_legacy, smem_size, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, X, Y, Z,
        dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz
    );
    { cudaError_t e=cudaGetLastError(); if(e!=cudaSuccess){fprintf(stderr,"LAUNCH FAIL [ricci]: %s (%s:%d)\n",cudaGetErrorString(e),__FILE__,__LINE__);MPI_Abort(MPI_COMM_WORLD,1);}}

    // 2. 物理演化引擎
    bssn_rhs_eval_kernel<ghost_width, GAUGE><<<grid_legacy, block_legacy, smem_size, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev,
        BHN, PBH[0], PBH[1], PBH[2], PBH[3], PBH[4], PBH[5], Mass[0], Mass[1],
        X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap, betax, betay, betaz, dtSfx, dtSfy, dtSfz,
        rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        chi_rhs, trK_rhs, gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs, betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs
    );
    { cudaError_t e=cudaGetLastError(); if(e!=cudaSuccess){fprintf(stderr,"LAUNCH FAIL [eval]: %s (%s:%d)\n",cudaGetErrorString(e),__FILE__,__LINE__);MPI_Abort(MPI_COMM_WORLD,1);}}

    // 3. 耗散与平流 (纯 L1 Cache 版本)
    // 强制 Warp 沿 X 轴对齐，彻底解决 75% 的 GM Uncoalesced
    dim3 block_adv(32, 8, 1); 
    dim3 grid_adv(
        (ex[0] + block_adv.x - 1) / block_adv.x,
        (ex[1] + block_adv.y - 1) / block_adv.y,
        (ex[2] + block_adv.z - 1) / block_adv.z
    );

    bssn_advection_dissipation_kernel<ghost_width><<<grid_adv, block_adv, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, eps, X, Y, Z,
        betax, betay, betaz, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, chi, trK, Gamx, Gamy, Gamz, Lap, dtSfx, dtSfy, dtSfz,
        gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        chi_rhs, trK_rhs, Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs, dtSfx_rhs, dtSfy_rhs, dtSfz_rhs
    );
    { cudaError_t e=cudaGetLastError(); if(e!=cudaSuccess){fprintf(stderr,"LAUNCH FAIL [advection]: %s (%s:%d)\n",cudaGetErrorString(e),__FILE__,__LINE__);MPI_Abort(MPI_COMM_WORLD,1);}}

    // 4. 诊断系统
    if (co == 0) {
        bssn_constraints_kernel<ghost_width><<<grid_legacy, block_legacy, smem_size, stream>>>(
            ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, X, Y, Z,
            chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
            Axx, Axy, Axz, Ayy, Ayz, Azz, Gamx, Gamy, Gamz, rho, Sx, Sy, Sz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
            ham_Res, movx_Res, movy_Res, movz_Res,
            Gmx_Res, Gmy_Res, Gmz_Res
        );
        { cudaError_t e=cudaGetLastError(); if(e!=cudaSuccess){fprintf(stderr,"LAUNCH FAIL [constraints]: %s (%s:%d)\n",cudaGetErrorString(e),__FILE__,__LINE__);MPI_Abort(MPI_COMM_WORLD,1);}}
    }
}

// K1-only: launch bssn_ricci_kernel, synchronize, then return.
// Used by tests that need to fill intermediate halos between K1 and K2.
void gpu_ricci_only_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
) {
    dim3 block(DIMX, DIMY, DIMZ);
    dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );
    constexpr int smem_size = 3 * TileParams<ghost_width>::SIZE * sizeof(double);
    static bool smem_attr_set = false;
    if (!smem_attr_set) {
        cudaFuncSetAttribute(bssn_ricci_kernel<ghost_width>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        smem_attr_set = true;
    }

    bssn_ricci_kernel<ghost_width><<<grid, block, smem_size, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, X, Y, Z,
        dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz
    );
    cudaStreamSynchronize(stream);
}

// K2+K3+K4: launch the rest of the pipeline after K1.
// Used by tests that fill intermediate halos between K1 and K2.
void gpu_rhs_rest_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
) {
    dim3 block(DIMX, DIMY, DIMZ);
    dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );
    constexpr int smem_size = 3 * TileParams<ghost_width>::SIZE * sizeof(double);
    static bool smem_attr_set = false;
    if (!smem_attr_set) {
        cudaFuncSetAttribute(bssn_rhs_eval_kernel<ghost_width, GAUGE>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        cudaFuncSetAttribute(bssn_constraints_kernel<ghost_width>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        smem_attr_set = true;
    }

    // K2: RHS eval
    bssn_rhs_eval_kernel<ghost_width, GAUGE><<<grid, block, smem_size, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev,
        BHN, PBH[0], PBH[1], PBH[2], PBH[3], PBH[4], PBH[5], Mass[0], Mass[1],
        X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap, betax, betay, betaz, dtSfx, dtSfy, dtSfz,
        rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        chi_rhs, trK_rhs, gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs, betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs
    );

    // K3: advection + dissipation
    dim3 block_adv(32, 8, 1);
    dim3 grid_adv(
        (ex[0] + block_adv.x - 1) / block_adv.x,
        (ex[1] + block_adv.y - 1) / block_adv.y,
        (ex[2] + block_adv.z - 1) / block_adv.z
    );
    bssn_advection_dissipation_kernel<ghost_width><<<grid_adv, block_adv, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, eps, X, Y, Z,
        betax, betay, betaz, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, chi, trK, Gamx, Gamy, Gamz, Lap, dtSfx, dtSfy, dtSfz,
        gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        chi_rhs, trK_rhs, Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs, dtSfx_rhs, dtSfy_rhs, dtSfz_rhs
    );

    // K4: constraints (if requested)
    if (co == 0) {
        bssn_constraints_kernel<ghost_width><<<grid, block, smem_size, stream>>>(
            ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, X, Y, Z,
            chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
            Axx, Axy, Axz, Ayy, Ayz, Azz, Gamx, Gamy, Gamz, rho, Sx, Sy, Sz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
            ham_Res, movx_Res, movy_Res, movz_Res,
            Gmx_Res, Gmy_Res, Gmz_Res
        );
    }
}

// Constraint-only: 只跑 K1 (ricci) + K4 (constraints)，跳过 K2+K3
void gpu_compute_constraints_only_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
) {
    dim3 block(DIMX, DIMY, DIMZ);
    dim3 grid(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );

    // Dynamic shared memory size for K1 and K4: 3 tiles of doubles
    constexpr int smem_size = 3 * TileParams<ghost_width>::SIZE * sizeof(double);
    static bool smem_attr_set = false;
    if (!smem_attr_set) {
        cudaFuncSetAttribute(bssn_ricci_kernel<ghost_width>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        cudaFuncSetAttribute(bssn_constraints_kernel<ghost_width>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        smem_attr_set = true;
    }

    // 1. 生成共形几何
    bssn_ricci_kernel<ghost_width><<<grid, block, smem_size, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, X, Y, Z,
        dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz
    );

    // 4. 诊断系统 (跳过 K2 rhs_eval 和 K3 advection_dissipation)
    bssn_constraints_kernel<ghost_width><<<grid, block, smem_size, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, Gamx, Gamy, Gamz, rho, Sx, Sy, Sz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        ham_Res, movx_Res, movy_Res, movz_Res,
        Gmx_Res, Gmy_Res, Gmz_Res
    );
}