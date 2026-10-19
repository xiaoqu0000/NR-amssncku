#include "macrodef.h"
#include "Z4c_rhs.h"
#include "bssn_constraint.h"
#include "bssn_constraint_gpu.cuh"
#include "rhs_smem_common.cuh"
#include "kodiss_gpu.cuh"
#include "lopsidediff_gpu.cuh"
#include "Z4c_rhs_gpu.cuh"
#include "gpu_info.h"
#include "gpu_manager.h"
#include "macrodef.fh"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math_constants.h>
#include <math.h>

#define IDX3D_PAD(i, j, k, pnx, pny) ((i) + (pnx) * ((j) + (pny) * (k)))

void gpu_ricci_gamma_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Gamx, const double* d_Gamy, const double* d_Gamz,
    double* d_Gamxxx, double* d_Gamxxy, double* d_Gamxxz,
    double* d_Gamxyy, double* d_Gamxyz, double* d_Gamxzz,
    double* d_Gamyxx, double* d_Gamyxy, double* d_Gamyxz,
    double* d_Gamyyy, double* d_Gamyyz, double* d_Gamyzz,
    double* d_Gamzxx, double* d_Gamzxy, double* d_Gamzxz,
    double* d_Gamzyy, double* d_Gamzyz, double* d_Gamzzz,
    double* d_Rxx, double* d_Rxy, double* d_Rxz,
    double* d_Ryy, double* d_Ryz, double* d_Rzz,
    int symmetry
);

// BH data globals (defined in transpbh.C)
extern int BHN;
extern double Mass[3];
extern double PBH[9];

namespace {

constexpr double ONE = 1.0;
constexpr double TWO = 2.0;
constexpr double HALF = 0.5;
constexpr double F1o3 = 1.0 / 3.0;
constexpr double FOUR_PI = 4.0 * M_PI;
constexpr double FF = 0.75;
constexpr double eta = 2.0;

inline int z4c_device_ghost_width(const int* ex, const int* pex)
{
    return (pex[1] - ex[1] - (GPU_DIM_Y - 1)) / 2;
}

inline int z4c_device_pad_x_left(int z4c_ghost_width)
{
    return ((z4c_ghost_width + GPU_ALIGNMENT_DOUBLE - 1) / GPU_ALIGNMENT_DOUBLE) *
           GPU_ALIGNMENT_DOUBLE;
}

inline int z4c_device_ghost_offset(const int* ex, const int* pex)
{
    const int z4c_ghost_width = z4c_device_ghost_width(ex, pex);
    return (z4c_ghost_width * pex[1] + z4c_ghost_width) * pex[0] +
           z4c_device_pad_x_left(z4c_ghost_width);
}

// Z4C owns these CUDA kernels. They are copied from the BSSN CUDA
// implementation where the active equations are identical, then launched by
// the Z4C path with Z4C-specific inputs and correction kernels.
template<int GW>
__launch_bounds__(256, 2)
__global__ void z4c_base_ricci_kernel(
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
__global__ void z4c_base_rhs_eval_kernel(
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
    const double F3o2 = 1.5; const double HALF = 0.5; const double EIGHT = 8.0; const double F16 = 16.0; const double PI = M_PI;

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
    } else {
        betax_rhs[idx] = 0.0;
        betay_rhs[idx] = 0.0;
        betaz_rhs[idx] = 0.0;
        dtSfx_rhs[idx] = 0.0;
        dtSfy_rhs[idx] = 0.0;
        dtSfz_rhs[idx] = 0.0;
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

template<int GW>
__device__ __forceinline__ double z4c_kodis_point(
    const int ex[3], const int pex[3], const double* f,
    const double* X, const double* Y, const double* Z,
    double sym1, double sym2, double sym3,
    int symmetry, double eps, int i, int j, int k)
{
    constexpr int NO_SYMM = 0;
    constexpr int EQ_SYMM = 1;
    constexpr int OCTANT = 2;
    constexpr double cof = double(1 << (2 * GW));
    constexpr double sign = ((GW - 1) % 2 == 0) ? 1.0 : -1.0;
    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];
    const int imax = ex[0] - 1;
    const int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;
    int imin = 0;
    int jmin = 0;
    int kmin = 0;
    if (symmetry > NO_SYMM && fabs(Z[0]) < dZ) kmin = -GW;
    if (symmetry > EQ_SYMM && fabs(X[0]) < dX) imin = -GW;
    if (symmetry > EQ_SYMM && fabs(Y[0]) < dY) jmin = -GW;
    if constexpr (GW == 3) {
        if (symmetry == OCTANT && fabs(X[0]) < dX) imin = -GW;
        if (symmetry == OCTANT && fabs(Y[0]) < dY) jmin = -GW;
    }
    if (i - GW < imin || i + GW > imax ||
        j - GW < jmin || j + GW > jmax ||
        k - GW < kmin || k + GW > kmax) {
        return 0.0;
    }
    const auto fh = [&](int ii, int jj, int kk) {
        return d_symmetry_bd_0b_pad(
            GW, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, sym1, sym2, sym3);
    };
    double sx = 0.0, sy = 0.0, sz = 0.0;
    if constexpr (GW == 2) {
        sx = (fh(i-2,j,k) + fh(i+2,j,k)) - 4.0 * (fh(i-1,j,k) + fh(i+1,j,k)) + 6.0 * fh(i,j,k);
        sy = (fh(i,j-2,k) + fh(i,j+2,k)) - 4.0 * (fh(i,j-1,k) + fh(i,j+1,k)) + 6.0 * fh(i,j,k);
        sz = (fh(i,j,k-2) + fh(i,j,k+2)) - 4.0 * (fh(i,j,k-1) + fh(i,j,k+1)) + 6.0 * fh(i,j,k);
    } else if constexpr (GW == 3) {
        sx = (fh(i-3,j,k) + fh(i+3,j,k)) - 6.0 * (fh(i-2,j,k) + fh(i+2,j,k)) + 15.0 * (fh(i-1,j,k) + fh(i+1,j,k)) - 20.0 * fh(i,j,k);
        sy = (fh(i,j-3,k) + fh(i,j+3,k)) - 6.0 * (fh(i,j-2,k) + fh(i,j+2,k)) + 15.0 * (fh(i,j-1,k) + fh(i,j+1,k)) - 20.0 * fh(i,j,k);
        sz = (fh(i,j,k-3) + fh(i,j,k+3)) - 6.0 * (fh(i,j,k-2) + fh(i,j,k+2)) + 15.0 * (fh(i,j,k-1) + fh(i,j,k+1)) - 20.0 * fh(i,j,k);
    } else if constexpr (GW == 4) {
        sx = (fh(i-4,j,k) + fh(i+4,j,k)) - 8.0 * (fh(i-3,j,k) + fh(i+3,j,k)) + 28.0 * (fh(i-2,j,k) + fh(i+2,j,k)) - 56.0 * (fh(i-1,j,k) + fh(i+1,j,k)) + 70.0 * fh(i,j,k);
        sy = (fh(i,j-4,k) + fh(i,j+4,k)) - 8.0 * (fh(i,j-3,k) + fh(i,j+3,k)) + 28.0 * (fh(i,j-2,k) + fh(i,j+2,k)) - 56.0 * (fh(i,j-1,k) + fh(i,j+1,k)) + 70.0 * fh(i,j,k);
        sz = (fh(i,j,k-4) + fh(i,j,k+4)) - 8.0 * (fh(i,j,k-3) + fh(i,j,k+3)) + 28.0 * (fh(i,j,k-2) + fh(i,j,k+2)) - 56.0 * (fh(i,j,k-1) + fh(i,j,k+1)) + 70.0 * fh(i,j,k);
    } else if constexpr (GW == 5) {
        sx = (fh(i-5,j,k) + fh(i+5,j,k)) - 10.0 * (fh(i-4,j,k) + fh(i+4,j,k)) + 45.0 * (fh(i-3,j,k) + fh(i+3,j,k)) - 120.0 * (fh(i-2,j,k) + fh(i+2,j,k)) + 210.0 * (fh(i-1,j,k) + fh(i+1,j,k)) - 252.0 * fh(i,j,k);
        sy = (fh(i,j-5,k) + fh(i,j+5,k)) - 10.0 * (fh(i,j-4,k) + fh(i,j+4,k)) + 45.0 * (fh(i,j-3,k) + fh(i,j+3,k)) - 120.0 * (fh(i,j-2,k) + fh(i,j+2,k)) + 210.0 * (fh(i,j-1,k) + fh(i,j+1,k)) - 252.0 * fh(i,j,k);
        sz = (fh(i,j,k-5) + fh(i,j,k+5)) - 10.0 * (fh(i,j,k-4) + fh(i,j,k+4)) + 45.0 * (fh(i,j,k-3) + fh(i,j,k+3)) - 120.0 * (fh(i,j,k-2) + fh(i,j,k+2)) + 210.0 * (fh(i,j,k-1) + fh(i,j,k+1)) - 252.0 * fh(i,j,k);
    }
    return sign * eps / cof * (sx / dX + sy / dY + sz / dZ);
}

// K3: Advection + KO Dissipation unified fused path for all GW
template<int GW>
__launch_bounds__(256, 2)
__global__ void z4c_base_advection_dissipation_kernel(
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
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;

    const int imax = ex0 - 1;
    const int jmax = ex1 - 1;
    const int kmax = ex2 - 1;

    if (i > imax || j > jmax || k > kmax) return;

    const int idx = i + pex0 * (j + pex1 * k);

    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    const double *f;
    double *f_rhs;

    #pragma unroll 1
    for (int step = 0; step < 24; ++step) {
        GET_SINGLE_STEP_POINTERS(step);
        double sym1 = ONE;
        double sym2 = ONE;
        double sym3 = ONE;
        switch (step) {
            case 1:
            case 7:
                sym1 = -ONE;
                sym2 = -ONE;
                break;
            case 2:
            case 8:
                sym1 = -ONE;
                sym3 = -ONE;
                break;
            case 4:
            case 10:
                sym2 = -ONE;
                sym3 = -ONE;
                break;
            case 14:
            case 18:
            case 21:
                sym1 = -ONE;
                break;
            case 15:
            case 19:
            case 22:
                sym2 = -ONE;
                break;
            case 16:
            case 20:
            case 23:
                sym3 = -ONE;
                break;
            default:
                break;
        }
        double res = d_lopsided_point<GW>(
            ex, pex, f, betax, betay, betaz, X, Y, Z,
            symmetry, sym1, sym2, sym3, i, j, k);
        if (eps > 0.0) {
            res += z4c_kodis_point<GW>(
                ex, pex, f, X, Y, Z, sym1, sym2, sym3,
                symmetry, eps, i, j, k);
        }
        f_rhs[idx] += res;
    }
    return;

}


template<int GW>
__launch_bounds__(256, 2)
__global__ void z4c_base_constraints_kernel(
    const int ex0, const int ex1, const int ex2,
    const int pex0, const int pex1, const int pex2,
    const int symmetry, const int lev,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi,
    const double* __restrict__ trK_hamiltonian,
    const double* __restrict__ trK_momentum,
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
    // 3. Batch 3: chi, trK_momentum (借用 load_3_fields，第三个参数传 trK 占位)
    // =========================================================================
    load_3_fields_to_smem<GW>(smem3, chi, trK_momentum, trK_momentum, pdims, block_i, block_j, block_k);
    __syncthreads();

    double chix, chiy, chiz, chixx, chixy, chixz, chiyy, chiyz, chizz;
    compute_all_derivs_smem<GW>(smem3[0], &chix, &chiy, &chiz, &chixx, &chixy, &chixz, &chiyy, &chiyz, &chizz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    double Kx, Ky, Kz;
    compute_first_derivs_smem<GW>(smem3[1], &Kx, &Ky, &Kz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    __syncthreads();

    double gxxx, gxxy, gxxz;
    double gxyx, gxyy, gxyz;
    double gxzx, gxzy, gxzz;
    load_3_fields_to_smem<GW>(smem3, dxx, gxy, gxz, pdims, block_i, block_j, block_k);
    __syncthreads();
    compute_first_derivs_smem<GW>(smem3[0], &gxxx, &gxxy, &gxxz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    compute_first_derivs_smem<GW>(smem3[1], &gxyx, &gxyy, &gxyz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    compute_first_derivs_smem<GW>(smem3[2], &gxzx, &gxzy, &gxzz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    __syncthreads();

    double gyyx, gyyy, gyyz;
    double gyzx, gyzy, gyzz;
    double gzzx, gzzy, gzzz;
    load_3_fields_to_smem<GW>(smem3, dyy, gyz, dzz, pdims, block_i, block_j, block_k);
    __syncthreads();
    compute_first_derivs_smem<GW>(smem3[0], &gyyx, &gyyy, &gyyz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    compute_first_derivs_smem<GW>(smem3[1], &gyzx, &gyzy, &gyzz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);
    compute_first_derivs_smem<GW>(smem3[2], &gzzx, &gzzy, &gzzz, dX, dY, dZ, i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    // =========================================================================
    // 4. 数学降维与物理组装 (完全保留原有数学逻辑)
    // =========================================================================
    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = IDX3D_PAD(i, j, k, pex0, pex1);

    const double TWO = 2.0; const double F2o3 = 2.0 / 3.0;
    const double EIGHT = 8.0; const double F16 = 16.0; const double PI = M_PI;

    double chin1 = chi[idx] + 1.0; double val_trK = trK_hamiltonian[idx];
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

    bssn_constraint_gpu::gamma_constraint_residual(
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        gxxx, gxyx, gxzx, gyyx, gyzx, gzzx,
        gxxy, gxyy, gxzy, gyyy, gyzy, gzzy,
        gxxz, gxyz, gxzz, gyyz, gyzz, gzzz,
        Gamx[idx], Gamy[idx], Gamz[idx], Gmx_Res[idx], Gmy_Res[idx], Gmz_Res[idx]);

    // Compute A^{ij} A_{ij} = sum_{i,j,k,l} gup^{ij} gup^{kl} A_{ik} A_{jl}
    // using the exact nested expression from Z4c_rhs.f90:1077-1102 to
    // preserve Fortran floating-point summation order.
    double trA2 =
        gupxx * (
            gupxx * l_Axx * l_Axx + gupyy * l_Axy * l_Axy + gupzz * l_Axz * l_Axz +
            TWO * (gupxy * l_Axx * l_Axy + gupxz * l_Axx * l_Axz + gupyz * l_Axy * l_Axz)
        ) +
        gupyy * (
            gupxx * l_Axy * l_Axy + gupyy * l_Ayy * l_Ayy + gupzz * l_Ayz * l_Ayz +
            TWO * (gupxy * l_Axy * l_Ayy + gupxz * l_Axy * l_Ayz + gupyz * l_Ayy * l_Ayz)
        ) +
        gupzz * (
            gupxx * l_Axz * l_Axz + gupyy * l_Ayz * l_Ayz + gupzz * l_Azz * l_Azz +
            TWO * (gupxy * l_Axz * l_Ayz + gupxz * l_Axz * l_Azz + gupyz * l_Ayz * l_Azz)
        ) +
        TWO * (
            gupxy * (
                gupxx * l_Axx * l_Axy + gupyy * l_Axy * l_Ayy + gupzz * l_Axz * l_Ayz +
                gupxy * (l_Axx * l_Ayy + l_Axy * l_Axy) +
                gupxz * (l_Axx * l_Ayz + l_Axz * l_Axy) +
                gupyz * (l_Axy * l_Ayz + l_Axz * l_Ayy)
            ) +
            gupxz * (
                gupxx * l_Axx * l_Axz + gupyy * l_Axy * l_Ayz + gupzz * l_Axz * l_Azz +
                gupxy * (l_Axx * l_Ayz + l_Axy * l_Axz) +
                gupxz * (l_Axx * l_Azz + l_Axz * l_Axz) +
                gupyz * (l_Axy * l_Azz + l_Axz * l_Ayz)
            ) +
            gupyz * (
                gupxx * l_Axy * l_Axz + gupyy * l_Ayy * l_Ayz + gupzz * l_Ayz * l_Azz +
                gupxy * (l_Axy * l_Ayz + l_Ayy * l_Axz) +
                gupxz * (l_Axy * l_Azz + l_Ayz * l_Axz) +
                gupyz * (l_Ayy * l_Azz + l_Ayz * l_Ayz)
            )
        );

    double conf_R = gupxx * Rxx_in[idx] + gupyy * Ryy_in[idx] + gupzz * Rzz_in[idx] + TWO * (gupxy * Rxy_in[idx] + gupxz * Rxz_in[idx] + gupyz * Ryz_in[idx]);
    ham_Res[idx] = chin1 * conf_R + F2o3 * val_trK * val_trK - trA2 - F16 * PI * rho[idx];

    bssn_constraint_gpu::covariant_momentum(
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        l_Axx, l_Axy, l_Axz, l_Ayy, l_Ayz, l_Azz,
        d_Axx_x, d_Axy_x, d_Axz_x, d_Ayy_x, d_Ayz_x, d_Azz_x,
        d_Axx_y, d_Axy_y, d_Axz_y, d_Ayy_y, d_Ayz_y, d_Azz_y,
        d_Axx_z, d_Axy_z, d_Axz_z, d_Ayy_z, d_Ayz_z, d_Azz_z,
        l_Gamxxx, l_Gamxxy, l_Gamxxz, l_Gamxyy, l_Gamxyz, l_Gamxzz,
        l_Gamyxx, l_Gamyxy, l_Gamyxz, l_Gamyyy, l_Gamyyz, l_Gamyzz,
        l_Gamzxx, l_Gamzxy, l_Gamzxz, l_Gamzyy, l_Gamzyz, l_Gamzzz,
        chix, chiy, chiz, inv_chin1, movx_Res[idx], movy_Res[idx], movz_Res[idx]);
    movx_Res[idx] -= F2o3 * Kx + EIGHT * PI * Sx[idx];
    movy_Res[idx] -= F2o3 * Ky + EIGHT * PI * Sy[idx];
    movz_Res[idx] -= F2o3 * Kz + EIGHT * PI * Sz[idx];
}

#ifndef Z4C_KAPPA1
#define Z4C_KAPPA1 0.02
#endif

#ifndef Z4C_KAPPA2
#define Z4C_KAPPA2 0.0
#endif

__global__ void z4c_build_trkd_kernel(
    int pex0, int pex1, int pex2,
    int pad_x_left, int z4c_ghost_width,
    const double* __restrict__ trK,
    const double* __restrict__ TZ,
    double* __restrict__ trKd,
    double* __restrict__ trKm)
{
    int pi = blockIdx.x * blockDim.x + threadIdx.x;
    int pj = blockIdx.y * blockDim.y + threadIdx.y;
    int pk = blockIdx.z * blockDim.z + threadIdx.z;
    if (pi >= pex0 || pj >= pex1 || pk >= pex2) return;

    int i = pi - pad_x_left;
    int j = pj - z4c_ghost_width;
    int k = pk - z4c_ghost_width;
    int idx = z4c_idx3d(i, j, k, pex0, pex1);
    trKd[idx] = trK[idx] + TWO * TZ[idx];
    trKm[idx] = trK[idx] + HALF * TZ[idx];
}

__global__ void z4c_lowerbound_chi_kernel(
    int total,
    double* __restrict__ chi, double chitiny)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= total) return;
    if (chi[tid] < chitiny) chi[tid] = chitiny;
}

__global__ void z4c_gamma_d_residual_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    const double* __restrict__ Gamxxx, const double* __restrict__ Gamxxy, const double* __restrict__ Gamxxz,
    const double* __restrict__ Gamxyy, const double* __restrict__ Gamxyz, const double* __restrict__ Gamxzz,
    const double* __restrict__ Gamyxx, const double* __restrict__ Gamyxy, const double* __restrict__ Gamyxz,
    const double* __restrict__ Gamyyy, const double* __restrict__ Gamyyz, const double* __restrict__ Gamyzz,
    const double* __restrict__ Gamzxx, const double* __restrict__ Gamzxy, const double* __restrict__ Gamzxz,
    const double* __restrict__ Gamzyy, const double* __restrict__ Gamzyz, const double* __restrict__ Gamzzz,
    double* __restrict__ Gmx_res, double* __restrict__ Gmy_res, double* __restrict__ Gmz_res)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int idx = z4c_idx3d(i, j, k, pex0, pex1);
    const double l_gxx = dxx[idx] + ONE;
    const double l_gxy = gxy[idx];
    const double l_gxz = gxz[idx];
    const double l_gyy = dyy[idx] + ONE;
    const double l_gyz = gyz[idx];
    const double l_gzz = dzz[idx] + ONE;
    const double det = l_gxx * l_gyy * l_gzz + TWO * l_gxy * l_gyz * l_gxz
                     - l_gxz * l_gyy * l_gxz - l_gxx * l_gyz * l_gyz - l_gzz * l_gxy * l_gxy;
    const double inv_det = ONE / det;
    const double gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) * inv_det;
    const double gupxy = (l_gyz * l_gxz - l_gxy * l_gzz) * inv_det;
    const double gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) * inv_det;
    const double gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) * inv_det;
    const double gupyz = (l_gxy * l_gxz - l_gxx * l_gyz) * inv_det;
    const double gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) * inv_det;
    const double Gamxa = gupxx * Gamxxx[idx] + gupyy * Gamxyy[idx] + gupzz * Gamxzz[idx]
        + TWO * (gupxy * Gamxxy[idx] + gupxz * Gamxxz[idx] + gupyz * Gamxyz[idx]);
    const double Gamya = gupxx * Gamyxx[idx] + gupyy * Gamyyy[idx] + gupzz * Gamyzz[idx]
        + TWO * (gupxy * Gamyxy[idx] + gupxz * Gamyxz[idx] + gupyz * Gamyyz[idx]);
    const double Gamza = gupxx * Gamzxx[idx] + gupyy * Gamzyy[idx] + gupzz * Gamzzz[idx]
        + TWO * (gupxy * Gamzxy[idx] + gupxz * Gamzxz[idx] + gupyz * Gamzyz[idx]);
    Gmx_res[idx] = Gamx[idx] - Gamxa;
    Gmy_res[idx] = Gamy[idx] - Gamya;
    Gmz_res[idx] = Gamz[idx] - Gamza;
}

__global__ void z4c_finalize_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* __restrict__ trK,
    const double* __restrict__ TZ,
    const double* __restrict__ Lap,
    const double* __restrict__ Hcon,
    const double* __restrict__ Gmx_d,
    const double* __restrict__ Gmy_d,
    const double* __restrict__ Gmz_d,
    double* __restrict__ trK_rhs,
    double* __restrict__ Gamx_rhs,
    double* __restrict__ Gamy_rhs,
    double* __restrict__ Gamz_rhs,
    double* __restrict__ Lap_rhs,
    double* __restrict__ TZ_rhs,
    int enable_damping)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int idx = z4c_idx3d(i, j, k, pex0, pex1);
    const double alpha = Lap[idx] + ONE;

    // Base RHS was intentionally evaluated with trKd = trK + 2*Theta.
    // Restore the Z4c lapse equation, which uses the unshifted trK.
    Lap_rhs[idx] += TWO * alpha * (TWO * TZ[idx]);

    // GPU path keeps the Hamiltonian constraint on device. Hcon is required
    // here and by the later Hamiltonian trace-correction launch.
    const double H = Hcon[idx];
    TZ_rhs[idx] = alpha * H * HALF;

    if (enable_damping) {
        const double kappa1 = Z4C_KAPPA1;
        const double kappa2 = Z4C_KAPPA2;
        // Note: TZ damping is applied in z4c_tz_advection_dissipation_kernel
        // BETWEEN the advection and the KO dissipation to match the CPU
        // Fortran update sequence (Z4c_rhs.f90 lines 1226 → 1229 → 1262).
        // Keeping it here breaks bit-exact agreement because (A + B - C + D) ≠
        // (A - C) + (B + D) in floating point.
        trK_rhs[idx] += alpha * kappa1 * (ONE - kappa2) * TZ[idx];
        Gamx_rhs[idx] -= TWO * alpha * kappa1 * Gmx_d[idx];
        Gamy_rhs[idx] -= TWO * alpha * kappa1 * Gmy_d[idx];
        Gamz_rhs[idx] -= TWO * alpha * kappa1 * Gmz_d[idx];
    }

}

__global__ void z4c_gamma_tz_gradient_correction_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    int symmetry, int lev,
    const double* __restrict__ X,
    const double* __restrict__ Y,
    const double* __restrict__ Z,
    const double* __restrict__ dxx,
    const double* __restrict__ gxy,
    const double* __restrict__ gxz,
    const double* __restrict__ dyy,
    const double* __restrict__ gyz,
    const double* __restrict__ dzz,
    const double* __restrict__ Lap,
    const double* __restrict__ TZ,
    double* __restrict__ Gamx_rhs,
    double* __restrict__ Gamy_rhs,
    double* __restrict__ Gamz_rhs,
    double* __restrict__ dtSfx_rhs,
    double* __restrict__ dtSfy_rhs,
    double* __restrict__ dtSfz_rhs)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;

    const int dims[3] = {ex0, ex1, ex2};
    const int pdims[3] = {pex0, pex1, pex2};
    const int block_i = blockIdx.x * blockDim.x;
    const int block_j = blockIdx.y * blockDim.y;
    const int block_k = blockIdx.z * blockDim.z;
    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    int imin = 0, jmin = 0, kmin = 0;
    const int imax = ex0 - 1;
    const int jmax = ex1 - 1;
    const int kmax = ex2 - 1;
    constexpr int ORD = ghost_width - 1;
    if (symmetry > 0 && fabs(Z[0]) < dZ) kmin = -ORD;
    if (symmetry > 1 && fabs(X[0]) < dX) imin = -ORD;
    if (symmetry > 1 && fabs(Y[0]) < dY) jmin = -ORD;

    __shared__ double smem[TileParams<ghost_width>::Z][TileParams<ghost_width>::Y][TileParams<ghost_width>::X];
    load_field_to_smem<ghost_width>(smem, TZ, dims, pdims, block_i, block_j, block_k, ONE, ONE, ONE);
    __syncthreads();

    double TZx = 0.0, TZy = 0.0, TZz = 0.0;
    compute_first_derivs_smem<ghost_width>(smem, &TZx, &TZy, &TZz, dX, dY, dZ,
                                           i, j, k, imin, jmin, kmin, imax, jmax, kmax);

    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int idx = z4c_idx3d(i, j, k, pex0, pex1);
    const double gxx = dxx[idx] + ONE;
    const double gyy = dyy[idx] + ONE;
    const double gzz = dzz[idx] + ONE;
    const double l_gxy = gxy[idx];
    const double l_gxz = gxz[idx];
    const double l_gyz = gyz[idx];
    const double det = gxx * gyy * gzz + TWO * l_gxy * l_gyz * l_gxz
                     - l_gxz * gyy * l_gxz - l_gxy * l_gxy * gzz - gxx * l_gyz * l_gyz;
    const double inv_det = ONE / det;
    const double gupxx = (gyy * gzz - l_gyz * l_gyz) * inv_det;
    const double gupxy = (l_gyz * l_gxz - l_gxy * gzz) * inv_det;
    const double gupxz = (l_gxy * l_gyz - gyy * l_gxz) * inv_det;
    const double gupyy = (gxx * gzz - l_gxz * l_gxz) * inv_det;
    const double gupyz = (l_gxy * l_gxz - gxx * l_gyz) * inv_det;
    const double gupzz = (gxx * gyy - l_gxy * l_gxy) * inv_det;
    const double alpha = Lap[idx] + ONE;

    // Base RHS was evaluated with d(trK + 2*Theta). Active Z4c uses the same
    // shifted trK in Aij/constraints, but only half of dTheta in Gamma.
    const double corr_x = TWO * alpha * (gupxx * TZx + gupxy * TZy + gupxz * TZz);
    const double corr_y = TWO * alpha * (gupxy * TZx + gupyy * TZy + gupyz * TZz);
    const double corr_z = TWO * alpha * (gupxz * TZx + gupyz * TZy + gupzz * TZz);
    Gamx_rhs[idx] += corr_x;
    Gamy_rhs[idx] += corr_y;
    Gamz_rhs[idx] += corr_z;
#if (GAUGE == 0)
    dtSfx_rhs[idx] += corr_x;
    dtSfy_rhs[idx] += corr_y;
    dtSfz_rhs[idx] += corr_z;
#endif

    (void)lev;
}

__global__ void z4c_aij_hamiltonian_trace_correction_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1,
    const double* __restrict__ dxx,
    const double* __restrict__ gxy,
    const double* __restrict__ gxz,
    const double* __restrict__ dyy,
    const double* __restrict__ gyz,
    const double* __restrict__ dzz,
    const double* __restrict__ Lap,
    const double* __restrict__ Hcon,
    double* __restrict__ Axx_rhs,
    double* __restrict__ Axy_rhs,
    double* __restrict__ Axz_rhs,
    double* __restrict__ Ayy_rhs,
    double* __restrict__ Ayz_rhs,
    double* __restrict__ Azz_rhs)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int idx = z4c_idx3d(i, j, k, pex0, pex1);
    // Active Z4c subtracts alpha*H/3 before the trace-free projection for Aij.
    // BSSN lacks this Z4 Hamiltonian trace term, so add the equivalent
    // post-projection correction: -alpha*H*g_ij/3.
    const double factor = -(Lap[idx] + ONE) * Hcon[idx] * F1o3;
    Axx_rhs[idx] += factor * (dxx[idx] + ONE);
    Axy_rhs[idx] += factor * gxy[idx];
    Axz_rhs[idx] += factor * gxz[idx];
    Ayy_rhs[idx] += factor * (dyy[idx] + ONE);
    Ayz_rhs[idx] += factor * gyz[idx];
    Azz_rhs[idx] += factor * (dzz[idx] + ONE);
}

__global__ void z4c_tz_advection_dissipation_kernel(
    int ex0, int ex1, int ex2, int pex0, int pex1, int pex2,
    int symmetry, double eps,
    const double* __restrict__ X,
    const double* __restrict__ Y,
    const double* __restrict__ Z,
    const double* __restrict__ betax,
    const double* __restrict__ betay,
    const double* __restrict__ betaz,
    const double* __restrict__ TZ,
    const double* __restrict__ Lap,
    double* __restrict__ TZ_rhs,
    int enable_damping)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int idx = z4c_idx3d(i, j, k, pex0, pex1);
    const int ex[3] = {ex0, ex1, ex2};
    const int pex[3] = {pex0, pex1, pex2};
    // Match CPU Fortran ordering exactly:
    //   TZ_rhs += advection  (line 1226 in Z4c_rhs.f90)
    //   TZ_rhs -= damping    (line 1229)
    //   TZ_rhs += KO         (line 1262)
    // Each line is a separate accumulating update so the FP rounding
    // matches the Fortran scalar evaluation order.
    const double tz_advection = d_lopsided_point<ghost_width>(
        ex, pex, TZ, betax, betay, betaz, X, Y, Z,
        symmetry, ONE, ONE, ONE, i, j, k);
    TZ_rhs[idx] += tz_advection;
    if (enable_damping) {
        const double alpha = Lap[idx] + ONE;
        const double kappa1 = Z4C_KAPPA1;
        const double kappa2 = Z4C_KAPPA2;
        TZ_rhs[idx] -= alpha * (TWO + kappa2) * kappa1 * TZ[idx];
    }
    if (eps > 0.0) {
        const double tz_kodis = z4c_kodis_point<ghost_width>(
            ex, pex, TZ, X, Y, Z, ONE, ONE, ONE,
            symmetry, eps, i, j, k);
        TZ_rhs[idx] += tz_kodis;
    }

}

void gpu_compute_rhs_z4c_common_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* constraint_chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz, double* TZ,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs, double* TZ_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* Hcon, double* Mxcon, double* Mycon, double* Mzcon,
    double* Gmxcon, double* Gmycon, double* Gmzcon,
    int symmetry, int lev, double eps, int enable_damping)
{
    dim3 block(8, 8, 4);
    dim3 grid((ex[0] + block.x - 1) / block.x,
              (ex[1] + block.y - 1) / block.y,
              (ex[2] + block.z - 1) / block.z);

    const size_t num_elements = static_cast<size_t>(pex[0]) * pex[1] * pex[2];
    const int z4c_ghost_width = z4c_device_ghost_width(ex, pex);
    const int z4c_pad_x_left = z4c_device_pad_x_left(z4c_ghost_width);
    const int ghost_offset = z4c_device_ghost_offset(ex, pex);

    dim3 grid_pad((pex[0] + block.x - 1) / block.x,
                  (pex[1] + block.y - 1) / block.y,
                  (pex[2] + block.z - 1) / block.z);
    GPUManager& gpu_manager = GPUManager::getInstance();
    double* trKd_storage = gpu_manager.allocate_device_memory<double>(num_elements, stream);
    double* trKm_storage = gpu_manager.allocate_device_memory<double>(num_elements, stream);
    double* Gmx_d_storage = gpu_manager.allocate_device_memory<double>(num_elements, stream);
    double* Gmy_d_storage = gpu_manager.allocate_device_memory<double>(num_elements, stream);
    double* Gmz_d_storage = gpu_manager.allocate_device_memory<double>(num_elements, stream);
    double* trKd = trKd_storage + ghost_offset;
    double* trKm = trKm_storage + ghost_offset;
    double* Gmx_d = Gmx_d_storage + ghost_offset;
    double* Gmy_d = Gmy_d_storage + ghost_offset;
    double* Gmz_d = Gmz_d_storage + ghost_offset;
    z4c_build_trkd_kernel<<<grid_pad, block, 0, stream>>>(
        pex[0], pex[1], pex[2], z4c_pad_x_left, z4c_ghost_width, trK, TZ, trKd, trKm);

    constexpr int smem_size = 3 * TileParams<ghost_width>::SIZE * sizeof(double);
    static bool smem_attr_set = false;
    if (!smem_attr_set) {
        cudaFuncSetAttribute(z4c_base_ricci_kernel<ghost_width>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        cudaFuncSetAttribute(z4c_base_rhs_eval_kernel<ghost_width, GAUGE>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        cudaFuncSetAttribute(z4c_base_constraints_kernel<ghost_width>,
                             cudaFuncAttributeMaxDynamicSharedMemorySize, smem_size);
        smem_attr_set = true;
    }

    z4c_base_ricci_kernel<ghost_width><<<grid, block, smem_size, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, X, Y, Z,
        dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz);

    z4c_base_rhs_eval_kernel<ghost_width, GAUGE><<<grid, block, smem_size, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev,
        BHN, PBH[0], PBH[1], PBH[2], PBH[3], PBH[4], PBH[5], Mass[0], Mass[1],
        X, Y, Z,
        chi, trKd, dxx, gxy, gxz, dyy, gyz, dzz,
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
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs);

    dim3 block_adv(32, 8, 1);
    dim3 grid_adv((ex[0] + block_adv.x - 1) / block_adv.x,
                  (ex[1] + block_adv.y - 1) / block_adv.y,
                  (ex[2] + block_adv.z - 1) / block_adv.z);
    z4c_base_advection_dissipation_kernel<ghost_width><<<grid_adv, block_adv, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, eps, X, Y, Z,
        betax, betay, betaz, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, chi, trK, Gamx, Gamy, Gamz, Lap, dtSfx, dtSfy, dtSfz,
        gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        chi_rhs, trK_rhs, Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs, dtSfx_rhs, dtSfy_rhs, dtSfz_rhs);

    z4c_gamma_d_residual_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1],
        dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Gmx_d, Gmy_d, Gmz_d);

    gpu_ricci_gamma_launch(
        stream, ex, pex, X, Y, Z, chi,
        dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, symmetry);

    z4c_base_constraints_kernel<ghost_width><<<grid, block, smem_size, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, X, Y, Z,
        chi, trKd, trKm, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, Gamx, Gamy, Gamz, rho, Sx, Sy, Sz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        Hcon, Mxcon, Mycon, Mzcon, Gmxcon, Gmycon, Gmzcon);

    z4c_gamma_tz_gradient_correction_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev,
        X, Y, Z, dxx, gxy, gxz, dyy, gyz, dzz, Lap, TZ,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, dtSfx_rhs, dtSfy_rhs, dtSfz_rhs);

    z4c_aij_hamiltonian_trace_correction_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1],
        dxx, gxy, gxz, dyy, gyz, dzz, Lap, Hcon,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs);

    z4c_finalize_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], trK, TZ, Lap, Hcon,
        Gmx_d, Gmy_d, Gmz_d,
        trK_rhs, Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs, TZ_rhs, enable_damping);

    z4c_tz_advection_dissipation_kernel<<<grid_adv, block_adv, 0, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        symmetry, eps, X, Y, Z, betax, betay, betaz, TZ, Lap, TZ_rhs,
        enable_damping);

    // Active Z4c_rhs.f90 ignores co and always returns constraints evaluated
    // with the unshifted trK after computing the Z4 RHS. The z4cnot wrapper
    // computes RHS with lowerbounded chi but recomputes final constraints with
    // the original chi, so keep the final-constraint chi pointer separate.
#if (ABV == 1)
    gpu_constraint_bssn_launch(
        stream, ex, pex, X, Y, Z, constraint_chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap, betax, betay, betaz,
        rho, Sx, Sy, Sz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        Hcon, Mxcon, Mycon, Mzcon, Gmxcon, Gmycon, Gmzcon,
        symmetry);
#else
    if (constraint_chi != chi) {
        gpu_ricci_gamma_launch(
            stream, ex, pex, X, Y, Z, constraint_chi,
            dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, symmetry);
    }
    z4c_base_constraints_kernel<ghost_width><<<grid, block, smem_size, stream>>>(
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2], symmetry, lev, X, Y, Z,
        constraint_chi, trK, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, Gamx, Gamy, Gamz, rho, Sx, Sy, Sz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        Hcon, Mxcon, Mycon, Mzcon, Gmxcon, Gmycon, Gmzcon);
#endif

    gpu_manager.free_device_memory(trKd_storage, num_elements, stream);
    gpu_manager.free_device_memory(trKm_storage, num_elements, stream);
    gpu_manager.free_device_memory(Gmx_d_storage, num_elements, stream);
    gpu_manager.free_device_memory(Gmy_d_storage, num_elements, stream);
    gpu_manager.free_device_memory(Gmz_d_storage, num_elements, stream);
}

} // namespace

void gpu_compute_rhs_z4c_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz, double* TZ,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs, double* TZ_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* Hcon, double* Mxcon, double* Mycon, double* Mzcon,
    double* Gmxcon, double* Gmycon, double* Gmzcon,
    int symmetry, int lev, double eps, int co)
{
    (void)co;
    gpu_compute_rhs_z4c_common_launch(
        stream, ex, pex, T, X, Y, Z, chi, chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz, Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap, betax, betay, betaz, dtSfx, dtSfy, dtSfz, TZ,
        chi_rhs, trK_rhs, gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs, betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs, TZ_rhs,
        rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, Hcon, Mxcon, Mycon, Mzcon, Gmxcon, Gmycon, Gmzcon,
        symmetry, lev, eps, 1);
}

void gpu_compute_rhs_z4cnot_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz, double* TZ,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs, double* TZ_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* Hcon, double* Mxcon, double* Mycon, double* Mzcon,
    double* Gmxcon, double* Gmycon, double* Gmzcon,
    int symmetry, int lev, double eps, int co, double chitiny)
{
    (void)co;
    const size_t num_elements = static_cast<size_t>(pex[0]) * pex[1] * pex[2];
    const size_t bytes = num_elements * sizeof(double);
    const int ghost_offset = z4c_device_ghost_offset(ex, pex);
    GPUManager& gpu_manager = GPUManager::getInstance();
    double* chi_storage = gpu_manager.allocate_device_memory<double>(num_elements, stream);
    cudaMemcpyAsync(chi_storage, chi - ghost_offset, bytes, cudaMemcpyDeviceToDevice, stream);

    double* chihere = chi_storage + ghost_offset;
    const int total = pex[0] * pex[1] * pex[2];
    dim3 lower_block(256);
    dim3 lower_grid((total + lower_block.x - 1) / lower_block.x);
    z4c_lowerbound_chi_kernel<<<lower_grid, lower_block, 0, stream>>>(
        total, chi_storage, chitiny);

    gpu_compute_rhs_z4c_common_launch(
        stream, ex, pex, T, X, Y, Z, chihere, chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz, Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap, betax, betay, betaz, dtSfx, dtSfy, dtSfz, TZ,
        chi_rhs, trK_rhs, gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs, betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs, TZ_rhs,
        rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, Hcon, Mxcon, Mycon, Mzcon, Gmxcon, Gmycon, Gmzcon,
        symmetry, lev, eps, 1);
    gpu_manager.free_device_memory(chi_storage, num_elements, stream);
}
