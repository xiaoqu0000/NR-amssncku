#include "getnpem2.h"

#include <cmath>
#include <cuda_runtime.h>
#include "macrodef.h"
#include "getnpem2_gpu.cuh"

// ---------------------------------------------------------
// Device point functions for Phi2 computation
// ---------------------------------------------------------
template<int TETRAD>
__device__ __forceinline__ void d_getnpem2_point(
    double X, double Y, double Z,
    double chi,
    double dxx, double gxy, double gxz, double dyy, double gyz, double dzz,
    double Ex, double Ey, double Ez,
    double Bx, double By, double Bz,
    double& Rphi2, double& Iphi2
)
{
    const double ONE = 1.0;

    double vx, vy, vz, ux, uy, uz, wx, wy, wz;
    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;

    d_setup_em_tetrad<TETRAD>(
        X, Y, Z,
        dxx, gxy, gxz, dyy, gyz, dzz,
        chi,
        vx, vy, vz, ux, uy, uz, wx, wy, wz,
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
    );

    double chipn1 = chi + ONE;

    // Physical metric (for lowering E)
    double l_gxx = dxx + ONE;
    double l_gyy = dyy + ONE;
    double l_gzz = dzz + ONE;
    double l_gxy = gxy;
    double l_gxz = gxz;
    double l_gyz = gyz;

    // Lower E with physical metric, scale by chipn1
    double HEx = (l_gxx*Ex + l_gxy*Ey + l_gxz*Ez) * chipn1;
    double HEy = (l_gxy*Ex + l_gyy*Ey + l_gyz*Ez) * chipn1;
    double HEz = (l_gxz*Ex + l_gyz*Ey + l_gzz*Ez) * chipn1;

    // f = chipn1^(-3/2) for EM field scaling
    double t = sqrt(ONE / chipn1);
    double f_em = t * t * t;

    // HB = f * (v x B)
    double HBx = (vy*Bz - vz*By) * f_em;
    double HBy = (vz*Bx - vx*Bz) * f_em;
    double HBz = (vx*By - vy*Bx) * f_em;

    // Difference: HE - HB
    double dx = HEx - HBx;
    double dy = HEy - HBy;
    double dz = HEz - HBz;

    if constexpr (TETRAD == 1) {
        // Sperhake convention: m = (w + i*u)/sqrt(2)
        // Real part: (HE-HB) dot w
        Rphi2 = dx*wx + dy*wy + dz*wz;
        // Imaginary part: -(HE-HB) dot u
        Iphi2 = -(dx*ux + dy*uy + dz*uz);
    } else {
        // Frans convention (TETRAD 0,2): m = (u - i*w)/sqrt(2)
        // Real part: (HE-HB) dot u
        Rphi2 = dx*ux + dy*uy + dz*uz;
        // Imaginary part: (HE-HB) dot w
        Iphi2 = dx*wx + dy*wy + dz*wz;
    }

    Rphi2 /= 2.0;
    Iphi2 /= 2.0;
}

// ---------------------------------------------------------
// getnpem2 kernel: compute Phi2 over 3D grid
// ---------------------------------------------------------
template<int TETRAD>
__global__ void getnpem2_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* X, const double* Y, const double* Z,
    const double* chi,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    const double* Ex, const double* Ey, const double* Ez,
    const double* Bx, const double* By, const double* Bz,
    double* Rphi2, double* Iphi2,
    int symmetry
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    size_t idx = i + j * pex0 + k * (size_t)pex0 * pex1;

    double l_X = X[i];
    double l_Y = Y[j];
    double l_Z = Z[k];

    double l_Rphi2, l_Iphi2;

    d_getnpem2_point<TETRAD>(
        l_X, l_Y, l_Z,
        chi[idx],
        dxx[idx], gxy[idx], gxz[idx],
        dyy[idx], gyz[idx], dzz[idx],
        Ex[idx], Ey[idx], Ez[idx],
        Bx[idx], By[idx], Bz[idx],
        l_Rphi2, l_Iphi2
    );

    Rphi2[idx] = l_Rphi2;
    Iphi2[idx] = l_Iphi2;
}

// ---------------------------------------------------------
// getnpem1 kernel: compute Phi1 over 3D grid
// ---------------------------------------------------------
template<int TETRAD>
__global__ void getnpem1_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1, int pex2,
    const double* X, const double* Y, const double* Z,
    const double* chi,
    const double* dxx, const double* gxy, const double* gxz,
    const double* dyy, const double* gyz, const double* dzz,
    const double* Ex, const double* Ey, const double* Ez,
    const double* Bx, const double* By, const double* Bz,
    double* Rphi1, double* Iphi1,
    int symmetry
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    size_t idx = i + j * pex0 + k * (size_t)pex0 * pex1;

    double l_X = X[i];
    double l_Y = Y[j];
    double l_Z = Z[k];

    double l_Rphi1, l_Iphi1;

    d_getnpem1_point<TETRAD>(
        l_X, l_Y, l_Z,
        chi[idx],
        dxx[idx], gxy[idx], gxz[idx],
        dyy[idx], gyz[idx], dzz[idx],
        Ex[idx], Ey[idx], Ez[idx],
        Bx[idx], By[idx], Bz[idx],
        l_Rphi1, l_Iphi1
    );

    Rphi1[idx] = l_Rphi1;
    Iphi1[idx] = l_Iphi1;
}

// ---------------------------------------------------------
// Host launch wrapper: getnpem2 (Phi2)
// Dispatches on compile-time tetradtype macro.
// ---------------------------------------------------------
void gpu_getnpem2_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Ex, const double* d_Ey, const double* d_Ez,
    const double* d_Bx, const double* d_By, const double* d_Bz,
    double* d_Rphi2, double* d_Iphi2,
    int symmetry
) {
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    if constexpr (tetradtype == 0) {
        getnpem2_kernel<0><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1], pex[2],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi2, d_Iphi2,
            symmetry
        );
    } else if constexpr (tetradtype == 1) {
        getnpem2_kernel<1><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1], pex[2],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi2, d_Iphi2,
            symmetry
        );
    } else {
        getnpem2_kernel<2><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1], pex[2],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi2, d_Iphi2,
            symmetry
        );
    }
}

// ---------------------------------------------------------
// Host launch wrapper: getnpem1 (Phi1)
// Dispatches on compile-time tetradtype macro.
// ---------------------------------------------------------
void gpu_getnpem1_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Ex, const double* d_Ey, const double* d_Ez,
    const double* d_Bx, const double* d_By, const double* d_Bz,
    double* d_Rphi1, double* d_Iphi1,
    int symmetry
) {
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    if constexpr (tetradtype == 0) {
        getnpem1_kernel<0><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1], pex[2],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi1, d_Iphi1,
            symmetry
        );
    } else if constexpr (tetradtype == 1) {
        getnpem1_kernel<1><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1], pex[2],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi1, d_Iphi1,
            symmetry
        );
    } else {
        getnpem1_kernel<2><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1], pex[2],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi1, d_Iphi1,
            symmetry
        );
    }
}

// ---------------------------------------------------------
// getnpem2_ss_kernel: Shell-patch Phi2 per-point kernel
// Loads X,Y,Z from 3D arrays (vs 1D in Cartesian getnpem2_kernel)
// ---------------------------------------------------------
template<int TETRAD>
__global__ void getnpem2_ss_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    const double* __restrict__ d_X, const double* __restrict__ d_Y, const double* __restrict__ d_Z,
    const double* __restrict__ d_chi,
    const double* __restrict__ d_dxx, const double* __restrict__ d_gxy, const double* __restrict__ d_gxz,
    const double* __restrict__ d_dyy, const double* __restrict__ d_gyz, const double* __restrict__ d_dzz,
    const double* __restrict__ d_Ex, const double* __restrict__ d_Ey, const double* __restrict__ d_Ez,
    const double* __restrict__ d_Bx, const double* __restrict__ d_By, const double* __restrict__ d_Bz,
    double* __restrict__ d_Rphi2, double* __restrict__ d_Iphi2
) {
    int gx = blockIdx.x * blockDim.x + threadIdx.x;
    int gy = blockIdx.y * blockDim.y + threadIdx.y;
    int gz = blockIdx.z * blockDim.z + threadIdx.z;
    if (gx >= ex0 || gy >= ex1 || gz >= ex2) return;

    size_t idx = (size_t)gx + (size_t)pex0 * ((size_t)gy + (size_t)pex1 * (size_t)gz);

    double Xc = d_X[idx];
    double Yc = d_Y[idx];
    double Zc = d_Z[idx];

    double chi = d_chi[idx];
    double gxx = d_dxx[idx];
    double gxy = d_gxy[idx];
    double gxz = d_gxz[idx];
    double gyy = d_dyy[idx];
    double gyz = d_gyz[idx];
    double gzz = d_dzz[idx];
    double Ex = d_Ex[idx], Ey = d_Ey[idx], Ez = d_Ez[idx];
    double Bx = d_Bx[idx], By = d_By[idx], Bz = d_Bz[idx];

    double Rphi2, Iphi2;
    d_getnpem2_point<TETRAD>(Xc, Yc, Zc, chi, gxx, gxy, gxz, gyy, gyz, gzz,
                              Ex, Ey, Ez, Bx, By, Bz, Rphi2, Iphi2);

    d_Rphi2[idx] = Rphi2;
    d_Iphi2[idx] = Iphi2;
}

// ---------------------------------------------------------
// getnpem1_ss_kernel: Shell-patch Phi1 per-point kernel
// Loads X,Y,Z from 3D arrays (vs 1D in Cartesian getnpem1_kernel)
// ---------------------------------------------------------
template<int TETRAD>
__global__ void getnpem1_ss_kernel(
    int ex0, int ex1, int ex2,
    int pex0, int pex1,
    const double* __restrict__ d_X, const double* __restrict__ d_Y, const double* __restrict__ d_Z,
    const double* __restrict__ d_chi,
    const double* __restrict__ d_dxx, const double* __restrict__ d_gxy, const double* __restrict__ d_gxz,
    const double* __restrict__ d_dyy, const double* __restrict__ d_gyz, const double* __restrict__ d_dzz,
    const double* __restrict__ d_Ex, const double* __restrict__ d_Ey, const double* __restrict__ d_Ez,
    const double* __restrict__ d_Bx, const double* __restrict__ d_By, const double* __restrict__ d_Bz,
    double* __restrict__ d_Rphi1, double* __restrict__ d_Iphi1
) {
    int gx = blockIdx.x * blockDim.x + threadIdx.x;
    int gy = blockIdx.y * blockDim.y + threadIdx.y;
    int gz = blockIdx.z * blockDim.z + threadIdx.z;
    if (gx >= ex0 || gy >= ex1 || gz >= ex2) return;

    size_t idx = (size_t)gx + (size_t)pex0 * ((size_t)gy + (size_t)pex1 * (size_t)gz);

    double Xc = d_X[idx];
    double Yc = d_Y[idx];
    double Zc = d_Z[idx];

    double chi = d_chi[idx];
    double gxx = d_dxx[idx];
    double gxy = d_gxy[idx];
    double gxz = d_gxz[idx];
    double gyy = d_dyy[idx];
    double gyz = d_gyz[idx];
    double gzz = d_dzz[idx];
    double Ex = d_Ex[idx], Ey = d_Ey[idx], Ez = d_Ez[idx];
    double Bx = d_Bx[idx], By = d_By[idx], Bz = d_Bz[idx];

    double Rphi1, Iphi1;
    d_getnpem1_point<TETRAD>(Xc, Yc, Zc, chi, gxx, gxy, gxz, gyy, gyz, gzz,
                              Ex, Ey, Ez, Bx, By, Bz, Rphi1, Iphi1);

    d_Rphi1[idx] = Rphi1;
    d_Iphi1[idx] = Iphi1;
}

// ---------------------------------------------------------
// Host launch wrapper: getnpem2_ss (Phi2, Shell-patch)
// Dispatches on compile-time tetradtype macro.
// ---------------------------------------------------------
void gpu_getnpem2_ss_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Ex, const double* d_Ey, const double* d_Ez,
    const double* d_Bx, const double* d_By, const double* d_Bz,
    double* d_Rphi2, double* d_Iphi2,
    int symmetry
) {
    (void)symmetry;
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    if constexpr (tetradtype == 0) {
        getnpem2_ss_kernel<0><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi2, d_Iphi2
        );
    } else if constexpr (tetradtype == 1) {
        getnpem2_ss_kernel<1><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi2, d_Iphi2
        );
    } else {
        getnpem2_ss_kernel<2><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi2, d_Iphi2
        );
    }
}

// ---------------------------------------------------------
// Host launch wrapper: getnpem1_ss (Phi1, Shell-patch)
// Dispatches on compile-time tetradtype macro.
// ---------------------------------------------------------
void gpu_getnpem1_ss_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Ex, const double* d_Ey, const double* d_Ez,
    const double* d_Bx, const double* d_By, const double* d_Bz,
    double* d_Rphi1, double* d_Iphi1,
    int symmetry
) {
    (void)symmetry;
    dim3 blockDim(8, 8, 4);
    dim3 gridDim(
        (ex[0] + blockDim.x - 1) / blockDim.x,
        (ex[1] + blockDim.y - 1) / blockDim.y,
        (ex[2] + blockDim.z - 1) / blockDim.z
    );

    if constexpr (tetradtype == 0) {
        getnpem1_ss_kernel<0><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi1, d_Iphi1
        );
    } else if constexpr (tetradtype == 1) {
        getnpem1_ss_kernel<1><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi1, d_Iphi1
        );
    } else {
        getnpem1_ss_kernel<2><<<gridDim, blockDim, 0, stream>>>(
            ex[0], ex[1], ex[2],
            pex[0], pex[1],
            d_X, d_Y, d_Z,
            d_chi,
            d_dxx, d_gxy, d_gxz, d_dyy, d_gyz, d_dzz,
            d_Ex, d_Ey, d_Ez, d_Bx, d_By, d_Bz,
            d_Rphi1, d_Iphi1
        );
    }
}
