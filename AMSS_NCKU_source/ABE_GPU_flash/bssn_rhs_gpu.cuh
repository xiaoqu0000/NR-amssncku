#ifndef BSSN_RHS_GPU_CUH
#define BSSN_RHS_GPU_CUH

#include "fmisc_gpu.cuh"
#include "macrodef.fh"

#if ghost_width == 2
static constexpr int kBssnCenteredBoundaryRadius = 1;
static constexpr int kBssnCenteredLoadRadius = 1;
static constexpr int kBssnTransportBoundaryRadius = 2;
static constexpr int kBssnTransportLoadRadius = 2;
#elif ghost_width == 3
static constexpr int kBssnCenteredBoundaryRadius = 2;
static constexpr int kBssnCenteredLoadRadius = 2;
static constexpr int kBssnTransportBoundaryRadius = 3;
static constexpr int kBssnTransportLoadRadius = 3;
#elif ghost_width == 4
static constexpr int kBssnCenteredBoundaryRadius = 3;
static constexpr int kBssnCenteredLoadRadius = 3;
static constexpr int kBssnTransportBoundaryRadius = 4;
static constexpr int kBssnTransportLoadRadius = 5;
#elif ghost_width == 5
static constexpr int kBssnCenteredBoundaryRadius = 4;
static constexpr int kBssnCenteredLoadRadius = 4;
static constexpr int kBssnTransportBoundaryRadius = 5;
static constexpr int kBssnTransportLoadRadius = 5;
#else
#error Unsupported ghost_width in bssn_rhs_gpu.cuh; expected ghost_width == 2, 3, 4, or 5.
#endif

template<typename Real, int SMEM_Z, int SMEM_Y, int SMEM_X>
__device__ __forceinline__ void load_field_to_smem(
    Real (&smem)[SMEM_Z][SMEM_Y][SMEM_X], const Real* f,
    const int ex[3], int block_i, int block_j, int block_k,
    double SYM1, double SYM2, double SYM3
) {
    int tid = threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;
    double SoA[3] = {SYM1, SYM2, SYM3};
    const int tile_x = blockDim.x + 2 * kBssnCenteredLoadRadius;
    const int tile_y = blockDim.y + 2 * kBssnCenteredLoadRadius;
    const int tile_z = blockDim.z + 2 * kBssnCenteredLoadRadius;
    const int tile_xy = tile_x * tile_y;
    const int tile_size = tile_xy * tile_z;

    for (int idx = tid; idx < tile_size; idx += 256) {
        int loc_k = idx / tile_xy;
        int rem = idx % tile_xy;
        int loc_j = rem / tile_x;
        int loc_i = rem % tile_x;

        int glob_i = block_i + loc_i - kBssnCenteredLoadRadius;
        int glob_j = block_j + loc_j - kBssnCenteredLoadRadius;
        int glob_k = block_k + loc_k - kBssnCenteredLoadRadius;

        smem[loc_k][loc_j][loc_i] = static_cast<Real>(d_symmetry_bd_0b_t(
            kBssnCenteredBoundaryRadius,
            ex[0], ex[1], ex[2], f, glob_i, glob_j, glob_k, SoA[0], SoA[1], SoA[2]
        ));
    }
}

template<typename Real, int SMEM_Z, int SMEM_Y, int SMEM_X>
__device__ __forceinline__ void compute_all_derivs_smem(
    const Real (&smem)[SMEM_Z][SMEM_Y][SMEM_X],
    Real* fx, Real* fy, Real* fz,
    Real* fxx, Real* fxy, Real* fxz,
    Real* fyy, Real* fyz, Real* fzz,
    Real dX, Real dY, Real dZ,
    int i, int j, int k,
    int imin, int jmin, int kmin,
    int imax, int jmax, int kmax
) {
    *fx = Real(0.0); *fy = Real(0.0); *fz = Real(0.0);
    *fxx = Real(0.0); *fyy = Real(0.0); *fzz = Real(0.0);
    *fxy = Real(0.0); *fxz = Real(0.0); *fyz = Real(0.0);

    if (i >= imax || j >= jmax || k >= kmax) return;

    const Real d2dx = Real(1.0) / Real(2.0) / dX;
    const Real d2dy = Real(1.0) / Real(2.0) / dY;
    const Real d2dz = Real(1.0) / Real(2.0) / dZ;

    const Real Sdxdx = Real(1.0) / (dX * dX);
    const Real Sdydy = Real(1.0) / (dY * dY);
    const Real Sdzdz = Real(1.0) / (dZ * dZ);
    const Real Sdxdy = Real(0.25) / (dX * dY);
    const Real Sdxdz = Real(0.25) / (dX * dZ);
    const Real Sdydz = Real(0.25) / (dY * dZ);

#if ghost_width == 3 || ghost_width == 4 || ghost_width == 5
    const Real d12dx = Real(1.0) / Real(12.0) / dX;
    const Real d12dy = Real(1.0) / Real(12.0) / dY;
    const Real d12dz = Real(1.0) / Real(12.0) / dZ;
    const Real Fdxdx = (Real(1.0) / Real(12.0)) / (dX * dX);
    const Real Fdydy = (Real(1.0) / Real(12.0)) / (dY * dY);
    const Real Fdzdz = (Real(1.0) / Real(12.0)) / (dZ * dZ);
    const Real Fdxdy = (Real(1.0) / Real(144.0)) / (dX * dY);
    const Real Fdxdz = (Real(1.0) / Real(144.0)) / (dX * dZ);
    const Real Fdydz = (Real(1.0) / Real(144.0)) / (dY * dZ);
#endif

#if ghost_width == 4 || ghost_width == 5
    const Real d60dx = Real(1.0) / Real(60.0) / dX;
    const Real d60dy = Real(1.0) / Real(60.0) / dY;
    const Real d60dz = Real(1.0) / Real(60.0) / dZ;
    const Real Xdxdx = (Real(1.0) / Real(180.0)) / (dX * dX);
    const Real Xdydy = (Real(1.0) / Real(180.0)) / (dY * dY);
    const Real Xdzdz = (Real(1.0) / Real(180.0)) / (dZ * dZ);
    const Real Xdxdy = (Real(1.0) / Real(3600.0)) / (dX * dY);
    const Real Xdxdz = (Real(1.0) / Real(3600.0)) / (dX * dZ);
    const Real Xdydz = (Real(1.0) / Real(3600.0)) / (dY * dZ);
#endif
#if ghost_width == 5
    const Real d840dx = Real(1.0) / Real(840.0) / dX;
    const Real d840dy = Real(1.0) / Real(840.0) / dY;
    const Real d840dz = Real(1.0) / Real(840.0) / dZ;
    const Real Edxdx = (Real(1.0) / Real(5040.0)) / (dX * dX);
    const Real Edydy = (Real(1.0) / Real(5040.0)) / (dY * dY);
    const Real Edzdz = (Real(1.0) / Real(5040.0)) / (dZ * dZ);
    const Real Edxdy = (Real(1.0) / Real(705600.0)) / (dX * dY);
    const Real Edxdz = (Real(1.0) / Real(705600.0)) / (dX * dZ);
    const Real Edydz = (Real(1.0) / Real(705600.0)) / (dY * dZ);
#endif

    int c_i = threadIdx.x + kBssnCenteredLoadRadius;
    int c_j = threadIdx.y + kBssnCenteredLoadRadius;
    int c_k = threadIdx.z + kBssnCenteredLoadRadius;

    auto fh = [&](int di, int dj, int dk) -> Real {
        return smem[c_k + dk][c_j + dj][c_i + di];
    };

    bool fit_x_2nd = (i + 1 <= imax && i - 1 >= imin);
    bool fit_y_2nd = (j + 1 <= jmax && j - 1 >= jmin);
    bool fit_z_2nd = (k + 1 <= kmax && k - 1 >= kmin);

#if ghost_width == 2
    if (fit_x_2nd) *fx = d2dx * (-fh(-1,0,0) + fh(1,0,0));
    if (fit_y_2nd) *fy = d2dy * (-fh(0,-1,0) + fh(0,1,0));
    if (fit_z_2nd) *fz = d2dz * (-fh(0,0,-1) + fh(0,0,1));

    if (fit_x_2nd) *fxx = Sdxdx * (fh(-1,0,0) - 2.0*fh(0,0,0) + fh(1,0,0));
    if (fit_y_2nd) *fyy = Sdydy * (fh(0,-1,0) - 2.0*fh(0,0,0) + fh(0,1,0));
    if (fit_z_2nd) *fzz = Sdzdz * (fh(0,0,-1) - 2.0*fh(0,0,0) + fh(0,0,1));

    if (fit_x_2nd && fit_y_2nd) *fxy = Sdxdy * (fh(-1,-1,0) - fh(1,-1,0) - fh(-1,1,0) + fh(1,1,0));
    if (fit_x_2nd && fit_z_2nd) *fxz = Sdxdz * (fh(-1,0,-1) - fh(1,0,-1) - fh(-1,0,1) + fh(1,0,1));
    if (fit_y_2nd && fit_z_2nd) *fyz = Sdydz * (fh(0,-1,-1) - fh(0,1,-1) - fh(0,-1,1) + fh(0,1,1));
#elif ghost_width == 3
    bool fit_2nd = fit_x_2nd && fit_y_2nd && fit_z_2nd;
    bool fit_4th = (i + 2 <= imax && i - 2 >= imin &&
                    j + 2 <= jmax && j - 2 >= jmin &&
                    k + 2 <= kmax && k - 2 >= kmin);
    if (fit_4th) {
        *fx = d12dx * (fh(-2,0,0) - 8.0*fh(-1,0,0) + 8.0*fh(1,0,0) - fh(2,0,0));
        *fy = d12dy * (fh(0,-2,0) - 8.0*fh(0,-1,0) + 8.0*fh(0,1,0) - fh(0,2,0));
        *fz = d12dz * (fh(0,0,-2) - 8.0*fh(0,0,-1) + 8.0*fh(0,0,1) - fh(0,0,2));

        *fxx = Fdxdx * (-fh(-2,0,0) + 16.0*fh(-1,0,0) - 30.0*fh(0,0,0) - fh(2,0,0) + 16.0*fh(1,0,0));
        *fyy = Fdydy * (-fh(0,-2,0) + 16.0*fh(0,-1,0) - 30.0*fh(0,0,0) - fh(0,2,0) + 16.0*fh(0,1,0));
        *fzz = Fdzdz * (-fh(0,0,-2) + 16.0*fh(0,0,-1) - 30.0*fh(0,0,0) - fh(0,0,2) + 16.0*fh(0,0,1));

        *fxy = Fdxdy * ((fh(-2,-2,0) - 8.0*fh(-1,-2,0) + 8.0*fh(1,-2,0) - fh(2,-2,0))
                 -8.0*(fh(-2,-1,0) - 8.0*fh(-1,-1,0) + 8.0*fh(1,-1,0) - fh(2,-1,0))
                 +8.0*(fh(-2, 1,0) - 8.0*fh(-1, 1,0) + 8.0*fh(1, 1,0) - fh(2, 1,0))
                     -(fh(-2, 2,0) - 8.0*fh(-1, 2,0) + 8.0*fh(1, 2,0) - fh(2, 2,0)));
        *fxz = Fdxdz * ((fh(-2,0,-2) - 8.0*fh(-1,0,-2) + 8.0*fh(1,0,-2) - fh(2,0,-2))
                 -8.0*(fh(-2,0,-1) - 8.0*fh(-1,0,-1) + 8.0*fh(1,0,-1) - fh(2,0,-1))
                 +8.0*(fh(-2,0, 1) - 8.0*fh(-1,0, 1) + 8.0*fh(1,0, 1) - fh(2,0, 1))
                     -(fh(-2,0, 2) - 8.0*fh(-1,0, 2) + 8.0*fh(1,0, 2) - fh(2,0, 2)));
        *fyz = Fdydz * ((fh(0,-2,-2) - 8.0*fh(0,-1,-2) + 8.0*fh(0,1,-2) - fh(0,2,-2))
                 -8.0*(fh(0,-2,-1) - 8.0*fh(0,-1,-1) + 8.0*fh(0,1,-1) - fh(0,2,-1))
                 +8.0*(fh(0,-2, 1) - 8.0*fh(0,-1, 1) + 8.0*fh(0,1, 1) - fh(0,2, 1))
                     -(fh(0,-2, 2) - 8.0*fh(0,-1, 2) + 8.0*fh(0,1, 2) - fh(0,2, 2)));
    } else if (fit_2nd) {
        *fx = d2dx * (-fh(-1,0,0) + fh(1,0,0));
        *fy = d2dy * (-fh(0,-1,0) + fh(0,1,0));
        *fz = d2dz * (-fh(0,0,-1) + fh(0,0,1));
        *fxx = Sdxdx * (fh(-1,0,0) - 2.0*fh(0,0,0) + fh(1,0,0));
        *fyy = Sdydy * (fh(0,-1,0) - 2.0*fh(0,0,0) + fh(0,1,0));
        *fzz = Sdzdz * (fh(0,0,-1) - 2.0*fh(0,0,0) + fh(0,0,1));
        *fxy = Sdxdy * (fh(-1,-1,0) - fh(1,-1,0) - fh(-1,1,0) + fh(1,1,0));
        *fxz = Sdxdz * (fh(-1,0,-1) - fh(1,0,-1) - fh(-1,0,1) + fh(1,0,1));
        *fyz = Sdydz * (fh(0,-1,-1) - fh(0,1,-1) - fh(0,-1,1) + fh(0,1,1));
    }
#elif ghost_width == 4
    bool fit_x_4th = (i + 2 <= imax && i - 2 >= imin);
    bool fit_y_4th = (j + 2 <= jmax && j - 2 >= jmin);
    bool fit_z_4th = (k + 2 <= kmax && k - 2 >= kmin);
    bool fit_x_6th = (i + 3 <= imax && i - 3 >= imin);
    bool fit_y_6th = (j + 3 <= jmax && j - 3 >= jmin);
    bool fit_z_6th = (k + 3 <= kmax && k - 3 >= kmin);

    auto dx6 = [&](int dj, int dk) { return -fh(-3,dj,dk) + 9.0*fh(-2,dj,dk) - 45.0*fh(-1,dj,dk) + 45.0*fh(1,dj,dk) - 9.0*fh(2,dj,dk) + fh(3,dj,dk); };
    auto dy6 = [&](int di, int dk) { return -fh(di,-3,dk) + 9.0*fh(di,-2,dk) - 45.0*fh(di,-1,dk) + 45.0*fh(di,1,dk) - 9.0*fh(di,2,dk) + fh(di,3,dk); };
    auto dz6 = [&](int di, int dj) { return -fh(di,dj,-3) + 9.0*fh(di,dj,-2) - 45.0*fh(di,dj,-1) + 45.0*fh(di,dj,1) - 9.0*fh(di,dj,2) + fh(di,dj,3); };
    auto dx4 = [&](int dj, int dk) { return fh(-2,dj,dk) - 8.0*fh(-1,dj,dk) + 8.0*fh(1,dj,dk) - fh(2,dj,dk); };
    auto dy4 = [&](int di, int dk) { return fh(di,-2,dk) - 8.0*fh(di,-1,dk) + 8.0*fh(di,1,dk) - fh(di,2,dk); };

    if (fit_x_6th) *fx = d60dx * dx6(0,0);
    else if (fit_x_4th) *fx = d12dx * dx4(0,0);
    else if (fit_x_2nd) *fx = d2dx * (-fh(-1,0,0) + fh(1,0,0));
    if (fit_y_6th) *fy = d60dy * dy6(0,0);
    else if (fit_y_4th) *fy = d12dy * dy4(0,0);
    else if (fit_y_2nd) *fy = d2dy * (-fh(0,-1,0) + fh(0,1,0));
    if (fit_z_6th) *fz = d60dz * dz6(0,0);
    else if (fit_z_4th) *fz = d12dz * (fh(0,0,-2) - 8.0*fh(0,0,-1) + 8.0*fh(0,0,1) - fh(0,0,2));
    else if (fit_z_2nd) *fz = d2dz * (-fh(0,0,-1) + fh(0,0,1));

    if (fit_x_6th) *fxx = Xdxdx * (2.0*fh(-3,0,0) - 27.0*fh(-2,0,0) + 270.0*fh(-1,0,0) - 490.0*fh(0,0,0) + 270.0*fh(1,0,0) - 27.0*fh(2,0,0) + 2.0*fh(3,0,0));
    else if (fit_x_4th) *fxx = Fdxdx * (-fh(-2,0,0) + 16.0*fh(-1,0,0) - 30.0*fh(0,0,0) - fh(2,0,0) + 16.0*fh(1,0,0));
    else if (fit_x_2nd) *fxx = Sdxdx * (fh(-1,0,0) - 2.0*fh(0,0,0) + fh(1,0,0));
    if (fit_y_6th) *fyy = Xdydy * (2.0*fh(0,-3,0) - 27.0*fh(0,-2,0) + 270.0*fh(0,-1,0) - 490.0*fh(0,0,0) + 270.0*fh(0,1,0) - 27.0*fh(0,2,0) + 2.0*fh(0,3,0));
    else if (fit_y_4th) *fyy = Fdydy * (-fh(0,-2,0) + 16.0*fh(0,-1,0) - 30.0*fh(0,0,0) - fh(0,2,0) + 16.0*fh(0,1,0));
    else if (fit_y_2nd) *fyy = Sdydy * (fh(0,-1,0) - 2.0*fh(0,0,0) + fh(0,1,0));
    if (fit_z_6th) *fzz = Xdzdz * (2.0*fh(0,0,-3) - 27.0*fh(0,0,-2) + 270.0*fh(0,0,-1) - 490.0*fh(0,0,0) + 270.0*fh(0,0,1) - 27.0*fh(0,0,2) + 2.0*fh(0,0,3));
    else if (fit_z_4th) *fzz = Fdzdz * (-fh(0,0,-2) + 16.0*fh(0,0,-1) - 30.0*fh(0,0,0) - fh(0,0,2) + 16.0*fh(0,0,1));
    else if (fit_z_2nd) *fzz = Sdzdz * (fh(0,0,-1) - 2.0*fh(0,0,0) + fh(0,0,1));

    if (fit_x_6th && fit_y_6th) *fxy = Xdxdy * (-dx6(-3,0) + 9.0*dx6(-2,0) - 45.0*dx6(-1,0) + 45.0*dx6(1,0) - 9.0*dx6(2,0) + dx6(3,0));
    else if (fit_x_4th && fit_y_4th) *fxy = Fdxdy * (dy4(-2,0) - 8.0*dy4(-1,0) + 8.0*dy4(1,0) - dy4(2,0));
    else if (fit_x_2nd && fit_y_2nd) *fxy = Sdxdy * (fh(-1,-1,0) - fh(1,-1,0) - fh(-1,1,0) + fh(1,1,0));

    if (fit_x_6th && fit_z_6th) *fxz = Xdxdz * (-dx6(0,-3) + 9.0*dx6(0,-2) - 45.0*dx6(0,-1) + 45.0*dx6(0,1) - 9.0*dx6(0,2) + dx6(0,3));
    else if (fit_x_4th && fit_z_4th) *fxz = Fdxdz * ((fh(-2,0,-2) - 8.0*fh(-1,0,-2) + 8.0*fh(1,0,-2) - fh(2,0,-2)) - 8.0*(fh(-2,0,-1) - 8.0*fh(-1,0,-1) + 8.0*fh(1,0,-1) - fh(2,0,-1)) + 8.0*(fh(-2,0,1) - 8.0*fh(-1,0,1) + 8.0*fh(1,0,1) - fh(2,0,1)) - (fh(-2,0,2) - 8.0*fh(-1,0,2) + 8.0*fh(1,0,2) - fh(2,0,2)));
    else if (fit_x_2nd && fit_z_2nd) *fxz = Sdxdz * (fh(-1,0,-1) - fh(1,0,-1) - fh(-1,0,1) + fh(1,0,1));

    if (fit_y_6th && fit_z_6th) *fyz = Xdydz * (-dy6(0,-3) + 9.0*dy6(0,-2) - 45.0*dy6(0,-1) + 45.0*dy6(0,1) - 9.0*dy6(0,2) + dy6(0,3));
    else if (fit_y_4th && fit_z_4th) *fyz = Fdydz * ((fh(0,-2,-2) - 8.0*fh(0,-1,-2) + 8.0*fh(0,1,-2) - fh(0,2,-2)) - 8.0*(fh(0,-2,-1) - 8.0*fh(0,-1,-1) + 8.0*fh(0,1,-1) - fh(0,2,-1)) + 8.0*(fh(0,-2,1) - 8.0*fh(0,-1,1) + 8.0*fh(0,1,1) - fh(0,2,1)) - (fh(0,-2,2) - 8.0*fh(0,-1,2) + 8.0*fh(0,1,2) - fh(0,2,2)));
    else if (fit_y_2nd && fit_z_2nd) *fyz = Sdydz * (fh(0,-1,-1) - fh(0,1,-1) - fh(0,-1,1) + fh(0,1,1));
#elif ghost_width == 5
    bool fit_x_4th = (i + 2 <= imax && i - 2 >= imin);
    bool fit_y_4th = (j + 2 <= jmax && j - 2 >= jmin);
    bool fit_z_4th = (k + 2 <= kmax && k - 2 >= kmin);
    bool fit_x_6th = (i + 3 <= imax && i - 3 >= imin);
    bool fit_y_6th = (j + 3 <= jmax && j - 3 >= jmin);
    bool fit_z_6th = (k + 3 <= kmax && k - 3 >= kmin);
    bool fit_x_8th = (i + 4 <= imax && i - 4 >= imin);
    bool fit_y_8th = (j + 4 <= jmax && j - 4 >= jmin);
    bool fit_z_8th = (k + 4 <= kmax && k - 4 >= kmin);

    auto dx8 = [&](int dj, int dk) { return 3.0*fh(-4,dj,dk) - 32.0*fh(-3,dj,dk) + 168.0*fh(-2,dj,dk) - 672.0*fh(-1,dj,dk) + 672.0*fh(1,dj,dk) - 168.0*fh(2,dj,dk) + 32.0*fh(3,dj,dk) - 3.0*fh(4,dj,dk); };
    auto dy8 = [&](int di, int dk) { return 3.0*fh(di,-4,dk) - 32.0*fh(di,-3,dk) + 168.0*fh(di,-2,dk) - 672.0*fh(di,-1,dk) + 672.0*fh(di,1,dk) - 168.0*fh(di,2,dk) + 32.0*fh(di,3,dk) - 3.0*fh(di,4,dk); };
    auto dz8 = [&](int di, int dj) { return 3.0*fh(di,dj,-4) - 32.0*fh(di,dj,-3) + 168.0*fh(di,dj,-2) - 672.0*fh(di,dj,-1) + 672.0*fh(di,dj,1) - 168.0*fh(di,dj,2) + 32.0*fh(di,dj,3) - 3.0*fh(di,dj,4); };
    auto dx6 = [&](int dj, int dk) { return -fh(-3,dj,dk) + 9.0*fh(-2,dj,dk) - 45.0*fh(-1,dj,dk) + 45.0*fh(1,dj,dk) - 9.0*fh(2,dj,dk) + fh(3,dj,dk); };
    auto dy6 = [&](int di, int dk) { return -fh(di,-3,dk) + 9.0*fh(di,-2,dk) - 45.0*fh(di,-1,dk) + 45.0*fh(di,1,dk) - 9.0*fh(di,2,dk) + fh(di,3,dk); };
    auto dz6 = [&](int di, int dj) { return -fh(di,dj,-3) + 9.0*fh(di,dj,-2) - 45.0*fh(di,dj,-1) + 45.0*fh(di,dj,1) - 9.0*fh(di,dj,2) + fh(di,dj,3); };
    auto dx4 = [&](int dj, int dk) { return fh(-2,dj,dk) - 8.0*fh(-1,dj,dk) + 8.0*fh(1,dj,dk) - fh(2,dj,dk); };
    auto dy4 = [&](int di, int dk) { return fh(di,-2,dk) - 8.0*fh(di,-1,dk) + 8.0*fh(di,1,dk) - fh(di,2,dk); };
    auto dz4 = [&](int di, int dj) { return fh(di,dj,-2) - 8.0*fh(di,dj,-1) + 8.0*fh(di,dj,1) - fh(di,dj,2); };

    if (fit_x_8th) *fx = d840dx * dx8(0,0);
    else if (fit_x_6th) *fx = d60dx * dx6(0,0);
    else if (fit_x_4th) *fx = d12dx * dx4(0,0);
    else if (fit_x_2nd) *fx = d2dx * (-fh(-1,0,0) + fh(1,0,0));
    if (fit_y_8th) *fy = d840dy * dy8(0,0);
    else if (fit_y_6th) *fy = d60dy * dy6(0,0);
    else if (fit_y_4th) *fy = d12dy * dy4(0,0);
    else if (fit_y_2nd) *fy = d2dy * (-fh(0,-1,0) + fh(0,1,0));
    if (fit_z_8th) *fz = d840dz * dz8(0,0);
    else if (fit_z_6th) *fz = d60dz * dz6(0,0);
    else if (fit_z_4th) *fz = d12dz * dz4(0,0);
    else if (fit_z_2nd) *fz = d2dz * (-fh(0,0,-1) + fh(0,0,1));

    if (fit_x_8th) *fxx = Edxdx * (-9.0*fh(-4,0,0) + 128.0*fh(-3,0,0) - 1008.0*fh(-2,0,0) + 8064.0*fh(-1,0,0) - 14350.0*fh(0,0,0) + 8064.0*fh(1,0,0) - 1008.0*fh(2,0,0) + 128.0*fh(3,0,0) - 9.0*fh(4,0,0));
    else if (fit_x_6th) *fxx = Xdxdx * (2.0*fh(-3,0,0) - 27.0*fh(-2,0,0) + 270.0*fh(-1,0,0) - 490.0*fh(0,0,0) + 270.0*fh(1,0,0) - 27.0*fh(2,0,0) + 2.0*fh(3,0,0));
    else if (fit_x_4th) *fxx = Fdxdx * (-fh(-2,0,0) + 16.0*fh(-1,0,0) - 30.0*fh(0,0,0) - fh(2,0,0) + 16.0*fh(1,0,0));
    else if (fit_x_2nd) *fxx = Sdxdx * (fh(-1,0,0) - 2.0*fh(0,0,0) + fh(1,0,0));
    if (fit_y_8th) *fyy = Edydy * (-9.0*fh(0,-4,0) + 128.0*fh(0,-3,0) - 1008.0*fh(0,-2,0) + 8064.0*fh(0,-1,0) - 14350.0*fh(0,0,0) + 8064.0*fh(0,1,0) - 1008.0*fh(0,2,0) + 128.0*fh(0,3,0) - 9.0*fh(0,4,0));
    else if (fit_y_6th) *fyy = Xdydy * (2.0*fh(0,-3,0) - 27.0*fh(0,-2,0) + 270.0*fh(0,-1,0) - 490.0*fh(0,0,0) + 270.0*fh(0,1,0) - 27.0*fh(0,2,0) + 2.0*fh(0,3,0));
    else if (fit_y_4th) *fyy = Fdydy * (-fh(0,-2,0) + 16.0*fh(0,-1,0) - 30.0*fh(0,0,0) - fh(0,2,0) + 16.0*fh(0,1,0));
    else if (fit_y_2nd) *fyy = Sdydy * (fh(0,-1,0) - 2.0*fh(0,0,0) + fh(0,1,0));
    if (fit_z_8th) *fzz = Edzdz * (-9.0*fh(0,0,-4) + 128.0*fh(0,0,-3) - 1008.0*fh(0,0,-2) + 8064.0*fh(0,0,-1) - 14350.0*fh(0,0,0) + 8064.0*fh(0,0,1) - 1008.0*fh(0,0,2) + 128.0*fh(0,0,3) - 9.0*fh(0,0,4));
    else if (fit_z_6th) *fzz = Xdzdz * (2.0*fh(0,0,-3) - 27.0*fh(0,0,-2) + 270.0*fh(0,0,-1) - 490.0*fh(0,0,0) + 270.0*fh(0,0,1) - 27.0*fh(0,0,2) + 2.0*fh(0,0,3));
    else if (fit_z_4th) *fzz = Fdzdz * (-fh(0,0,-2) + 16.0*fh(0,0,-1) - 30.0*fh(0,0,0) - fh(0,0,2) + 16.0*fh(0,0,1));
    else if (fit_z_2nd) *fzz = Sdzdz * (fh(0,0,-1) - 2.0*fh(0,0,0) + fh(0,0,1));

    if (fit_x_8th && fit_y_8th) *fxy = Edxdy * (3.0*dx8(-4,0) - 32.0*dx8(-3,0) + 168.0*dx8(-2,0) - 672.0*dx8(-1,0) + 672.0*dx8(1,0) - 168.0*dx8(2,0) + 32.0*dx8(3,0) - 3.0*dx8(4,0));
    else if (fit_x_6th && fit_y_6th) *fxy = Xdxdy * (-dx6(-3,0) + 9.0*dx6(-2,0) - 45.0*dx6(-1,0) + 45.0*dx6(1,0) - 9.0*dx6(2,0) + dx6(3,0));
    else if (fit_x_4th && fit_y_4th) *fxy = Fdxdy * (dy4(-2,0) - 8.0*dy4(-1,0) + 8.0*dy4(1,0) - dy4(2,0));
    else if (fit_x_2nd && fit_y_2nd) *fxy = Sdxdy * (fh(-1,-1,0) - fh(1,-1,0) - fh(-1,1,0) + fh(1,1,0));

    if (fit_x_8th && fit_z_8th) *fxz = Edxdz * (3.0*dx8(0,-4) - 32.0*dx8(0,-3) + 168.0*dx8(0,-2) - 672.0*dx8(0,-1) + 672.0*dx8(0,1) - 168.0*dx8(0,2) + 32.0*dx8(0,3) - 3.0*dx8(0,4));
    else if (fit_x_6th && fit_z_6th) *fxz = Xdxdz * (-dx6(0,-3) + 9.0*dx6(0,-2) - 45.0*dx6(0,-1) + 45.0*dx6(0,1) - 9.0*dx6(0,2) + dx6(0,3));
    else if (fit_x_4th && fit_z_4th) *fxz = Fdxdz * (dz4(-2,0) - 8.0*dz4(-1,0) + 8.0*dz4(1,0) - dz4(2,0));
    else if (fit_x_2nd && fit_z_2nd) *fxz = Sdxdz * (fh(-1,0,-1) - fh(1,0,-1) - fh(-1,0,1) + fh(1,0,1));

    if (fit_y_8th && fit_z_8th) *fyz = Edydz * (3.0*dy8(0,-4) - 32.0*dy8(0,-3) + 168.0*dy8(0,-2) - 672.0*dy8(0,-1) + 672.0*dy8(0,1) - 168.0*dy8(0,2) + 32.0*dy8(0,3) - 3.0*dy8(0,4));
    else if (fit_y_6th && fit_z_6th) *fyz = Xdydz * (-dy6(0,-3) + 9.0*dy6(0,-2) - 45.0*dy6(0,-1) + 45.0*dy6(0,1) - 9.0*dy6(0,2) + dy6(0,3));
    else if (fit_y_4th && fit_z_4th) *fyz = Fdydz * (dz4(0,-2) - 8.0*dz4(0,-1) + 8.0*dz4(0,1) - dz4(0,2));
    else if (fit_y_2nd && fit_z_2nd) *fyz = Sdydz * (fh(0,-1,-1) - fh(0,1,-1) - fh(0,-1,1) + fh(0,1,1));
#endif
}

template<typename Real, int SMEM_Z, int SMEM_Y, int SMEM_X>
__device__ __forceinline__ void load_field_to_smem_rad3(
    Real (&smem)[SMEM_Z][SMEM_Y][SMEM_X], const Real* f,
    const int ex[3], int block_i, int block_j, int block_k,
    double SYM1, double SYM2, double SYM3
) {
    int tid = threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;
    double SoA[3] = {SYM1, SYM2, SYM3};
    const int tile_x = blockDim.x + 2 * kBssnTransportLoadRadius;
    const int tile_y = blockDim.y + 2 * kBssnTransportLoadRadius;
    const int tile_z = blockDim.z + 2 * kBssnTransportLoadRadius;
    const int tile_xy = tile_x * tile_y;
    const int tile_size = tile_xy * tile_z;

    for (int idx = tid; idx < tile_size; idx += 256) {
        int loc_k = idx / tile_xy;
        int rem = idx % tile_xy;
        int loc_j = rem / tile_x;
        int loc_i = rem % tile_x;

        int glob_i = block_i + loc_i - kBssnTransportLoadRadius;
        int glob_j = block_j + loc_j - kBssnTransportLoadRadius;
        int glob_k = block_k + loc_k - kBssnTransportLoadRadius;

        smem[loc_k][loc_j][loc_i] = static_cast<Real>(d_symmetry_bd_0b_t(
            kBssnTransportBoundaryRadius,
            ex[0], ex[1], ex[2], f, glob_i, glob_j, glob_k, SoA[0], SoA[1], SoA[2]
        ));
    }
}

template<typename Real, int SMEM_Z, int SMEM_Y, int SMEM_X>
__device__ __forceinline__ void compute_first_derivs_smem(
    const Real (&smem)[SMEM_Z][SMEM_Y][SMEM_X],
    Real* fx, Real* fy, Real* fz,
    Real dX, Real dY, Real dZ,
    int i, int j, int k,
    int imin, int jmin, int kmin,
    int imax, int jmax, int kmax
) {
    *fx = Real(0.0); *fy = Real(0.0); *fz = Real(0.0);

    if (i >= imax || j >= jmax || k >= kmax) return;

    const Real d2dx = Real(1.0) / Real(2.0) / dX;
    const Real d2dy = Real(1.0) / Real(2.0) / dY;
    const Real d2dz = Real(1.0) / Real(2.0) / dZ;
#if ghost_width == 3 || ghost_width == 4 || ghost_width == 5
    const Real d12dx = Real(1.0) / Real(12.0) / dX;
    const Real d12dy = Real(1.0) / Real(12.0) / dY;
    const Real d12dz = Real(1.0) / Real(12.0) / dZ;
#endif
#if ghost_width == 4 || ghost_width == 5
    const Real d60dx = Real(1.0) / Real(60.0) / dX;
    const Real d60dy = Real(1.0) / Real(60.0) / dY;
    const Real d60dz = Real(1.0) / Real(60.0) / dZ;
#endif
#if ghost_width == 5
    const Real d840dx = Real(1.0) / Real(840.0) / dX;
    const Real d840dy = Real(1.0) / Real(840.0) / dY;
    const Real d840dz = Real(1.0) / Real(840.0) / dZ;
#endif

    int c_i = threadIdx.x + kBssnCenteredLoadRadius;
    int c_j = threadIdx.y + kBssnCenteredLoadRadius;
    int c_k = threadIdx.z + kBssnCenteredLoadRadius;

    auto fh = [&](int di, int dj, int dk) -> Real {
        return smem[c_k + dk][c_j + dj][c_i + di];
    };

    bool fit_x_2nd = (i + 1 <= imax && i - 1 >= imin);
    bool fit_y_2nd = (j + 1 <= jmax && j - 1 >= jmin);
    bool fit_z_2nd = (k + 1 <= kmax && k - 1 >= kmin);

#if ghost_width == 2
    if (fit_x_2nd) *fx = d2dx * (-fh(-1,0,0) + fh(1,0,0));
    if (fit_y_2nd) *fy = d2dy * (-fh(0,-1,0) + fh(0,1,0));
    if (fit_z_2nd) *fz = d2dz * (-fh(0,0,-1) + fh(0,0,1));
#elif ghost_width == 3
    bool fit_2nd = fit_x_2nd && fit_y_2nd && fit_z_2nd;
    bool fit_4th = (i + 2 <= imax && i - 2 >= imin &&
                    j + 2 <= jmax && j - 2 >= jmin &&
                    k + 2 <= kmax && k - 2 >= kmin);
    if (fit_4th) {
        *fx = d12dx * (fh(-2,0,0) - 8.0*fh(-1,0,0) + 8.0*fh(1,0,0) - fh(2,0,0));
        *fy = d12dy * (fh(0,-2,0) - 8.0*fh(0,-1,0) + 8.0*fh(0,1,0) - fh(0,2,0));
        *fz = d12dz * (fh(0,0,-2) - 8.0*fh(0,0,-1) + 8.0*fh(0,0,1) - fh(0,0,2));
    } else if (fit_2nd) {
        *fx = d2dx * (-fh(-1,0,0) + fh(1,0,0));
        *fy = d2dy * (-fh(0,-1,0) + fh(0,1,0));
        *fz = d2dz * (-fh(0,0,-1) + fh(0,0,1));
    }
#elif ghost_width == 4
    bool fit_x_4th = (i + 2 <= imax && i - 2 >= imin);
    bool fit_y_4th = (j + 2 <= jmax && j - 2 >= jmin);
    bool fit_z_4th = (k + 2 <= kmax && k - 2 >= kmin);
    bool fit_x_6th = (i + 3 <= imax && i - 3 >= imin);
    bool fit_y_6th = (j + 3 <= jmax && j - 3 >= jmin);
    bool fit_z_6th = (k + 3 <= kmax && k - 3 >= kmin);

    if (fit_x_6th) *fx = d60dx * (-fh(-3,0,0) + 9.0*fh(-2,0,0) - 45.0*fh(-1,0,0) + 45.0*fh(1,0,0) - 9.0*fh(2,0,0) + fh(3,0,0));
    else if (fit_x_4th) *fx = d12dx * (fh(-2,0,0) - 8.0*fh(-1,0,0) + 8.0*fh(1,0,0) - fh(2,0,0));
    else if (fit_x_2nd) *fx = d2dx * (-fh(-1,0,0) + fh(1,0,0));

    if (fit_y_6th) *fy = d60dy * (-fh(0,-3,0) + 9.0*fh(0,-2,0) - 45.0*fh(0,-1,0) + 45.0*fh(0,1,0) - 9.0*fh(0,2,0) + fh(0,3,0));
    else if (fit_y_4th) *fy = d12dy * (fh(0,-2,0) - 8.0*fh(0,-1,0) + 8.0*fh(0,1,0) - fh(0,2,0));
    else if (fit_y_2nd) *fy = d2dy * (-fh(0,-1,0) + fh(0,1,0));

    if (fit_z_6th) *fz = d60dz * (-fh(0,0,-3) + 9.0*fh(0,0,-2) - 45.0*fh(0,0,-1) + 45.0*fh(0,0,1) - 9.0*fh(0,0,2) + fh(0,0,3));
    else if (fit_z_4th) *fz = d12dz * (fh(0,0,-2) - 8.0*fh(0,0,-1) + 8.0*fh(0,0,1) - fh(0,0,2));
    else if (fit_z_2nd) *fz = d2dz * (-fh(0,0,-1) + fh(0,0,1));
#elif ghost_width == 5
    bool fit_x_4th = (i + 2 <= imax && i - 2 >= imin);
    bool fit_y_4th = (j + 2 <= jmax && j - 2 >= jmin);
    bool fit_z_4th = (k + 2 <= kmax && k - 2 >= kmin);
    bool fit_x_6th = (i + 3 <= imax && i - 3 >= imin);
    bool fit_y_6th = (j + 3 <= jmax && j - 3 >= jmin);
    bool fit_z_6th = (k + 3 <= kmax && k - 3 >= kmin);
    bool fit_x_8th = (i + 4 <= imax && i - 4 >= imin);
    bool fit_y_8th = (j + 4 <= jmax && j - 4 >= jmin);
    bool fit_z_8th = (k + 4 <= kmax && k - 4 >= kmin);

    if (fit_x_8th) *fx = d840dx * (3.0*fh(-4,0,0) - 32.0*fh(-3,0,0) + 168.0*fh(-2,0,0) - 672.0*fh(-1,0,0) + 672.0*fh(1,0,0) - 168.0*fh(2,0,0) + 32.0*fh(3,0,0) - 3.0*fh(4,0,0));
    else if (fit_x_6th) *fx = d60dx * (-fh(-3,0,0) + 9.0*fh(-2,0,0) - 45.0*fh(-1,0,0) + 45.0*fh(1,0,0) - 9.0*fh(2,0,0) + fh(3,0,0));
    else if (fit_x_4th) *fx = d12dx * (fh(-2,0,0) - 8.0*fh(-1,0,0) + 8.0*fh(1,0,0) - fh(2,0,0));
    else if (fit_x_2nd) *fx = d2dx * (-fh(-1,0,0) + fh(1,0,0));

    if (fit_y_8th) *fy = d840dy * (3.0*fh(0,-4,0) - 32.0*fh(0,-3,0) + 168.0*fh(0,-2,0) - 672.0*fh(0,-1,0) + 672.0*fh(0,1,0) - 168.0*fh(0,2,0) + 32.0*fh(0,3,0) - 3.0*fh(0,4,0));
    else if (fit_y_6th) *fy = d60dy * (-fh(0,-3,0) + 9.0*fh(0,-2,0) - 45.0*fh(0,-1,0) + 45.0*fh(0,1,0) - 9.0*fh(0,2,0) + fh(0,3,0));
    else if (fit_y_4th) *fy = d12dy * (fh(0,-2,0) - 8.0*fh(0,-1,0) + 8.0*fh(0,1,0) - fh(0,2,0));
    else if (fit_y_2nd) *fy = d2dy * (-fh(0,-1,0) + fh(0,1,0));

    if (fit_z_8th) *fz = d840dz * (3.0*fh(0,0,-4) - 32.0*fh(0,0,-3) + 168.0*fh(0,0,-2) - 672.0*fh(0,0,-1) + 672.0*fh(0,0,1) - 168.0*fh(0,0,2) + 32.0*fh(0,0,3) - 3.0*fh(0,0,4));
    else if (fit_z_6th) *fz = d60dz * (-fh(0,0,-3) + 9.0*fh(0,0,-2) - 45.0*fh(0,0,-1) + 45.0*fh(0,0,1) - 9.0*fh(0,0,2) + fh(0,0,3));
    else if (fit_z_4th) *fz = d12dz * (fh(0,0,-2) - 8.0*fh(0,0,-1) + 8.0*fh(0,0,1) - fh(0,0,2));
    else if (fit_z_2nd) *fz = d2dz * (-fh(0,0,-1) + fh(0,0,1));
#endif
}

void gpu_compute_constraint_fr_launch(
    cudaStream_t &stream,
    int* ex,
    double* d_chi,
    double* d_trK,
    double* d_rho,
    double* d_Sphi,
    double* d_dxx,
    double* d_gxy,
    double* d_gxz,
    double* d_dyy,
    double* d_gyz,
    double* d_dzz,
    double* d_Axx,
    double* d_Axy,
    double* d_Axz,
    double* d_Ayy,
    double* d_Ayz,
    double* d_Azz,
    double* d_Rxx,
    double* d_Rxy,
    double* d_Rxz,
    double* d_Ryy,
    double* d_Ryz,
    double* d_Rzz,
    double* d_Sxx,
    double* d_Sxy,
    double* d_Sxz,
    double* d_Syy,
    double* d_Syz,
    double* d_Szz,
    double* d_Cons_fR,
    double a2
);

#endif // BSSN_RHS_GPU_CUH
