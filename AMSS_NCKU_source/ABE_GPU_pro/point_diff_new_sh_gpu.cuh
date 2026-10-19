#ifndef POINT_DIFF_NEW_SH_GPU_CUH
#define POINT_DIFF_NEW_SH_GPU_CUH

#include "diff_new_sh_gpu.cuh"
#include "macrodef.fh"

// ============================================================
// Single-direction first derivatives (templated on GW)
// Dedicated implementations that compute only one direction.
// ============================================================

// --- df/dx ---
template<int GW>
__device__ __forceinline__ void d_point_fdx_sh_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fx,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW - 1;
    const double dX = X[1] - X[0];
    const int imax = ex[0] - 1;

    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SoA0, SoA1);
    };

    int imin = 0;
    if (symmetry == 2) { // OCTANT
        if (fabs(X[0]) < dX) imin = -ORD;
    }

    *fx = 0.0;

    do {
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin) {
                const double d840dx = 1.0/840.0/dX;
                *fx = d840dx*(  3*fh(i-4,j,k) - 32*fh(i-3,j,k) + 168*fh(i-2,j,k) - 672*fh(i-1,j,k)
                              + 672*fh(i+1,j,k) - 168*fh(i+2,j,k) + 32*fh(i+3,j,k) - 3*fh(i+4,j,k));
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (i+3 <= imax && i-3 >= imin) {
                const double d60dx = 1.0/60.0/dX;
                *fx = d60dx*(-fh(i-3,j,k) + 9*fh(i-2,j,k) - 45*fh(i-1,j,k)
                             + 45*fh(i+1,j,k) - 9*fh(i+2,j,k) + fh(i+3,j,k));
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (i+2 <= imax && i-2 >= imin) {
                const double d12dx = 1.0/12.0/dX;
                *fx = d12dx*(fh(i-2,j,k) - 8*fh(i-1,j,k) + 8*fh(i+1,j,k) - fh(i+2,j,k));
                break;
            }
        }
        if (i+1 <= imax && i-1 >= imin) {
            const double d2dx = 0.5/dX;
            *fx = d2dx*(-fh(i-1,j,k) + fh(i+1,j,k));
        }
    } while(false);

    (void)onoff;
}

// --- df/dy ---
template<int GW>
__device__ __forceinline__ void d_point_fdy_sh_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fy,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW - 1;
    const double dY = Y[1] - Y[0];
    int jmax = ex[1] - 1;

    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SoA0, SoA1);
    };

    int jmin = 0;
    if (symmetry == 2) { // OCTANT
        if (fabs(Y[0]) < dY) jmin = -ORD;
    } else if (symmetry == 1) { // EQ_SYMM
        if ((sst == 2 || sst == 4) && fabs(Y[0]) < dY) jmin = -ORD;
        if ((sst == 3 || sst == 5) && fabs(Y[ex[1]-1]) < dY) jmax = ex[1] - 1 + ORD;
    }

    *fy = 0.0;

    do {
        if constexpr (GW >= 5) {
            if (j+4 <= jmax && j-4 >= jmin) {
                const double d840dy = 1.0/840.0/dY;
                *fy = d840dy*(  3*fh(i,j-4,k) - 32*fh(i,j-3,k) + 168*fh(i,j-2,k) - 672*fh(i,j-1,k)
                              + 672*fh(i,j+1,k) - 168*fh(i,j+2,k) + 32*fh(i,j+3,k) - 3*fh(i,j+4,k));
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (j+3 <= jmax && j-3 >= jmin) {
                const double d60dy = 1.0/60.0/dY;
                *fy = d60dy*(-fh(i,j-3,k) + 9*fh(i,j-2,k) - 45*fh(i,j-1,k)
                             + 45*fh(i,j+1,k) - 9*fh(i,j+2,k) + fh(i,j+3,k));
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (j+2 <= jmax && j-2 >= jmin) {
                const double d12dy = 1.0/12.0/dY;
                *fy = d12dy*(fh(i,j-2,k) - 8*fh(i,j-1,k) + 8*fh(i,j+1,k) - fh(i,j+2,k));
                break;
            }
        }
        if (j+1 <= jmax && j-1 >= jmin) {
            const double d2dy = 0.5/dY;
            *fy = d2dy*(-fh(i,j-1,k) + fh(i,j+1,k));
        }
    } while(false);

    (void)onoff;
}

// --- df/dz ---
template<int GW>
__device__ __forceinline__ void d_point_fdz_sh_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW - 1;
    const double dZ = Z[1] - Z[0];
    const int kmax = ex[2] - 1;

    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SoA0, SoA1);
    };

    // z has no symmetry: kmin=0 always

    *fz = 0.0;

    do {
        if constexpr (GW >= 5) {
            if (k+4 <= kmax && k-4 >= 0) {
                const double d840dz = 1.0/840.0/dZ;
                *fz = d840dz*(  3*fh(i,j,k-4) - 32*fh(i,j,k-3) + 168*fh(i,j,k-2) - 672*fh(i,j,k-1)
                              + 672*fh(i,j,k+1) - 168*fh(i,j,k+2) + 32*fh(i,j,k+3) - 3*fh(i,j,k+4));
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (k+3 <= kmax && k-3 >= 0) {
                const double d60dz = 1.0/60.0/dZ;
                *fz = d60dz*(-fh(i,j,k-3) + 9*fh(i,j,k-2) - 45*fh(i,j,k-1)
                             + 45*fh(i,j,k+1) - 9*fh(i,j,k+2) + fh(i,j,k+3));
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (k+2 <= kmax && k-2 >= 0) {
                const double d12dz = 1.0/12.0/dZ;
                *fz = d12dz*(fh(i,j,k-2) - 8*fh(i,j,k-1) + 8*fh(i,j,k+1) - fh(i,j,k+2));
                break;
            }
        }
        if (k+1 <= kmax && k-1 >= 0) {
            const double d2dz = 0.5/dZ;
            *fz = d2dz*(-fh(i,j,k-1) + fh(i,j,k+1));
        }
    } while(false);

    (void)onoff;
}

// ============================================================
// Single-direction second derivatives (templated on GW)
// Dedicated implementations that compute only one component.
// ============================================================

// --- d^2f/dx^2 ---
template<int GW>
__device__ __forceinline__ void d_point_fddxx_sh_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fxx,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW - 1;
    const double dX = X[1] - X[0];
    const int imax = ex[0] - 1;

    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SoA0, SoA1);
    };

    int imin = 0;
    if (symmetry == 2) { // OCTANT
        if (fabs(X[0]) < dX) imin = -ORD;
    }

    *fxx = 0.0;

    do {
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin) {
                const double Edxdx = 1.0/5040.0/(dX*dX);
                *fxx = Edxdx*(- 9*fh(i-4,j,k)+128*fh(i-3,j,k)-1008*fh(i-2,j,k)+8064*fh(i-1,j,k)-14350*fh(i,j,k)
                              +8064*fh(i+1,j,k)-1008*fh(i+2,j,k)+128*fh(i+3,j,k)- 9*fh(i+4,j,k));
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (i+3 <= imax && i-3 >= imin) {
                const double Xdxdx = 1.0/180.0/(dX*dX);
                *fxx = Xdxdx*(2*fh(i-3,j,k)-27*fh(i-2,j,k)+270*fh(i-1,j,k)-490*fh(i,j,k)
                              +270*fh(i+1,j,k)-27*fh(i+2,j,k)+2*fh(i+3,j,k));
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (i+2 <= imax && i-2 >= imin) {
                const double Fdxdx = 1.0/12.0/(dX*dX);
                *fxx = Fdxdx*(-fh(i-2,j,k)+16*fh(i-1,j,k)-30*fh(i,j,k)
                              -fh(i+2,j,k)+16*fh(i+1,j,k));
                break;
            }
        }
        if (i+1 <= imax && i-1 >= imin) {
            const double Sdxdx = 1.0/(dX*dX);
            *fxx = Sdxdx*(fh(i-1,j,k) - 2*fh(i,j,k) + fh(i+1,j,k));
        }
    } while(false);

    (void)onoff;
}

// --- d^2f/dy^2 ---
template<int GW>
__device__ __forceinline__ void d_point_fddyy_sh_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fyy,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW - 1;
    const double dY = Y[1] - Y[0];
    int jmax = ex[1] - 1;

    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SoA0, SoA1);
    };

    int jmin = 0;
    if (symmetry == 2) { // OCTANT
        if (fabs(Y[0]) < dY) jmin = -ORD;
    } else if (symmetry == 1) { // EQ_SYMM
        if ((sst == 2 || sst == 4) && fabs(Y[0]) < dY) jmin = -ORD;
        if ((sst == 3 || sst == 5) && fabs(Y[ex[1]-1]) < dY) jmax = ex[1] - 1 + ORD;
    }

    *fyy = 0.0;

    do {
        if constexpr (GW >= 5) {
            if (j+4 <= jmax && j-4 >= jmin) {
                const double Edydy = 1.0/5040.0/(dY*dY);
                *fyy = Edydy*(- 9*fh(i,j-4,k)+128*fh(i,j-3,k)-1008*fh(i,j-2,k)+8064*fh(i,j-1,k)-14350*fh(i,j,k)
                              +8064*fh(i,j+1,k)-1008*fh(i,j+2,k)+128*fh(i,j+3,k)- 9*fh(i,j+4,k));
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (j+3 <= jmax && j-3 >= jmin) {
                const double Xdydy = 1.0/180.0/(dY*dY);
                *fyy = Xdydy*(2*fh(i,j-3,k)-27*fh(i,j-2,k)+270*fh(i,j-1,k)-490*fh(i,j,k)
                              +270*fh(i,j+1,k)-27*fh(i,j+2,k)+2*fh(i,j+3,k));
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (j+2 <= jmax && j-2 >= jmin) {
                const double Fdydy = 1.0/12.0/(dY*dY);
                *fyy = Fdydy*(-fh(i,j-2,k)+16*fh(i,j-1,k)-30*fh(i,j,k)
                              -fh(i,j+2,k)+16*fh(i,j+1,k));
                break;
            }
        }
        if (j+1 <= jmax && j-1 >= jmin) {
            const double Sdydy = 1.0/(dY*dY);
            *fyy = Sdydy*(fh(i,j-1,k) - 2*fh(i,j,k) + fh(i,j+1,k));
        }
    } while(false);

    (void)onoff;
}

// --- d^2f/dz^2 ---
template<int GW>
__device__ __forceinline__ void d_point_fddzz_sh_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW - 1;
    const double dZ = Z[1] - Z[0];
    const int kmax = ex[2] - 1;

    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SoA0, SoA1);
    };

    // z has no symmetry: kmin=0 always

    *fzz = 0.0;

    do {
        if constexpr (GW >= 5) {
            if (k+4 <= kmax && k-4 >= 0) {
                const double Edzdz = 1.0/5040.0/(dZ*dZ);
                *fzz = Edzdz*(- 9*fh(i,j,k-4)+128*fh(i,j,k-3)-1008*fh(i,j,k-2)+8064*fh(i,j,k-1)-14350*fh(i,j,k)
                              +8064*fh(i,j,k+1)-1008*fh(i,j,k+2)+128*fh(i,j,k+3)- 9*fh(i,j,k+4));
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (k+3 <= kmax && k-3 >= 0) {
                const double Xdzdz = 1.0/180.0/(dZ*dZ);
                *fzz = Xdzdz*(2*fh(i,j,k-3)-27*fh(i,j,k-2)+270*fh(i,j,k-1)-490*fh(i,j,k)
                              +270*fh(i,j,k+1)-27*fh(i,j,k+2)+2*fh(i,j,k+3));
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (k+2 <= kmax && k-2 >= 0) {
                const double Fdzdz = 1.0/12.0/(dZ*dZ);
                *fzz = Fdzdz*(-fh(i,j,k-2)+16*fh(i,j,k-1)-30*fh(i,j,k)
                              -fh(i,j,k+2)+16*fh(i,j,k+1));
                break;
            }
        }
        if (k+1 <= kmax && k-1 >= 0) {
            const double Sdzdz = 1.0/(dZ*dZ);
            *fzz = Sdzdz*(fh(i,j,k-1) - 2*fh(i,j,k) + fh(i,j,k+1));
        }
    } while(false);

    (void)onoff;
}

// --- d^2f/dxdy ---
template<int GW>
__device__ __forceinline__ void d_point_fddxy_sh_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fxy,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW - 1;
    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const int imax = ex[0] - 1;
    int jmax = ex[1] - 1;

    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SoA0, SoA1);
    };

    int imin = 0, jmin = 0;
    if (symmetry == 2) { // OCTANT
        if (fabs(X[0]) < dX) imin = -ORD;
        if (fabs(Y[0]) < dY) jmin = -ORD;
    } else if (symmetry == 1) { // EQ_SYMM
        if ((sst == 2 || sst == 4) && fabs(Y[0]) < dY) jmin = -ORD;
        if ((sst == 3 || sst == 5) && fabs(Y[ex[1]-1]) < dY) jmax = ex[1] - 1 + ORD;
    }

    *fxy = 0.0;

    do {
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin && j+4 <= jmax && j-4 >= jmin) {
                double c[8]; int o[8]; double fd1_denom_;
                diff_new_detail::get_fd1_8(c, o, fd1_denom_);
                const double d2 = fd1_denom_ * fd1_denom_;
                *fxy = diff_new_detail::cross_deriv<8, 0, 1>(fh, i, j, k, c, o, d2, 1.0/(dX*dY));
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (i+3 <= imax && i-3 >= imin && j+3 <= jmax && j-3 >= jmin) {
                double c[6]; int o[6]; double fd1_denom_;
                diff_new_detail::get_fd1_6(c, o, fd1_denom_);
                const double d2 = fd1_denom_ * fd1_denom_;
                *fxy = diff_new_detail::cross_deriv<6, 0, 1>(fh, i, j, k, c, o, d2, 1.0/(dX*dY));
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (i+2 <= imax && i-2 >= imin && j+2 <= jmax && j-2 >= jmin) {
                const double Fdxdy = 1.0/144.0/(dX*dY);
                *fxy = Fdxdy*(   (fh(i-2,j-2,k)- 8*fh(i-1,j-2,k)+ 8*fh(i+1,j-2,k)-fh(i+2,j-2,k))
                            - 8*(fh(i-2,j-1,k)- 8*fh(i-1,j-1,k)+ 8*fh(i+1,j-1,k)-fh(i+2,j-1,k))
                            + 8*(fh(i-2,j+1,k)- 8*fh(i-1,j+1,k)+ 8*fh(i+1,j+1,k)-fh(i+2,j+1,k))
                            -   (fh(i-2,j+2,k)- 8*fh(i-1,j+2,k)+ 8*fh(i+1,j+2,k)-fh(i+2,j+2,k)));
                break;
            }
        }
        if (i+1 <= imax && i-1 >= imin && j+1 <= jmax && j-1 >= jmin) {
            const double Sdxdy = 0.25/(dX*dY);
            *fxy = Sdxdy*(fh(i-1,j-1,k) - fh(i+1,j-1,k) - fh(i-1,j+1,k) + fh(i+1,j+1,k));
        }
    } while(false);

    (void)onoff;
}

// --- d^2f/dxdz ---
template<int GW>
__device__ __forceinline__ void d_point_fddxz_sh_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fxz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW - 1;
    const double dX = X[1] - X[0];
    const double dZ = Z[1] - Z[0];
    const int imax = ex[0] - 1;
    const int kmax = ex[2] - 1;

    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SoA0, SoA1);
    };

    int imin = 0;
    // z has no symmetry: kmin=0 always
    if (symmetry == 2) { // OCTANT
        if (fabs(X[0]) < dX) imin = -ORD;
    }

    *fxz = 0.0;

    do {
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin && k+4 <= kmax && k-4 >= 0) {
                double c[8]; int o[8]; double fd1_denom_;
                diff_new_detail::get_fd1_8(c, o, fd1_denom_);
                const double d2 = fd1_denom_ * fd1_denom_;
                *fxz = diff_new_detail::cross_deriv<8, 0, 2>(fh, i, j, k, c, o, d2, 1.0/(dX*dZ));
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (i+3 <= imax && i-3 >= imin && k+3 <= kmax && k-3 >= 0) {
                double c[6]; int o[6]; double fd1_denom_;
                diff_new_detail::get_fd1_6(c, o, fd1_denom_);
                const double d2 = fd1_denom_ * fd1_denom_;
                *fxz = diff_new_detail::cross_deriv<6, 0, 2>(fh, i, j, k, c, o, d2, 1.0/(dX*dZ));
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (i+2 <= imax && i-2 >= imin && k+2 <= kmax && k-2 >= 0) {
                const double Fdxdz = 1.0/144.0/(dX*dZ);
                *fxz = Fdxdz*(   (fh(i-2,j,k-2)- 8*fh(i-1,j,k-2)+ 8*fh(i+1,j,k-2)-fh(i+2,j,k-2))
                            - 8*(fh(i-2,j,k-1)- 8*fh(i-1,j,k-1)+ 8*fh(i+1,j,k-1)-fh(i+2,j,k-1))
                            + 8*(fh(i-2,j,k+1)- 8*fh(i-1,j,k+1)+ 8*fh(i+1,j,k+1)-fh(i+2,j,k+1))
                            -   (fh(i-2,j,k+2)- 8*fh(i-1,j,k+2)+ 8*fh(i+1,j,k+2)-fh(i+2,j,k+2)));
                break;
            }
        }
        if (i+1 <= imax && i-1 >= imin && k+1 <= kmax && k-1 >= 0) {
            const double Sdxdz = 0.25/(dX*dZ);
            *fxz = Sdxdz*(fh(i-1,j,k-1) - fh(i+1,j,k-1) - fh(i-1,j,k+1) + fh(i+1,j,k+1));
        }
    } while(false);

    (void)onoff;
}

// --- d^2f/dydz ---
template<int GW>
__device__ __forceinline__ void d_point_fddyz_sh_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fyz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW - 1;
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];
    int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;

    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD, ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            f, ii, jj, kk, SoA0, SoA1);
    };

    int jmin = 0;
    // z has no symmetry: kmin=0 always
    if (symmetry == 2) { // OCTANT
        if (fabs(Y[0]) < dY) jmin = -ORD;
    } else if (symmetry == 1) { // EQ_SYMM
        if ((sst == 2 || sst == 4) && fabs(Y[0]) < dY) jmin = -ORD;
        if ((sst == 3 || sst == 5) && fabs(Y[ex[1]-1]) < dY) jmax = ex[1] - 1 + ORD;
    }

    *fyz = 0.0;

    do {
        if constexpr (GW >= 5) {
            if (j+4 <= jmax && j-4 >= jmin && k+4 <= kmax && k-4 >= 0) {
                double c[8]; int o[8]; double fd1_denom_;
                diff_new_detail::get_fd1_8(c, o, fd1_denom_);
                const double d2 = fd1_denom_ * fd1_denom_;
                *fyz = diff_new_detail::cross_deriv<8, 1, 2>(fh, i, j, k, c, o, d2, 1.0/(dY*dZ));
                break;
            }
        }
        if constexpr (GW >= 4) {
            if (j+3 <= jmax && j-3 >= jmin && k+3 <= kmax && k-3 >= 0) {
                double c[6]; int o[6]; double fd1_denom_;
                diff_new_detail::get_fd1_6(c, o, fd1_denom_);
                const double d2 = fd1_denom_ * fd1_denom_;
                *fyz = diff_new_detail::cross_deriv<6, 1, 2>(fh, i, j, k, c, o, d2, 1.0/(dY*dZ));
                break;
            }
        }
        if constexpr (GW >= 3) {
            if (j+2 <= jmax && j-2 >= jmin && k+2 <= kmax && k-2 >= 0) {
                const double Fdydz = 1.0/144.0/(dY*dZ);
                *fyz = Fdydz*(   (fh(i,j-2,k-2)- 8*fh(i,j-1,k-2)+ 8*fh(i,j+1,k-2)-fh(i,j+2,k-2))
                            - 8*(fh(i,j-2,k-1)- 8*fh(i,j-1,k-1)+ 8*fh(i,j+1,k-1)-fh(i,j+2,k-1))
                            + 8*(fh(i,j-2,k+1)- 8*fh(i,j-1,k+1)+ 8*fh(i,j+1,k+1)-fh(i,j+2,k+1))
                            -   (fh(i,j-2,k+2)- 8*fh(i,j-1,k+2)+ 8*fh(i,j+1,k+2)-fh(i,j+2,k+2)));
                break;
            }
        }
        if (j+1 <= jmax && j-1 >= jmin && k+1 <= kmax && k-1 >= 0) {
            const double Sdydz = 0.25/(dY*dZ);
            *fyz = Sdydz*(fh(i,j-1,k-1) - fh(i,j+1,k-1) - fh(i,j-1,k+1) + fh(i,j+1,k+1));
        }
    } while(false);

    (void)onoff;
}

// ============================================================
// Non-templated convenience wrappers for single-direction
// ============================================================

__device__ __forceinline__ void d_point_fdx_sh(
    const int ex[3], const int pex[3], const double* f,
    double* fx,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    d_point_fdx_sh_gw<ghost_width>(ex, pex, f, fx,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);
}

__device__ __forceinline__ void d_point_fdy_sh(
    const int ex[3], const int pex[3], const double* f,
    double* fy,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    d_point_fdy_sh_gw<ghost_width>(ex, pex, f, fy,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);
}

__device__ __forceinline__ void d_point_fdz_sh(
    const int ex[3], const int pex[3], const double* f,
    double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    d_point_fdz_sh_gw<ghost_width>(ex, pex, f, fz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);
}

__device__ __forceinline__ void d_point_fddxx_sh(
    const int ex[3], const int pex[3], const double* f,
    double* fxx,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    d_point_fddxx_sh_gw<ghost_width>(ex, pex, f, fxx,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);
}

__device__ __forceinline__ void d_point_fddyy_sh(
    const int ex[3], const int pex[3], const double* f,
    double* fyy,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    d_point_fddyy_sh_gw<ghost_width>(ex, pex, f, fyy,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);
}

__device__ __forceinline__ void d_point_fddzz_sh(
    const int ex[3], const int pex[3], const double* f,
    double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    d_point_fddzz_sh_gw<ghost_width>(ex, pex, f, fzz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);
}

__device__ __forceinline__ void d_point_fddxy_sh(
    const int ex[3], const int pex[3], const double* f,
    double* fxy,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    d_point_fddxy_sh_gw<ghost_width>(ex, pex, f, fxy,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);
}

__device__ __forceinline__ void d_point_fddxz_sh(
    const int ex[3], const int pex[3], const double* f,
    double* fxz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    d_point_fddxz_sh_gw<ghost_width>(ex, pex, f, fxz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);
}

__device__ __forceinline__ void d_point_fddyz_sh(
    const int ex[3], const int pex[3], const double* f,
    double* fyz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    d_point_fddyz_sh_gw<ghost_width>(ex, pex, f, fyz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);
}

// ============================================================
// Coordinate-transformed first derivatives (templated on GW)
// ============================================================

template<int GW>
__device__ __forceinline__ void d_point_fderivs_shc_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fx, double* fy, double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    const double* drhodx,  const double* drhody,  const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx,    const double* dRdy,    const double* dRdz,
    int i, int j, int k)
{
    double gx, gy, gz;
    d_fderivs_sh_point_gw<GW>(ex, pex, f, &gx, &gy, &gz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);

    int idx = (k * pex[1] + j) * pex[0] + i;

    *fx = dRdx[idx]*gz + drhodx[idx]*gx + dsigmadx[idx]*gy;
    *fy = dRdy[idx]*gz + drhody[idx]*gx + dsigmady[idx]*gy;
    *fz = dRdz[idx]*gz + drhodz[idx]*gx + dsigmadz[idx]*gy;
}

// ============================================================
// Coordinate-transformed second derivatives (templated on GW)
// ============================================================

template<int GW>
__device__ __forceinline__ void d_point_fdderivs_shc_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    const double* drhodx,  const double* drhody,  const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx,    const double* dRdy,    const double* dRdz,
    const double* drhodxx, const double* drhodxy, const double* drhodxz,
    const double* drhodyy, const double* drhodyz, const double* drhodzz,
    const double* dsigmadxx, const double* dsigmadxy, const double* dsigmadxz,
    const double* dsigmadyy, const double* dsigmadyz, const double* dsigmadzz,
    const double* dRdxx,   const double* dRdxy,   const double* dRdxz,
    const double* dRdyy,   const double* dRdyz,   const double* dRdzz,
    int i, int j, int k)
{
    double gx, gy, gz;
    d_fderivs_sh_point_gw<GW>(ex, pex, f, &gx, &gy, &gz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);

    double gxx, gxy, gxz, gyy, gyz, gzz;
    d_fdderivs_sh_point_gw<GW>(ex, pex, f, &gxx, &gxy, &gxz, &gyy, &gyz, &gzz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst, i, j, k);

    int idx = (k * pex[1] + j) * pex[0] + i;

    *fxx = dRdxx[idx]*gz + drhodxx[idx]*gx + dsigmadxx[idx]*gy
         + dRdx[idx]*dRdx[idx]*gzz + drhodx[idx]*drhodx[idx]*gxx + dsigmadx[idx]*dsigmadx[idx]*gyy
         + 2.0*(dRdx[idx]*drhodx[idx]*gxz + dRdx[idx]*dsigmadx[idx]*gyz + drhodx[idx]*dsigmadx[idx]*gxy);

    *fyy = dRdyy[idx]*gz + drhodyy[idx]*gx + dsigmadyy[idx]*gy
         + dRdy[idx]*dRdy[idx]*gzz + drhody[idx]*drhody[idx]*gxx + dsigmady[idx]*dsigmady[idx]*gyy
         + 2.0*(dRdy[idx]*drhody[idx]*gxz + dRdy[idx]*dsigmady[idx]*gyz + drhody[idx]*dsigmady[idx]*gxy);

    *fzz = dRdzz[idx]*gz + drhodzz[idx]*gx + dsigmadzz[idx]*gy
         + dRdz[idx]*dRdz[idx]*gzz + drhodz[idx]*drhodz[idx]*gxx + dsigmadz[idx]*dsigmadz[idx]*gyy
         + 2.0*(dRdz[idx]*drhodz[idx]*gxz + dRdz[idx]*dsigmadz[idx]*gyz + drhodz[idx]*dsigmadz[idx]*gxy);

    *fxy = dRdxy[idx]*gz + drhodxy[idx]*gx + dsigmadxy[idx]*gy
         + dRdx[idx]*drhody[idx]*gxz + dRdx[idx]*dsigmady[idx]*gyz + drhodx[idx]*dsigmady[idx]*gxy
         + dRdy[idx]*drhodx[idx]*gxz + dRdy[idx]*dsigmadx[idx]*gyz + drhody[idx]*dsigmadx[idx]*gxy
         + dRdx[idx]*dRdy[idx]*gzz + drhodx[idx]*drhody[idx]*gxx + dsigmadx[idx]*dsigmady[idx]*gyy;

    *fxz = dRdxz[idx]*gz + drhodxz[idx]*gx + dsigmadxz[idx]*gy
         + dRdx[idx]*drhodz[idx]*gxz + dRdx[idx]*dsigmadz[idx]*gyz + drhodx[idx]*dsigmadz[idx]*gxy
         + dRdz[idx]*drhodx[idx]*gxz + dRdz[idx]*dsigmadx[idx]*gyz + drhodz[idx]*dsigmadx[idx]*gxy
         + dRdx[idx]*dRdz[idx]*gzz + drhodx[idx]*drhodz[idx]*gxx + dsigmadx[idx]*dsigmadz[idx]*gyy;

    *fyz = dRdyz[idx]*gz + drhodyz[idx]*gx + dsigmadyz[idx]*gy
         + dRdz[idx]*drhody[idx]*gxz + dRdz[idx]*dsigmady[idx]*gyz + drhodz[idx]*dsigmady[idx]*gxy
         + dRdy[idx]*drhodz[idx]*gxz + dRdy[idx]*dsigmadz[idx]*gyz + drhody[idx]*dsigmadz[idx]*gxy
         + dRdz[idx]*dRdy[idx]*gzz + drhodz[idx]*drhody[idx]*gxx + dsigmadz[idx]*dsigmady[idx]*gyy;
}

// ============================================================
// Non-templated convenience wrappers for coordinate-transformed
// ============================================================

__device__ __forceinline__ void d_point_fderivs_shc(
    const int ex[3], const int pex[3], const double* f,
    double* fx, double* fy, double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    const double* drhodx,  const double* drhody,  const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx,    const double* dRdy,    const double* dRdz,
    int i, int j, int k)
{
    d_point_fderivs_shc_gw<ghost_width>(ex, pex, f, fx, fy, fz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        i, j, k);
}

__device__ __forceinline__ void d_point_fdderivs_shc(
    const int ex[3], const int pex[3], const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    const double* drhodx,  const double* drhody,  const double* drhodz,
    const double* dsigmadx, const double* dsigmady, const double* dsigmadz,
    const double* dRdx,    const double* dRdy,    const double* dRdz,
    const double* drhodxx, const double* drhodxy, const double* drhodxz,
    const double* drhodyy, const double* drhodyz, const double* drhodzz,
    const double* dsigmadxx, const double* dsigmadxy, const double* dsigmadxz,
    const double* dsigmadyy, const double* dsigmadyz, const double* dsigmadzz,
    const double* dRdxx,   const double* dRdxy,   const double* dRdxz,
    const double* dRdyy,   const double* dRdyz,   const double* dRdzz,
    int i, int j, int k)
{
    d_point_fdderivs_shc_gw<ghost_width>(ex, pex, f, fxx, fxy, fxz, fyy, fyz, fzz,
        X, Y, Z, SYM1, SYM2, SYM3, symmetry, onoff, sst,
        drhodx, drhody, drhodz, dsigmadx, dsigmady, dsigmadz, dRdx, dRdy, dRdz,
        drhodxx, drhodxy, drhodxz, drhodyy, drhodyz, drhodzz,
        dsigmadxx, dsigmadxy, dsigmadxz, dsigmadyy, dsigmadyz, dsigmadzz,
        dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz,
        i, j, k);
}

#endif // POINT_DIFF_NEW_SH_GPU_CUH
