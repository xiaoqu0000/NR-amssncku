#ifndef DIFF_NEW_SH_GPU_CUH
#define DIFF_NEW_SH_GPU_CUH

#include "diff_new_gpu.cuh"
#include "fmisc_gpu.cuh"
#include "macrodef.fh"

#include <cmath>

#ifndef IDX3D_PAD
#define IDX3D_PAD(i, j, k, pnx, pny) ((i) + (pnx) * ((j) + (pny) * (k)))
#endif

// ============================================================
// Shell boundary accessor with upper+lower reflection
// (only x and y have ghost zones; z/radial has none)
// ============================================================
__device__ __forceinline__ double d_symmetry_stbd_0b_pad(
    int ord,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    const double* __restrict__ func,
    int i, int j, int k,
    double SoA0, double SoA1)
{
    // z (radial): no ghost zone
    if (k < 0 || k >= extc2) return 0.0;

    double factor = 1.0;

    // x lower boundary
    if (i < 0) {
        if (i < -ord) return 0.0;
        factor *= SoA0;
        i = -i - 1;
    }
    // x upper boundary
    if (i >= extc0) {
        if (i >= extc0 + ord) return 0.0;
        factor *= SoA0;
        i = 2 * extc0 - 1 - i;
    }

    // y lower boundary
    if (j < 0) {
        if (j < -ord) return 0.0;
        factor *= SoA1;
        j = -j - 1;
    }
    // y upper boundary
    if (j >= extc1) {
        if (j >= extc1 + ord) return 0.0;
        factor *= SoA1;
        j = 2 * extc1 - 1 - j;
    }

    // final bounds check
    if (i < 0 || i >= extc0 || j < 0 || j >= extc1) return 0.0;

    return func[(k * pextc1 + j) * pextc0 + i] * factor;
}

// ============================================================
// SoA pair selection helper for shell coordinates
// ============================================================
__device__ __forceinline__ void d_sh_get_SoA(
    int sst,
    double SYM1, double SYM2, double SYM3,
    double& SoA0, double& SoA1)
{
    if (sst == 0 || sst == 1) {
        SoA0 = SYM1;
        SoA1 = SYM2;
    } else if (sst == 2 || sst == 3) {
        SoA0 = SYM2;
        SoA1 = SYM3;
    } else { // sst == 4 || sst == 5
        SoA0 = SYM1;
        SoA1 = SYM3;
    }
}

// ============================================================
// First derivatives — shell coordinates, template on GW
// ============================================================
template<int GW>
__device__ __forceinline__ void d_fderivs_sh_point_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fx, double* fy, double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW - 1;

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    const int imax = ex[0] - 1;
    int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;

    // select SoA pair
    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD,
            ex[0], ex[1], ex[2],
            pex[0], pex[1], pex[2],
            f, ii, jj, kk,
            SoA0, SoA1
        );
    };

    // compute boundary limits (0-based, inclusive max index)
    const int NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2;
    (void)NO_SYMM;

    int imin = 0, jmin = 0;
    // kmin always 0, kmax always ex[2]-1 (no symmetry in z/radial)

    if (symmetry == OCTANT) {
        if (fabs(X[0]) < dX) imin = -ORD;
        if (fabs(Y[0]) < dY) jmin = -ORD;
    } else if (symmetry == EQ_SYMM) {
        if ((sst == 2 || sst == 4) && fabs(Y[0]) < dY) jmin = -ORD;
        if ((sst == 3 || sst == 5) && fabs(Y[ex[1]-1]) < dY) jmax = ex[1] - 1 + ORD;
    }

    *fx = 0.0; *fy = 0.0; *fz = 0.0;

    // ---------- x direction ----------
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

    // ---------- y direction ----------
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

    // ---------- z direction (kmin=0, kmax=ex[2], no symmetry) ----------
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
// Second derivatives — shell coordinates, template on GW
// ============================================================
template<int GW>
__device__ __forceinline__ void d_fdderivs_sh_point_gw(
    const int ex[3], const int pex[3], const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    constexpr int ORD = GW - 1;

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    const int imax = ex[0] - 1;
    int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;

    // select SoA pair
    double SoA0, SoA1;
    d_sh_get_SoA(sst, SYM1, SYM2, SYM3, SoA0, SoA1);

    const auto fh = [&](int ii, int jj, int kk) -> double {
        return d_symmetry_stbd_0b_pad(
            ORD,
            ex[0], ex[1], ex[2],
            pex[0], pex[1], pex[2],
            f, ii, jj, kk,
            SoA0, SoA1
        );
    };

    // compute boundary limits (0-based, inclusive max index)
    const int NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2;
    (void)NO_SYMM;

    int imin = 0, jmin = 0;
    // kmin always 0, kmax always ex[2]-1 (no symmetry in z/radial)

    if (symmetry == OCTANT) {
        if (fabs(X[0]) < dX) imin = -ORD;
        if (fabs(Y[0]) < dY) jmin = -ORD;
    } else if (symmetry == EQ_SYMM) {
        if ((sst == 2 || sst == 4) && fabs(Y[0]) < dY) jmin = -ORD;
        if ((sst == 3 || sst == 5) && fabs(Y[ex[1]-1]) < dY) jmax = ex[1] - 1 + ORD;
    }

    *fxx = 0.0; *fyy = 0.0; *fzz = 0.0;
    *fxy = 0.0; *fxz = 0.0; *fyz = 0.0;

    // ---------- fxx ----------
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

    // ---------- fyy ----------
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

    // ---------- fzz (kmin=0, kmax=ex[2]) ----------
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

    // ---------- fxy (cross: x*y) ----------
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

    // ---------- fxz (cross: x*z) ----------
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

    // ---------- fyz (cross: y*z) ----------
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
// Backward-compatible wrappers (default GW from macrodef.fh)
// ============================================================
__device__ __forceinline__ void d_fderivs_sh_point(
    const int ex[3], const int pex[3], const double* f,
    double* fx, double* fy, double* fz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    d_fderivs_sh_point_gw<ghost_width>(
        ex, pex, f, fx, fy, fz,
        X, Y, Z, SYM1, SYM2, SYM3,
        symmetry, onoff, sst, i, j, k
    );
}

__device__ __forceinline__ void d_fdderivs_sh_point(
    const int ex[3], const int pex[3], const double* f,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    const double* X, const double* Y, const double* Z,
    double SYM1, double SYM2, double SYM3,
    int symmetry, int onoff, int sst,
    int i, int j, int k)
{
    d_fdderivs_sh_point_gw<ghost_width>(
        ex, pex, f, fxx, fxy, fxz, fyy, fyz, fzz,
        X, Y, Z, SYM1, SYM2, SYM3,
        symmetry, onoff, sst, i, j, k
    );
}

#endif // DIFF_NEW_SH_GPU_CUH
