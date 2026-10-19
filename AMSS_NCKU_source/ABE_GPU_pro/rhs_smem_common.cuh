#ifndef RHS_SMEM_COMMON_CUH
#define RHS_SMEM_COMMON_CUH

// Shared-memory + stencil helpers shared by every RHS GPU kernel
// (BSSN / BSSN-EM / BSSN-EScalar / Z4c).  Nothing here is BSSN-specific —
// only TileParams<GW>, load_*_to_smem, compute_*_derivs_smem and the tile
// block geometry constants.

#include "fmisc_gpu.cuh"
#include <cuda_pipeline.h>

constexpr int DIMX = 8;
constexpr int DIMY = 8;
constexpr int DIMZ = 4;
constexpr int BLOCK_SIZE = DIMX * DIMY * DIMZ;

// ============================================================
// TileParams<GW>: compile-time tile geometry parameterized on
// ghost width (valid values 2-5, corresponding to 2nd-8th order).
// ============================================================
template<int GW>
struct TileParams {
    static_assert(GW >= 2 && GW <= 5, "GW must be between 2 and 5");
    static constexpr int HALO = GW - 1;
    static constexpr int X = DIMX + 2 * HALO;
    static constexpr int Y = DIMY + 2 * HALO;
    static constexpr int Z = DIMZ + 2 * HALO;
    static constexpr int XY = X * Y;
    static constexpr int SIZE = X * Y * Z;
};

// ============================================================
// Shared-memory load helpers (templated on GW)
// ============================================================

template<int GW>
__device__ __forceinline__ void load_3_fields_to_smem(
    double smem[3][TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X],
    const double* fA, const double* fB, const double* fC,
    const int pex[3], int block_i, int block_j, int block_k
) {
    constexpr int HALO = GW - 1;
    constexpr int TX = TileParams<GW>::X;
    constexpr int TXY = TileParams<GW>::XY;
    constexpr int TSZ = TileParams<GW>::SIZE;

    int tid = threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;
    for (int idx = tid; idx < TSZ; idx += BLOCK_SIZE) {
        int loc_k = idx / TXY;
        int rem   = idx % TXY;
        int loc_j = rem / TX;
        int loc_i = rem % TX;

        int glob_i = block_i + loc_i - HALO;
        int glob_j = block_j + loc_j - HALO;
        int glob_k = block_k + loc_k - HALO;

        int gm_idx = (glob_k * pex[1] + glob_j) * pex[0] + glob_i;
        smem[0][loc_k][loc_j][loc_i] = fA[gm_idx];
        smem[1][loc_k][loc_j][loc_i] = fB[gm_idx];
        smem[2][loc_k][loc_j][loc_i] = fC[gm_idx];
    }
}

template<int GW>
__device__ __forceinline__ void load_field_to_smem(
    double smem[TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X],
    const double* f,
    const int ex[3], const int pex[3],
    int block_i, int block_j, int block_k,
    double SYM1, double SYM2, double SYM3
) {
    constexpr int HALO = GW - 1;
    constexpr int TX = TileParams<GW>::X;
    constexpr int TXY = TileParams<GW>::XY;
    constexpr int TSZ = TileParams<GW>::SIZE;

    int tid = threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;
    for (int idx = tid; idx < TSZ; idx += BLOCK_SIZE) {
        int loc_k = idx / TXY;
        int rem   = idx % TXY;
        int loc_j = rem / TX;
        int loc_i = rem % TX;

        int glob_i = block_i + loc_i - HALO;
        int glob_j = block_j + loc_j - HALO;
        int glob_k = block_k + loc_k - HALO;

        smem[loc_k][loc_j][loc_i] = f[(glob_k * pex[1] + glob_j) * pex[0] + glob_i];
    }
}

// ============================================================
// compute_all_derivs_smem<GW>
//   Computes all first and second derivatives from shared memory.
//   GW=2: per-direction boundary check (matches Fortran GW=2 per-direction code).
//   GW=3: combined boundary check (matches Fortran GW=3 "for bam comparison").
//   GW >= 4: per-direction boundary check (matches Fortran per-direction code)
//            where each axis independently selects its stencil order.
//   GW=3 produces EXACTLY the same code as the old non-templated
//   version (zero regression).
// ============================================================
template<int GW>
__device__ __forceinline__ void compute_all_derivs_smem(
    const double smem[TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X],
    double* fx, double* fy, double* fz,
    double* fxx, double* fxy, double* fxz,
    double* fyy, double* fyz, double* fzz,
    double dX, double dY, double dZ,
    int i, int j, int k,
    int imin, int jmin, int kmin,
    int imax, int jmax, int kmax
) {
    constexpr int HALO = GW - 1;

    *fx = 0.0; *fy = 0.0; *fz = 0.0;
    *fxx = 0.0; *fyy = 0.0; *fzz = 0.0;
    *fxy = 0.0; *fxz = 0.0; *fyz = 0.0;

    if (i >= imax || j >= jmax || k >= kmax) return;

    int c_i = threadIdx.x + HALO;
    int c_j = threadIdx.y + HALO;
    int c_k = threadIdx.z + HALO;

    auto fh = [&](int di, int dj, int dk) -> double {
        return smem[c_k + dk][c_j + dj][c_i + di];
    };

    if constexpr (GW == 2) {
        // ===== Per-direction check (matches Fortran GW=2 per-direction code) =====
        bool fit_x = (i + 1 <= imax && i - 1 >= imin);
        bool fit_y = (j + 1 <= jmax && j - 1 >= jmin);
        bool fit_z = (k + 1 <= kmax && k - 1 >= kmin);

        const double d2dx = 1.0 / 2.0 / dX;
        const double d2dy = 1.0 / 2.0 / dY;
        const double d2dz = 1.0 / 2.0 / dZ;

        const double Sdxdx = 1.0 / (dX * dX);
        const double Sdydy = 1.0 / (dY * dY);
        const double Sdzdz = 1.0 / (dZ * dZ);
        const double Sdxdy = 0.25 / (dX * dY);
        const double Sdxdz = 0.25 / (dX * dZ);
        const double Sdydz = 0.25 / (dY * dZ);

        if (fit_x) {
            *fx = d2dx * (-fh(-1,0,0) + fh(1,0,0));
            *fxx = Sdxdx * (fh(-1,0,0) - 2.0*fh(0,0,0) + fh(1,0,0));
        }
        if (fit_y) {
            *fy = d2dy * (-fh(0,-1,0) + fh(0,1,0));
            *fyy = Sdydy * (fh(0,-1,0) - 2.0*fh(0,0,0) + fh(0,1,0));
        }
        if (fit_z) {
            *fz = d2dz * (-fh(0,0,-1) + fh(0,0,1));
            *fzz = Sdzdz * (fh(0,0,-1) - 2.0*fh(0,0,0) + fh(0,0,1));
        }
        if (fit_x && fit_y) {
            *fxy = Sdxdy * (fh(-1,-1,0) - fh(1,-1,0) - fh(-1,1,0) + fh(1,1,0));
        }
        if (fit_x && fit_z) {
            *fxz = Sdxdz * (fh(-1,0,-1) - fh(1,0,-1) - fh(-1,0,1) + fh(1,0,1));
        }
        if (fit_y && fit_z) {
            *fyz = Sdydz * (fh(0,-1,-1) - fh(0,1,-1) - fh(0,-1,1) + fh(0,1,1));
        }
    } else if constexpr (GW == 3) {
        // ===== Combined check (matches Fortran GW=3 "for bam comparison") =====

        // 4th order
        {
            bool fit_4th = (i + 2 <= imax && i - 2 >= imin && j + 2 <= jmax && j - 2 >= jmin && k + 2 <= kmax && k - 2 >= kmin);
            if (fit_4th) {
                const double d12dx = 1.0 / 12.0 / dX;
                const double d12dy = 1.0 / 12.0 / dY;
                const double d12dz = 1.0 / 12.0 / dZ;

                const double Fdxdx = (1.0 / 12.0) / (dX * dX);
                const double Fdydy = (1.0 / 12.0) / (dY * dY);
                const double Fdzdz = (1.0 / 12.0) / (dZ * dZ);
                const double Fdxdy = (1.0 / 144.0) / (dX * dY);
                const double Fdxdz = (1.0 / 144.0) / (dX * dZ);
                const double Fdydz = (1.0 / 144.0) / (dY * dZ);

                *fx = d12dx * (fh(-2,0,0) - 8.0*fh(-1,0,0) + 8.0*fh(1,0,0) - fh(2,0,0));
                *fy = d12dy * (fh(0,-2,0) - 8.0*fh(0,-1,0) + 8.0*fh(0,1,0) - fh(0,2,0));
                *fz = d12dz * (fh(0,0,-2) - 8.0*fh(0,0,-1) + 8.0*fh(0,0,1) - fh(0,0,2));

                *fxx = Fdxdx * (-fh(-2,0,0) + 16.0*fh(-1,0,0) - 30.0*fh(0,0,0) - fh(2,0,0) + 16.0*fh(1,0,0));
                *fyy = Fdydy * (-fh(0,-2,0) + 16.0*fh(0,-1,0) - 30.0*fh(0,0,0) - fh(0,2,0) + 16.0*fh(0,1,0));
                *fzz = Fdzdz * (-fh(0,0,-2) + 16.0*fh(0,0,-1) - 30.0*fh(0,0,0) - fh(0,0,2) + 16.0*fh(0,0,1));

                *fxy = Fdxdy * (    (fh(-2,-2,0) - 8.0*fh(-1,-2,0) + 8.0*fh(1,-2,0) - fh(2,-2,0))
                                -8.0*(fh(-2,-1,0) - 8.0*fh(-1,-1,0) + 8.0*fh(1,-1,0) - fh(2,-1,0))
                                +8.0*(fh(-2, 1,0) - 8.0*fh(-1, 1,0) + 8.0*fh(1, 1,0) - fh(2, 1,0))
                                -   (fh(-2, 2,0) - 8.0*fh(-1, 2,0) + 8.0*fh(1, 2,0) - fh(2, 2,0)) );
                *fxz = Fdxdz * (    (fh(-2,0,-2) - 8.0*fh(-1,0,-2) + 8.0*fh(1,0,-2) - fh(2,0,-2))
                                -8.0*(fh(-2,0,-1) - 8.0*fh(-1,0,-1) + 8.0*fh(1,0,-1) - fh(2,0,-1))
                                +8.0*(fh(-2,0, 1) - 8.0*fh(-1,0, 1) + 8.0*fh(1,0, 1) - fh(2,0, 1))
                                -   (fh(-2,0, 2) - 8.0*fh(-1,0, 2) + 8.0*fh(1,0, 2) - fh(2,0, 2)) );
                *fyz = Fdydz * (    (fh(0,-2,-2) - 8.0*fh(0,-1,-2) + 8.0*fh(0,1,-2) - fh(0,2,-2))
                                -8.0*(fh(0,-2,-1) - 8.0*fh(0,-1,-1) + 8.0*fh(0,1,-1) - fh(0,2,-1))
                                +8.0*(fh(0,-2, 1) - 8.0*fh(0,-1, 1) + 8.0*fh(0,1, 1) - fh(0,2, 1))
                                -   (fh(0,-2, 2) - 8.0*fh(0,-1, 2) + 8.0*fh(0,1, 2) - fh(0,2, 2)) );

                return;
            }
        }

        // 2nd order fallback
        {
            bool fit_2nd = (i + 1 <= imax && i - 1 >= imin && j + 1 <= jmax && j - 1 >= jmin && k + 1 <= kmax && k - 1 >= kmin);

            const double d2dx = 1.0 / 2.0 / dX;
            const double d2dy = 1.0 / 2.0 / dY;
            const double d2dz = 1.0 / 2.0 / dZ;

            const double Sdxdx = 1.0 / (dX * dX);
            const double Sdydy = 1.0 / (dY * dY);
            const double Sdzdz = 1.0 / (dZ * dZ);
            const double Sdxdy = 0.25 / (dX * dY);
            const double Sdxdz = 0.25 / (dX * dZ);
            const double Sdydz = 0.25 / (dY * dZ);

            if (fit_2nd) {
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
        }
    } else {
        // ===== Per-direction check (matches Fortran GW>=4 per-direction code) =====

        // Determine per-axis stencil order: 0=none, 1=2nd, 2=4th, 3=6th, 4=8th
        int ox = 0, oy = 0, oz = 0;
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin) ox = 4;
            if (j+4 <= jmax && j-4 >= jmin) oy = 4;
            if (k+4 <= kmax && k-4 >= kmin) oz = 4;
        }
        if (ox == 0 && i+3 <= imax && i-3 >= imin) ox = 3;
        if (oy == 0 && j+3 <= jmax && j-3 >= jmin) oy = 3;
        if (oz == 0 && k+3 <= kmax && k-3 >= kmin) oz = 3;
        if (ox == 0 && i+2 <= imax && i-2 >= imin) ox = 2;
        if (oy == 0 && j+2 <= jmax && j-2 >= jmin) oy = 2;
        if (oz == 0 && k+2 <= kmax && k-2 >= kmin) oz = 2;
        if (ox == 0 && i+1 <= imax && i-1 >= imin) ox = 1;
        if (oy == 0 && j+1 <= jmax && j-1 >= jmin) oy = 1;
        if (oz == 0 && k+1 <= kmax && k-1 >= kmin) oz = 1;

        // For each derivative, compute using the appropriate order for its axis.
        // if constexpr eliminates dead branches for lower GW.

        // ---- fx (order = ox) ----
        if (ox == 1) {
            *fx = (1.0 / 2.0 / dX) * (-fh(-1,0,0) + fh(1,0,0));
        } else if (ox == 2) {
            *fx = (1.0 / 12.0 / dX) * (fh(-2,0,0) - 8.0*fh(-1,0,0) + 8.0*fh(1,0,0) - fh(2,0,0));
        }
        if constexpr (GW >= 4) {
            if (ox == 3) {
                *fx = (1.0 / 60.0 / dX) * (-fh(-3,0,0) + 9.0*fh(-2,0,0) - 45.0*fh(-1,0,0)
                       + 45.0*fh(1,0,0) - 9.0*fh(2,0,0) + fh(3,0,0));
            }
        }
        if constexpr (GW >= 5) {
            if (ox == 4) {
                *fx = (1.0 / 840.0 / dX) * (3.0*fh(-4,0,0) - 32.0*fh(-3,0,0) + 168.0*fh(-2,0,0) - 672.0*fh(-1,0,0)
                       + 672.0*fh(1,0,0) - 168.0*fh(2,0,0) + 32.0*fh(3,0,0) - 3.0*fh(4,0,0));
            }
        }

        // ---- fy (order = oy) ----
        if (oy == 1) {
            *fy = (1.0 / 2.0 / dY) * (-fh(0,-1,0) + fh(0,1,0));
        } else if (oy == 2) {
            *fy = (1.0 / 12.0 / dY) * (fh(0,-2,0) - 8.0*fh(0,-1,0) + 8.0*fh(0,1,0) - fh(0,2,0));
        }
        if constexpr (GW >= 4) {
            if (oy == 3) {
                *fy = (1.0 / 60.0 / dY) * (-fh(0,-3,0) + 9.0*fh(0,-2,0) - 45.0*fh(0,-1,0)
                       + 45.0*fh(0,1,0) - 9.0*fh(0,2,0) + fh(0,3,0));
            }
        }
        if constexpr (GW >= 5) {
            if (oy == 4) {
                *fy = (1.0 / 840.0 / dY) * (3.0*fh(0,-4,0) - 32.0*fh(0,-3,0) + 168.0*fh(0,-2,0) - 672.0*fh(0,-1,0)
                       + 672.0*fh(0,1,0) - 168.0*fh(0,2,0) + 32.0*fh(0,3,0) - 3.0*fh(0,4,0));
            }
        }

        // ---- fz (order = oz) ----
        if (oz == 1) {
            *fz = (1.0 / 2.0 / dZ) * (-fh(0,0,-1) + fh(0,0,1));
        } else if (oz == 2) {
            *fz = (1.0 / 12.0 / dZ) * (fh(0,0,-2) - 8.0*fh(0,0,-1) + 8.0*fh(0,0,1) - fh(0,0,2));
        }
        if constexpr (GW >= 4) {
            if (oz == 3) {
                *fz = (1.0 / 60.0 / dZ) * (-fh(0,0,-3) + 9.0*fh(0,0,-2) - 45.0*fh(0,0,-1)
                       + 45.0*fh(0,0,1) - 9.0*fh(0,0,2) + fh(0,0,3));
            }
        }
        if constexpr (GW >= 5) {
            if (oz == 4) {
                *fz = (1.0 / 840.0 / dZ) * (3.0*fh(0,0,-4) - 32.0*fh(0,0,-3) + 168.0*fh(0,0,-2) - 672.0*fh(0,0,-1)
                       + 672.0*fh(0,0,1) - 168.0*fh(0,0,2) + 32.0*fh(0,0,3) - 3.0*fh(0,0,4));
            }
        }

        // ---- fxx (order = ox) ----
        if (ox == 1) {
            *fxx = (1.0 / (dX*dX)) * (fh(-1,0,0) - 2.0*fh(0,0,0) + fh(1,0,0));
        } else if (ox == 2) {
            *fxx = ((1.0/12.0)/(dX*dX)) * (-fh(-2,0,0) + 16.0*fh(-1,0,0) - 30.0*fh(0,0,0) - fh(2,0,0) + 16.0*fh(1,0,0));
        }
        if constexpr (GW >= 4) {
            if (ox == 3) {
                *fxx = (1.0/180.0/(dX*dX)) * ( 2.0*fh(-3,0,0) - 27.0*fh(-2,0,0) + 270.0*fh(-1,0,0) - 490.0*fh(0,0,0)
                        + 270.0*fh(1,0,0) - 27.0*fh(2,0,0) + 2.0*fh(3,0,0));
            }
        }
        if constexpr (GW >= 5) {
            if (ox == 4) {
                *fxx = (1.0/5040.0/(dX*dX)) * (- 9.0*fh(-4,0,0) + 128.0*fh(-3,0,0) - 1008.0*fh(-2,0,0) + 8064.0*fh(-1,0,0) - 14350.0*fh(0,0,0)
                        + 8064.0*fh(1,0,0) - 1008.0*fh(2,0,0) + 128.0*fh(3,0,0) - 9.0*fh(4,0,0));
            }
        }

        // ---- fyy (order = oy) ----
        if (oy == 1) {
            *fyy = (1.0 / (dY*dY)) * (fh(0,-1,0) - 2.0*fh(0,0,0) + fh(0,1,0));
        } else if (oy == 2) {
            *fyy = ((1.0/12.0)/(dY*dY)) * (-fh(0,-2,0) + 16.0*fh(0,-1,0) - 30.0*fh(0,0,0) - fh(0,2,0) + 16.0*fh(0,1,0));
        }
        if constexpr (GW >= 4) {
            if (oy == 3) {
                *fyy = (1.0/180.0/(dY*dY)) * ( 2.0*fh(0,-3,0) - 27.0*fh(0,-2,0) + 270.0*fh(0,-1,0) - 490.0*fh(0,0,0)
                        + 270.0*fh(0,1,0) - 27.0*fh(0,2,0) + 2.0*fh(0,3,0));
            }
        }
        if constexpr (GW >= 5) {
            if (oy == 4) {
                *fyy = (1.0/5040.0/(dY*dY)) * (- 9.0*fh(0,-4,0) + 128.0*fh(0,-3,0) - 1008.0*fh(0,-2,0) + 8064.0*fh(0,-1,0) - 14350.0*fh(0,0,0)
                        + 8064.0*fh(0,1,0) - 1008.0*fh(0,2,0) + 128.0*fh(0,3,0) - 9.0*fh(0,4,0));
            }
        }

        // ---- fzz (order = oz) ----
        if (oz == 1) {
            *fzz = (1.0 / (dZ*dZ)) * (fh(0,0,-1) - 2.0*fh(0,0,0) + fh(0,0,1));
        } else if (oz == 2) {
            *fzz = ((1.0/12.0)/(dZ*dZ)) * (-fh(0,0,-2) + 16.0*fh(0,0,-1) - 30.0*fh(0,0,0) - fh(0,0,2) + 16.0*fh(0,0,1));
        }
        if constexpr (GW >= 4) {
            if (oz == 3) {
                *fzz = (1.0/180.0/(dZ*dZ)) * ( 2.0*fh(0,0,-3) - 27.0*fh(0,0,-2) + 270.0*fh(0,0,-1) - 490.0*fh(0,0,0)
                        + 270.0*fh(0,0,1) - 27.0*fh(0,0,2) + 2.0*fh(0,0,3));
            }
        }
        if constexpr (GW >= 5) {
            if (oz == 4) {
                *fzz = (1.0/5040.0/(dZ*dZ)) * (- 9.0*fh(0,0,-4) + 128.0*fh(0,0,-3) - 1008.0*fh(0,0,-2) + 8064.0*fh(0,0,-1) - 14350.0*fh(0,0,0)
                        + 8064.0*fh(0,0,1) - 1008.0*fh(0,0,2) + 128.0*fh(0,0,3) - 9.0*fh(0,0,4));
            }
        }

        // ---- fxy (order = min(ox, oy)) ----
        int oxy = ox < oy ? ox : oy;
        if (oxy == 1) {
            *fxy = (0.25/(dX*dY)) * (fh(-1,-1,0) - fh(1,-1,0) - fh(-1,1,0) + fh(1,1,0));
        } else if (oxy == 2) {
            *fxy = ((1.0/144.0)/(dX*dY)) * (    (fh(-2,-2,0) - 8.0*fh(-1,-2,0) + 8.0*fh(1,-2,0) - fh(2,-2,0))
                            -8.0*(fh(-2,-1,0) - 8.0*fh(-1,-1,0) + 8.0*fh(1,-1,0) - fh(2,-1,0))
                            +8.0*(fh(-2, 1,0) - 8.0*fh(-1, 1,0) + 8.0*fh(1, 1,0) - fh(2, 1,0))
                            -   (fh(-2, 2,0) - 8.0*fh(-1, 2,0) + 8.0*fh(1, 2,0) - fh(2, 2,0)) );
        }
        if constexpr (GW >= 4) {
            if (oxy == 3) {
                *fxy = (1.0/3600.0/(dX*dY)) * (
                    // 6th order outer product: c_x = {-1,9,-45,45,-9,1} x c_y = {-1,9,-45,45,-9,1}
                      -1.0*(-fh(-3,-3,0) +  9.0*fh(-2,-3,0) - 45.0*fh(-1,-3,0) + 45.0*fh(1,-3,0) - 9.0*fh(2,-3,0) + fh(3,-3,0))
                    +  9.0*(-fh(-3,-2,0) +  9.0*fh(-2,-2,0) - 45.0*fh(-1,-2,0) + 45.0*fh(1,-2,0) - 9.0*fh(2,-2,0) + fh(3,-2,0))
                    - 45.0*(-fh(-3,-1,0) +  9.0*fh(-2,-1,0) - 45.0*fh(-1,-1,0) + 45.0*fh(1,-1,0) - 9.0*fh(2,-1,0) + fh(3,-1,0))
                    + 45.0*(-fh(-3, 1,0) +  9.0*fh(-2, 1,0) - 45.0*fh(-1, 1,0) + 45.0*fh(1, 1,0) - 9.0*fh(2, 1,0) + fh(3, 1,0))
                    -  9.0*(-fh(-3, 2,0) +  9.0*fh(-2, 2,0) - 45.0*fh(-1, 2,0) + 45.0*fh(1, 2,0) - 9.0*fh(2, 2,0) + fh(3, 2,0))
                    +     1.0*(-fh(-3, 3,0) +  9.0*fh(-2, 3,0) - 45.0*fh(-1, 3,0) + 45.0*fh(1, 3,0) - 9.0*fh(2, 3,0) + fh(3, 3,0))
                );
            }
        }
        if constexpr (GW >= 5) {
            if (oxy == 4) {
                *fxy = (1.0/705600.0/(dX*dY)) * (
                    // 8th order outer product: c_x = {3,-32,168,-672,672,-168,32,-3} x c_y = same
                      3.0*(  3.0*fh(-4,-4,0) - 32.0*fh(-3,-4,0) + 168.0*fh(-2,-4,0) - 672.0*fh(-1,-4,0)
                          + 672.0*fh( 1,-4,0) - 168.0*fh( 2,-4,0) +  32.0*fh( 3,-4,0) -   3.0*fh(4,-4,0))
                    - 32.0*(  3.0*fh(-4,-3,0) - 32.0*fh(-3,-3,0) + 168.0*fh(-2,-3,0) - 672.0*fh(-1,-3,0)
                          + 672.0*fh( 1,-3,0) - 168.0*fh( 2,-3,0) +  32.0*fh( 3,-3,0) -   3.0*fh(4,-3,0))
                    +168.0*(  3.0*fh(-4,-2,0) - 32.0*fh(-3,-2,0) + 168.0*fh(-2,-2,0) - 672.0*fh(-1,-2,0)
                          + 672.0*fh( 1,-2,0) - 168.0*fh( 2,-2,0) +  32.0*fh( 3,-2,0) -   3.0*fh(4,-2,0))
                    -672.0*(  3.0*fh(-4,-1,0) - 32.0*fh(-3,-1,0) + 168.0*fh(-2,-1,0) - 672.0*fh(-1,-1,0)
                          + 672.0*fh( 1,-1,0) - 168.0*fh( 2,-1,0) +  32.0*fh( 3,-1,0) -   3.0*fh(4,-1,0))
                    +672.0*(  3.0*fh(-4, 1,0) - 32.0*fh(-3, 1,0) + 168.0*fh(-2, 1,0) - 672.0*fh(-1, 1,0)
                          + 672.0*fh( 1, 1,0) - 168.0*fh( 2, 1,0) +  32.0*fh( 3, 1,0) -   3.0*fh(4, 1,0))
                    -168.0*(  3.0*fh(-4, 2,0) - 32.0*fh(-3, 2,0) + 168.0*fh(-2, 2,0) - 672.0*fh(-1, 2,0)
                          + 672.0*fh( 1, 2,0) - 168.0*fh( 2, 2,0) +  32.0*fh( 3, 2,0) -   3.0*fh(4, 2,0))
                    + 32.0*(  3.0*fh(-4, 3,0) - 32.0*fh(-3, 3,0) + 168.0*fh(-2, 3,0) - 672.0*fh(-1, 3,0)
                          + 672.0*fh( 1, 3,0) - 168.0*fh( 2, 3,0) +  32.0*fh( 3, 3,0) -   3.0*fh(4, 3,0))
                    -  3.0*(  3.0*fh(-4, 4,0) - 32.0*fh(-3, 4,0) + 168.0*fh(-2, 4,0) - 672.0*fh(-1, 4,0)
                          + 672.0*fh( 1, 4,0) - 168.0*fh( 2, 4,0) +  32.0*fh( 3, 4,0) -   3.0*fh(4, 4,0))
                );
            }
        }

        // ---- fxz (order = min(ox, oz)) ----
        int oxz = ox < oz ? ox : oz;
        if (oxz == 1) {
            *fxz = (0.25/(dX*dZ)) * (fh(-1,0,-1) - fh(1,0,-1) - fh(-1,0,1) + fh(1,0,1));
        } else if (oxz == 2) {
            *fxz = ((1.0/144.0)/(dX*dZ)) * (    (fh(-2,0,-2) - 8.0*fh(-1,0,-2) + 8.0*fh(1,0,-2) - fh(2,0,-2))
                            -8.0*(fh(-2,0,-1) - 8.0*fh(-1,0,-1) + 8.0*fh(1,0,-1) - fh(2,0,-1))
                            +8.0*(fh(-2,0, 1) - 8.0*fh(-1,0, 1) + 8.0*fh(1,0, 1) - fh(2,0, 1))
                            -   (fh(-2,0, 2) - 8.0*fh(-1,0, 2) + 8.0*fh(1,0, 2) - fh(2,0, 2)) );
        }
        if constexpr (GW >= 4) {
            if (oxz == 3) {
                *fxz = (1.0/3600.0/(dX*dZ)) * (
                      -1.0*(-fh(-3,0,-3) +  9.0*fh(-2,0,-3) - 45.0*fh(-1,0,-3) + 45.0*fh(1,0,-3) - 9.0*fh(2,0,-3) + fh(3,0,-3))
                    +  9.0*(-fh(-3,0,-2) +  9.0*fh(-2,0,-2) - 45.0*fh(-1,0,-2) + 45.0*fh(1,0,-2) - 9.0*fh(2,0,-2) + fh(3,0,-2))
                    - 45.0*(-fh(-3,0,-1) +  9.0*fh(-2,0,-1) - 45.0*fh(-1,0,-1) + 45.0*fh(1,0,-1) - 9.0*fh(2,0,-1) + fh(3,0,-1))
                    + 45.0*(-fh(-3,0, 1) +  9.0*fh(-2,0, 1) - 45.0*fh(-1,0, 1) + 45.0*fh(1,0, 1) - 9.0*fh(2,0, 1) + fh(3,0, 1))
                    -  9.0*(-fh(-3,0, 2) +  9.0*fh(-2,0, 2) - 45.0*fh(-1,0, 2) + 45.0*fh(1,0, 2) - 9.0*fh(2,0, 2) + fh(3,0, 2))
                    +     1.0*(-fh(-3,0, 3) +  9.0*fh(-2,0, 3) - 45.0*fh(-1,0, 3) + 45.0*fh(1,0, 3) - 9.0*fh(2,0, 3) + fh(3,0, 3))
                );
            }
        }
        if constexpr (GW >= 5) {
            if (oxz == 4) {
                *fxz = (1.0/705600.0/(dX*dZ)) * (
                      3.0*(  3.0*fh(-4,0,-4) - 32.0*fh(-3,0,-4) + 168.0*fh(-2,0,-4) - 672.0*fh(-1,0,-4)
                          + 672.0*fh( 1,0,-4) - 168.0*fh( 2,0,-4) +  32.0*fh( 3,0,-4) -   3.0*fh(4,0,-4))
                    - 32.0*(  3.0*fh(-4,0,-3) - 32.0*fh(-3,0,-3) + 168.0*fh(-2,0,-3) - 672.0*fh(-1,0,-3)
                          + 672.0*fh( 1,0,-3) - 168.0*fh( 2,0,-3) +  32.0*fh( 3,0,-3) -   3.0*fh(4,0,-3))
                    +168.0*(  3.0*fh(-4,0,-2) - 32.0*fh(-3,0,-2) + 168.0*fh(-2,0,-2) - 672.0*fh(-1,0,-2)
                          + 672.0*fh( 1,0,-2) - 168.0*fh( 2,0,-2) +  32.0*fh( 3,0,-2) -   3.0*fh(4,0,-2))
                    -672.0*(  3.0*fh(-4,0,-1) - 32.0*fh(-3,0,-1) + 168.0*fh(-2,0,-1) - 672.0*fh(-1,0,-1)
                          + 672.0*fh( 1,0,-1) - 168.0*fh( 2,0,-1) +  32.0*fh( 3,0,-1) -   3.0*fh(4,0,-1))
                    +672.0*(  3.0*fh(-4,0, 1) - 32.0*fh(-3,0, 1) + 168.0*fh(-2,0, 1) - 672.0*fh(-1,0, 1)
                          + 672.0*fh( 1,0, 1) - 168.0*fh( 2,0, 1) +  32.0*fh( 3,0, 1) -   3.0*fh(4,0, 1))
                    -168.0*(  3.0*fh(-4,0, 2) - 32.0*fh(-3,0, 2) + 168.0*fh(-2,0, 2) - 672.0*fh(-1,0, 2)
                          + 672.0*fh( 1,0, 2) - 168.0*fh( 2,0, 2) +  32.0*fh( 3,0, 2) -   3.0*fh(4,0, 2))
                    + 32.0*(  3.0*fh(-4,0, 3) - 32.0*fh(-3,0, 3) + 168.0*fh(-2,0, 3) - 672.0*fh(-1,0, 3)
                          + 672.0*fh( 1,0, 3) - 168.0*fh( 2,0, 3) +  32.0*fh( 3,0, 3) -   3.0*fh(4,0, 3))
                    -  3.0*(  3.0*fh(-4,0, 4) - 32.0*fh(-3,0, 4) + 168.0*fh(-2,0, 4) - 672.0*fh(-1,0, 4)
                          + 672.0*fh( 1,0, 4) - 168.0*fh( 2,0, 4) +  32.0*fh( 3,0, 4) -   3.0*fh(4,0, 4))
                );
            }
        }

        // ---- fyz (order = min(oy, oz)) ----
        int oyz = oy < oz ? oy : oz;
        if (oyz == 1) {
            *fyz = (0.25/(dY*dZ)) * (fh(0,-1,-1) - fh(0,1,-1) - fh(0,-1,1) + fh(0,1,1));
        } else if (oyz == 2) {
            *fyz = ((1.0/144.0)/(dY*dZ)) * (    (fh(0,-2,-2) - 8.0*fh(0,-1,-2) + 8.0*fh(0,1,-2) - fh(0,2,-2))
                            -8.0*(fh(0,-2,-1) - 8.0*fh(0,-1,-1) + 8.0*fh(0,1,-1) - fh(0,2,-1))
                            +8.0*(fh(0,-2, 1) - 8.0*fh(0,-1, 1) + 8.0*fh(0,1, 1) - fh(0,2, 1))
                            -   (fh(0,-2, 2) - 8.0*fh(0,-1, 2) + 8.0*fh(0,1, 2) - fh(0,2, 2)) );
        }
        if constexpr (GW >= 4) {
            if (oyz == 3) {
                *fyz = (1.0/3600.0/(dY*dZ)) * (
                      -1.0*(-fh(0,-3,-3) +  9.0*fh(0,-2,-3) - 45.0*fh(0,-1,-3) + 45.0*fh(0,1,-3) - 9.0*fh(0,2,-3) + fh(0,3,-3))
                    +  9.0*(-fh(0,-3,-2) +  9.0*fh(0,-2,-2) - 45.0*fh(0,-1,-2) + 45.0*fh(0,1,-2) - 9.0*fh(0,2,-2) + fh(0,3,-2))
                    - 45.0*(-fh(0,-3,-1) +  9.0*fh(0,-2,-1) - 45.0*fh(0,-1,-1) + 45.0*fh(0,1,-1) - 9.0*fh(0,2,-1) + fh(0,3,-1))
                    + 45.0*(-fh(0,-3, 1) +  9.0*fh(0,-2, 1) - 45.0*fh(0,-1, 1) + 45.0*fh(0,1, 1) - 9.0*fh(0,2, 1) + fh(0,3, 1))
                    -  9.0*(-fh(0,-3, 2) +  9.0*fh(0,-2, 2) - 45.0*fh(0,-1, 2) + 45.0*fh(0,1, 2) - 9.0*fh(0,2, 2) + fh(0,3, 2))
                    +     1.0*(-fh(0,-3, 3) +  9.0*fh(0,-2, 3) - 45.0*fh(0,-1, 3) + 45.0*fh(0,1, 3) - 9.0*fh(0,2, 3) + fh(0,3, 3))
                );
            }
        }
        if constexpr (GW >= 5) {
            if (oyz == 4) {
                *fyz = (1.0/705600.0/(dY*dZ)) * (
                      3.0*(  3.0*fh(0,-4,-4) - 32.0*fh(0,-3,-4) + 168.0*fh(0,-2,-4) - 672.0*fh(0,-1,-4)
                          + 672.0*fh(0, 1,-4) - 168.0*fh(0, 2,-4) +  32.0*fh(0, 3,-4) -   3.0*fh(0,4,-4))
                    - 32.0*(  3.0*fh(0,-4,-3) - 32.0*fh(0,-3,-3) + 168.0*fh(0,-2,-3) - 672.0*fh(0,-1,-3)
                          + 672.0*fh(0, 1,-3) - 168.0*fh(0, 2,-3) +  32.0*fh(0, 3,-3) -   3.0*fh(0,4,-3))
                    +168.0*(  3.0*fh(0,-4,-2) - 32.0*fh(0,-3,-2) + 168.0*fh(0,-2,-2) - 672.0*fh(0,-1,-2)
                          + 672.0*fh(0, 1,-2) - 168.0*fh(0, 2,-2) +  32.0*fh(0, 3,-2) -   3.0*fh(0,4,-2))
                    -672.0*(  3.0*fh(0,-4,-1) - 32.0*fh(0,-3,-1) + 168.0*fh(0,-2,-1) - 672.0*fh(0,-1,-1)
                          + 672.0*fh(0, 1,-1) - 168.0*fh(0, 2,-1) +  32.0*fh(0, 3,-1) -   3.0*fh(0,4,-1))
                    +672.0*(  3.0*fh(0,-4, 1) - 32.0*fh(0,-3, 1) + 168.0*fh(0,-2, 1) - 672.0*fh(0,-1, 1)
                          + 672.0*fh(0, 1, 1) - 168.0*fh(0, 2, 1) +  32.0*fh(0, 3, 1) -   3.0*fh(0,4, 1))
                    -168.0*(  3.0*fh(0,-4, 2) - 32.0*fh(0,-3, 2) + 168.0*fh(0,-2, 2) - 672.0*fh(0,-1, 2)
                          + 672.0*fh(0, 1, 2) - 168.0*fh(0, 2, 2) +  32.0*fh(0, 3, 2) -   3.0*fh(0,4, 2))
                    + 32.0*(  3.0*fh(0,-4, 3) - 32.0*fh(0,-3, 3) + 168.0*fh(0,-2, 3) - 672.0*fh(0,-1, 3)
                          + 672.0*fh(0, 1, 3) - 168.0*fh(0, 2, 3) +  32.0*fh(0, 3, 3) -   3.0*fh(0,4, 3))
                    -  3.0*(  3.0*fh(0,-4, 4) - 32.0*fh(0,-3, 4) + 168.0*fh(0,-2, 4) - 672.0*fh(0,-1, 4)
                          + 672.0*fh(0, 1, 4) - 168.0*fh(0, 2, 4) +  32.0*fh(0, 3, 4) -   3.0*fh(0,4, 4))
                );
            }
        }
    }
}

// ============================================================
// compute_first_derivs_smem<GW>
//   Computes only first derivatives (fx, fy, fz) from shared memory.
//   GW=2: per-direction boundary check (matches Fortran GW=2 per-direction code).
//   GW=3: combined boundary check (matches Fortran GW=3 "for bam comparison").
//   GW >= 4: per-direction boundary check (matches Fortran per-direction code)
//            where each axis independently selects its stencil order.
// ============================================================
template<int GW>
__device__ __forceinline__ void compute_first_derivs_smem(
    const double smem[TileParams<GW>::Z][TileParams<GW>::Y][TileParams<GW>::X],
    double* fx, double* fy, double* fz,
    double dX, double dY, double dZ,
    int i, int j, int k,
    int imin, int jmin, int kmin,
    int imax, int jmax, int kmax
) {
    constexpr int HALO = GW - 1;

    *fx = 0.0; *fy = 0.0; *fz = 0.0;

    if (i >= imax || j >= jmax || k >= kmax) return;

    int c_i = threadIdx.x + HALO;
    int c_j = threadIdx.y + HALO;
    int c_k = threadIdx.z + HALO;

    auto fh = [&](int di, int dj, int dk) -> double {
        return smem[c_k + dk][c_j + dj][c_i + di];
    };

    if constexpr (GW == 2) {
        // ===== Per-direction check (matches Fortran GW=2 per-direction code) =====
        bool fit_x = (i + 1 <= imax && i - 1 >= imin);
        bool fit_y = (j + 1 <= jmax && j - 1 >= jmin);
        bool fit_z = (k + 1 <= kmax && k - 1 >= kmin);

        const double d2dx = 1.0 / 2.0 / dX;
        const double d2dy = 1.0 / 2.0 / dY;
        const double d2dz = 1.0 / 2.0 / dZ;

        if (fit_x) *fx = d2dx * (-fh(-1,0,0) + fh(1,0,0));
        if (fit_y) *fy = d2dy * (-fh(0,-1,0) + fh(0,1,0));
        if (fit_z) *fz = d2dz * (-fh(0,0,-1) + fh(0,0,1));
    } else if constexpr (GW == 3) {
        // ===== Combined check (matches Fortran GW=3 "for bam comparison") =====

        // 4th order
        {
            bool fit_4th = (i + 2 <= imax && i - 2 >= imin && j + 2 <= jmax && j - 2 >= jmin && k + 2 <= kmax && k - 2 >= kmin);
            if (fit_4th) {
                const double d12dx = 1.0 / 12.0 / dX;
                const double d12dy = 1.0 / 12.0 / dY;
                const double d12dz = 1.0 / 12.0 / dZ;

                *fx = d12dx * (fh(-2,0,0) - 8.0*fh(-1,0,0) + 8.0*fh(1,0,0) - fh(2,0,0));
                *fy = d12dy * (fh(0,-2,0) - 8.0*fh(0,-1,0) + 8.0*fh(0,1,0) - fh(0,2,0));
                *fz = d12dz * (fh(0,0,-2) - 8.0*fh(0,0,-1) + 8.0*fh(0,0,1) - fh(0,0,2));
                return;
            }
        }

        // 2nd order fallback
        {
            bool fit_2nd = (i + 1 <= imax && i - 1 >= imin && j + 1 <= jmax && j - 1 >= jmin && k + 1 <= kmax && k - 1 >= kmin);

            const double d2dx = 1.0 / 2.0 / dX;
            const double d2dy = 1.0 / 2.0 / dY;
            const double d2dz = 1.0 / 2.0 / dZ;

            if (fit_2nd) {
                *fx = d2dx * (-fh(-1,0,0) + fh(1,0,0));
                *fy = d2dy * (-fh(0,-1,0) + fh(0,1,0));
                *fz = d2dz * (-fh(0,0,-1) + fh(0,0,1));
            }
        }
    } else {
        // ===== Per-direction check (matches Fortran GW>=4 per-direction code) =====

        // Determine per-axis stencil order: 0=none, 1=2nd, 2=4th, 3=6th, 4=8th
        int ox = 0, oy = 0, oz = 0;
        if constexpr (GW >= 5) {
            if (i+4 <= imax && i-4 >= imin) ox = 4;
            if (j+4 <= jmax && j-4 >= jmin) oy = 4;
            if (k+4 <= kmax && k-4 >= kmin) oz = 4;
        }
        if (ox == 0 && i+3 <= imax && i-3 >= imin) ox = 3;
        if (oy == 0 && j+3 <= jmax && j-3 >= jmin) oy = 3;
        if (oz == 0 && k+3 <= kmax && k-3 >= kmin) oz = 3;
        if (ox == 0 && i+2 <= imax && i-2 >= imin) ox = 2;
        if (oy == 0 && j+2 <= jmax && j-2 >= jmin) oy = 2;
        if (oz == 0 && k+2 <= kmax && k-2 >= kmin) oz = 2;
        if (ox == 0 && i+1 <= imax && i-1 >= imin) ox = 1;
        if (oy == 0 && j+1 <= jmax && j-1 >= jmin) oy = 1;
        if (oz == 0 && k+1 <= kmax && k-1 >= kmin) oz = 1;

        // ---- fx (order = ox) ----
        if (ox == 1) {
            *fx = (1.0 / 2.0 / dX) * (-fh(-1,0,0) + fh(1,0,0));
        } else if (ox == 2) {
            *fx = (1.0 / 12.0 / dX) * (fh(-2,0,0) - 8.0*fh(-1,0,0) + 8.0*fh(1,0,0) - fh(2,0,0));
        }
        if constexpr (GW >= 4) {
            if (ox == 3) {
                *fx = (1.0 / 60.0 / dX) * (-fh(-3,0,0) + 9.0*fh(-2,0,0) - 45.0*fh(-1,0,0)
                       + 45.0*fh(1,0,0) - 9.0*fh(2,0,0) + fh(3,0,0));
            }
        }
        if constexpr (GW >= 5) {
            if (ox == 4) {
                *fx = (1.0 / 840.0 / dX) * (3.0*fh(-4,0,0) - 32.0*fh(-3,0,0) + 168.0*fh(-2,0,0) - 672.0*fh(-1,0,0)
                       + 672.0*fh(1,0,0) - 168.0*fh(2,0,0) + 32.0*fh(3,0,0) - 3.0*fh(4,0,0));
            }
        }

        // ---- fy (order = oy) ----
        if (oy == 1) {
            *fy = (1.0 / 2.0 / dY) * (-fh(0,-1,0) + fh(0,1,0));
        } else if (oy == 2) {
            *fy = (1.0 / 12.0 / dY) * (fh(0,-2,0) - 8.0*fh(0,-1,0) + 8.0*fh(0,1,0) - fh(0,2,0));
        }
        if constexpr (GW >= 4) {
            if (oy == 3) {
                *fy = (1.0 / 60.0 / dY) * (-fh(0,-3,0) + 9.0*fh(0,-2,0) - 45.0*fh(0,-1,0)
                       + 45.0*fh(0,1,0) - 9.0*fh(0,2,0) + fh(0,3,0));
            }
        }
        if constexpr (GW >= 5) {
            if (oy == 4) {
                *fy = (1.0 / 840.0 / dY) * (3.0*fh(0,-4,0) - 32.0*fh(0,-3,0) + 168.0*fh(0,-2,0) - 672.0*fh(0,-1,0)
                       + 672.0*fh(0,1,0) - 168.0*fh(0,2,0) + 32.0*fh(0,3,0) - 3.0*fh(0,4,0));
            }
        }

        // ---- fz (order = oz) ----
        if (oz == 1) {
            *fz = (1.0 / 2.0 / dZ) * (-fh(0,0,-1) + fh(0,0,1));
        } else if (oz == 2) {
            *fz = (1.0 / 12.0 / dZ) * (fh(0,0,-2) - 8.0*fh(0,0,-1) + 8.0*fh(0,0,1) - fh(0,0,2));
        }
        if constexpr (GW >= 4) {
            if (oz == 3) {
                *fz = (1.0 / 60.0 / dZ) * (-fh(0,0,-3) + 9.0*fh(0,0,-2) - 45.0*fh(0,0,-1)
                       + 45.0*fh(0,0,1) - 9.0*fh(0,0,2) + fh(0,0,3));
            }
        }
        if constexpr (GW >= 5) {
            if (oz == 4) {
                *fz = (1.0 / 840.0 / dZ) * (3.0*fh(0,0,-4) - 32.0*fh(0,0,-3) + 168.0*fh(0,0,-2) - 672.0*fh(0,0,-1)
                       + 672.0*fh(0,0,1) - 168.0*fh(0,0,2) + 32.0*fh(0,0,3) - 3.0*fh(0,0,4));
            }
        }
    }
}

#endif // RHS_SMEM_COMMON_CUH
