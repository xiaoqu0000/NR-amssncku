
#ifndef BSSN_GPU_H_
#define BSSN_GPU_H_
#include "bssn_macro.h"
#include "macrodef.fh"

#define DEVICE_ID 0
// #define DEVICE_ID_BY_MPI_RANK
#define GRID_DIM 256
#define BLOCK_DIM 128

#define _FH2_(i, j, k) fh[(i) + (j) * _1D_SIZE[2] + (k) * _2D_SIZE[2]]
#define _FH3_(i, j, k) fh[(i) + (j) * _1D_SIZE[3] + (k) * _2D_SIZE[3]]
#define pow2(x) ((x) * (x))
#define TimeBetween(a, b) ((b.tv_sec - a.tv_sec) + (b.tv_usec - a.tv_usec) / 1000000.0f)
#define M_ metac.
#define Mh_ meta->
#define Ms_ metassc.
#define Msh_ metass->

// #define TIMING

#define RHS_SS_PARA int calledby, int mpi_rank, int *ex, double &T, double *crho, double *sigma, double *R, double *X, double *Y, double *Z, double *drhodx, double *drhody, double *drhodz, double *dsigmadx, double *dsigmady, double *dsigmadz, double *dRdx, double *dRdy, double *dRdz, double *drhodxx, double *drhodxy, double *drhodxz, double *drhodyy, double *drhodyz, double *drhodzz, double *dsigmadxx, double *dsigmadxy, double *dsigmadxz, double *dsigmadyy, double *dsigmadyz, double *dsigmadzz, double *dRdxx, double *dRdxy, double *dRdxz, double *dRdyy, double *dRdyz, double *dRdzz, double *chi, double *trK, double *dxx, double *gxy, double *gxz, double *dyy, double *gyz, double *dzz, double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz, double *Gamx, double *Gamy, double *Gamz, double *Lap, double *betax, double *betay, double *betaz, double *dtSfx, double *dtSfy, double *dtSfz, double *chi_rhs, double *trK_rhs, double *gxx_rhs, double *gxy_rhs, double *gxz_rhs, double *gyy_rhs, double *gyz_rhs, double *gzz_rhs, double *Axx_rhs, double *Axy_rhs, double *Axz_rhs, double *Ayy_rhs, double *Ayz_rhs, double *Azz_rhs, double *Gamx_rhs, double *Gamy_rhs, double *Gamz_rhs, double *Lap_rhs, double *betax_rhs, double *betay_rhs, double *betaz_rhs, double *dtSfx_rhs, double *dtSfy_rhs, double *dtSfz_rhs, double *rho, double *Sx, double *Sy, double *Sz, double *Sxx, double *Sxy, double *Sxz, double *Syy, double *Syz, double *Szz, double *Gamxxx, double *Gamxxy, double *Gamxxz, double *Gamxyy, double *Gamxyz, double *Gamxzz, double *Gamyxx, double *Gamyxy, double *Gamyxz, double *Gamyyy, double *Gamyyz, double *Gamyzz, double *Gamzxx, double *Gamzxy, double *Gamzxz, double *Gamzyy, double *Gamzyz, double *Gamzzz, double *Rxx, double *Rxy, double *Rxz, double *Ryy, double *Ryz, double *Rzz, double *ham_Res, double *movx_Res, double *movy_Res, double *movz_Res, double *Gmx_Res, double *Gmy_Res, double *Gmz_Res, int &Symmetry, int &Lev, double &eps, int &sst, int &co

/**  main function */
int gpu_rhs(int calledby, int mpi_rank, int *ex, double &T,
            double *X, double *Y, double *Z,

            double *chi, double *trK,

            double *dxx, double *gxy, double *gxz, double *dyy, double *gyz, double *dzz,

            double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz,

            double *Gamx, double *Gamy, double *Gamz,

            double *Lap, double *betax, double *betay, double *betaz,

            double *dtSfx, double *dtSfy, double *dtSfz,

            double *chi_rhs, double *trK_rhs,

            double *gxx_rhs, double *gxy_rhs, double *gxz_rhs, double *gyy_rhs, double *gyz_rhs, double *gzz_rhs,

            double *Axx_rhs, double *Axy_rhs, double *Axz_rhs, double *Ayy_rhs, double *Ayz_rhs, double *Azz_rhs,

            double *Gamx_rhs, double *Gamy_rhs, double *Gamz_rhs,

            double *Lap_rhs, double *betax_rhs, double *betay_rhs, double *betaz_rhs,

            double *dtSfx_rhs, double *dtSfy_rhs, double *dtSfz_rhs,

            double *rho, double *Sx, double *Sy, double *Sz, double *Sxx,
            double *Sxy, double *Sxz, double *Syy, double *Syz, double *Szz,

            double *Gamxxx, double *Gamxxy, double *Gamxxz, double *Gamxyy, double *Gamxyz, double *Gamxzz,

            double *Gamyxx, double *Gamyxy, double *Gamyxz, double *Gamyyy, double *Gamyyz, double *Gamyzz,

            double *Gamzxx, double *Gamzxy, double *Gamzxz, double *Gamzyy, double *Gamzyz, double *Gamzzz,

            double *Rxx, double *Rxy, double *Rxz, double *Ryy, double *Ryz, double *Rzz,

            double *ham_Res, double *movx_Res, double *movy_Res, double *movz_Res,
            double *Gmx_Res, double *Gmy_Res, double *Gmz_Res,
            int &Symmetry, int &Lev, double &eps, int &co);

int gpu_rhs_ss(RHS_SS_PARA);

/** Init GPU side data in GPUMeta. */
// void init_fluid_meta_gpu(GPUMeta *gpu_meta);

#endif
