
#ifndef RICCI_GAMMA_H
#define RICCI_GAMMA_H

#ifdef fortran1
#define f_ricci_gamma ricci_gamma
#define f_ricci_gamma_ss ricci_gamma_ss
#endif
#ifdef fortran2
#define f_ricci_gamma RICCI_GAMMA
#define f_ricci_gamma_ss RICCI_GAMMA_SS
#endif
#ifdef fortran3
#define f_ricci_gamma ricci_gamma_
#define f_ricci_gamma_ss ricci_gamma_ss_
#endif
extern "C"
{
        void f_ricci_gamma(int *, double *, double *, double *,
                           double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           int &);
}

extern "C"
{
        void f_ricci_gamma_ss(int *, double *, double *, double *, double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              int &, int &, int &);
}
#include "macrodef.h"
#ifdef USE_GPU
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

void gpu_ricci_gamma_ss_launch(
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
    int symmetry, int sst
);
#endif /* USE_GPU */

#endif /* RICCI_GAMMA_H */
