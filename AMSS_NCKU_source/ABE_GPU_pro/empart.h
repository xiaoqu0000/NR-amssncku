
#ifndef EMPART_H
#define EMPART_H

#ifdef fortran1
#define f_compute_rhs_empart compute_rhs_empart
#define f_compute_rhs_empart_ss compute_rhs_empart_ss
#endif
#ifdef fortran2
#define f_compute_rhs_empart COMPUTE_RHS_EMPART
#define f_compute_rhs_empart_ss COMPUTE_RHS_EMPART_SS
#endif
#ifdef fortran3
#define f_compute_rhs_empart compute_rhs_empart_
#define f_compute_rhs_empart_ss compute_rhs_empart_ss_
#endif

extern "C"
{
    int f_compute_rhs_empart(int *, double *, double *, double *,
                             double *, double *, double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                             int &, int &, double &);
}

extern "C"
{
    int f_compute_rhs_empart_ss(int *, double *, double *, double *, double *, double *, double *,
                                double *, double *, double *,
                                double *, double *, double *,
                                double *, double *, double *,
                                double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                                int &, int &, double &, int &);
}
#ifdef USE_GPU

void gpu_apply_symmetry_bd_empart_launch(
    cudaStream_t stream,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    double* Ex,  double* Ey,  double* Ez,
    double* Bx,  double* By,  double* Bz,
    double* Kpsi, double* Kphi
);

void gpu_compute_rhs_empart_launch(
    cudaStream_t &stream,
    int* ex, int* pex,
    double* X, double* Y, double* Z,
    double* chi, double* Jx, double* Jy, double* Jz, double* qchar,
    double* dxx, double* dxy, double* dxz, double* dyy, double* dyz, double* dzz,
    double* Lap, double* betax, double* betay, double* betaz, double* trK,
    double* Ex, double* Ey, double* Ez, double* Bx, double* By, double* Bz,
    double* Kpsi, double* Kphi,
    double* Ex_rhs, double* Ey_rhs, double* Ez_rhs,
    double* Bx_rhs, double* By_rhs, double* Bz_rhs,
    double* Kpsi_rhs, double* Kphi_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz,
    double* Syy, double* Syz, double* Szz,
    int symmetry, int lev, double eps
);

void gpu_compute_rhs_empart_ss_launch(
    cudaStream_t &stream,
    int* ex, int* pex,
    double* crho, double* sigma, double* R,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* chi, double* Jx, double* Jy, double* Jz, double* qchar,
    double* dxx, double* dxy, double* dxz, double* dyy, double* dyz, double* dzz,
    double* Lap, double* betax, double* betay, double* betaz, double* trK,
    double* Ex, double* Ey, double* Ez, double* Bx, double* By, double* Bz,
    double* Kpsi, double* Kphi,
    double* Ex_rhs, double* Ey_rhs, double* Ez_rhs,
    double* Bx_rhs, double* By_rhs, double* Bz_rhs,
    double* Kpsi_rhs, double* Kphi_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz,
    double* Syy, double* Syz, double* Szz,
    int symmetry, int lev, double eps, int sst
);
#endif /* USE_GPU */

#endif /* EMPART_H */
